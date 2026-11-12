function reconstruction = M87_EHT_Headless_GPU_Pipeline(CFG)
% =========================================================================
% HEADLESS GPU-ACCELERATED VLBI IMAGE RECONSTRUCTION  (v2)
% EHT M87* -- Stokes I imaging + visibility telemetry + polarimetry
%
% New in v2 (all additions are optional / fail-soft; the original Stokes I
% imaging path is unchanged):
%   - Reader also extracts: observation time, baseline (station pair) and
%     the AIPS AN antenna table (station names).
%   - Polarimetric files: if RR/LL/RL/LR planes exist, Stokes I is formed as
%     (RR+LL)/2 and a complex P = Q + iU image is gridded and CLEANed
%     (complex Hogbom) -> Q, U, |P|, debiased P, fractional pol, EVPA.
%   - Visibility telemetry for the dashboard:
%       visSample  : random sample of amplitude/phase/SNR/baseline/time and
%                    CLEAN-model amplitude at the same uv points
%       visRadial  : weighted-binned amplitude vs |uv| with model overlay
%       closure    : closure phases vs time for all available triangles
%   - Exports: stationNames, blNames, dataInfo.
%
% POLARIMETRY NOTES
%   Q + iU is taken as RL (and Q - iU as LR), so P(+u) = RL(u) and
%   P(-u) = conj(LR(u)). Use the D-term-calibrated files for polarimetric
%   imaging; the base "*_hops_ALMArot.uvfits" release has no D-term
%   calibration applied. The EVPA sign convention follows the file's
%   RL/LR definition; verify against a published EVPA before quoting.
% =========================================================================

    if nargin < 1
        % Default Configuration (Linux Paths)
        % Polarimetric release (2023-D01-01). For the 2018 Stokes I run use
        % '.../2024-D01-01-main' and the M87_EHT_FAST_RECONSTRUCTION.mat output.
        CFG.dataRoot = '/home/hellblazer/Documents/MATLAB/VLBI/2023-D01-01-main';
        % Release layout: hops_data/AprilXX/*.uvfits and casa_data/AprilXX/*.uvfits
        % Pick ONE variant per run (the folders are searched recursively):
        %   HOPS, D-term calibrated : hops_data + '*hops_zbl-dtcal+selfcal.uvfits'
        %   CASA, D-term calibrated : casa_data + '*casa_EVPA_rotation+dcal*.uvfits'
        %   HOPS, NO D-terms (avoid for pol) : '*hops_ALMArot.uvfits'
        CFG.uvfitsDir = fullfile(CFG.dataRoot, 'hops_data');
        CFG.fileFilter = '*hops_zbl-dtcal+selfcal.uvfits';
        CFG.Nuv = 2048;
        CFG.du_Mlambda = 20;
        CFG.kernelHalfWidth = 3;
        CFG.kernelSigma = 1.2;
        CFG.cleanGain = 0.08;
        CFG.cleanMaxIter = 5000;
        CFG.cleanSigma = 3.0;
        CFG.cleanRadius_mas = 0.20;
        CFG.cleanPatchRadius = 512;
        CFG.positiveClean = true;
        CFG.maxModelFluxJy = 10;
        CFG.outputMAT = '/home/hellblazer/Documents/MATLAB/VLBI/M87_EHT_POL_RECONSTRUCTION.mat';
    end

    % ---- v2 options (defaults only fill what is missing) ------------------
    CFG = cfgDefault(CFG, 'polEnable',          true);
    CFG = cfgDefault(CFG, 'polCleanMaxIter',    3000);
    CFG = cfgDefault(CFG, 'polCleanSigma',      4.0);
    CFG = cfgDefault(CFG, 'nRadialBins',        60);
    CFG = cfgDefault(CFG, 'visSampleN',         30000);
    CFG = cfgDefault(CFG, 'minTriangleSamples', 50);
    CFG = cfgDefault(CFG, 'maxClosurePoints',   60000);
    CFG = cfgDefault(CFG, 'fileFilter',         '*.uvfits'); % e.g. '*hops*dterm*.uvfits'
    CFG = cfgDefault(CFG, 'excludeRegex',       '');         % e.g. 'casa|nodterm' (case-sensitive)
    CFG = cfgDefault(CFG, 'evpaOffsetDeg',      0);          % add 90 if your EVPA map is rotated by 90 deg

    % ---- v2.1 low-load options ----------------------------------------------
    CFG = cfgDefault(CFG, 'computeMode',        'HYBRID');   % HYBRID = GPU FFTs + CPU CLEAN | CPU = no GPU at all
    CFG = cfgDefault(CFG, 'maxFiles',           Inf);        % set 2 for a quick, cool test run
    CFG = cfgDefault(CFG, 'thermalGuard',       true);       % pause/abort when too hot (Linux sensors)
    CFG = cfgDefault(CFG, 'maxTempC',           85);
    CFG = cfgDefault(CFG, 'maxCoolWaitSec',     300);
    CFG = cfgDefault(CFG, 'stagePauseSec',      2);          % short cool-down between stages

    totalTimer = tic;

    % ---- share the load: limit CPU threads, GPU only for FFTs ----------------
    try
        nCores = feature('numcores');
    catch
        nCores = 4;
    end
    CFG = cfgDefault(CFG, 'cpuThreads', max(1, floor(nCores / 2)));
    try
        wState = warning('off', 'all');
        maxNumCompThreads(CFG.cpuThreads);
        warning(wState);
    catch
    end

    useGPU = false;
    if ~strcmpi(CFG.computeMode, 'CPU')
        try
            useGPU = canUseGPU();
        catch
            useGPU = false;
        end
        if ~useGPU
            warning('GPU not detected. Falling back to CPU compute.');
        end
    end

    if useGPU
        computeLabel = sprintf('GPU FFT + CPU CLEAN, %d CPU threads', CFG.cpuThreads);
    else
        computeLabel = sprintf('CPU only, %d threads', CFG.cpuThreads);
    end
    fprintf('Compute mode: %s\n', computeLabel);
    fprintf('\n======================================================\n');
    fprintf('INITIATING HEADLESS GPU RECONSTRUCTION PIPELINE (v2)\n');
    fprintf('======================================================\n');

    %% ============================================================
    % STEP 1 & 2 — FIND AND READ UVFITS (CPU)
    % ============================================================
    fprintf('[1/9] Fetching and parsing UVFITS metadata...\n');
    ioTimer = tic;

    files = dir(fullfile(CFG.uvfitsDir, '**', CFG.fileFilter));
    files = files(~[files.isdir]);
    if isempty(files)
        files = dir(fullfile(CFG.uvfitsDir, '**', '*.UVFITS'));
        files = files(~[files.isdir]);
    end

    if ~isempty(CFG.excludeRegex) && ~isempty(files)
        drop = ~cellfun(@isempty, regexp({files.name}, CFG.excludeRegex, 'once'));
        files = files(~drop);
    end

    if isempty(files), error('No matching .uvfits files found in: %s', CFG.uvfitsDir); end

    fullNames = cellfun(@(f, n) fullfile(f, n), {files.folder}, {files.name}, ...
        'UniformOutput', false);
    [~, ord] = sort(fullNames);
    files = files(ord);
    files = files(1:min(numel(files), CFG.maxFiles));
    nFiles = numel(files);

    fprintf('      -> Merging %d file(s):\n', nFiles);
    for kf = 1:nFiles
        fprintf('         %s\n', files(kf).name);
    end

    U = []; V = []; W = []; VIS = []; WT = []; FREQ = [];
    TIME = []; S1 = []; S2 = []; FID = [];
    RL = []; LR = []; WRL = []; WLR = [];
    stationNames = {};

    for k = 1:nFiles
        filename = fullfile(files(k).folder, files(k).name);
        D = readRawUVFITS(filename);
        n = numel(D.visibility);

        [g1, g2, stationNames] = mapStations(D, stationNames);

        U = [U; D.uSec(:)]; %#ok<AGROW>
        V = [V; D.vSec(:)]; %#ok<AGROW>
        W = [W; D.wSec(:)]; %#ok<AGROW>
        VIS = [VIS; D.visibility(:)]; %#ok<AGROW>
        WT = [WT; D.weight(:)]; %#ok<AGROW>
        FREQ = [FREQ; repmat(D.frequencyHz, n, 1)]; %#ok<AGROW>
        TIME = [TIME; D.time(:)]; %#ok<AGROW>
        S1 = [S1; g1(:)]; %#ok<AGROW>
        S2 = [S2; g2(:)]; %#ok<AGROW>
        FID = [FID; k * ones(n, 1)]; %#ok<AGROW>

        if D.hasPol
            RL  = [RL;  D.RL(:)];  %#ok<AGROW>
            LR  = [LR;  D.LR(:)];  %#ok<AGROW>
            WRL = [WRL; D.WRL(:)]; %#ok<AGROW>
            WLR = [WLR; D.WLR(:)]; %#ok<AGROW>
        else
            RL  = [RL;  complex(zeros(n, 1))]; %#ok<AGROW>
            LR  = [LR;  complex(zeros(n, 1))]; %#ok<AGROW>
            WRL = [WRL; zeros(n, 1)]; %#ok<AGROW>
            WLR = [WLR; zeros(n, 1)]; %#ok<AGROW>
        end
    end

    valid = isfinite(U) & isfinite(V) & isfinite(W) & isfinite(real(VIS)) & ...
            isfinite(imag(VIS)) & isfinite(WT) & WT > 0 & isfinite(FREQ) & FREQ > 1e9;

    U = U(valid); V = V(valid); W = W(valid); VIS = VIS(valid); WT = WT(valid); FREQ = FREQ(valid);
    TIME = TIME(valid); S1 = S1(valid); S2 = S2(valid); FID = FID(valid);
    RL = RL(valid); LR = LR(valid); WRL = WRL(valid); WLR = WLR(valid);
    nMeasured = numel(VIS);

    hasBL  = all(S1 > 0 & S2 > 0) && numel(stationNames) >= 2;
    anyPol = CFG.polEnable && any(WRL > 0 & WLR > 0);

    uML = U .* FREQ / 1e6;
    vML = V .* FREQ / 1e6;
    fprintf('      -> Loaded %d valid visibility samples from %d files. (%.2f s)\n', ...
        nMeasured, nFiles, toc(ioTimer));
    fprintf('      -> Stations: %s\n', strjoin(stationNames, ', '));
    if anyPol
        fprintf('      -> Polarimetric RL/LR products detected.\n');
    end

    %% ============================================================
    % STEP 3 — HERMITIAN COMPLETION & UI PREP (CPU)
    % ============================================================
    fprintf('[2/9] Enforcing Hermitian symmetry...\n');
    uH = [uML; -uML];
    vH = [vML; -vML];
    visH = [VIS; conj(VIS)];
    wtH = [WT; WT];

    % Save a 25k point downsampled subset strictly for the UI Dashboard
    plotIdx = randperm(numel(uH), min(numel(uH), 25000));
    uvPlotU = uH(plotIdx);
    uvPlotV = vH(plotIdx);
    uvPlotVis = abs(visH(plotIdx));

    %% ============================================================
    % STEP 4 — UV GRIDDING (CPU)
    % ============================================================
    fprintf('[3/9] Convolutional UV Gridding...\n');
    gridTimer = tic;

    N = CFG.Nuv;
    du = CFG.du_Mlambda;
    c = N/2 + 1;

    [uvNum, uvWt] = gridVisibilities(uH, vH, visH, wtH, N, du, ...
        CFG.kernelHalfWidth, CFG.kernelSigma, @() thermalGuard(CFG, 'grid', true));

    totalWeight = sum(uvWt(:));
    fprintf('      -> Gridding completed. (%.2f s)\n', toc(gridTimer));

    % Coordinate Setup
    fovRad = 1 / (du * 1e6);
    fovMas = fovRad * (180/pi) * 3600 * 1000;
    pixelMas = fovMas / N;
    x = ((1:N) - c) * pixelMas;
    [XX, YY] = meshgrid(x, x);

    %% ============================================================
    % STEP 5 — DIRTY IMAGE & PSF (FFTs on GPU if available, results kept on CPU)
    % ============================================================
    thermalGuard(CFG, 'before FFT', false);
    fprintf('[4/9] FFTs (%s)...\n', computeLabel);
    gpuTimer = tic;

    dirtyComplex = fftImage(uvNum, useGPU);
    beamComplex  = fftImage(uvWt, useGPU);
    clear uvNum uvWt

    dirty = real(dirtyComplex) / totalWeight;
    beamPeak = real(beamComplex(c, c));
    psf = real(beamComplex / beamPeak);
    % exact point symmetry about the centre pixel (c,c) for an even-sized grid
    psf = 0.5 * (psf + circshift(rot90(psf, 2), [1 1]));
    psf = psf / max(abs(psf(:)));
    clear dirtyComplex beamComplex

    rXY = hypot(XX, YY);
    noiseMask = rXY > 0.30 & rXY < 1.0;
    dirtyRMS = robustRMS(dirty(noiseMask));

    psf_cpu = psf;
    beam = fitBeam(psf_cpu, XX, YY);
    fprintf('      -> Synthesized Beam: %.4f x %.4f mas at %.1f deg\n', beam.major, beam.minor, beam.pa);
    stageBreak(CFG);

    %% ============================================================
    % STEP 6 — HÖGBOM CLEAN ON A CENTRAL WINDOW (CPU)
    %   Peak search + component subtraction only touch the small window
    %   that contains the CLEAN mask. The exact full-field residual is then
    %   rebuilt once with an FFT convolution (no truncated PSF patch).
    % ============================================================
    thermalGuard(CFG, 'before CLEAN', false);
    fprintf('[5/9] Running windowed Hogbom CLEAN (CPU)...\n');
    cleanTimer = tic;

    Rpix = ceil(CFG.cleanRadius_mas / pixelMas);
    Wc   = min(Rpix + 4, floor(N/2) - 4);
    rows = (c - Wc):(c + Wc);
    Mw   = numel(rows);
    cleanMaskW = rXY(rows, rows) <= CFG.cleanRadius_mas;

    Pw = 2 * Wc;                                   % PSF half-size covering all in-window offsets
    if (c - Pw < 1) || (c + Pw > N)
        error('CLEAN window too large for the grid. Reduce CFG.cleanRadius_mas.');
    end
    psfPatch = psf(c-Pw:c+Pw, c-Pw:c+Pw);

    resW   = dirty(rows, rows);
    modelW = zeros(Mw, Mw);
    threshold = CFG.cleanSigma * dirtyRMS;

    modelFlux = 0;
    nIter = 0;
    cleanHistory = zeros(CFG.cleanMaxIter, 2);

    for iter = 1:CFG.cleanMaxIter
        search = resW;
        if CFG.positiveClean
            search(~cleanMaskW) = -Inf;
            [peak, ind] = max(search(:));
        else
            search(~cleanMaskW) = 0;
            [~, ind] = max(abs(search(:)));
            peak = search(ind);
        end

        if ~isfinite(peak) || abs(peak) <= threshold
            break;
        end

        [py, px] = ind2sub([Mw, Mw], ind);
        component = CFG.cleanGain * peak;

        modelW(py, px) = modelW(py, px) + component;
        modelFlux = modelFlux + component;

        sy = (Pw + 1 + (1 - py)):(Pw + 1 + (Mw - py));
        sx = (Pw + 1 + (1 - px)):(Pw + 1 + (Mw - px));
        resW = resW - component * psfPatch(sy, sx);

        nIter = iter;
        cleanHistory(iter, :) = [iter, peak];

        if modelFlux > CFG.maxModelFluxJy
            fprintf('      -> Clean halted: Flux limit reached.\n');
            break;
        end
    end
    cleanHistory = cleanHistory(1:nIter, :);

    model = zeros(N, N);
    model(rows, rows) = modelW;
    residual = dirty - convolveFFT(model, psf, useGPU);

    fprintf('      -> Completed %d iterations. Clean Flux: %.4f Jy (%.2f s)\n', ...
        nIter, modelFlux, toc(cleanTimer));
    stageBreak(CFG);

    %% ============================================================
    % STEP 7 — RESTORATION & METRICS
    %   Small Gaussian kernel applied only to the crop that holds the model.
    % ============================================================
    thermalGuard(CFG, 'before restore', false);
    fprintf('[6/9] Restoring Image...\n');

    kh = ceil(2.5 * max(beam.major, beam.minor) / pixelMas);
    kx = (-kh:kh) * pixelMas;
    [KX, KY] = meshgrid(kx, kx);
    restoreKernel = gaussianBeam(KX, KY, beam.major, beam.minor, beam.pa);
    restoreKernel = restoreKernel / max(restoreKernel(:));

    cr = max(1, c - Wc - kh):min(N, c + Wc + kh);       % crop that contains model + kernel tail

    restoredModel = zeros(N, N);
    restoredModel(cr, cr) = conv2(model(cr, cr), restoreKernel, 'same');

    beamRadius = 1.5 * max(beam.major, beam.minor);
    beamAreaMask = rXY <= beamRadius & psf_cpu > 0;
    dirtyBeamArea = sum(psf_cpu(beamAreaMask)) * pixelMas^2;
    cleanBeamArea = pi / (4 * log(2)) * beam.major * beam.minor;

    residualScale = cleanBeamArea / dirtyBeamArea;
    if ~isfinite(residualScale) || residualScale < 0.1 || residualScale > 10
        residualScale = 1;
    end

    restored_cpu = restoredModel + residual * residualScale;
    restored_cpu(isnan(restored_cpu)) = 0;           % NaN sanitisation for plotting

    restoredPeak = max(restored_cpu(:));
    restoredRMS = robustRMS(restored_cpu(noiseMask));
    restoredDR = restoredPeak / max(restoredRMS, eps);

    fprintf('      -> Compute Time (FFT + CLEAN + restore): %.2f s\n', toc(gpuTimer));
    stageBreak(CFG);

    %% ============================================================
    % STEP 8 — POLARIMETRIC IMAGING (optional, same low-load scheme)
    % ============================================================
    pol = [];
    if anyPol
        thermalGuard(CFG, 'before polarimetry', false);
        fprintf('[7/9] Polarimetric imaging (Q, U) ...\n');
        polTimer = tic;
        try
            visP = [RL; conj(LR)];
            wP   = [WRL; WLR];
            okP  = isfinite(real(visP)) & isfinite(imag(visP)) & isfinite(wP) & wP > 0;

            [uvNumP, uvWtP] = gridVisibilities(uH(okP), vH(okP), visP(okP), wP(okP), ...
                N, du, CFG.kernelHalfWidth, CFG.kernelSigma, ...
                @() thermalGuard(CFG, 'pol grid', true));
            totWP = sum(uvWtP(:));

            dirtyP = fftImage(uvNumP, useGPU) / totWP;            % complex: Q + iU
            beamP  = fftImage(uvWtP, useGPU);
            clear uvNumP uvWtP

            psfP = real(beamP / real(beamP(c, c)));
            psfP = 0.5 * (psfP + circshift(rot90(psfP, 2), [1 1]));
            psfP = psfP / max(abs(psfP(:)));
            clear beamP

            sigmaP = robustRMS([real(dirtyP(noiseMask)); imag(dirtyP(noiseMask))]);
            thrP   = CFG.polCleanSigma * sigmaP;

            psfPatchP = psfP(c-Pw:c+Pw, c-Pw:c+Pw);
            resP    = dirtyP(rows, rows);
            modelPW = complex(zeros(Mw, Mw));

            nIterP = 0;
            for iter = 1:CFG.polCleanMaxIter
                search = abs(resP);
                search(~cleanMaskW) = 0;
                [peakP, indP] = max(search(:));

                if ~isfinite(peakP) || peakP <= thrP
                    break;
                end

                [py, px] = ind2sub([Mw, Mw], indP);
                compP = CFG.cleanGain * resP(py, px);

                modelPW(py, px) = modelPW(py, px) + compP;

                sy = (Pw + 1 + (1 - py)):(Pw + 1 + (Mw - py));
                sx = (Pw + 1 + (1 - px)):(Pw + 1 + (Mw - px));
                resP = resP - compP * psfPatchP(sy, sx);

                nIterP = iter;
            end

            modelP = complex(zeros(N, N));
            modelP(rows, rows) = modelPW;

            residualP = dirtyP - convolveFFT(modelP, psfP, useGPU);

            restoredP = complex(zeros(N, N));
            restoredP(cr, cr) = conv2(modelP(cr, cr), restoreKernel, 'same');
            restoredP = restoredP + residualP * residualScale;
            restoredP(~isfinite(restoredP)) = 0;

            Qimg = real(restoredP);
            Uimg = imag(restoredP);
            sigQU = robustRMS([Qimg(noiseMask); Uimg(noiseMask)]);

            Pabs = hypot(Qimg, Uimg);
            Pdeb = sqrt(max(Pabs.^2 - sigQU^2, 0));

            frac = Pdeb ./ max(restored_cpu, eps);
            frac(restored_cpu < 5 * restoredRMS) = NaN;
            frac(frac > 1) = NaN;

            evpa = 0.5 * atan2d(Uimg, Qimg) + CFG.evpaOffsetDeg;
            evpa = mod(evpa + 90, 180) - 90;          % wrap to [-90, 90)
            evpa(Pdeb < 3 * sigQU) = NaN;

            region = restored_cpu >= 0.05 * restoredPeak;
            sumI = sum(restored_cpu(region));
            netP = sum(restoredP(region));

            pol = struct();
            pol.Q = Qimg;
            pol.U = Uimg;
            pol.P = Pabs;
            pol.Pdeb = Pdeb;
            pol.frac = frac;
            pol.evpa = evpa;
            pol.sigma = sigQU;
            pol.peak = max(Pdeb(:));
            pol.dynamicRange = max(Pdeb(:)) / max(sigQU, eps);
            pol.mNet = abs(netP) / max(sumI, eps);
            pol.mAvg = sum(Pdeb(region)) / max(sumI, eps);
            pol.evpaMean = mod(0.5 * angle(netP) * 180/pi + CFG.evpaOffsetDeg + 90, 180) - 90;
            pol.cleanIterations = nIterP;
            pol.nVis = nnz(okP);

            fprintf('      -> Pol CLEAN %d iter | |P|peak %.4g Jy/beam | m_net %.1f %% | (%.2f s)\n', ...
                nIterP, pol.peak, 100*pol.mNet, toc(polTimer));
        catch ME
            warning('Polarimetric imaging failed: %s', ME.message);
            pol = [];
        end
    else
        fprintf('[7/9] No RL/LR products found -> skipping polarimetry.\n');
    end
    stageBreak(CFG);

    %% ============================================================
    % STEP 9 — VISIBILITY-DOMAIN TELEMETRY
    % ============================================================
    fprintf('[8/9] Building visibility telemetry...\n');
    telTimer = tic;

    rMeas = hypot(uML, vML) / 1e3;          % |B| in Glambda
    ampMeas = abs(VIS);

    % CLEAN-model visibilities interpolated at the measured uv points
    model_cpu = gather(model);
    modelGrid = fftshift(fft2(ifftshift(model_cpu)));
    mAmp = modelAmpAtUV(modelGrid, uML, vML, du, N);

    % ---- radial (|B|) binned amplitude with model overlay ----------------
    nb = CFG.nRadialBins;
    edges = linspace(0, max(rMeas) * 1.001, nb + 1);
    bi = discretize(rMeas, edges);
    okb = ~isnan(bi);
    sw = accumarray(bi(okb), WT(okb), [nb 1], @sum, 0);
    sa = accumarray(bi(okb), WT(okb) .* ampMeas(okb), [nb 1], @sum, 0);
    meanA = sa ./ sw;

    dev = ampMeas(okb) - meanA(bi(okb));
    sv = accumarray(bi(okb), WT(okb) .* dev.^2, [nb 1], @sum, 0);
    stdA = sqrt(sv ./ sw);

    okm = okb & isfinite(mAmp);
    swm = accumarray(bi(okm), WT(okm), [nb 1], @sum, 0);
    sm  = accumarray(bi(okm), WT(okm) .* mAmp(okm), [nb 1], @sum, 0);
    meanM = sm ./ swm;

    cnt = accumarray(bi(okb), 1, [nb 1], @sum, 0);

    visRadial = struct();
    visRadial.r     = 0.5 * (edges(1:end-1) + edges(2:end))';
    visRadial.amp   = meanA;
    visRadial.std   = stdA;
    visRadial.err   = 1 ./ sqrt(sw);
    visRadial.model = meanM;
    visRadial.n     = cnt;

    % ---- baselines -------------------------------------------------------
    if hasBL
        [blList, ~, blIdx] = unique([min(S1, S2), max(S1, S2)], 'rows');
        blNames = cell(size(blList, 1), 1);
        for q = 1:size(blList, 1)
            blNames{q} = [stationNames{blList(q, 1)} '-' stationNames{blList(q, 2)}];
        end
    else
        blIdx = ones(nMeasured, 1);
        blNames = {'ALL'};
    end

    % ---- random sample of raw visibilities --------------------------------
    sIdx = randperm(nMeasured, min(nMeasured, CFG.visSampleN));
    visSample = struct();
    visSample.u        = uML(sIdx) / 1e3;
    visSample.v        = vML(sIdx) / 1e3;
    visSample.r        = rMeas(sIdx);
    visSample.amp      = ampMeas(sIdx);
    visSample.phase    = angle(VIS(sIdx)) * 180/pi;
    visSample.snr      = ampMeas(sIdx) .* sqrt(WT(sIdx));
    visSample.modelAmp = mAmp(sIdx);
    visSample.bl       = blIdx(sIdx);
    visSample.t        = TIME(sIdx);
    visSample.file     = FID(sIdx);

    % ---- closure phases ----------------------------------------------------
    closure = [];
    if hasBL && numel(stationNames) >= 3
        try
            closure = computeClosurePhases(TIME, S1, S2, VIS, FID, ...
                stationNames, CFG.minTriangleSamples, CFG.maxClosurePoints);
            fprintf('      -> Closure phases: %d points, %d triangles\n', ...
                numel(closure.val), numel(closure.triNames));
        catch ME
            warning('Closure-phase computation failed: %s', ME.message);
            closure = [];
        end
    end

    fprintf('      -> Telemetry ready. (%.2f s)\n', toc(telTimer));

    %% ============================================================
    % PACKAGING OUTPUT
    % ============================================================
    fprintf('[9/9] Compiling Data Structure & Saving...\n');

    reconstruction = struct();
    reconstruction.config = CFG;
    reconstruction.dirtyImage = gather(dirty);
    reconstruction.psf = psf_cpu;
    reconstruction.dirtyRMS = dirtyRMS;
    reconstruction.beam = beam;
    reconstruction.cleanModel = model_cpu;
    reconstruction.cleanResidual = gather(residual);
    reconstruction.cleanIterations = nIter;
    reconstruction.cleanFlux = modelFlux;
    reconstruction.cleanHistory = cleanHistory;
    reconstruction.restoredImage = restored_cpu;
    reconstruction.restoredPeak = restoredPeak;
    reconstruction.restoredRMS = restoredRMS;
    reconstruction.restoredDynamicRange = restoredDR;

    % UV Telemetry for UI Dashboard (original fields)
    reconstruction.uvPlotU = uvPlotU;
    reconstruction.uvPlotV = uvPlotV;
    reconstruction.uvPlotVis = uvPlotVis;

    % v2 telemetry
    reconstruction.stationNames = stationNames;
    reconstruction.blNames = blNames;
    reconstruction.visSample = visSample;
    reconstruction.visRadial = visRadial;
    reconstruction.closure = closure;
    reconstruction.dataInfo = struct('nFiles', nFiles, 'nVis', nMeasured, ...
        'hasBaselines', hasBL, 'hasPol', ~isempty(pol));
    if ~isempty(pol)
        reconstruction.pol = pol;
    end

    save(CFG.outputMAT, 'reconstruction', '-v7.3');

    fprintf('\n======================================================\n');
    fprintf('PIPELINE COMPLETE (Total Time: %.2f s)\n', toc(totalTimer));
    fprintf('Peak Intensity: %.6f Jy/B | Dynamic Range: %.1f\n', restoredPeak, restoredDR);
    fprintf('Data saved to: %s\n', CFG.outputMAT);
    fprintf('======================================================\n');
end

%% ========================================================================
% CONFIG / GRIDDING / TELEMETRY HELPERS
% =========================================================================
function CFG = cfgDefault(CFG, name, value)
    if ~isfield(CFG, name) || isempty(CFG.(name))
        CFG.(name) = value;
    end
end

function [uvNum, uvWt] = gridVisibilities(uH, vH, visH, wtH, N, du, kh, ks, guard)
    if nargin < 9, guard = []; end
    uvNum = complex(zeros(N, N));
    uvWt = zeros(N, N);
    c = N/2 + 1;

    gu = round(uH/du);
    gv = round(vH/du);
    fu = uH/du - gu;
    fv = vH/du - gv;

    lim = N/2 - kh - 2;
    validGrid = gu >= -lim & gu <= lim & gv >= -lim & gv <= lim;

    gu = gu(validGrid); gv = gv(validGrid); fu = fu(validGrid); fv = fv(validGrid);
    vv = visH(validGrid); ww = wtH(validGrid);

    for oy = -kh:kh
        if ~isempty(guard), guard(); end
        dy = oy - fv;
        wy = exp(-0.5 * (dy / ks).^2);
        for ox = -kh:kh
            dx = ox - fu;
            wx = exp(-0.5 * (dx / ks).^2);
            kw = wx .* wy;
            idx = (gv + oy + c) + (gu + ox + c - 1) * N;

            aW = accumarray(idx, ww .* kw, [N*N, 1], @sum, 0);
            aV = accumarray(idx, ww .* vv .* kw, [N*N, 1], @sum, 0);

            uvWt = uvWt + reshape(aW, [N, N]);
            uvNum = uvNum + reshape(aV, [N, N]);
        end
    end
end

function a = modelAmpAtUV(Vgrid, u, v, du, N)
    % Bilinear sampling of the gridded CLEAN-model transform at (u,v) [Mlambda].
    c = N/2 + 1;
    xq = u / du + c;
    yq = v / du + c;
    re = interp2(real(Vgrid), xq, yq, 'linear', NaN);
    im = interp2(imag(Vgrid), xq, yq, 'linear', NaN);
    a = hypot(re, im);
end

function cp = computeClosurePhases(TIME, S1, S2, VIS, FID, stNames, minSamples, maxPts)
    nSta = numel(stNames);
    triAll = nchoosek(1:nSta, 3);

    tOut = cell(0, 1); vOut = cell(0, 1); iOut = cell(0, 1);

    for f = unique(FID(:))'
        m = FID == f;
        t = TIME(m); a = S1(m); b = S2(m); z = VIS(m);

        tk = round((t - min(t)) * 86400);       % 1 s time keys
        [~, ~, ti] = unique(tk);
        nT = max(ti);

        Vc = complex(zeros(nT, nSta, nSta));
        Vc(sub2ind(size(Vc), ti, a, b)) = z;
        Vc(sub2ind(size(Vc), ti, b, a)) = conj(z);

        tJD = accumarray(ti, t, [nT 1], @mean);

        for q = 1:size(triAll, 1)
            i1 = triAll(q, 1); i2 = triAll(q, 2); i3 = triAll(q, 3);
            vab = Vc(:, i1, i2);
            vbc = Vc(:, i2, i3);
            vca = Vc(:, i3, i1);
            ok = vab ~= 0 & vbc ~= 0 & vca ~= 0;
            if ~any(ok), continue; end

            tOut{end+1, 1} = tJD(ok); %#ok<AGROW>
            vOut{end+1, 1} = angle(vab(ok) .* vbc(ok) .* vca(ok)) * 180/pi; %#ok<AGROW>
            iOut{end+1, 1} = q * ones(nnz(ok), 1); %#ok<AGROW>
        end
    end

    if isempty(tOut)
        error('No closed triangles found.');
    end

    t = vertcat(tOut{:});
    val = vertcat(vOut{:});
    tri = vertcat(iOut{:});

    counts = accumarray(tri, 1, [size(triAll, 1) 1]);
    keep = find(counts >= minSamples);
    if isempty(keep)
        error('No triangle reached %d samples.', minSamples);
    end

    sel = ismember(tri, keep);
    t = t(sel); val = val(sel); tri = tri(sel);

    [~, ~, triNew] = unique(tri);

    if numel(val) > maxPts
        r = randperm(numel(val), maxPts);
        t = t(r); val = val(r); triNew = triNew(r);
    end

    names = cell(numel(keep), 1);
    for q = 1:numel(keep)
        tt = triAll(keep(q), :);
        names{q} = strjoin(stNames(tt), '-');
    end

    cp = struct('t', t, 'val', val, 'tri', triNew, 'triNames', {names}, ...
                'triList', triAll(keep, :));
end

function [g1, g2, names] = mapStations(D, names)
    g1 = zeros(numel(D.ant1), 1);
    g2 = g1;

    ids = unique([D.ant1(:); D.ant2(:)]);
    ids = ids(ids > 0);
    lut = zeros(max([ids; 1]), 1);

    for q = 1:numel(ids)
        id = ids(q);
        nm = sprintf('A%d', id);
        j = find(D.antNum == id, 1);
        if ~isempty(j)
            nm = stationLabel(D.antName{j});
        end
        gi = find(strcmp(names, nm), 1);
        if isempty(gi)
            names{end+1} = nm; %#ok<AGROW>
            gi = numel(names);
        end
        lut(id) = gi;
    end

    m1 = D.ant1 > 0; g1(m1) = lut(D.ant1(m1));
    m2 = D.ant2 > 0; g2(m2) = lut(D.ant2(m2));
end

function s = stationLabel(code)
    code = upper(strtrim(code));
    switch code
        case 'AA', s = 'ALMA';
        case 'AP', s = 'APEX';
        case 'AZ', s = 'SMT';
        case 'JC', s = 'JCMT';
        case 'LM', s = 'LMT';
        case 'PV', s = 'IRAM30';
        case 'SM', s = 'SMA';
        case 'SP', s = 'SPT';
        case 'GL', s = 'GLT';
        case 'KT', s = 'KP';
        case 'KP', s = 'KP';
        case 'NN', s = 'NOEMA';
        case 'PB', s = 'PdB';
        otherwise, s = code;
    end
    if isempty(s), s = 'UNK'; end
end

%% ========================================================================
% NON-UI COMPUTATIONAL HELPERS
% =========================================================================
function D = readRawUVFITS(filename)
    D = struct('uSec',[],'vSec',[],'wSec',[],'visibility',[],'weight',[],'time',[], ...
               'frequencyHz',NaN,'ant1',[],'ant2',[],'antNum',[],'antName',{{}}, ...
               'hasPol',false,'RL',[],'LR',[],'WRL',[],'WLR',[]);
    fid = fopen(filename, 'r', 'ieee-be');
    if fid < 0, error('Could not open UVFITS file: %s', filename); end
    cleanup = onCleanup(@()fclose(fid));

    header = containers.Map('KeyType','char','ValueType','char');
    cardCount = 0;
    while true
        bytes = fread(fid, 80, '*uint8');
        if numel(bytes) ~= 80, break; end
        card = char(bytes(:).');
        cardCount = cardCount + 1;
        key = upper(strtrim(card(1:8)));
        if strcmp(key, 'END'), break; end
        if numel(card) >= 10, header(key) = strtrim(card(11:end)); end
    end

    headerBytes = ceil((cardCount * 80) / 2880) * 2880;
    fseek(fid, headerBytes, 'bof');

    bitpix = getHeaderNumber(header, 'BITPIX', -32);
    naxis = round(getHeaderNumber(header, 'NAXIS', 0));
    naxisValues = zeros(1, naxis);
    for k = 1:naxis, naxisValues(k) = getHeaderNumber(header, sprintf('NAXIS%d', k), 0); end

    pcount = round(getHeaderNumber(header, 'PCOUNT', 0));
    gcount = round(getHeaderNumber(header, 'GCOUNT', 1));
    arrayDims = naxisValues(2:end);
    groupLength = pcount + prod(arrayDims);

    precision = 'single';
    if bitpix == -64, precision = 'double'; end
    raw = fread(fid, groupLength * gcount, ['*' precision]);
    raw = reshape(raw, groupLength, gcount);

    % byte offset just after the primary (random-groups) HDU
    bytesPer = abs(bitpix) / 8;
    primaryEnd = headerBytes + ceil(groupLength * gcount * bytesPer / 2880) * 2880;

    ptype = cell(pcount, 1);
    for k = 1:pcount
        if isKey(header, sprintf('PTYPE%d', k))
            ptype{k} = normalizeFITSString(header(sprintf('PTYPE%d', k)));
        else
            ptype{k} = '';
        end
    end

    iu = findUVParameter(ptype, 'UU');
    iv = findUVParameter(ptype, 'VV');
    iw = findUVParameter(ptype, 'WW');
    u = applyPScale(header, double(raw(iu, :))', iu);
    v = applyPScale(header, double(raw(iv, :))', iv);
    w = applyPScale(header, double(raw(iw, :))', iw);

    dataRaw = raw(pcount + 1:end, :);
    data = reshape(dataRaw, [arrayDims, gcount]);
    axis3 = find(arrayDims == 3, 1);

    order = 1:(numel(arrayDims) + 1);
    tmp = order(1); order(1) = order(axis3); order(axis3) = tmp;
    data = permute(data, order);
    data = reshape(data, 3, [], gcount);

    selectedPlane = 1;
    for p = 1:size(data, 2)
        if any(double(squeeze(data(3, p, :))) > 0)
            selectedPlane = p;
            break;
        end
    end

    bscale = getHeaderNumber(header, 'BSCALE', 1);
    bzero = getHeaderNumber(header, 'BZERO', 0);

    [zI, wI] = planeVis(data, selectedPlane, bscale, bzero);

    % ---- optional: observation time ---------------------------------------
    tJD = zeros(gcount, 1);
    try
        idt = find(strcmp(ptype, 'DATE') | strcmp(ptype, '_DATE'));
        for q = idt(:)'
            tJD = tJD + applyPScale(header, double(raw(q, :))', q);
        end
    catch
        tJD = zeros(gcount, 1);
    end

    % ---- optional: baseline -> antenna numbers ----------------------------
    a1 = zeros(gcount, 1); a2 = zeros(gcount, 1);
    try
        ib = findUVParameter(ptype, 'BASELINE');
        if ~isempty(ib)
            bl = round(applyPScale(header, double(raw(ib, :))', ib));
            a1 = floor(bl / 256);
            a2 = mod(bl, 256);
        end
    catch
    end

    % ---- optional: polarisation planes (RR, LL, RL, LR) --------------------
    zRL = complex(nan(gcount, 1), nan(gcount, 1)); wRL = zeros(gcount, 1);
    zLR = zRL; wLR = zeros(gcount, 1);
    hasPol = false;
    try
        if axis3 == 1 && numel(arrayDims) >= 2
            namesArr = cell(1, numel(arrayDims));
            for kk = 1:numel(arrayDims)
                key = sprintf('CTYPE%d', kk + 1);
                if isKey(header, key)
                    namesArr{kk} = normalizeFITSString(header(key));
                else
                    namesArr{kk} = '';
                end
            end
            sk = find(strcmp(namesArr, 'STOKES'), 1);

            if ~isempty(sk) && arrayDims(sk) >= 4
                nStk = arrayDims(sk);
                crval = getHeaderNumber(header, sprintf('CRVAL%d', sk + 1), 1);
                crpix = getHeaderNumber(header, sprintf('CRPIX%d', sk + 1), 1);
                cdelt = getHeaderNumber(header, sprintf('CDELT%d', sk + 1), 1);
                codes = round(crval + ((1:nStk) - crpix) * cdelt);

                restDims = arrayDims(2:end);
                subs = cell(1, numel(restDims));
                [subs{:}] = ind2sub(restDims, selectedPlane);
                sPos = sk - 1;

                getPlane = @(code) planeIndexForCode(codes, code, subs, sPos, restDims);

                iRR = getPlane(-1); iLL = getPlane(-2);
                iRL = getPlane(-3); iLR = getPlane(-4);

                if ~isempty(iRR) && ~isempty(iLL) && ~isempty(iRL) && ~isempty(iLR)
                    [zRR, wRR] = planeVis(data, iRR, bscale, bzero);
                    [zLL, wLL] = planeVis(data, iLL, bscale, bzero);
                    [zRL, wRL] = planeVis(data, iRL, bscale, bzero);
                    [zLR, wLR] = planeVis(data, iLR, bscale, bzero);

                    % Stokes I = (RR + LL) / 2
                    ok2 = wRR > 0 & wLL > 0;
                    wI = zeros(gcount, 1);
                    wI(ok2) = 4 ./ (1 ./ wRR(ok2) + 1 ./ wLL(ok2));
                    zI = (zRR + zLL) / 2;
                    zI(~ok2) = complex(NaN, NaN);

                    hasPol = true;

                else
                    % Fallback: files already stored as Stokes I, Q, U (codes 1,2,3)
                    iI = getPlane(1); iQ = getPlane(2); iU = getPlane(3);

                    if ~isempty(iI) && ~isempty(iQ) && ~isempty(iU)
                        [zIs, wIs] = planeVis(data, iI, bscale, bzero);
                        [zQ,  wQ ] = planeVis(data, iQ, bscale, bzero);
                        [zU,  wU ] = planeVis(data, iU, bscale, bzero);

                        okq = wQ > 0 & wU > 0;
                        wPq = zeros(gcount, 1);
                        wPq(okq) = 1 ./ (1 ./ wQ(okq) + 1 ./ wU(okq));

                        zI = zIs;
                        wI = wIs;

                        % Map to the RL/LR convention used downstream:
                        % RL = Q + iU, LR = Q - iU
                        zRL = zQ + 1i * zU;
                        zLR = zQ - 1i * zU;
                        zRL(~okq) = complex(NaN, NaN);
                        zLR(~okq) = complex(NaN, NaN);
                        wRL = wPq;
                        wLR = wPq;

                        hasPol = true;
                    end
                end
            end
        end
    catch
        hasPol = false;
    end

    frequencyHz = findFrequency(header);
    if ~isfinite(frequencyHz) || frequencyHz <= 1e9, frequencyHz = frequencyFromFilename(filename); end

    valid = isfinite(u) & isfinite(v) & isfinite(w) & isfinite(real(zI)) & ...
            isfinite(imag(zI)) & isfinite(wI) & wI > 0;

    D.uSec = u(valid); D.vSec = v(valid); D.wSec = w(valid);
    D.visibility = zI(valid); D.weight = wI(valid); D.frequencyHz = frequencyHz;
    D.time = tJD(valid);
    D.ant1 = a1(valid); D.ant2 = a2(valid);

    if hasPol
        D.hasPol = true;
        D.RL = zRL(valid); D.LR = zLR(valid);
        D.WRL = wRL(valid); D.WLR = wLR(valid);
        D.RL(~isfinite(D.RL)) = 0; D.LR(~isfinite(D.LR)) = 0;
        D.WRL(~isfinite(D.WRL)) = 0; D.WLR(~isfinite(D.WLR)) = 0;
    end

    % ---- optional: AIPS AN antenna table ----------------------------------
    try
        [D.antNum, D.antName] = readANTable(fid, primaryEnd);
    catch
        D.antNum = []; D.antName = {};
    end
end

function [z, wt] = planeVis(data, planeIdx, bscale, bzero)
    re = double(squeeze(data(1, planeIdx, :)));
    im = double(squeeze(data(2, planeIdx, :)));
    wt = double(squeeze(data(3, planeIdx, :)));
    re = re(:) * bscale + bzero;
    im = im(:) * bscale + bzero;
    wt = wt(:) * bscale + bzero;
    z = complex(re, im);
end

function lin = planeIndexForCode(codes, code, subs, sPos, restDims)
    s = find(codes == code, 1);
    if isempty(s)
        lin = [];
        return;
    end
    subs{sPos} = s;
    lin = sub2ind(restDims, subs{:});
end

function [num, name] = readANTable(fid, startByte)
    num = []; name = {};
    pos = startByte;

    for hdu = 1:20 %#ok<NASGU>
        if fseek(fid, pos, 'bof') ~= 0, return; end

        hdr = containers.Map('KeyType','char','ValueType','char');
        nCards = 0;
        while true
            bytes = fread(fid, 80, '*uint8');
            if numel(bytes) ~= 80, return; end
            card = char(bytes(:).');
            nCards = nCards + 1;
            key = upper(strtrim(card(1:8)));
            if strcmp(key, 'END'), break; end
            if numel(card) >= 10, hdr(key) = strtrim(card(11:end)); end
        end

        hdrBytes = ceil(nCards * 80 / 2880) * 2880;
        dataStart = pos + hdrBytes;

        naxis1 = round(getHeaderNumber(hdr, 'NAXIS1', 0));
        naxis2 = round(getHeaderNumber(hdr, 'NAXIS2', 0));
        pc = round(getHeaderNumber(hdr, 'PCOUNT', 0));
        dataBytes = naxis1 * naxis2 + pc;

        extname = '';
        if isKey(hdr, 'EXTNAME'), extname = normalizeFITSString(hdr('EXTNAME')); end

        if contains(extname, 'AN')
            nF = round(getHeaderNumber(hdr, 'TFIELDS', 0));
            off = 0;
            nameOff = -1; nameLen = 0; staOff = -1; staBytes = 0;

            for j = 1:nF
                ttype = ''; tform = '';
                if isKey(hdr, sprintf('TTYPE%d', j)), ttype = normalizeFITSString(hdr(sprintf('TTYPE%d', j))); end
                if isKey(hdr, sprintf('TFORM%d', j)), tform = normalizeFITSString(hdr(sprintf('TFORM%d', j))); end

                tok = regexp(tform, '^(\d*)([A-Z])', 'tokens', 'once');
                if isempty(tok), return; end
                rep = str2double(tok{1}); if isnan(rep), rep = 1; end
                switch tok{2}
                    case {'A', 'B', 'L'}, sz = 1;
                    case 'I', sz = 2;
                    case {'J', 'E'}, sz = 4;
                    case {'K', 'D'}, sz = 8;
                    otherwise, return;
                end

                if strcmp(ttype, 'ANNAME')
                    nameOff = off; nameLen = rep;
                elseif strcmp(ttype, 'NOSTA')
                    staOff = off; staBytes = rep * sz;
                end
                off = off + rep * sz;
            end

            if nameOff < 0, return; end

            fseek(fid, dataStart, 'bof');
            rows = fread(fid, [naxis1, naxis2], '*uint8');

            name = cell(naxis2, 1);
            for r = 1:naxis2
                name{r} = strtrim(char(rows(nameOff + 1:nameOff + nameLen, r).'));
            end

            if staOff >= 0 && staBytes >= 1 && staBytes <= 4
                num = zeros(naxis2, 1);
                for b = 1:staBytes
                    num = num * 256 + double(rows(staOff + b, :)).';
                end
            else
                num = (1:naxis2).';
            end
            return;
        end

        pos = dataStart + ceil(dataBytes / 2880) * 2880;
    end
end

function s = normalizeFITSString(s)
    s = strtrim(s);
    if numel(s) >= 2 && s(1) == '''' && s(end) == ''''
        s = s(2:end-1);
    end
    s = upper(strtrim(s));
end

function idx = findUVParameter(ptype, target)
    idx = []; target = upper(strtrim(target));
    for k = 1:numel(ptype)
        name = normalizeFITSString(ptype{k});
        if strcmp(name, target) || startsWith(name, [target '-']) || startsWith(name, [target '_'])
            idx = k; return;
        end
    end
end

function x = applyPScale(header, x, index)
    pscal = getHeaderNumber(header, sprintf('PSCAL%d', index), 1);
    pzero = getHeaderNumber(header, sprintf('PZERO%d', index), 0);
    x = x * pscal + pzero;
end

function val = getHeaderNumber(header, key, def)
    if ~isKey(header, key), val = def; return; end
    s = strtrim(header(key));
    if isempty(s) || s(1) == '''', val = def; return; end
    slash = find(s == '/', 1);
    if ~isempty(slash), s = s(1:slash-1); end
    val = str2double(strtrim(s));
    if ~isfinite(val), val = def; end
end

function freq = findFrequency(header)
    freq = NaN;
    naxis = round(getHeaderNumber(header, 'NAXIS', 0));
    for k = 1:naxis
        key = sprintf('CTYPE%d', k);
        if ~isKey(header, key), continue; end
        ctype = normalizeFITSString(header(key));
        if contains(ctype, 'FREQ')
            crval = getHeaderNumber(header, sprintf('CRVAL%d', k), NaN);
            crpix = getHeaderNumber(header, sprintf('CRPIX%d', k), 1);
            cdelt = getHeaderNumber(header, sprintf('CDELT%d', k), 0);
            freq = crval + (1 - crpix) * cdelt; return;
        end
    end
    if isKey(header, 'RESTFRQ'), freq = getHeaderNumber(header, 'RESTFRQ', NaN); end
end

function freq = frequencyFromFilename(name)
    freq = NaN; name = lower(name);
    if contains(name, '_lo_'), freq = 227.1e9;
    elseif contains(name, '_hi_'), freq = 229.1e9;
    elseif contains(name, '_b1_'), freq = 227.1e9;
    elseif contains(name, '_b2_'), freq = 229.1e9;
    elseif contains(name, '_b3_'), freq = 241.1e9;
    elseif contains(name, '_b4_'), freq = 243.1e9;
    end
end

function s = robustRMS(x)
    x = x(isfinite(x));
    if isempty(x), s = NaN; return; end
    m = median(x); s = 1.4826 * median(abs(x - m));
    if ~isfinite(s) || s <= 0, s = std(x); end
    if ~isfinite(s) || s <= 0, s = eps; end
end

function beam = fitBeam(psf, X, Y)
    r = hypot(X, Y); mask = r <= 0.05 & psf >= 0.20;
    if nnz(mask) < 20
        beam.major = 0.020; beam.minor = 0.015; beam.pa = 0; return;
    end
    w = psf(mask); xx = X(mask); yy = Y(mask); w = w / sum(w);
    mx = sum(w .* xx); my = sum(w .* yy); dx = xx - mx; dy = yy - my;
    C = [sum(w .* dx.^2), sum(w .* dx .* dy); sum(w .* dx .* dy), sum(w .* dy.^2)];
    [V, L] = eig(C); [lambda, ord] = sort(diag(L), 'descend'); V = V(:, ord);
    sigmaMajor = sqrt(max(lambda(1), eps)); sigmaMinor = sqrt(max(lambda(2), eps));
    beam.major = 2 * sqrt(2 * log(2)) * sigmaMajor;
    beam.minor = 2 * sqrt(2 * log(2)) * sigmaMinor;
    beam.pa = atan2d(V(2, 1), V(1, 1));
    beam.major = max(min(beam.major, 0.08), 0.005);
    beam.minor = max(min(beam.minor, 0.08), 0.005);
    if beam.minor > beam.major
        tmp = beam.major; beam.major = beam.minor; beam.minor = tmp;
        beam.pa = beam.pa + 90;
    end
end

function G = gaussianBeam(X, Y, major, minor, pa)
    theta = deg2rad(pa);
    xr = X * cos(theta) + Y * sin(theta); yr = -X * sin(theta) + Y * cos(theta);
    G = exp(-4 * log(2) * ((xr / major).^2 + (yr / minor).^2));
end

%% ========================================================================
% LOW-LOAD COMPUTE HELPERS (v2.1)
% =========================================================================
function img = fftImage(A, useGPU)
    % Centred inverse FFT, scaled by N^2. GPU FFT when requested; result on CPU.
    N = size(A, 1);
    try
        if useGPU
            img = gather(fftshift(ifft2(ifftshift(gpuArray(A)))) * N^2);
        else
            img = fftshift(ifft2(ifftshift(A))) * N^2;
        end
    catch
        img = fftshift(ifft2(ifftshift(A))) * N^2;
    end
end

function out = convolveFFT(A, K, useGPU)
    % Circular convolution of two centred N x N arrays via FFT (exact PSF,
    % no patch truncation). Complex in -> complex out; real in -> real out.
    isCplx = ~isreal(A);
    try
        if useGPU
            F = fft2(ifftshift(gpuArray(A))) .* fft2(ifftshift(gpuArray(K)));
            out = gather(fftshift(ifft2(F)));
        else
            F = fft2(ifftshift(A)) .* fft2(ifftshift(K));
            out = fftshift(ifft2(F));
        end
    catch
        F = fft2(ifftshift(A)) .* fft2(ifftshift(K));
        out = fftshift(ifft2(F));
    end
    if ~isCplx
        out = real(out);
    end
end

function t = readMaxTempC()
    % Hottest ACPI thermal zone on Linux (deg C); NaN if unavailable.
    t = NaN;
    try
        d = dir('/sys/class/thermal/thermal_zone*/temp');
        v = [];
        for q = 1:numel(d)
            fid = fopen(fullfile(d(q).folder, d(q).name), 'r');
            if fid > 0
                x = fscanf(fid, '%f');
                fclose(fid);
                if ~isempty(x), v(end+1) = x(1) / 1000; end %#ok<AGROW>
            end
        end
        if ~isempty(v), t = max(v); end
    catch
    end
end

function thermalGuard(CFG, label, quiet)
    % Waits for the machine to cool, and aborts if it stays dangerously hot.
    if ~CFG.thermalGuard
        return;
    end

    t = readMaxTempC();
    if isnan(t)
        return;
    end

    if ~quiet
        fprintf('      -> [%s] max temperature %.0f C\n', label, t);
    end

    waited = 0;
    while t > CFG.maxTempC && waited < CFG.maxCoolWaitSec
        fprintf('      -> Cooling pause: %.0f C > %d C (waited %d s)\n', t, CFG.maxTempC, waited);
        pause(10);
        waited = waited + 10;
        t = readMaxTempC();
    end

    if t > CFG.maxTempC + 5
        error('Aborting to protect the laptop: %.0f C is too hot. Let it cool down and rerun.', t);
    end
end

function stageBreak(CFG)
    if CFG.stagePauseSec > 0
        pause(CFG.stagePauseSec);
    end
end