function M87_EHT_Master_Explorer_Dashboard_v3

clc;
close all;

CFG.dataRoot = '/home/hellblazer/Documents/MATLAB/VLBI/2023-D01-01-main';
CFG.uvfitsDir = fullfile(CFG.dataRoot,'uvfits');
CFG.outputMAT = ...
    '/home/hellblazer/Documents/MATLAB/VLBI/M87_EHT_POL_RECONSTRUCTION.mat';

CFG.du_Mlambda     = 20;     % fallback if MAT file has no config
CFG.fovMode        = 'AUTO'; % AUTO | FULL | MANUAL
CFG.displayFOV_uas = 400;    % full FOV used only in MANUAL mode
CFG.raEastLeft     = false;  % true -> RA increases to the left
CFG.autoLoad       = true;
CFG.maxUVPoints    = 150000;

%% ========================================================================
% SPECIAL CHARACTERS
% =========================================================================

MU    = char(181);
TIMES = char(215);
DEG   = char(176);
LAM   = char(955);
SIG   = char(963);

%% ========================================================================
% PALETTE
% =========================================================================

BG_COLOR       = [0.030 0.038 0.052];
CARD_BG        = [0.060 0.075 0.100];
CARD_BG_ALT    = [0.045 0.058 0.078];
CTRL_BG        = [0.100 0.130 0.180];
BORDER_COL     = [0.16 0.22 0.29];

TEXT_PRIMARY   = [0.92 0.96 1.00];
TEXT_SECONDARY = [0.72 0.79 0.88];
TEXT_MUTED     = [0.48 0.57 0.68];

ACCENT_BLUE    = [0.20 0.55 0.95];
ACCENT_CYAN    = [0.05 0.82 0.95];
ACCENT_GREEN   = [0.18 0.82 0.45];
ACCENT_ORANGE  = [0.96 0.62 0.16];
ACCENT_RED     = [0.92 0.25 0.25];
ACCENT_MAGENTA = [0.85 0.40 0.95];
ACCENT_YELLOW  = [0.95 0.88 0.30];

SERIES = [ACCENT_CYAN; ACCENT_ORANGE; ACCENT_GREEN; ACCENT_MAGENTA; ...
          ACCENT_YELLOW; ACCENT_BLUE; ACCENT_RED; 0.80 0.82 0.88];

FONT_NAME      = 'DejaVu Sans';
MONO_FONT      = 'DejaVu Sans Mono';

%% ========================================================================
% SHARED STATE
% =========================================================================

S.recon   = [];
S.x       = [];
S.pixUas  = NaN;
S.du      = CFG.du_Mlambda;
S.fovHalf = CFG.displayFOV_uas/2;
S.beam    = struct('valid',false);
S.uv      = struct('u',[],'v',[],'bmax',NaN);
S.ring    = struct('valid',false);
S.view    = struct();
S.page    = 1;
S.tmove   = tic;

pageOfKey = struct('dirty',1,'model',1,'resid',1,'restored',1,'psf',1, ...
                   'polI',4,'polP',4,'polM',4,'polQ',4,'polU',4);

idleFooter = '  CURSOR  --    |    hover over any image panel for coordinates and pixel values';

%% ========================================================================
% MAIN FIGURE
% =========================================================================

scr  = get(groot,'ScreenSize');
figW = min(1850, scr(3)-40);
figH = min(1040, scr(4)-80);

fig = uifigure( ...
    'Name','M87* EHT | VLBI Scientific Telemetry Explorer v3', ...
    'Position',[20 40 figW figH], ...
    'Color',BG_COLOR);

main = uigridlayout(fig,[5 1]);
main.RowHeight   = {56,62,30,'1x',24};
main.ColumnWidth = {'1x'};
main.BackgroundColor = BG_COLOR;
main.Padding = [10 10 10 8];
main.RowSpacing = 8;

%% ========================================================================
% HEADER
% =========================================================================

headerPanel = makePanel(main,CARD_BG);
headerPanel.Layout.Row = 1;

hg = uigridlayout(headerPanel,[1 12]);
hg.ColumnWidth = {'1x',105,95,150,82,62,78,92,92,88,128,100};
hg.BackgroundColor = CARD_BG;
hg.Padding = [12 8 12 8];
hg.ColumnSpacing = 8;

titleLabel = uilabel(hg, ...
    'Text','  M87* EHT  |  VLBI SCIENTIFIC TELEMETRY EXPLORER', ...
    'FontName',FONT_NAME,'FontSize',14,'FontWeight','bold', ...
    'FontColor',TEXT_PRIMARY); %#ok<NASGU>

gpuOnline = false;
gpuName   = 'N/A';
try
    gpuOnline = canUseGPU();
    if gpuOnline
        gd = gpuDevice;
        gpuName = gd.Name;
    end
catch
end

gpuStatus = uilabel(hg, ...
    'FontName',FONT_NAME,'FontSize',10,'FontWeight','bold', ...
    'HorizontalAlignment','center');
if gpuOnline
    gpuStatus.Text = [char(9679) ' GPU ONLINE'];
    gpuStatus.FontColor = ACCENT_GREEN;
else
    gpuStatus.Text = [char(9679) ' GPU OFFLINE'];
    gpuStatus.FontColor = ACCENT_RED;
end

statusLabel = uilabel(hg, ...
    'Text',[char(9679) ' IDLE'], ...
    'FontName',FONT_NAME,'FontSize',10,'FontWeight','bold', ...
    'FontColor',TEXT_MUTED,'HorizontalAlignment','center');

datasetLabel = uilabel(hg, ...
    'Text','DATA: NONE', ...
    'FontName',MONO_FONT,'FontSize',9,'FontWeight','bold', ...
    'FontColor',TEXT_MUTED,'HorizontalAlignment','center');

ddFovMode = uidropdown(hg, ...
    'Items',{'AUTO','FULL','MANUAL'},'Value',CFG.fovMode, ...
    'FontName',FONT_NAME,'FontSize',9, ...
    'BackgroundColor',CTRL_BG,'FontColor',TEXT_PRIMARY, ...
    'ValueChangedFcn',@onControlChanged);

lblFov = uilabel(hg, ...
    'Text',['FOV (' MU 'as)'], ...
    'FontName',FONT_NAME,'FontSize',9, ...
    'FontColor',TEXT_SECONDARY,'HorizontalAlignment','right'); %#ok<NASGU>

spnFOV = uispinner(hg, ...
    'Limits',[20 20000],'Value',CFG.displayFOV_uas,'Step',20, ...
    'ValueDisplayFormat','%.0f', ...
    'FontName',MONO_FONT,'FontSize',10, ...
    'BackgroundColor',CTRL_BG,'FontColor',TEXT_PRIMARY, ...
    'ValueChangedFcn',@onControlChanged);

ddCmap = uidropdown(hg, ...
    'Items',{'EHT','HOT','TURBO','PARULA','GRAY'},'Value','EHT', ...
    'FontName',FONT_NAME,'FontSize',9, ...
    'BackgroundColor',CTRL_BG,'FontColor',TEXT_PRIMARY, ...
    'ValueChangedFcn',@onControlChanged);

ddStretch = uidropdown(hg, ...
    'Items',{'LINEAR','SQRT','LOG','ASINH'},'Value','LINEAR', ...
    'FontName',FONT_NAME,'FontSize',9, ...
    'BackgroundColor',CTRL_BG,'FontColor',TEXT_PRIMARY, ...
    'ValueChangedFcn',@onControlChanged);

chkSmooth = uicheckbox(hg, ...
    'Text','SMOOTH','Value',false, ...
    'FontName',FONT_NAME,'FontSize',9,'FontColor',TEXT_SECONDARY, ...
    'ValueChangedFcn',@onControlChanged);

loadButton = uibutton(hg,'push', ...
    'Text','LOAD .MAT DATA', ...
    'FontName',FONT_NAME,'FontSize',10,'FontWeight','bold', ...
    'BackgroundColor',[0.10 0.30 0.50],'FontColor',[1 1 1], ...
    'ButtonPushedFcn',@loadSavedData);

exportButton = uibutton(hg,'push', ...
    'Text','EXPORT PNG', ...
    'FontName',FONT_NAME,'FontSize',10,'FontWeight','bold', ...
    'BackgroundColor',CTRL_BG,'FontColor',TEXT_SECONDARY, ...
    'ButtonPushedFcn',@exportFigure);

%% ========================================================================
% KPI STRIP
% =========================================================================

kpiGrid = uigridlayout(main,[1 10]);
kpiGrid.Layout.Row = 2;
kpiGrid.ColumnWidth = repmat({'1x'},1,10);
kpiGrid.BackgroundColor = BG_COLOR;
kpiGrid.ColumnSpacing = 8;
kpiGrid.Padding = [0 0 0 0];

kpiTitles = {'PEAK','NOISE RMS','DYNAMIC RANGE','CLEAN FLUX','SYNTH. BEAM', ...
             'UV SAMPLES','MAX BASELINE','CLEAN ITER.','RING DIAMETER','NET POLARIZATION'};
kpiAccent = {ACCENT_ORANGE,ACCENT_CYAN,ACCENT_GREEN,ACCENT_BLUE,ACCENT_CYAN, ...
             ACCENT_BLUE,ACCENT_ORANGE,ACCENT_GREEN,ACCENT_MAGENTA,ACCENT_YELLOW};

KPI.val = cell(1,10);
KPI.sub = cell(1,10);
for iKpi = 1:10
    [KPI.val{iKpi},KPI.sub{iKpi}] = ...
        makeKPI(kpiGrid,iKpi,kpiTitles{iKpi},kpiAccent{iKpi});
end

%% ========================================================================
% TAB BAR + PAGES
% =========================================================================

tabGrid = uigridlayout(main,[1 5]);
tabGrid.Layout.Row = 3;
tabGrid.ColumnWidth = {130,130,130,150,'1x'};
tabGrid.BackgroundColor = BG_COLOR;
tabGrid.Padding = [0 0 0 0];
tabGrid.ColumnSpacing = 6;

tabNames = {'IMAGING','VISIBILITY','ANALYSIS','POLARIMETRY'};
tabBtn = cell(1,4);
for iTab = 1:4
    tabBtn{iTab} = uibutton(tabGrid,'push', ...
        'Text',tabNames{iTab}, ...
        'FontName',FONT_NAME,'FontSize',10,'FontWeight','bold', ...
        'BackgroundColor',CTRL_BG,'FontColor',TEXT_SECONDARY, ...
        'ButtonPushedFcn',@(~,~) showPage(iTab));
    tabBtn{iTab}.Layout.Column = iTab;
end

pgImg = makePage(2,4);
pgVis = makePage(2,2);
pgAna = makePage(2,3);
pgPol = makePage(2,3);

%% ---- IMAGING page ---------------------------------------------------------

mapCard = makePanel(pgImg,CARD_BG);
mapCard.Layout.Row = 1;
mapCard.Layout.Column = 1;

mapGrid = uigridlayout(mapCard,[2 1]);
mapGrid.RowHeight = {22,'1x'};
mapGrid.ColumnWidth = {'1x'};
mapGrid.BackgroundColor = CARD_BG;
mapGrid.Padding = [7 6 7 6];
mapGrid.RowSpacing = 3;

mapLabel = uilabel(mapGrid, ...
    'Text','  ARRAY GEOMETRY  |  EHT STATIONS', ...
    'FontName',FONT_NAME,'FontSize',9,'FontWeight','bold', ...
    'FontColor',TEXT_SECONDARY);
mapLabel.Layout.Row = 1;

axTelescopes = geoaxes(mapGrid);
axTelescopes.Layout.Row = 2;
try
    geobasemap(axTelescopes,'darkwater');
catch
    try
        geobasemap(axTelescopes,'streets');
    catch
    end
end

PN.dirty    = createCard(pgImg,1,2,'Dirty image',true);
PN.model    = createCard(pgImg,1,3,'CLEAN model',true);
PN.resid    = createCard(pgImg,1,4,'CLEAN residual',true);
PN.restored = createCard(pgImg,2,1,'Restored M87* image',true);
PN.psf      = createCard(pgImg,2,2,'Synthesized beam / PSF',true);
PN.uv       = createCard(pgImg,2,3,'Baseline coverage  |  UV plane',true);

convCard = makePanel(pgImg,CARD_BG);
convCard.Layout.Row = 2;
convCard.Layout.Column = 4;

convGrid = uigridlayout(convCard,[3 1]);
convGrid.RowHeight = {22,'1.2x','1x'};
convGrid.ColumnWidth = {'1x'};
convGrid.BackgroundColor = CARD_BG;
convGrid.Padding = [7 6 7 6];
convGrid.RowSpacing = 3;

convLabel = uilabel(convGrid, ...
    'Text','  CLEAN CONVERGENCE  |  SYSTEM CONSOLE', ...
    'FontName',FONT_NAME,'FontSize',9,'FontWeight','bold', ...
    'FontColor',TEXT_SECONDARY);
convLabel.Layout.Row = 1;

axConv = uiaxes(convGrid);
axConv.Layout.Row = 2;
formatAxes(axConv);

infoText = uitextarea(convGrid, ...
    'BackgroundColor',[0.012 0.017 0.025], ...
    'FontColor',ACCENT_CYAN, ...
    'FontName',MONO_FONT,'FontSize',8, ...
    'Editable','off');
infoText.Layout.Row = 3;
infoText.Value = {''};

%% ---- VISIBILITY page ------------------------------------------------------

PN.ampuv   = createCard(pgVis,1,1,['Visibility amplitude vs |B|  |  data, bins, CLEAN model'],false);
PN.cphase  = createCard(pgVis,1,2,'Closure phase vs time (top triangles)',false);
PN.amptime = createCard(pgVis,2,1,'Amplitude vs time (top baselines)',false);
PN.cphist  = createCard(pgVis,2,2,'Closure phase distribution (top triangle)',false);

%% ---- ANALYSIS page --------------------------------------------------------

PN.radial  = createCard(pgAna,1,1,'Ring radial profile',false);
PN.azimuth = createCard(pgAna,1,2,'Ring azimuthal profile',false);
PN.xsec    = createCard(pgAna,1,3,'Peak cross-sections',false);
PN.reshist = createCard(pgAna,2,1,'Residual histogram vs Gaussian',false);
PN.psfrad  = createCard(pgAna,2,2,'PSF sidelobe envelope',false);
PN.fourier = createCard(pgAna,2,3,'Fourier plane  |  restored image + uv sampling',true);

%% ---- POLARIMETRY page -----------------------------------------------------

PN.polI    = createCard(pgPol,1,1,'Stokes I + EVPA ticks',true);
PN.polP    = createCard(pgPol,1,2,'Polarized intensity |P| (debiased)',true);
PN.polM    = createCard(pgPol,1,3,'Fractional polarization',true);
PN.polQ    = createCard(pgPol,2,1,'Stokes Q',true);
PN.polU    = createCard(pgPol,2,2,'Stokes U',true);
PN.polEvpa = createCard(pgPol,2,3,'EVPA vs position angle',true);

%% ---- Footer ---------------------------------------------------------------

footer = uilabel(main, ...
    'Text',idleFooter, ...
    'FontName',MONO_FONT,'FontSize',9, ...
    'FontColor',TEXT_MUTED,'BackgroundColor',CARD_BG_ALT);
footer.Layout.Row = 5;

%% ========================================================================
% INITIALISE
% =========================================================================

allKeys = fieldnames(PN);
for iPh = 1:numel(allKeys)
    showPlaceholder(PN.(allKeys{iPh}).ax);
end
showPlaceholder(axConv);

drawTelescopeMap(axTelescopes);

fig.WindowButtonMotionFcn = @onMouseMove;

if strcmp(ddFovMode.Value,'MANUAL')
    spnFOV.Enable = 'on';
else
    spnFOV.Enable = 'off';
end

showPage(1);

logMsg('INFO','M87* EHT MASTER EXPLORER v3 ONLINE.');
logMsg('INFO',sprintf('MATLAB %s | %s',version('-release'),computer('arch')));
logMsg('INFO',sprintf('GPU: %s',gpuName));
logMsg('INFO','AWAITING RECONSTRUCTION DATA.');

if CFG.autoLoad && isfile(CFG.outputMAT)
    loadReconstruction(CFG.outputMAT);
end

%% ========================================================================
% PAGES
% =========================================================================

    function g = makePage(nr,nc)

        g = uigridlayout(main,[nr nc]);
        g.Layout.Row = 4;
        g.RowHeight = repmat({'1x'},1,nr);
        g.ColumnWidth = repmat({'1x'},1,nc);
        g.BackgroundColor = BG_COLOR;
        g.Padding = [0 0 0 0];
        g.RowSpacing = 8;
        g.ColumnSpacing = 8;

    end

    function showPage(idx)

        pgList = {pgImg,pgVis,pgAna,pgPol};

        for q = 1:4

            if q == idx
                pgList{q}.Visible = 'on';
                tabBtn{q}.BackgroundColor = [0.10 0.30 0.50];
                tabBtn{q}.FontColor = [1 1 1];
            else
                pgList{q}.Visible = 'off';
                tabBtn{q}.BackgroundColor = CTRL_BG;
                tabBtn{q}.FontColor = TEXT_SECONDARY;
            end

        end

        S.page = idx;

    end

%% ========================================================================
% DATA LOADING
% =========================================================================

    function loadSavedData(~,~)

        startPath = CFG.dataRoot;
        if isfile(CFG.outputMAT)
            startPath = CFG.outputMAT;
        end

        [file,path] = uigetfile('*.mat','Select Reconstruction File',startPath);

        if isequal(file,0)
            return;
        end

        loadReconstruction(fullfile(path,file));

    end

    function loadReconstruction(fullpath)

        [~,fname,fext] = fileparts(fullpath);

        try
            setStatus('LOADING',ACCENT_ORANGE);
            datasetLabel.Text = 'DATA: LOADING';
            drawnow;

            logMsg('INFO',['LOADING FILE: ' fname fext]);

            data = load(fullpath);

            if ~isfield(data,'reconstruction')
                error('Missing "reconstruction" struct in MAT file.');
            end

            recon = data.reconstruction;

            required = {'dirtyImage','cleanModel','cleanResidual', ...
                        'restoredImage','psf'};
            missing = required(~isfield(recon,required));
            if ~isempty(missing)
                error('Missing fields: %s',strjoin(missing,', '));
            end

            % ---- Image coordinates (microarcseconds) ----------------------

            N = size(recon.restoredImage,1);

            du = CFG.du_Mlambda;
            if isfield(recon,'config') && isfield(recon.config,'du_Mlambda')
                du = recon.config.du_Mlambda;
            end

            fovRad = 1/(du*1e6);
            fovUas = fovRad*(180/pi)*3600*1e6;
            pixUas = fovUas/N;

            S.x      = ((1:N)-(N/2+1))*pixUas;
            S.pixUas = pixUas;
            S.du     = du;
            S.recon  = recon;

            % ---- Beam: prefer the pipeline's restoring beam -----------------

            S.beam = beamFromRecon(recon);
            if ~S.beam.valid
                S.beam = fitBeam(recon.psf,S.x);
            end

            % ---- UV data ----------------------------------------------------

            S.uv = struct('u',[],'v',[],'bmax',NaN);
            if isfield(recon,'uvPlotU') && ~isempty(recon.uvPlotU)
                u = double(recon.uvPlotU(:))/1e3;   % Mlambda -> Glambda
                v = double(recon.uvPlotV(:))/1e3;
                ok = isfinite(u) & isfinite(v);
                S.uv.u = u(ok);
                S.uv.v = v(ok);
                if ~isempty(S.uv.u)
                    S.uv.bmax = max(hypot(S.uv.u,S.uv.v));
                end
            end

            logMsg('INFO',sprintf('GRID %d x %d | PIXEL %.3f %sas | FULL FIELD %.2f mas', ...
                N,N,pixUas,MU,fovUas/1e3));

            if S.beam.valid
                logMsg('INFO',sprintf('BEAM: %.1f x %.1f %sas, PA %+.1f%s', ...
                    S.beam.maj,S.beam.min,MU,S.beam.pa,DEG));
            else
                logMsg('WARN','BEAM FIT FAILED.');
            end

            logMsg('INFO',sprintf('TELEMETRY | VIS:%s  CLOSURE:%s  POL:%s', ...
                yesno(hasField(recon,'visSample')), ...
                yesno(hasField(recon,'closure')), ...
                yesno(hasField(recon,'pol'))));

            % ---- Ring analysis, then render ---------------------------------

            S.ring = analyzeRing();

            renderImagePanels();
            renderLinePanels();
            renderKPIs();

            datasetLabel.Text = ['DATA: ' upper(fname)];

            if isfield(recon,'restoredDynamicRange')
                logMsg('OK',sprintf('RESTORED DYNAMIC RANGE: %.1f', ...
                    recon.restoredDynamicRange));
            end

            logMsg('OK','RENDER COMPLETE.');
            setStatus('ONLINE',ACCENT_GREEN);
            drawnow;

        catch ME

            logMsg('ERR',ME.message);

            uialert(fig,ME.message,'Loading Error','Icon','error');

            setStatus('ERROR',ACCENT_RED);
            datasetLabel.Text = 'DATA: ERROR';

        end

    end

    function s = yesno(tf)

        if tf
            s = 'YES';
        else
            s = 'NO';
        end

    end

%% ========================================================================
% RENDER ORCHESTRATION
% =========================================================================

    function onControlChanged(~,~)

        if isempty(S.recon)
            return;
        end

        try
            renderImagePanels();
        catch ME
            logMsg('ERR',ME.message);
        end

    end

    function renderLinePanels()

        safeRun(@renderUV,'UV');
        safeRun(@renderConvergence,'CONVERGENCE');

        safeRun(@renderAmpUV,'AMP-UV');
        safeRun(@renderClosure,'CLOSURE');
        safeRun(@renderAmpTime,'AMP-TIME');
        safeRun(@renderCPHist,'CP-HIST');

        safeRun(@renderRadial,'RADIAL');
        safeRun(@renderAzimuth,'AZIMUTH');
        safeRun(@renderXsec,'XSEC');
        safeRun(@renderResHist,'RES-HIST');
        safeRun(@renderPsfRadial,'PSF-RADIAL');
        safeRun(@renderFourier,'FOURIER');

        safeRun(@renderPolEvpa,'EVPA');

    end

    function safeRun(fn,label)

        try
            fn();
        catch ME
            logMsg('ERR',[label ': ' ME.message]);
        end

    end

%% ========================================================================
% IMAGE PANELS
% =========================================================================

    function renderImagePanels()

        recon = S.recon;

        % ---- Field of view (half-width, microarcsec) ----------------------

        hmax = max(abs(S.x));

        switch ddFovMode.Value
            case 'AUTO'
                halfEmis = autoHalfFov();
                halfPSF  = autoHalfFovPSF();
                spnFOV.Enable = 'off';
                spnFOV.Value = min(max(round(2*halfEmis),20),20000);
                logMsg('INFO',sprintf('AUTO FOV: %.0f %sas (emission) | %.0f %sas (PSF)', ...
                    2*halfEmis,MU,2*halfPSF,MU));
            case 'FULL'
                halfEmis = hmax;
                halfPSF  = hmax;
                spnFOV.Enable = 'off';
                spnFOV.Value = min(max(round(2*hmax),20),20000);
            otherwise
                halfEmis = min(spnFOV.Value/2,hmax);
                halfPSF  = halfEmis;
                spnFOV.Enable = 'on';
        end

        S.fovHalf = halfEmis;

        cmap    = getCmap(ddCmap.Value);
        divmap  = divergingMap();
        stretch = ddStretch.Value;
        smooth  = chkSmooth.Value;

        S.view = struct();

        dirtyRMS    = getf(recon,'dirtyRMS',NaN);
        restoredRMS = getf(recon,'restoredRMS',NaN);

        renderImage('dirty',recon.dirtyImage, ...
            mkOpts('unit','Jy/beam','beam',true),cmap,stretch,smooth);

        renderImage('model',recon.cleanModel, ...
            mkOpts('unit','Jy','sigma',dirtyRMS),cmap,stretch,smooth);

        renderImage('resid',recon.cleanResidual, ...
            mkOpts('unit','Jy/beam','symmetric',true),divmap,stretch,smooth);

        renderImage('restored',recon.restoredImage, ...
            mkOpts('unit','Jy/beam','sigma',restoredRMS,'peak',true, ...
                   'beam',true,'scalebar',true),cmap,stretch,smooth);

        S.fovHalf = halfPSF;
        renderImage('psf',recon.psf, ...
            mkOpts('unit','','symmetric',true,'halfmax',true),divmap,stretch,smooth);
        S.fovHalf = halfEmis;

        % ---- Polarimetry maps --------------------------------------------

        if hasField(recon,'pol')

            pol = recon.pol;

            renderImage('polI',recon.restoredImage, ...
                mkOpts('unit','Jy/beam','beam',true,'scalebar',true,'ticks',true), ...
                cmap,stretch,smooth);

            renderImage('polP',pol.Pdeb, ...
                mkOpts('unit','Jy/beam','beam',true),cmap,stretch,smooth);

            renderImage('polM',pol.frac, ...
                mkOpts('unit',''),cmap,stretch,smooth);

            renderImage('polQ',pol.Q, ...
                mkOpts('unit','Jy/beam','symmetric',true),divmap,stretch,smooth);

            renderImage('polU',pol.U, ...
                mkOpts('unit','Jy/beam','symmetric',true),divmap,stretch,smooth);

        else

            polKeys = {'polI','polP','polM','polQ','polU'};
            for q = 1:numel(polKeys)
                showPlaceholder(PN.(polKeys{q}).ax,'NO POLARIMETRIC DATA');
                PN.(polKeys{q}).meta.Text = '';
            end

        end

    end

    function o = mkOpts(varargin)

        o = struct('unit','','symmetric',false,'sigma',NaN, ...
                   'peak',false,'beam',false,'scalebar',false, ...
                   'halfmax',false,'ticks',false);

        for q = 1:2:numel(varargin)
            o.(varargin{q}) = varargin{q+1};
        end

    end

    function renderImage(key,img,o,cmap,stretch,smooth)

        P  = PN.(key);
        ax = P.ax;

        cla(ax);

        xv  = S.x;
        idx = find(abs(xv) <= S.fovHalf);

        if numel(idx) < 4
            idx = 1:numel(xv);
        end

        xc  = xv(idx);
        sub = double(img(idx,idx));

        fin  = isfinite(sub);
        vals = sub(fin);

        if isempty(vals)
            showPlaceholder(ax,'NO VALID DATA');
            P.meta.Text = '';
            return;
        end

        % ---- Normalise + stretch ----------------------------------------

        if o.symmetric

            vmax = max(abs(vals));
            if vmax <= 0, vmax = 1; end

            tn = sub/vmax;
            Z  = (sign(tn).*applyStretch(abs(tn),stretch) + 1)/2;

            tt   = linspace(-1,1,7);
            pos  = (sign(tt).*applyStretch(abs(tt),stretch) + 1)/2;
            labv = tt*vmax;

        else

            vmin = min(vals);
            vmax = max(vals);
            if vmax <= vmin, vmax = vmin + 1; end

            t = min(max((sub-vmin)/(vmax-vmin),0),1);
            Z = applyStretch(t,stretch);

            tt   = linspace(0,1,6);
            pos  = applyStretch(tt,stretch);
            labv = vmin + tt*(vmax-vmin);

        end

        % ---- Display image (NaN -> transparent) --------------------------

        xd = xc;
        Zf = Z;
        Zf(~fin) = 0;
        alphaMask = double(fin);

        if smooth && numel(xc) >= 4
            Zd = min(max(interp2(Zf,3,'cubic'),0),1);
            alphaMask = double(interp2(double(fin),3,'nearest') > 0.5);
            xd = linspace(xc(1),xc(end),size(Zd,2));
        else
            Zd = Zf;
        end

        imagesc(ax,xd,xd,Zd,'AlphaData',alphaMask);

        ax.YDir = 'normal';
        axis(ax,'image');

        if CFG.raEastLeft
            ax.XDir = 'reverse';
        else
            ax.XDir = 'normal';
        end

        colormap(ax,cmap);
        clim(ax,[0 1]);

        hold(ax,'on');

        % ---- Reference crosshair at phase centre -------------------------

        plot(ax,[xc(1) xc(end)],[0 0],':','Color',[0.45 0.52 0.62],'LineWidth',0.5);
        plot(ax,[0 0],[xc(1) xc(end)],':','Color',[0.45 0.52 0.62],'LineWidth',0.5);

        % ---- Contours ----------------------------------------------------

        [XX,YY] = meshgrid(xc,xc);

        if isfinite(o.sigma) && o.sigma > 0
            lev = sigmaLevels(o.sigma,max(vals));
            if ~isempty(lev)
                if numel(lev) == 1, lev = [lev lev]; end
                contour(ax,XX,YY,sub,lev,'LineColor',ACCENT_CYAN,'LineWidth',0.7);
            end
        end

        if o.halfmax
            contour(ax,XX,YY,sub,[0.5 0.5],'LineColor',[1 1 1],'LineWidth',1.0);
        end

        % ---- Peak marker -------------------------------------------------

        if o.peak
            [~,ip] = max(img(:));
            [py,px] = ind2sub(size(img),ip);
            plot(ax,xv(px),xv(py),'+','Color',[1 1 1],'MarkerSize',9,'LineWidth',1.4);
        end

        % ---- Overlays ------------------------------------------------------

        if o.scalebar
            addScaleBar(ax,xc);
        end

        if o.beam
            addBeam(ax,xc);
        end

        if o.ticks
            addEvpaTicks(ax,xc,idx);
        end

        hold(ax,'off');

        xlim(ax,[xc(1) xc(end)]);
        ylim(ax,[xc(1) xc(end)]);

        xlabel(ax,['RA offset (' MU 'as)'],'Color',TEXT_SECONDARY,'FontSize',8);
        ylabel(ax,['Dec offset (' MU 'as)'],'Color',TEXT_SECONDARY,'FontSize',8);

        % ---- Colorbar ----------------------------------------------------

        unitTxt = o.unit;
        if isempty(unitTxt), unitTxt = 'a.u.'; end

        cb = P.cb;
        cb.Limits = [0 1];
        cb.Ticks = pos;
        cb.TickLabels = arrayfun(@(q) sprintf('%.3g',q),labv,'UniformOutput',false);
        cb.Label.String = sprintf('%s  [%s]',unitTxt,stretch);
        cb.Label.Color = TEXT_MUTED;
        cb.Label.FontSize = 7;

        % ---- Meta line ---------------------------------------------------

        P.meta.Text = sprintf('MAX %s  |  STD %s', ...
            fmtSI(max(vals),o.unit),fmtSI(std(vals),o.unit));

        % ---- Store for cursor readout -----------------------------------

        S.view.(key) = struct('ax',ax,'x',xc,'img',sub,'unit',o.unit, ...
                              'page',pageOfKey.(key));

    end

    function addScaleBar(ax,xc)

        span = xc(end)-xc(1);

        nice = [1 2 5 10 20 50 100 200 500 1000 2000 5000];
        L = nice(find(nice <= 0.25*span,1,'last'));

        if isempty(L)
            return;
        end

        if CFG.raEastLeft
            sgn = -1;  x0 = xc(end);
        else
            sgn = 1;   x0 = xc(1);
        end

        xs = x0 + sgn*0.06*span;
        xe = xs + sgn*L;
        ys = xc(1) + 0.07*span;

        plot(ax,[xs xe],[ys ys],'-','Color',[1 1 1],'LineWidth',2);

        text(ax,(xs+xe)/2,ys+0.015*span,sprintf('%g %sas',L,MU), ...
            'Color',[1 1 1],'FontName',MONO_FONT,'FontSize',8, ...
            'HorizontalAlignment','center','VerticalAlignment','bottom');

    end

    function addBeam(ax,xc)

        B = S.beam;

        if ~B.valid
            return;
        end

        span = xc(end)-xc(1);

        if 2*B.a > 0.30*span
            return;
        end

        if CFG.raEastLeft
            sgn = -1;  xr = xc(1);
        else
            sgn = 1;   xr = xc(end);
        end

        cx = xr - sgn*(0.08*span + B.a);
        cy = xc(1) + 0.08*span + B.a;

        t = linspace(0,2*pi,72);
        E = B.V * [B.a*cos(t); B.b*sin(t)];

        patch(ax,cx+E(1,:),cy+E(2,:),[1 1 1], ...
            'FaceAlpha',0.18,'EdgeColor',[1 1 1],'LineWidth',1);

        text(ax,cx,cy+B.a+0.015*span,'BEAM', ...
            'Color',[1 1 1],'FontName',MONO_FONT,'FontSize',7, ...
            'HorizontalAlignment','center','VerticalAlignment','bottom');

    end

    function addEvpaTicks(ax,xc,idx)

        pol = S.recon.pol;

        ev = pol.evpa(idx,idx);
        Pd = pol.Pdeb(idx,idx);
        Ii = double(S.recon.restoredImage(idx,idx));

        step = max(1,round(numel(xc)/28));
        ii = 1:step:numel(xc);

        ev = ev(ii,ii);
        Pd = Pd(ii,ii);
        Ii = Ii(ii,ii);

        [XS,YS] = meshgrid(xc(ii),xc(ii));

        ok = isfinite(ev) & Pd > 0 & Ii > 0.05*max(Ii(:));

        if ~any(ok(:))
            return;
        end

        pmax = max(Pd(ok));
        Lmax = 0.9*step*abs(xc(2)-xc(1));

        len = Lmax*(0.35 + 0.65*Pd(ok)/pmax);
        th  = deg2rad(ev(ok));

        dx = 0.5*len.*sin(th);     % EVPA: east of north
        dy = 0.5*len.*cos(th);

        xs = XS(ok);
        ys = YS(ok);

        X = [xs-dx, xs+dx, nan(size(xs))].';
        Y = [ys-dy, ys+dy, nan(size(xs))].';

        plot(ax,X(:),Y(:),'-','Color',[0.45 1.00 0.55],'LineWidth',1.2);

    end

%% ------------------------------------------------------------------------
% FIELD OF VIEW
% -------------------------------------------------------------------------

    function h = autoHalfFov()

        xv   = S.x;
        pix  = S.pixUas;
        hmax = max(abs(xv));
        h    = NaN;

        try
            M = double(S.recon.cleanModel);
            M(~isfinite(M)) = 0;
            M = max(M,0);

            [ri,ci] = find(M > 0);

            if numel(ri) >= 3
                w = M(sub2ind(size(M),ri,ci));
                r = hypot(xv(ci(:)),xv(ri(:)));
                [r,o] = sort(r(:));
                cw = cumsum(w(o));
                if cw(end) > 0
                    k = find(cw >= 0.98*cw(end),1);
                    h = r(k);
                end
            end
        catch
        end

        if ~isfinite(h) || h <= 0
            try
                R = double(S.recon.restoredImage);
                pk = max(R(:));
                [ri,ci] = find(R >= 0.05*pk);
                if ~isempty(ri)
                    r = sort(hypot(xv(ci(:)),xv(ri(:))));
                    h = r(ceil(0.98*numel(r)));
                end
            catch
            end
        end

        if ~isfinite(h) || h <= 0
            h = hmax;
        end

        h = 1.5*h;

        if S.beam.valid
            hmin = 5*S.beam.maj;
        else
            hmin = 16*pix;
        end

        h = min(max([h hmin 8*pix]),hmax);

    end

    function h = autoHalfFovPSF()

        hmax = max(abs(S.x));

        if S.beam.valid
            h = max(8*S.beam.maj,12*S.pixUas);
        else
            h = autoHalfFov();
        end

        h = min(h,hmax);

    end

%% ------------------------------------------------------------------------
% UV PLANE + CONVERGENCE
% -------------------------------------------------------------------------

    function renderUV()

        P  = PN.uv;
        ax = P.ax;

        resetAxes(ax);

        u = S.uv.u;
        v = S.uv.v;

        if isempty(u)
            showPlaceholder(ax,'NO UV DATA');
            P.meta.Text = '';
            return;
        end

        step = max(1,ceil(numel(u)/CFG.maxUVPoints));
        u = u(1:step:end);
        v = v(1:step:end);

        r = hypot(u,v);

        scatter(ax,[u; -u],[v; -v],3,[r; r],'filled','MarkerFaceAlpha',0.65);

        axis(ax,'equal');

        bmax = S.uv.bmax;
        lim  = 1.08*bmax;

        xlim(ax,[-lim lim]);
        ylim(ax,[-lim lim]);

        colormap(ax,getCmap('TURBO'));
        clim(ax,[0 bmax]);

        hold(ax,'on');

        steps = [0.5 1 2 5 10 20];
        st = steps(find(bmax./steps <= 5,1,'first'));
        if isempty(st), st = steps(end); end

        th = linspace(0,2*pi,200);

        for rr = st:st:bmax

            plot(ax,rr*cos(th),rr*sin(th),':', ...
                'Color',[0.50 0.58 0.68],'LineWidth',0.6);

            resUas = 206.264806/rr;   % lambda/B in microarcsec, rr in Glambda

            text(ax,0.707*rr,0.707*rr, ...
                sprintf('%g G%s | %.0f %sas',rr,LAM,resUas,MU), ...
                'Color',TEXT_MUTED,'FontName',MONO_FONT,'FontSize',7, ...
                'VerticalAlignment','bottom');

        end

        plot(ax,[-lim lim],[0 0],'-','Color',[0.28 0.34 0.42],'LineWidth',0.5);
        plot(ax,[0 0],[-lim lim],'-','Color',[0.28 0.34 0.42],'LineWidth',0.5);

        hold(ax,'off');

        xlabel(ax,['u (G' LAM ')'],'Color',TEXT_SECONDARY,'FontSize',8);
        ylabel(ax,['v (G' LAM ')'],'Color',TEXT_SECONDARY,'FontSize',8);

        cb = P.cb;
        cb.Limits = [0 bmax];
        cb.TicksMode = 'auto';
        cb.TickLabelsMode = 'auto';
        cb.Label.String = ['|B|  (G' LAM ')'];
        cb.Label.Color = TEXT_MUTED;
        cb.Label.FontSize = 7;

        P.meta.Text = sprintf('N=%d (%s2 conj.)  |  |B|max %.2f G%s', ...
            numel(S.uv.u),TIMES,bmax,LAM);

    end

    function renderConvergence()

        ax = axConv;
        resetAxes(ax);

        hist = getf(S.recon,'cleanHistory',[]);

        if isempty(hist) || size(hist,2) < 2
            showPlaceholder(ax,'NO CLEAN HISTORY');
            return;
        end

        it = double(hist(:,1));
        pk = double(hist(:,2));

        ok = isfinite(it) & isfinite(pk) & pk > 0;
        it = it(ok);
        pk = pk(ok);

        if isempty(it)
            showPlaceholder(ax,'NO VALID HISTORY');
            return;
        end

        semilogy(ax,it,pk,'-','Color',ACCENT_CYAN,'LineWidth',1.4);

        hold(ax,'on');

        plot(ax,it(end),pk(end),'o', ...
            'MarkerSize',5,'MarkerFaceColor',ACCENT_ORANGE, ...
            'MarkerEdgeColor',[1 1 1],'LineWidth',0.6);

        yl = [min(pk)*0.7, max(pk)*1.4];

        dirtyRMS = getf(S.recon,'dirtyRMS',NaN);
        if isfinite(dirtyRMS) && dirtyRMS > 0
            thr = 3*dirtyRMS;
            plot(ax,[it(1) max(it(end),it(1)+1)],[thr thr],'--', ...
                'Color',ACCENT_ORANGE,'LineWidth',0.9);
            text(ax,it(1),thr,['  3' SIG '_{dirty}'], ...
                'Color',ACCENT_ORANGE,'FontSize',7,'FontName',FONT_NAME, ...
                'VerticalAlignment','bottom','Interpreter','tex');
            yl = [min(yl(1),thr*0.7), max(yl(2),thr*1.4)];
        end

        hold(ax,'off');

        ax.YMinorGrid = 'on';
        ax.XMinorGrid = 'on';

        xlim(ax,[it(1) max(it(end),it(1)+1)]);
        ylim(ax,yl);

        xlabel(ax,'Iteration','Color',TEXT_SECONDARY,'FontSize',8);
        ylabel(ax,'Peak residual (Jy)','Color',TEXT_SECONDARY,'FontSize',8);

        title(ax,sprintf('FINAL PEAK %s  |  %d ITER', ...
            fmtSI(pk(end),'Jy'),round(it(end))), ...
            'Color',TEXT_MUTED,'FontName',FONT_NAME,'FontSize',8, ...
            'FontWeight','normal');

    end

%% ========================================================================
% VISIBILITY PAGE
% =========================================================================

    function renderAmpUV()

        P  = PN.ampuv;
        ax = P.ax;

        resetAxes(ax);

        if ~hasField(S.recon,'visSample') || ~hasField(S.recon,'visRadial')
            showPlaceholder(ax,'NO VISIBILITY TELEMETRY');
            P.meta.Text = '';
            return;
        end

        vs = S.recon.visSample;
        vr = S.recon.visRadial;

        ok = vs.amp > 0 & isfinite(vs.r);

        ax.YScale = 'log';
        hold(ax,'on');

        scatter(ax,vs.r(ok),vs.amp(ok),3,ACCENT_BLUE,'filled', ...
            'MarkerFaceAlpha',0.22,'DisplayName','samples');

        nz = vr.n > 0 & isfinite(vr.amp) & vr.amp > 0;
        lo = min(vr.std(nz),0.9*vr.amp(nz));

        errorbar(ax,vr.r(nz),vr.amp(nz),lo,vr.std(nz),'o', ...
            'Color',ACCENT_ORANGE,'MarkerSize',3,'LineWidth',0.9, ...
            'CapSize',2,'DisplayName',['binned ' char(177) '1' SIG]);

        nm = nz & isfinite(vr.model) & vr.model > 0;
        plot(ax,vr.r(nm),vr.model(nm),'-','Color',ACCENT_GREEN, ...
            'LineWidth',1.6,'DisplayName','CLEAN model');

        hold(ax,'off');

        s = sort(vs.amp(ok));
        lowV = s(max(1,round(0.005*numel(s))));

        xlim(ax,[0 max(vr.r)*1.02]);
        ylim(ax,[max(lowV*0.7,1e-4), max(vs.amp(ok))*1.6]);

        xlabel(ax,['|B| (G' LAM ')'],'Color',TEXT_SECONDARY,'FontSize',8);
        ylabel(ax,'Amplitude (Jy)','Color',TEXT_SECONDARY,'FontSize',8);

        styleLegend(legend(ax,'Location','northeast'));

        rr = vs.amp(ok)./vs.modelAmp(ok);
        rr = rr(isfinite(rr) & rr > 0);

        P.meta.Text = sprintf('%d samples  |  data/model median %.2f', ...
            nnz(ok),median(rr));

    end

    function renderClosure()

        P  = PN.cphase;
        ax = P.ax;

        resetAxes(ax);

        if ~hasField(S.recon,'closure')
            showPlaceholder(ax,'NO CLOSURE DATA');
            P.meta.Text = '';
            return;
        end

        cl  = S.recon.closure;
        hrs = mod(cl.t + 0.5,1)*24;      % UT hours from Julian date

        top = topGroups(cl.tri,6);

        hold(ax,'on');

        for q = 1:numel(top)
            m = cl.tri == top(q);
            scatter(ax,hrs(m),cl.val(m),4,SERIES(q,:),'filled', ...
                'MarkerFaceAlpha',0.55,'DisplayName',cl.triNames{top(q)});
        end

        plot(ax,[0 24],[0 0],'-','Color',[0.40 0.46 0.55], ...
            'LineWidth',0.6,'HandleVisibility','off');

        hold(ax,'off');

        x0 = floor(min(hrs));
        x1 = ceil(max(hrs));
        if x1 <= x0, x1 = x0 + 1; end

        xlim(ax,[x0 x1]);
        ylim(ax,[-185 185]);
        ax.YTick = -180:90:180;

        xlabel(ax,'UT (h)','Color',TEXT_SECONDARY,'FontSize',8);
        ylabel(ax,['Closure phase (' DEG ')'],'Color',TEXT_SECONDARY,'FontSize',8);

        styleLegend(legend(ax,'Location','northeast'));

        P.meta.Text = sprintf('%d pts  |  %d triangles',numel(cl.val),numel(cl.triNames));

    end

    function renderAmpTime()

        P  = PN.amptime;
        ax = P.ax;

        resetAxes(ax);

        if ~hasField(S.recon,'visSample') || ~hasField(S.recon,'blNames')
            showPlaceholder(ax,'NO VISIBILITY TELEMETRY');
            P.meta.Text = '';
            return;
        end

        vs = S.recon.visSample;
        bn = S.recon.blNames;

        ok = vs.amp > 0 & isfinite(vs.t) & vs.t > 0;

        if ~any(ok)
            showPlaceholder(ax,'NO TIMESTAMPS IN FILE');
            P.meta.Text = '';
            return;
        end

        hrs = mod(vs.t + 0.5,1)*24;
        top = topGroups(vs.bl(ok),8);

        ax.YScale = 'log';
        hold(ax,'on');

        for q = 1:numel(top)
            m = ok & vs.bl == top(q);
            scatter(ax,hrs(m),vs.amp(m),4,SERIES(q,:),'filled', ...
                'MarkerFaceAlpha',0.55,'DisplayName',bn{top(q)});
        end

        hold(ax,'off');

        x0 = floor(min(hrs(ok)));
        x1 = ceil(max(hrs(ok)));
        if x1 <= x0, x1 = x0 + 1; end

        xlim(ax,[x0 x1]);

        xlabel(ax,'UT (h)','Color',TEXT_SECONDARY,'FontSize',8);
        ylabel(ax,'Amplitude (Jy)','Color',TEXT_SECONDARY,'FontSize',8);

        styleLegend(legend(ax,'Location','northeast'));

        P.meta.Text = sprintf('%d baselines shown of %d',numel(top),numel(bn));

    end

    function renderCPHist()

        P  = PN.cphist;
        ax = P.ax;

        resetAxes(ax);

        if ~hasField(S.recon,'closure')
            showPlaceholder(ax,'NO CLOSURE DATA');
            P.meta.Text = '';
            return;
        end

        cl = S.recon.closure;

        t1   = topGroups(cl.tri,1);
        vals = cl.val(cl.tri == t1);

        histogram(ax,vals,36,'BinLimits',[-180 180], ...
            'FaceColor',ACCENT_CYAN,'EdgeColor','none','FaceAlpha',0.80);

        z    = mean(exp(1i*deg2rad(vals)));
        mu   = rad2deg(angle(z));
        Rlen = abs(z);
        sc   = rad2deg(sqrt(-2*log(max(Rlen,eps))));

        xlim(ax,[-180 180]);
        ax.XTick = -180:60:180;

        xlabel(ax,['Closure phase (' DEG ')'],'Color',TEXT_SECONDARY,'FontSize',8);
        ylabel(ax,'Count','Color',TEXT_SECONDARY,'FontSize',8);

        P.meta.Text = sprintf('%s | mean %+.1f%s | circ. %s %.1f%s | N=%d', ...
            cl.triNames{t1},mu,DEG,SIG,sc,DEG,numel(vals));

    end

%% ========================================================================
% ANALYSIS PAGE
% =========================================================================

    function A = analyzeRing()

        A = struct('valid',false,'cx',0,'cy',0,'r0',NaN,'fwhm',NaN, ...
                   'lo',NaN,'hi',NaN,'rb',[],'prof',[],'ps',[], ...
                   'pa',[],'azi',[],'asym',NaN,'beta1',NaN,'paPeak',NaN);

        try

            R = double(S.recon.restoredImage);
            R(~isfinite(R)) = 0;

            xv  = S.x;
            pix = S.pixUas;

            rWin = autoHalfFov();
            idx  = find(abs(xv) <= rWin);

            if numel(idx) < 8
                return;
            end

            xc  = xv(idx);
            sub = R(idx,idx);

            [XX,YY] = meshgrid(xc,xc);

            w = max(sub - 0.10*max(sub(:)),0);
            if sum(w(:)) <= 0
                return;
            end

            A.cx = sum(w(:).*XX(:))/sum(w(:));
            A.cy = sum(w(:).*YY(:))/sum(w(:));

            rr = hypot(XX-A.cx,YY-A.cy);

            nb    = max(8,floor(rWin/pix));
            edges = (0:nb)*pix;

            bi = discretize(rr(:),edges);
            ok = ~isnan(bi);

            prof = accumarray(bi(ok),sub(ok),[nb 1],@mean,NaN);
            rb   = 0.5*(edges(1:end-1)+edges(2:end)).';

            ps = movmean(prof,3,'omitnan');

            [pk,ip] = max(ps);

            A.rb   = rb;
            A.prof = prof;
            A.ps   = ps;

            if ip <= 1 || ip >= nb
                A.valid = true;       % profile valid, ring not resolved
                return;
            end

            y1 = ps(ip-1); y2 = ps(ip); y3 = ps(ip+1);
            den = y1 - 2*y2 + y3;
            dlt = 0;
            if den ~= 0
                dlt = 0.5*(y1 - y3)/den;
            end
            A.r0 = rb(ip) + dlt*pix;

            half = 0.5*pk;

            lo = rb(1);
            for q = ip-1:-1:1
                if ps(q) < half
                    lo = rb(q) + (half-ps(q))/(ps(q+1)-ps(q))*pix;
                    break;
                end
            end

            hi = rb(end);
            for q = ip+1:nb
                if ps(q) < half
                    hi = rb(q-1) + (half-ps(q-1))/(ps(q)-ps(q-1))*pix;
                    break;
                end
            end

            A.lo = lo;
            A.hi = hi;
            A.fwhm = hi - lo;

            % ---- azimuthal profile in the ring annulus ---------------------

            ann = abs(rr - A.r0) <= max(A.fwhm/2,1.5*pix);
            pa  = mod(atan2d(XX-A.cx,YY-A.cy),360);

            eAz = 0:10:360;
            bj  = discretize(pa(ann),eAz);
            vv  = sub(ann);
            okj = ~isnan(bj);

            A.azi = accumarray(bj(okj),vv(okj),[36 1],@mean,NaN);
            A.pa  = 0.5*(eAz(1:end-1)+eAz(2:end)).';

            mx = max(A.azi);
            mn = min(A.azi);
            A.asym = (mx-mn)/(mx+mn);

            zc = sum(vv.*exp(1i*deg2rad(pa(ann))));
            A.beta1  = abs(zc)/sum(vv);
            A.paPeak = mod(rad2deg(angle(zc)),360);

            A.valid = true;

        catch
        end

    end

    function renderRadial()

        P  = PN.radial;
        ax = P.ax;

        resetAxes(ax);

        A = S.ring;

        if ~A.valid || isempty(A.rb)
            showPlaceholder(ax,'NO RING PROFILE');
            P.meta.Text = '';
            return;
        end

        yl = [min(0,min(A.ps)), max(A.ps)*1.15];

        hold(ax,'on');

        plot(ax,A.rb,A.prof,'.','Color',[0.50 0.60 0.72],'MarkerSize',6);
        plot(ax,A.rb,A.ps,'-','Color',ACCENT_CYAN,'LineWidth',1.5);

        if isfinite(A.r0)

            patch(ax,[A.lo A.hi A.hi A.lo],[yl(1) yl(1) yl(2) yl(2)], ...
                ACCENT_ORANGE,'FaceAlpha',0.10,'EdgeColor','none');

            plot(ax,[A.r0 A.r0],yl,'--','Color',ACCENT_ORANGE,'LineWidth',1.0);

        end

        if S.beam.valid
            by = 0.92*yl(2);
            plot(ax,[0 S.beam.maj],[by by],'-','Color',[1 1 1],'LineWidth',2);
            text(ax,S.beam.maj/2,by,'BEAM','Color',[1 1 1],'FontSize',7, ...
                'FontName',MONO_FONT,'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom');
        end

        hold(ax,'off');

        xlim(ax,[0 max(A.rb)]);
        ylim(ax,yl);

        xlabel(ax,['Radius from ring centre (' MU 'as)'],'Color',TEXT_SECONDARY,'FontSize',8);
        ylabel(ax,'Azimuthal mean I (Jy/beam)','Color',TEXT_SECONDARY,'FontSize',8);

        if isfinite(A.r0)
            P.meta.Text = sprintf('D = %.1f %sas  |  FWHM %.1f %sas (beam-convolved)', ...
                2*A.r0,MU,A.fwhm,MU);
        else
            P.meta.Text = 'RING NOT RESOLVED';
        end

    end

    function renderAzimuth()

        P  = PN.azimuth;
        ax = P.ax;

        resetAxes(ax);

        A = S.ring;

        if ~A.valid || ~isfinite(A.r0) || isempty(A.azi)
            showPlaceholder(ax,'RING NOT RESOLVED');
            P.meta.Text = '';
            return;
        end

        hold(ax,'on');

        plot(ax,[0 360],[mean(A.azi,'omitnan') mean(A.azi,'omitnan')],':', ...
            'Color',[0.55 0.62 0.72],'LineWidth',0.8);
        plot(ax,A.pa,A.azi,'-o','Color',ACCENT_MAGENTA,'LineWidth',1.4, ...
            'MarkerSize',3,'MarkerFaceColor',ACCENT_MAGENTA);

        hold(ax,'off');

        xlim(ax,[0 360]);
        ax.XTick = 0:60:360;

        xlabel(ax,['Position angle, E of N (' DEG ')'],'Color',TEXT_SECONDARY,'FontSize',8);
        ylabel(ax,'Mean I in ring (Jy/beam)','Color',TEXT_SECONDARY,'FontSize',8);

        P.meta.Text = sprintf('A = %.2f  |  %s1 = %.2f  |  peak PA %.0f%s', ...
            A.asym,char(946),A.beta1,A.paPeak,DEG);

    end

    function renderXsec()

        P  = PN.xsec;
        ax = P.ax;

        resetAxes(ax);

        R = double(S.recon.restoredImage);
        R(~isfinite(R)) = 0;

        [pk,ip] = max(R(:));
        [py,px] = ind2sub(size(R),ip);

        hw  = autoHalfFov();
        sel = abs(S.x) <= hw;
        xc  = S.x(sel);

        hx = R(py,sel);
        vy = R(sel,px).';

        hold(ax,'on');

        plot(ax,xc,hx,'-','Color',ACCENT_ORANGE,'LineWidth',1.4,'DisplayName','RA cut (row)');
        plot(ax,xc,vy,'-','Color',ACCENT_CYAN,'LineWidth',1.4,'DisplayName','Dec cut (column)');

        if S.beam.valid
            bm = S.beam.maj;
            plot(ax,S.x(px)+[-0.5 0.5]*bm,[0.5 0.5]*pk,'-','Color',[1 1 1], ...
                'LineWidth',2,'HandleVisibility','off');
            text(ax,S.x(px),0.5*pk,'  beam FWHM','Color',[1 1 1],'FontSize',7, ...
                'FontName',MONO_FONT,'VerticalAlignment','bottom');
        end

        hold(ax,'off');

        xlim(ax,[xc(1) xc(end)]);

        xlabel(ax,['Offset (' MU 'as)'],'Color',TEXT_SECONDARY,'FontSize',8);
        ylabel(ax,'I (Jy/beam)','Color',TEXT_SECONDARY,'FontSize',8);

        styleLegend(legend(ax,'Location','northeast'));

        P.meta.Text = sprintf('PEAK %s @ (%.1f, %.1f) %sas', ...
            fmtSI(pk,'Jy/beam'),S.x(px),S.x(py),MU);

    end

    function renderResHist()

        P  = PN.reshist;
        ax = P.ax;

        resetAxes(ax);

        res = double(S.recon.cleanResidual);

        xv = S.x(:);
        R2 = xv.^2 + xv.'.^2;

        mask = R2 > 300^2 & R2 < 1000^2;

        if nnz(mask) < 1000
            mask = true(size(res));
        end

        vals = res(mask);
        vals = vals(isfinite(vals));

        step = max(1,ceil(numel(vals)/500000));
        vals = vals(1:step:end);

        mu  = median(vals);
        sg  = 1.4826*median(abs(vals-mu));
        if ~isfinite(sg) || sg <= 0, sg = std(vals); end

        ax.YScale = 'log';
        hold(ax,'on');

        histogram(ax,vals,120,'Normalization','pdf', ...
            'BinLimits',mu+[-6 6]*sg, ...
            'FaceColor',ACCENT_BLUE,'EdgeColor','none','FaceAlpha',0.75);

        xg = linspace(mu-6*sg,mu+6*sg,240);
        g  = exp(-0.5*((xg-mu)/sg).^2)/(sg*sqrt(2*pi));

        plot(ax,xg,g,'-','Color',ACCENT_ORANGE,'LineWidth',1.5);

        hold(ax,'off');

        pkg = 1/(sg*sqrt(2*pi));
        xlim(ax,mu+[-6 6]*sg);
        ylim(ax,[pkg*1e-4 pkg*2]);

        z  = (vals-mean(vals))/std(vals);
        sk = mean(z.^3);
        ku = mean(z.^4) - 3;

        xlabel(ax,'Residual (Jy/beam)','Color',TEXT_SECONDARY,'FontSize',8);
        ylabel(ax,'PDF','Color',TEXT_SECONDARY,'FontSize',8);

        P.meta.Text = sprintf('%s_{rob} %s | skew %.2f | ex.kurt %.2f', ...
            SIG,fmtSI(sg,'Jy/beam'),sk,ku);

    end

    function renderPsfRadial()

        P  = PN.psfrad;
        ax = P.ax;

        resetAxes(ax);

        psf = double(S.recon.psf);

        hw  = min(1500,max(abs(S.x)));
        sel = abs(S.x) <= hw;
        xc  = S.x(sel);
        sub = abs(psf(sel,sel));

        [XX,YY] = meshgrid(xc,xc);
        r = hypot(XX,YY);

        pix = S.pixUas;
        nb  = floor(hw/pix);

        edges = (0:nb)*pix;
        bi = discretize(r(:),edges);
        ok = ~isnan(bi);

        env = accumarray(bi(ok),sub(ok),[nb 1],@max,NaN);
        rc  = 0.5*(edges(1:end-1)+edges(2:end)).';

        dB = 20*log10(max(env,1e-4));

        hold(ax,'on');

        plot(ax,rc,dB,'-','Color',ACCENT_CYAN,'LineWidth',1.4);

        slTxt = '';

        if S.beam.valid

            r15 = 1.5*S.beam.maj;
            plot(ax,[r15 r15],[-80 5],'--','Color',ACCENT_ORANGE,'LineWidth',0.9);

            far = rc > r15 & isfinite(dB);
            if any(far)
                slTxt = sprintf('PEAK SIDELOBE %.1f dB (r > 1.5 FWHM)',max(dB(far)));
            end

        end

        hold(ax,'off');

        xlim(ax,[0 hw]);
        ylim(ax,[max(-70,min(dB)-3) 3]);

        xlabel(ax,['Radius (' MU 'as)'],'Color',TEXT_SECONDARY,'FontSize',8);
        ylabel(ax,'|PSF| envelope (dB)','Color',TEXT_SECONDARY,'FontSize',8);

        P.meta.Text = slTxt;

    end

    function renderFourier()

        P  = PN.fourier;
        ax = P.ax;

        cla(ax);

        R = double(S.recon.restoredImage);
        R(~isfinite(R)) = 0;

        N = size(R,1);
        c = N/2 + 1;

        Vf = fftshift(fft2(ifftshift(R)));

        duG = S.du/1e3;                    % uv cell size in Glambda

        bmax = S.uv.bmax;
        if ~isfinite(bmax), bmax = 10; end
        lim = 1.15*bmax;

        k  = ((1:N)-c)*duG;
        ik = abs(k) <= lim;
        kc = k(ik);

        A = abs(Vf(ik,ik));
        A = log10(max(A,max(A(:))*1e-5));

        imagesc(ax,kc,kc,A);

        ax.YDir = 'normal';
        ax.XDir = 'normal';
        axis(ax,'image');

        colormap(ax,getCmap('TURBO'));
        clim(ax,[min(A(:)) max(A(:))]);

        hold(ax,'on');

        if ~isempty(S.uv.u)
            step = max(1,ceil(numel(S.uv.u)/20000));
            uu = S.uv.u(1:step:end);
            vv = S.uv.v(1:step:end);
            scatter(ax,[uu; -uu],[vv; -vv],1.5,[1 1 1],'filled', ...
                'MarkerFaceAlpha',0.18);
        end

        hold(ax,'off');

        xlim(ax,[kc(1) kc(end)]);
        ylim(ax,[kc(1) kc(end)]);

        xlabel(ax,['u (G' LAM ')'],'Color',TEXT_SECONDARY,'FontSize',8);
        ylabel(ax,['v (G' LAM ')'],'Color',TEXT_SECONDARY,'FontSize',8);

        cb = P.cb;
        cb.TicksMode = 'auto';
        cb.TickLabelsMode = 'auto';
        cb.Label.String = 'log10 |I~(u,v)|  (Jy)';
        cb.Label.Color = TEXT_MUTED;
        cb.Label.FontSize = 7;

        P.meta.Text = sprintf('|u|,|v| <= %.1f G%s  |  uv samples overlaid',lim,LAM);

    end

%% ========================================================================
% POLARIMETRY PAGE (scatter panel)
% =========================================================================

    function renderPolEvpa()

        P  = PN.polEvpa;
        ax = P.ax;

        resetAxes(ax);

        if ~hasField(S.recon,'pol')
            showPlaceholder(ax,'NO POLARIMETRIC DATA');
            P.meta.Text = '';
            return;
        end

        pol = S.recon.pol;

        hw  = autoHalfFov();
        sel = abs(S.x) <= hw;
        xc  = S.x(sel);

        ev = pol.evpa(sel,sel);
        Pd = pol.Pdeb(sel,sel);

        [XX,YY] = meshgrid(xc,xc);

        A = S.ring;
        cx = 0; cy = 0;
        if A.valid
            cx = A.cx;
            cy = A.cy;
        end

        pa = mod(atan2d(XX-cx,YY-cy),360);

        m = isfinite(ev) & Pd > 0;

        if ~any(m(:))
            showPlaceholder(ax,'NO SIGNIFICANT POLARIZATION');
            P.meta.Text = '';
            return;
        end

        paV = pa(m);
        evV = ev(m);
        pdV = Pd(m);

        step = max(1,ceil(numel(paV)/6000));
        paV = paV(1:step:end);
        evV = evV(1:step:end);
        pdV = pdV(1:step:end);

        scatter(ax,paV,evV,9,pdV,'filled','MarkerFaceAlpha',0.8);

        colormap(ax,getCmap('TURBO'));
        clim(ax,[0 max(pdV)]);

        xlim(ax,[0 360]);
        ylim(ax,[-90 90]);
        ax.XTick = 0:60:360;
        ax.YTick = -90:45:90;

        xlabel(ax,['Position angle, E of N (' DEG ')'],'Color',TEXT_SECONDARY,'FontSize',8);
        ylabel(ax,['EVPA (' DEG ')'],'Color',TEXT_SECONDARY,'FontSize',8);

        cb = P.cb;
        cb.TicksMode = 'auto';
        cb.TickLabelsMode = 'auto';
        cb.Label.String = '|P|  (Jy/beam)';
        cb.Label.Color = TEXT_MUTED;
        cb.Label.FontSize = 7;

        P.meta.Text = sprintf('EVPA mean %+.0f%s  |  m_{net} %.1f %%', ...
            pol.evpaMean,DEG,100*pol.mNet);

    end

%% ========================================================================
% KPI STRIP
% =========================================================================

    function renderKPIs()

        recon = S.recon;

        restored = double(recon.restoredImage);

        peak = getf(recon,'restoredPeak',max(restored(:)));
        rmsv = getf(recon,'restoredRMS',NaN);
        dr   = getf(recon,'restoredDynamicRange',peak/rmsv);
        flux = sum(double(recon.cleanModel(:)),'omitnan');

        KPI.val{1}.Text = fmtSI(peak,'Jy/beam');
        KPI.sub{1}.Text = 'RESTORED IMAGE';

        KPI.val{2}.Text = fmtSI(rmsv,'Jy/beam');
        KPI.sub{2}.Text = 'OFF-SOURCE ESTIMATE';

        KPI.val{3}.Text = sprintf('%.1f',dr);
        KPI.sub{3}.Text = 'PEAK / RMS';

        KPI.val{4}.Text = fmtSI(flux,'Jy');
        KPI.sub{4}.Text = 'SUM OF COMPONENTS';

        if S.beam.valid
            KPI.val{5}.Text = sprintf('%.1f%s%.1f %sas', ...
                S.beam.maj,TIMES,S.beam.min,MU);
            KPI.sub{5}.Text = sprintf('PA %+.0f%s  |  FWHM',S.beam.pa,DEG);
        else
            KPI.val{5}.Text = '--';
            KPI.sub{5}.Text = 'FIT FAILED';
        end

        nUV = numel(S.uv.u);
        if nUV > 0
            KPI.val{6}.Text = sprintf('%d',nUV);
            KPI.sub{6}.Text = 'PLOTTED POINTS';
            KPI.val{7}.Text = sprintf('%.2f G%s',S.uv.bmax,LAM);
            KPI.sub{7}.Text = sprintf('%s/B = %.1f %sas',LAM,206.264806/S.uv.bmax,MU);
        else
            KPI.val{6}.Text = '--';
            KPI.sub{6}.Text = 'NO UV DATA';
            KPI.val{7}.Text = '--';
            KPI.sub{7}.Text = '';
        end

        hist = getf(recon,'cleanHistory',[]);
        if ~isempty(hist) && size(hist,2) >= 2
            KPI.val{8}.Text = sprintf('%d',round(max(hist(:,1))));
            KPI.sub{8}.Text = ['FINAL ' fmtSI(hist(end,2),'Jy')];
        else
            KPI.val{8}.Text = '--';
            KPI.sub{8}.Text = 'NO HISTORY';
        end

        A = S.ring;
        if A.valid && isfinite(A.r0)
            KPI.val{9}.Text = sprintf('%.1f %sas',2*A.r0,MU);
            KPI.sub{9}.Text = sprintf('FWHM %.1f %sas | BEAM-CONVOLVED',A.fwhm,MU);
        else
            KPI.val{9}.Text = '--';
            KPI.sub{9}.Text = 'NOT RESOLVED';
        end

        if hasField(recon,'pol')
            KPI.val{10}.Text = sprintf('%.1f %%',100*recon.pol.mNet);
            KPI.sub{10}.Text = sprintf('EVPA %+.0f%s',recon.pol.evpaMean,DEG);
        else
            KPI.val{10}.Text = '--';
            KPI.sub{10}.Text = 'NO POL DATA';
        end

    end

    function [valLbl,subLbl] = makeKPI(parent,col,titleTxt,accent)

        card = makePanel(parent,CARD_BG);
        card.Layout.Column = col;

        g = uigridlayout(card,[3 1]);
        g.RowHeight = {14,'1x',14};
        g.ColumnWidth = {'1x'};
        g.BackgroundColor = CARD_BG;
        g.Padding = [10 4 10 4];
        g.RowSpacing = 0;

        t = uilabel(g,'Text',titleTxt, ...
            'FontName',FONT_NAME,'FontSize',8,'FontWeight','bold', ...
            'FontColor',TEXT_MUTED);
        t.Layout.Row = 1;

        valLbl = uilabel(g,'Text','--', ...
            'FontName',MONO_FONT,'FontSize',14,'FontWeight','bold', ...
            'FontColor',accent);
        valLbl.Layout.Row = 2;

        subLbl = uilabel(g,'Text','AWAITING DATA', ...
            'FontName',FONT_NAME,'FontSize',7, ...
            'FontColor',TEXT_MUTED);
        subLbl.Layout.Row = 3;

    end

%% ========================================================================
% CURSOR READOUT
% =========================================================================

    function onMouseMove(~,~)

        if isempty(S.recon) || toc(S.tmove) < 0.05
            return;
        end

        S.tmove = tic;

        try

            fp   = fig.CurrentPoint;
            keys = fieldnames(S.view);

            for ii = 1:numel(keys)

                V = S.view.(keys{ii});

                if V.page ~= S.page
                    continue;
                end

                pp = getpixelposition(V.ax,true);

                inside = fp(1) >= pp(1) && fp(1) <= pp(1)+pp(3) && ...
                         fp(2) >= pp(2) && fp(2) <= pp(2)+pp(4);

                if ~inside || numel(V.x) < 2
                    continue;
                end

                cp = V.ax.CurrentPoint;
                xp = cp(1,1);
                yp = cp(1,2);

                dxp = V.x(2)-V.x(1);
                n   = numel(V.x);

                ix = round((xp-V.x(1))/dxp) + 1;
                iy = round((yp-V.x(1))/dxp) + 1;

                if ix >= 1 && ix <= n && iy >= 1 && iy <= n

                    footer.Text = sprintf( ...
                        '  CURSOR [%-8s]   RA %+8.1f %sas   DEC %+8.1f %sas   |   VALUE %s   |   r = %.1f %sas', ...
                        upper(keys{ii}),xp,MU,yp,MU, ...
                        fmtSI(V.img(iy,ix),V.unit),hypot(xp,yp),MU);

                    return;

                end

            end

            footer.Text = idleFooter;

        catch
        end

    end

%% ========================================================================
% EXPORT
% =========================================================================

    function exportFigure(~,~)

        [f,p] = uiputfile({'*.png','PNG image'}, ...
            'Export dashboard','M87_EHT_dashboard.png');

        if isequal(f,0)
            return;
        end

        try
            exportapp(fig,fullfile(p,f));
            logMsg('OK',['EXPORTED: ' f]);
        catch ME
            logMsg('ERR',['EXPORT FAILED: ' ME.message]);
        end

    end

%% ========================================================================
% BEAM
% =========================================================================

    function B = beamFromRecon(recon)

        B = struct('valid',false);

        if ~isfield(recon,'beam') || ~isstruct(recon.beam) || ...
                ~isfield(recon.beam,'major')
            return;
        end

        maj = double(recon.beam.major)*1000;     % mas -> microarcsec
        mn  = double(recon.beam.minor)*1000;
        th  = deg2rad(double(recon.beam.pa));

        if ~(isfinite(maj) && isfinite(mn) && maj > 0 && mn > 0 && isfinite(th))
            return;
        end

        v1 = [cos(th); sin(th)];
        v2 = [-sin(th); cos(th)];

        pa = atan2d(v1(1),v1(2));
        if pa > 90,   pa = pa - 180; end
        if pa <= -90, pa = pa + 180; end

        B = struct('valid',true,'maj',maj,'min',mn,'pa',pa, ...
                   'a',maj/2,'b',mn/2,'V',[v1 v2]);

    end

    function B = fitBeam(psf,x)

        B = struct('valid',false,'maj',NaN,'min',NaN,'pa',NaN, ...
                   'a',NaN,'b',NaN,'V',eye(2));

        try

            psf = double(psf);

            [~,ip] = max(psf(:));
            [r0,c0] = ind2sub(size(psf),ip);

            dx = abs(x(2)-x(1));
            hw = max(4,round(60/dx));

            r = max(1,r0-hw):min(size(psf,1),r0+hw);
            c = max(1,c0-hw):min(size(psf,2),c0+hw);

            sub = psf(r,c);
            [CX,CY] = meshgrid(x(c),x(r));

            M = sub >= 0.5*max(sub(:));

            seed = false(size(M));
            seed(r0-r(1)+1,c0-c(1)+1) = true;

            for it = 1:(2*hw)
                nxt = (conv2(double(seed),ones(3),'same') > 0) & M;
                if isequal(nxt,seed), break; end
                seed = nxt;
            end

            if nnz(seed) < 5
                return;
            end

            X = CX(seed);
            Y = CY(seed);

            X = X - mean(X);
            Y = Y - mean(Y);

            C = [X Y].'*[X Y]/numel(X);

            [V,D] = eig(C);
            [lam,order] = sort(diag(D),'descend');
            V = V(:,order);

            lam = max(lam - dx^2/12,0);

            B.maj = 4*sqrt(lam(1));
            B.min = 4*sqrt(lam(2));
            B.a = B.maj/2;
            B.b = B.min/2;
            B.V = V;

            pa = atan2d(V(1,1),V(2,1));
            if pa > 90,   pa = pa - 180; end
            if pa <= -90, pa = pa + 180; end
            B.pa = pa;

            B.valid = isfinite(B.maj) && isfinite(B.min) && B.maj > 0;

        catch
        end

    end

%% ========================================================================
% TELESCOPE MAP
% =========================================================================

    function drawTelescopeMap(gx)

        cla(gx);

        stations = {
            'ALMA', -23.02, -67.75;
            'APEX', -23.00, -67.76;
            'JCMT',  19.82,-155.47;
            'SMA',   19.82,-155.48;
            'LMT',   18.98, -97.31;
            'SMT',   32.70,-109.88;
            'IRAM',  37.07,  -3.39;
            'PdB',   44.54,   5.91;
            'SPT',  -84.90,   0.00
        };

        names = stations(:,1);
        lats  = cell2mat(stations(:,2));
        lons  = cell2mat(stations(:,3));

        geoscatter(gx,lats,lons,46,'filled', ...
            'MarkerFaceColor',ACCENT_CYAN, ...
            'MarkerEdgeColor',[1 1 1],'LineWidth',0.6);

        hold(gx,'on');

        for ii = 2:numel(lats)
            geoplot(gx,[lats(1) lats(ii)],[lons(1) lons(ii)],'-', ...
                'Color',[0.15 0.40 0.68],'LineWidth',0.8);
        end

        try
            for ii = 1:numel(lats)
                text(gx,lats(ii),lons(ii),['  ' names{ii}], ...
                    'Color',TEXT_PRIMARY,'FontName',FONT_NAME,'FontSize',8, ...
                    'FontWeight','bold');
            end
        catch
        end

        hold(gx,'off');

        try
            geolimits(gx, ...
                [max(-80,min(lats)-10), min(80,max(lats)+10)], ...
                [min(lons)-15, max(lons)+15]);
        catch
        end

    end

%% ========================================================================
% UI HELPERS
% =========================================================================

    function P = createCard(parentGrid,rowIdx,colIdx,titleTxt,hasCbar)

        card = makePanel(parentGrid,CARD_BG);
        card.Layout.Row = rowIdx;
        card.Layout.Column = colIdx;

        cg = uigridlayout(card,[2 1]);
        cg.RowHeight = {22,'1x'};
        cg.ColumnWidth = {'1x'};
        cg.BackgroundColor = CARD_BG;
        cg.Padding = [7 6 7 6];
        cg.RowSpacing = 3;

        hh = uigridlayout(cg,[1 2]);
        hh.Layout.Row = 1;
        hh.ColumnWidth = {'1x',300};
        hh.BackgroundColor = CARD_BG;
        hh.Padding = [0 0 0 0];
        hh.ColumnSpacing = 4;

        lbl = uilabel(hh,'Text',['  ' upper(titleTxt)], ...
            'FontName',FONT_NAME,'FontSize',9,'FontWeight','bold', ...
            'FontColor',TEXT_SECONDARY);
        lbl.Layout.Column = 1;

        meta = uilabel(hh,'Text','', ...
            'FontName',MONO_FONT,'FontSize',7, ...
            'FontColor',TEXT_MUTED,'HorizontalAlignment','right');
        meta.Layout.Column = 2;

        ax = uiaxes(cg);
        ax.Layout.Row = 2;
        formatAxes(ax);

        P.ax   = ax;
        P.meta = meta;
        P.cb   = [];

        if hasCbar
            cb = colorbar(ax,'eastoutside');
            cb.Color = TEXT_MUTED;
            cb.FontSize = 7;
            cb.FontName = FONT_NAME;
            cb.LineWidth = 0.6;
            P.cb = cb;
        end

    end

    function p = makePanel(parent,bg)

        p = uipanel(parent,'BackgroundColor',bg,'BorderType','line');

        try
            p.BorderColor = BORDER_COL;
        catch
            try
                p.HighlightColor = BORDER_COL;
            catch
            end
        end

    end

    function formatAxes(ax)

        ax.Color = [0.020 0.027 0.038];

        ax.XColor = TEXT_MUTED;
        ax.YColor = TEXT_MUTED;

        ax.FontName = FONT_NAME;
        ax.FontSize = 8;

        ax.Box = 'on';
        ax.LineWidth = 0.8;
        ax.TickDir = 'in';
        ax.Layer = 'top';
        ax.TickLength = [0.012 0.025];
        ax.XMinorTick = 'on';
        ax.YMinorTick = 'on';

        ax.GridColor = [0.35 0.42 0.52];
        ax.GridAlpha = 0.18;
        ax.GridLineStyle = ':';

        grid(ax,'on');

    end

    function resetAxes(ax)

        cla(ax);

        try
            legend(ax,'off');
        catch
        end

        ax.XScale = 'linear';
        ax.YScale = 'linear';
        ax.XDir = 'normal';
        ax.YDir = 'normal';
        ax.XLimMode = 'auto';
        ax.YLimMode = 'auto';
        ax.XTickMode = 'auto';
        ax.YTickMode = 'auto';
        ax.DataAspectRatioMode = 'auto';

        formatAxes(ax);

    end

    function styleLegend(lg)

        try
            lg.TextColor = TEXT_SECONDARY;
            lg.Color = CARD_BG;
            lg.EdgeColor = BORDER_COL;
            lg.FontSize = 7;
            lg.Box = 'on';
        catch
        end

    end

    function showPlaceholder(ax,msg)

        if nargin < 2
            msg = 'AWAITING DATA';
        end

        cla(ax);

        try
            legend(ax,'off');
        catch
        end

        ax.XScale = 'linear';
        ax.YScale = 'linear';

        formatAxes(ax);

        ax.XLim = [0 1];
        ax.YLim = [0 1];

        text(ax,0.5,0.5,msg,'Units','normalized', ...
            'HorizontalAlignment','center', ...
            'Color',TEXT_MUTED,'FontName',MONO_FONT,'FontSize',9);

    end

    function setStatus(txt,col)

        statusLabel.Text = [char(9679) ' ' txt];
        statusLabel.FontColor = col;

    end

    function logMsg(level,txt)

        stamp = datestr(now,'HH:MM:SS.FFF');
        line  = sprintf('[%s] %-4s %s',stamp,level,txt);

        cur = infoText.Value;

        if numel(cur) == 1 && isempty(cur{1})
            cur = {};
        end

        cur = [cur(:); {line}];

        if numel(cur) > 200
            cur = cur(end-199:end);
        end

        infoText.Value = cur;

        try
            scroll(infoText,'bottom');
        catch
        end

        drawnow limitrate;

    end

%% ========================================================================
% NUMERIC / DISPLAY HELPERS
% =========================================================================

    function v = getf(s,name,default)

        if isfield(s,name) && ~isempty(s.(name))
            v = s.(name);
        else
            v = default;
        end

    end

    function tf = hasField(s,name)

        tf = isstruct(s) && isfield(s,name) && ~isempty(s.(name));

    end

    function top = topGroups(g,n)

        cnt = accumarray(g(:),1);
        [~,ord] = sort(cnt,'descend');
        top = ord(1:min(n,numel(ord)));
        top = top(cnt(top) > 0);

    end

    function y = applyStretch(t,name)

        switch name
            case 'SQRT'
                y = sqrt(t);
            case 'LOG'
                y = log10(1 + 1000*t)/log10(1001);
            case 'ASINH'
                y = asinh(10*t)/asinh(10);
            otherwise
                y = t;
        end

    end

    function lev = sigmaLevels(sigma,peak)

        % Standard radio-astronomy contours: 3*sigma * 2^n, below the peak
        lev = [];
        L = 3*sigma;

        while L < peak && numel(lev) < 12
            lev(end+1) = L; %#ok<AGROW>
            L = 2*L;
        end

    end

    function s = fmtSI(val,unit)

        if ~isfinite(val)
            s = '--';
            return;
        end

        if isempty(unit)
            s = sprintf('%.3g',val);
            return;
        end

        a = abs(val);

        if a == 0
            s = sprintf('0 %s',unit);
        elseif a >= 1
            s = sprintf('%.3g %s',val,unit);
        elseif a >= 1e-3
            s = sprintf('%.3g m%s',val*1e3,unit);
        elseif a >= 1e-6
            s = sprintf('%.3g %s%s',val*1e6,MU,unit);
        else
            s = sprintf('%.3g n%s',val*1e9,unit);
        end

        s = strtrim(s);

    end

%% ========================================================================
% COLORMAPS
% =========================================================================

    function cmap = getCmap(name)

        switch name
            case 'HOT'
                cmap = hot(256);
            case 'TURBO'
                try
                    cmap = turbo(256);
                catch
                    cmap = parula(256);
                end
            case 'PARULA'
                cmap = parula(256);
            case 'GRAY'
                cmap = gray(256);
            otherwise
                cmap = ehtMap();
        end

    end

    function cmap = ehtMap()

        stops = [
            0.00  0.00 0.00 0.00
            0.12  0.12 0.12 0.00
            0.28  0.45 0.00 0.00
            0.45  0.82 0.05 0.00
            0.62  1.00 0.35 0.00
            0.80  1.00 0.78 0.20
            1.00  1.00 1.00 1.00
        ];

        cmap = stopsToMap(stops);

    end

    function cmap = divergingMap()

        stops = [
            0.00  0.25 0.65 1.00
            0.25  0.08 0.28 0.55
            0.50  0.03 0.04 0.06
            0.75  0.60 0.15 0.10
            1.00  1.00 0.60 0.25
        ];

        cmap = stopsToMap(stops);

    end

    function cmap = stopsToMap(stops)

        xx = linspace(0,1,256);
        cmap = zeros(256,3);

        for ch = 1:3
            cmap(:,ch) = interp1(stops(:,1),stops(:,ch+1),xx,'linear');
        end

    end

end
