# M87* EHT — GPU-Accelerated VLBI Image Reconstruction & Telemetry Explorer

A pure-MATLAB, end-to-end Very Long Baseline Interferometry (VLBI) imaging pipeline and interactive scientific telemetry suite designed for Event Horizon Telescope (EHT) observations of **M87\***.

The project consists of two core modules:
1. **`M87_EHT_Headless_GPU_Pipeline.m` (v2.1):** A headless, hybrid CPU/GPU synthesis imaging pipeline that natively parses raw binary `.uvfits` radio interferometry files, performs convolutional $uv$-gridding, executes windowed Högbom CLEAN (Stokes $I$ and full-Stokes linear polarimetry $Q, U$), and extracts visibility-domain and closure-phase telemetry.
2. **`M87_EHT_Master_Explorer_Dashboard_v3.m` (v3):** A multi-tab `uifigure` scientific visualization dashboard for inspecting reconstructed images, black hole shadow ring geometry, synthesized beams, closure phases, Fourier-plane sampling, and linear polarization (EVPA) maps.

---

## Repository Structure

```text
├── M87_EHT_Headless_GPU_Pipeline.m          # Headless UVFITS reader, hybrid GPU/CPU CLEAN, & polarimetry pipeline
├── M87_EHT_Master_Explorer_Dashboard_v3.m   # Interactive 4-tab GUI dashboard for post-reconstruction analysis
└── README.md                                # Project documentation
```

---

## Key Capabilities

### 1. Headless Reconstruction Pipeline (`M87_EHT_Headless_GPU_Pipeline`)
* **Zero-Dependency Native UVFITS Parser:** Reads IEEE big-endian FITS headers, random-group $(u, v, w)$ coordinates, visibility weights, observation Julian dates (`DATE` / `_DATE`), baseline antenna pairs, and AIPS `AN` station tables without requiring external CFITSIO or Python bindings.
* **Automatic Station & Frequency Resolution:** Maps standard EHT two-letter station codes (`AA` $\rightarrow$ ALMA, `AP` $\rightarrow$ APEX, `AZ` $\rightarrow$ SMT, `JC` $\rightarrow$ JCMT, `LM` $\rightarrow$ LMT, `PV` $\rightarrow$ IRAM30, `SM` $\rightarrow$ SMA, `SP` $\rightarrow$ SPT, `GL` $\rightarrow$ GLT, `KT`/`KP` $\rightarrow$ KP, `NN` $\rightarrow$ NOEMA, `PB` $\rightarrow$ PdB) and resolves observing frequencies from FITS headers or EHT band filename tags (`_lo_`, `_hi_`, `_b1_`–`_b4_`).
* **Hybrid GPU/CPU Execution:** Offloads large $2048 \times 2048$ 2D inverse FFTs and exact full-field FFT residual convolutions to CUDA GPUs (`gpuArray`) while running cache-friendly windowed Högbom CLEAN on the CPU—avoiding PSF patch truncation errors while keeping VRAM usage low.
* **Full-Stokes Polarimetric Imaging:** Automatically detects `RR/LL/RL/LR` (or pre-formed `I/Q/U`) polarization planes, forms Stokes $I = (RR + LL)/2$, runs complex Högbom CLEAN on $P = Q + iU \approx RL$, applies Ricean noise debiasing ($P_{\text{deb}} = \sqrt{\max(\vert{}P\vert{}^2 - \sigma_{QU}^2, 0)}$), and computes fractional linear polarization $m$ and Electric Vector Position Angle ($\text{EVPA}$).
* **Visibility & Closure-Phase Telemetry:** Computes weighted radial visibility amplitude profiles against the CLEAN model, downsampled visibility scatter sets, and closed-triangle closure phases across all valid station triads.
* **Thermal & System Load Guard (v2.1):** Monitors Linux ACPI thermal zones (`/sys/class/thermal/thermal_zone*/temp`), caps BLAS thread counts, and inserts configurable cool-down pauses during sustained runs on laptops.

### 2. Interactive Telemetry Dashboard (`M87_EHT_Master_Explorer_Dashboard_v3`)
* **10 Live KPI Cards:** Peak brightness, off-source noise RMS, dynamic range, total CLEAN flux, synthesized beam dimensions & position angle, plotted $uv$ sample count, maximum baseline $\vert{}B\vert{}_{\max}$ & diffraction limit $\lambda/B$, CLEAN iteration count, automated shadow ring diameter & FWHM, and net linear polarization.
* **4 Dedicated Analysis Views:**
  * **IMAGING:** Geographic EHT array baseline map (`geoaxes`), Dirty Image, CLEAN Model, CLEAN Residual, Restored M87* Image (with $3\sigma \cdot 2^n$ contours, scale bar, and restoring beam ellipse), Synthesized Beam / PSF, $uv$-plane coverage with angular-resolution rings, and CLEAN convergence curve with live system console logs.
  * **VISIBILITY:** Visibility amplitude vs. baseline length $\vert{}B\vert{}$ (raw samples, weighted radial bins $\pm 1\sigma$, and CLEAN model overlay), closure phase vs. UT time, baseline amplitude vs. UT time, and circular-statistics closure-phase histogram.
  * **ANALYSIS:** Automated ring-centering and radial brightness profile (measuring shadow ring diameter $D$ and beam-convolved FWHM), azimuthal brightness modulation (asymmetry $A$ and $m=1$ mode $\beta_1$), RA/Dec peak cross-sections, residual PDF vs. Gaussian noise model (skewness & excess kurtosis), PSF sidelobe envelope in dB, and 2D Fourier transform $\vert{}\tilde{I}(u,v)\vert{}$ overlaid with $uv$ sampling.
  * **POLARIMETRY:** Stokes $I$ overlaid with EVPA polarization vectors, debiased polarized intensity $\vert{}P\vert{}$, fractional polarization $m$, Stokes $Q$, Stokes $U$, and EVPA vs. azimuthal position angle.
* **Interactive Controls & Cursor Readout:** Real-time hover readout displaying RA/Dec offset ($\mu\text{as}$), radial distance, and pixel intensity in SI units, alongside live FOV modes (`AUTO` / `FULL` / `MANUAL`), colormaps (`EHT`, `HOT`, `TURBO`, `PARULA`, `GRAY`), intensity stretches (`LINEAR`, `SQRT`, `LOG`, `ASINH`), and high-DPI PNG export.

---

## Requirements

* **MATLAB:** R2021a or newer recommended (uses `uifigure`, `uigridlayout`, `geoaxes`, and `exportapp`).
* **Core MATLAB Only:** All FITS parsing, $uv$-gridding, Högbom CLEAN, beam fitting, and polarimetric math use built-in MATLAB functions—no third-party astronomy toolboxes required.
* **Optional Toolboxes:**
  * *Parallel Computing Toolbox* — enables GPU acceleration (`computeMode = 'HYBRID'`). Automatically falls back to multi-threaded CPU execution if no compatible GPU is detected.
  * *Mapping Toolbox* — enhances the geographic basemap (`darkwater` / `streets`) in the **Array Geometry** panel.

---

## Dataset Setup (EHT Public `.uvfits` Releases)

This pipeline is structured around official Event Horizon Telescope (EHT) public data releases (e.g., 2017 M87* Polarimetric Release `2023-D01-01` or 2018 M87* Release `2024-D01-01`):

1. Download the EHT `.uvfits` release archive and extract it locally.
2. Choose **one** calibration pipeline variant per run so redundant datasets are not merged together:
   * **HOPS (D-term calibrated — recommended for polarimetry):** `hops_data` with filter `'*hops_zbl-dtcal+selfcal.uvfits'`
   * **CASA (D-term calibrated):** `casa_data` with filter `'*casa_EVPA_rotation+dcal*.uvfits'`
   * **HOPS (No D-terms — Stokes $I$ only):** `'*hops_ALMArot.uvfits'`

> **Polarimetry Calibration Note:** The pipeline defines complex linear polarization as $P = Q + iU \approx RL$ (with $P(-u) = LR^*(u)$). Always use D-term-calibrated files when running with `CFG.polEnable = true`. Verify the EVPA sign convention against published EHT calibration notes for your specific dataset release before quoting absolute EVPA angles.

---

## Quick Start

### Step 1: Run the Headless Reconstruction Pipeline
You can run the pipeline using its default configuration or pass a custom `CFG` struct in MATLAB:

```matlab
CFG = struct();
CFG.dataRoot   = '/path/to/EHT/2023-D01-01-main';
CFG.uvfitsDir  = fullfile(CFG.dataRoot, 'hops_data');
CFG.fileFilter = '*hops_zbl-dtcal+selfcal.uvfits';
CFG.outputMAT  = '/path/to/output/M87_EHT_POL_RECONSTRUCTION.mat';

% Optional: Force CPU-only mode or limit file count for a fast test run
CFG.computeMode = 'HYBRID';   % 'HYBRID' (GPU FFT + CPU CLEAN) or 'CPU'
CFG.maxFiles    = Inf;        % Set to 2 for a quick single-band test

reconstruction = M87_EHT_Headless_GPU_Pipeline(CFG);
```

### Step 2: Launch the Interactive Explorer Dashboard
Update `CFG.outputMAT` at the top of `M87_EHT_Master_Explorer_Dashboard_v3.m` (or click **LOAD .MAT DATA** inside the GUI) and run:

```matlab
M87_EHT_Master_Explorer_Dashboard_v3
```

---

## Configuration Reference (`CFG`)

| Parameter | Default | Description |
| :--- | :--- | :--- |
| `dataRoot` | `.../2023-D01-01-main` | Root directory of the EHT data release. |
| `uvfitsDir` | `fullfile(dataRoot, 'hops_data')` | Directory searched recursively for `.uvfits` files. |
| `fileFilter` | `'*hops_zbl-dtcal+selfcal.uvfits'` | Glob pattern selecting a single calibration variant. |
| `excludeRegex` | `''` | Optional case-sensitive regex to exclude matching filenames. |
| `outputMAT` | `.../M87_EHT_POL_RECONSTRUCTION.mat` | Output `-v7.3` `.mat` file path containing the `reconstruction` struct. |
| `Nuv` | `2048` | Square $uv$-grid and image dimension ($N \times N$ pixels). |
| `du_Mlambda` | `20` | $uv$-plane cell spacing in $\text{M}\lambda$ (sets total field of view $\text{FOV} = 1/\Delta u \approx 10.31\text{ mas}$). |
| `kernelHalfWidth` | `3` | Half-width (in grid cells) of the Gaussian convolutional gridding kernel. |
| `kernelSigma` | `1.2` | Standard deviation (in grid cells) of the gridding kernel. |
| `cleanGain` | `0.08` | Högbom CLEAN loop gain $\gamma$. |
| `cleanMaxIter` | `5000` | Maximum CLEAN iterations for Stokes $I$. |
| `cleanSigma` | `3.0` | Stopping threshold multiplier ($\text{threshold} = \text{cleanSigma} \times \sigma_{\text{dirty}}$). |
| `cleanRadius_mas` | `0.20` | Circular CLEAN search window radius in milliarcseconds ($\text{mas}$). |
| `positiveClean` | `true` | Restricts Stokes $I$ CLEAN components to positive flux. |
| `maxModelFluxJy` | `10` | Safety cap on cumulative cleaned flux ($\text{Jy}$). |
| `polEnable` | `true` | Enables full-Stokes $Q, U$ complex CLEAN when `RL/LR` planes exist. |
| `polCleanMaxIter` | `3000` | Maximum complex CLEAN iterations for $P = Q + iU$. |
| `polCleanSigma` | `4.0` | Polarimetric CLEAN stopping threshold multiplier. |
| `evpaOffsetDeg` | `0` | Optional constant EVPA rotation offset in degrees (e.g., `90`). |
| `computeMode` | `'HYBRID'` | `'HYBRID'` (GPU FFTs + CPU CLEAN) or `'CPU'` (no GPU usage). |
| `thermalGuard` | `true` | Pauses execution if Linux ACPI thermal sensors exceed `maxTempC` ($85^\circ\text{C}$). |

---

## Platform & Usage Notes

* **Windows & macOS Paths:** The default paths in both scripts use Linux directories (`/home/hellblazer/...`). On Windows or macOS, pass a custom `CFG` struct with your local paths to `M87_EHT_Headless_GPU_Pipeline(CFG)` and update lines 4–7 of `M87_EHT_Master_Explorer_Dashboard_v3.m` (or load any `.mat` file interactively via the **LOAD .MAT DATA** button).
* **Large `.mat` Output Files:** Because `Nuv = 2048` stores multiple full-resolution double-precision image planes (`dirtyImage`, `cleanModel`, `cleanResidual`, `restoredImage`, `psf`, and `pol.*`), generated `.mat` files are saved in HDF5 `-v7.3` format and typically exceed GitHub's 100 MB standard file limit. Keep `.mat` outputs excluded via `.gitignore` (`*.mat`) or distribute sample reconstructions via **GitHub Releases**.
