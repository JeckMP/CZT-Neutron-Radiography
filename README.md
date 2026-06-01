# CZT Data Processing

MATLAB GUI for processing CZT detector calibration and neutron spectroscopy data.

Developed as part of the Ph.D. research of Edcer Laguda at McMaster University.

## Overview

This software provides a graphical user interface (GUI) for processing data acquired from a pixelated Cadmium Zinc Telluride (CZT) detector system. The workflow includes background correction, pixel gain correction, energy calibration, energy resolution analysis, and neutron spectroscopy processing.

## Features

- Background subtraction
- Dual-background averaging
- Pixel gain correction
- Energy calibration
- Energy spectrum generation
- FWHM and energy resolution analysis
- Calibration curve generation
- Neutron exposure processing
- Data export and visualization

## Processing Workflow

1. Load background acquisition(s)
2. Load calibration source data
3. Perform background subtraction
4. Apply pixel gain correction
5. Generate energy calibration
6. Calculate FWHM and energy resolution
7. Process neutron exposure data
8. Export processed results

## Supported Calibration Sources

- Co-57
- Am-241
- Eu-152

## Requirements

- MATLAB R2020b or newer
- Signal Processing Toolbox (recommended)

## Input Data

The GUI expects MATLAB `.mat` files containing:

```matlab
LMmap_noEcorrect
LMmap_Ecorrect
```

where:

- `LMmap_noEcorrect` = detector counts before energy correction
- `LMmap_Ecorrect` = detector counts after energy correction

## Physics

The software implements the following processing methods:

### Background Subtraction

Ccorr = Craw − Cbackground

### Energy Calibration

E = m × Channel + b

### Energy Resolution

Resolution (%) = 100 × FWHM / Epeak

## Citation

If this software is used in academic work, please cite:

Laguda, E.

*Development of a CZT-Based Neutron Imaging Detector for Neutron Radiography and Tomography.*

Ph.D. Thesis, McMaster University, 2026.

### Research Supervision

This work was conducted under the supervision of:

- Dr. Soo Hyun Byun, McMaster University
- Dr. Troy Farncombe, McMaster University

The author gratefully acknowledges the guidance, mentorship, and support provided by Dr. Byun and Dr. Farncombe throughout the development of this research and the associated simulation framework.
