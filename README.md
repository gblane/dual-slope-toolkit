# dual-slope-toolkit

Tools for the design, optimization, and spatial mapping of Dual-Slope (DS) Near-Infrared Spectroscopy (NIRS) data.

## Overview

Dual-Slope (DS) is a self-calibrating measurement technique that minimizes the impact of optode-coupling errors and surface-layer heterogeneity. This toolkit provides the core algorithms for discovering valid DS pairings, optimizing array layouts, and reconstructing 2D spatial maps of optical properties and hemodynamics.

## Repository Structure

### Source Code (`src/`)
- **`physics/`**: Core DS and conversion tools.
  - `DSdmua.m`: Standalone DS absorption change calculator.
  - `mua2OandD.m`: Hemodynamic converter (HbO, HbR).
  - `DPF_DSF_calc.m`: Pathlength and slope factor calculator.
  - `GammaDelta.m`: Array resolution and localization metrics.
- **`geometry/`**: Array discovery and layout optimization.
  - `DSdisc.m`: Automated pairing discovery.
  - `findDS_SDinds.m` / `findSS_SDinds.m`: Index mapping helpers.
  - `makeArraySPs.m`: Subplot layout optimizer.
- **`imaging/`**: Spatial mapping and reconstruction.
  - `arrayAbsMap.m`: 2D absolute property mapping.
  - `arrayData2dmua.m`: Absorption and hemodynamic change mapping.
  - `arrayRecon.m`: Regularized image reconstruction.
  - `arrayNoiseMap.m`: Spatial noise mapping.
  - `plotVectorizedMap.m`: 2D visualization utility.
- **`io/`**: Data parsing and signal cleaning.
  - `parseArrayData.m`: Raw ISS data parser.
  - `rmBadChans.m`: Noise-based channel filtering.
  - `calAndAddDatatypes.m`: Calibration and data-type propagation.

### Examples (`examples/`)
- **`standalone/`**: Direct hemodynamics calculation (Blaney et al., *J. Biophotonics* 2019).
- **`array_design/`**: Pairing discovery and array evaluation (Blaney et al., *Rev. Sci. Instrum.* 2020).
- **`neurow_imaging/`**: Step-by-step cerebral imaging workflow (Blaney et al., *Neurophotonics* 2023).

### Shared Data (`data/`)
- Consolidated `.mat` files for examples and model parameters.

## Author
Developed by Giles Blaney, Ph.D.

---
*Documentation written by Gemini CLI.*
