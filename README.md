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
- **`standalone/`**: Direct hemodynamics calculation.
- **`array_design/`**: Pairing discovery and array evaluation.
- **`neurow_imaging/`**: Step-by-step cerebral imaging workflow.

### Shared Data (`data/`)
- Consolidated `.mat` files for examples and model parameters.

## Citations

If you use this toolkit in your research, please cite the following publications:

1.  **Dual-Slope Foundations:** Blaney, G., Sassaroli, A., Pham, T., Fernandez, C., & Fantini, S. (2019). Phase dual-slopes in frequency-domain near-infrared spectroscopy for enhanced sensitivity to brain tissue: First applications to human subjects. *Journal of Biophotonics*, 12(11), e201960018. [https://doi.org/10.1002/jbio.201960018](https://doi.org/10.1002/jbio.201960018)
2.  **Array Design:** Blaney, G., Sassaroli, A., & Fantini, S. (2020). Design of a source-detector array for dual-slope diffuse optical imaging. *Review of Scientific Instruments*, 91(11), 114102. [https://doi.org/10.1063/5.0015512](https://doi.org/10.1063/5.0015512)
3.  **Cerebral Imaging:** Blaney, G., Fernandez, C., Sassaroli, A., & Fantini, S. (2023). Dual-slope imaging of cerebral hemodynamics with frequency-domain near-infrared spectroscopy. *Neurophotonics*, 10(1), 013508. [http://doi.org/10.1117/1.NPh.10.1.013508](http://doi.org/10.1117/1.NPh.10.1.013508)

## Author
Developed by Giles Blaney, Ph.D.

---
*This repository is a reorganized and documented version of a personal codebase, performed by Gemini CLI.*
