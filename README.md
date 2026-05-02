# dual-slope-toolkit

Tools for the design, optimization, and spatial mapping of Dual-Slope (DS) Near-Infrared Spectroscopy (NIRS) data.

## Overview

Dual-Slope (DS) is a self-calibrating measurement technique that minimizes the impact of optode-coupling errors and surface-layer heterogeneity. This toolkit provides the core algorithms for discovering valid DS pairings, optimizing array layouts, and reconstructing 2D spatial maps of optical properties and hemodynamics.

## Contents

### Array Discovery & Design (`DSarrays/`)
- **`DSdisc.m`**: Core function to discover valid Single-Distance (SD), Single-Slope (SS), and Dual-Slope (DS) pairings from physical coordinates.
- **`findGoodOpts.m`**: Identifies usable optode pairs based on signal quality and noise criteria.
- **`findDS_SDinds.m` / `findSS_SDinds.m`**: Index mapping helpers for data propagation.
- **`GammaDelta.m`**: Calculates resolution (Gamma) and localization error (Delta) metrics for array evaluation.

### Spatial Mapping & Reconstruction (`DSarrays/`)
- **`arrayAbsMap.m`**: Generates continuous 2D maps of absolute optical properties using Gaussian smoothing.
- **`arrayData2dmua.m`**: Calculates absorption changes ($\Delta\mu_a$) and hemodynamics (HbO, HbR) from processed data.
- **`arrayRecon.m`**: Performs regularized image reconstruction for NIRS array data.
- **`arrayNoiseMap.m`**: Interpolates discrete noise metrics onto a continuous 2D spatial map.

### Core DS Calculations (`abs_multiDist/`)
- **`DPF_DSF_calc.m`**: Primary calculator for Differential Pathlength Factors (DPF) and Differential Slope Factors (DSF), including two-layer support.

### Standalone DS Tools (Root)
- **`DSdmua.m`**: Standalone function for calculating Dual-Slope absorption changes from 4-channel intensity or phase data.
- **`mua2OandD.m`**: Converts absorption changes into oxygenated (HbO) and deoxygenated (HbR) hemoglobin concentrations.

### Examples (`examples/`)
- **`example_hemodynamics.m`**: Demonstrates calculating hemodynamics from raw data using standalone DS functions.
- **`example_array_discovery.m`**: Shows how to use `DSdisc.m` to identify valid optode pairings for a custom array layout.
- **`example_array_evaluation.m`**: Demonstrates array performance evaluation using the `GammaDelta.m` metric calculator.

### Visualization Utilities (`DSarrays/`)
- **`plotVectorizedMap.m`**: Flexible 2D visualization for discrete spatial data, including complex value support.
- **`makeArraySPs.m`**: Automatically determines optimal subplot layouts based on physical optode arrangements.

## Author
Developed by Giles Blaney, Ph.D.

---
*Documentation written by Gemini CLI.*
