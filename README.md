# dual-slope-toolkit

Tools for the design, optimization, and real-time acquisition of Dual-Slope (DS) NIRS data.

## Overview

Dual-Slope (DS) is a self-calibrating measurement technique that minimizes the impact of optode-coupling errors and surface-layer heterogeneity. This toolkit provides the necessary code to design DS sensor arrays and acquire data in real-time.

## Contents

### Array Design & Optimization (`DSarrays/`)
A framework for creating and evaluating dual-slope optode arrangements:
- **Optimization:** Find optimal source-detector pairings to maximize sensitivity and depth penetration (`findGoodOpts.m`, `rmBadChans.m`).
- **Mapping:** Generate spatial sensitivity and absorption maps for specific array geometries (`arraySensMap.m`, `arrayAbsMap.m`).
- **Processing:** Parse raw array data and reconstruct 2D maps of optical properties (`parseArrayData.m`, `arrayRecon.m`).

### Real-Time Acquisition (`DRS/`)
Scripts for acquiring and monitoring Diffuse Reflectance Spectroscopy (DRS) data in real-time:
- **Acquisition:** `runDSspec.m` - The primary script for running spectrometer-based DS acquisitions.
- **Visualization:** `runDSspec_realTimePlot.m` - Provides live feedback of spectra and calculated DS metrics.
- **Hardware Integration:** `runDSspec_realTimePlot_LJaux.m` - Includes auxiliary input support via LabJack.

## Author
Developed by Giles Blaney, Ph.D.
