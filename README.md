# SAIM: Systemic Attractor Instability Metric (v1.0)

This repository contains the official Python implementation of the **SAIM Analysis Pipeline v1.0**, as described in the paper:

> **"Informational Perturbation Resolves Precision Collapse and Restores Adaptive Neural Dynamics"**
> *Takafumi Shiga* (2026)

## Overview
This software processes raw EEG and fNIRS signals (specifically from the Muse S Gen 2 Athena device) to compute the **Free Energy Proxy (F)** and associated metrics based on the Free Energy Principle.

It features:
- **Dual-wavelength fNIRS analysis** (730nm/850nm) with ambient light correction.
- **Precision-Weighted Prediction Error (PE)** calculation.
- **Systemic Phase Transition** detection.

## Requirements
- Python 3.8+
- See `requirements.txt` for dependencies.

## Usage
1. Place your raw CSV data files in the same directory or a data folder.
2. Run the pipeline:
   ```bash
   python SAIM_pipeline_v1.0.py