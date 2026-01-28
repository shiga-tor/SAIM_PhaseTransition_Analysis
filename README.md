# SAIM: Systemic Attractor Instability Metric (v1.0)

This repository contains the official Python implementation of the **SAIM Analysis Pipeline v1.0** and the **Generative Model Simulation**, as described in the manuscript:

> **"Informational Perturbation Resolves Precision Collapse and Restores Adaptive Neural Dynamics"**
> *Takafumi Shiga* (Submitted to PNAS, 2026)

---

## ⚠️ Important Note for Reviewers

**The analysis logic and preprocessing pipeline presented in this repository were fixed prior to the commencement of the study (pre-registered) and have NOT been modified based on post-hoc observations.**

While we acknowledge the existence of advanced artifact removal techniques (e.g., ASR, CCA, or spectral parameterization), we prioritized **transparency, reproducibility, and minimal signal distortion** for this mechanistic investigation. The validity of this "fixed pipeline" approach is supported by rigorous experimental constraints and control analyses (detailed in the manuscript's *Supplementary Information, Text S2 & S3*).

---

## Overview

This repository provides two key components of the study:

1. **Empirical Analysis Pipeline (`SAIM_pipeline_v1.0.py`)**:
Processes raw EEG and fNIRS signals (Muse S Gen 2 Athena) to compute the **Free Energy Proxy (F)**, **Hemodynamic Coupling (HEMO)**, and associated metrics. It implements the automated detection of **PPEN** (Proprioceptive Prediction Error Neglect), including the "Paradoxical Rigidity" subtype (Rule C).
* *Corresponds to Supplementary Text S9.*


2. **Generative Model Simulation (`simulation.py`)**:
Reproduces the computational dynamics of the "Frozen Attractor" and the phase transition mechanism induced by Specific Informational Perturbation (SIP).
* *Corresponds to Supplementary Text S8.*



## Methodological Rationale

To ensure physiological interpretability given the portable hardware constraints, this pipeline adopts a **Proxy-Based Approach**:

* **Gamma Band (PE Proxy):** Gamma power () is treated as a macro-scale proxy for prediction error signaling. Myogenic contamination is mitigated via strict experimental protocols rather than aggressive filtering that might distort phase dynamics.
* **HEMO (Stability Proxy):** Optical signals are processed to quantify **Neurovascular Coupling Stability** (signal-to-noise consistency) based on raw photodiode variance, avoiding assumptions about optical pathlength required for absolute hemoglobin quantification.

For a detailed technical validation, please refer to our **Technical Report**:
📄 **[Muse_Technical_Validation.pdf](https://www.google.com/search?q=./docs/Muse_Technical_Validation.pdf)**

## Requirements

* Python 3.8+
* Dependencies:
* `pandas>=1.3.0`
* `numpy>=1.21.0`
* `matplotlib>=3.4.0`
* `seaborn>=0.11.0`
* `scipy>=1.7.0`
* `scikit-learn>=0.24.0`



## Directory Structure

```
.
├── data/                   # Place raw CSV files here
├── docs/                   # Documentation & Validation Reports
│   └── Muse_Technical_Validation.pdf  <-- Deep Research Report
├── output/                 # Generated CSVs and PNGs
├── SAIM_pipeline_v1.0.py   # Main Analysis Script
├── simulation.py           # Generative Model Script
├── LICENSE
└── README.md

```

## Usage

### 1. Empirical Data Analysis (S9)

To analyze raw physiological data:

1. Place your raw CSV data files (Muse S format) in the `data/` directory.
2. Run the pipeline:
```bash
python SAIM_pipeline_v1.0.py

```


3. **Outputs**:
* `*_TimeSeries.csv`: Frame-by-frame metric calculations.
* `*_Overall_Stats.csv`: Session-level statistics (Pre vs Post).
* `*_Continuous_Dynamics.png`: Visualization of the phase transition with interpolation.
* `*_OmniPanel.png`: The 20-panel spectral profiling plot (Figure 2).



### 2. Computational Simulation (S8)

To reproduce the theoretical bistable dynamics and SIP mechanism (Figure 1):

1. Run the simulation script:
```bash
python simulation.py

```


2. **Outputs**:
* Displays/Saves the "Phase Transition via SIP" plot, illustrating the transition from the Frozen Attractor to the Adaptive State.



## License

This project is licensed under the MIT License - see the [LICENSE](https://www.google.com/search?q=LICENSE) file for details.

---

Copyright (c) 2026 TIC-DO Institute.
