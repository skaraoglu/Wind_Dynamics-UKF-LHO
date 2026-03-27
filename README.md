# Wind Dynamics Experiment

[![R](https://img.shields.io/badge/R-≥4.2-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-≥3.9-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-Academic_Use-lightgrey)]()
[![Status](https://img.shields.io/badge/Status-Pilot_Complete-brightgreen)]()

## Overview

We apply a UKF-based linearized harmonic oscillator (LHO) estimation pipeline to daily wind speed observations from two meteorological stations in the Sürmeli Depression of eastern Türkiye, spanning nine years (2014-2022,
3287 daily records). The analysis framework evolved through a principled sequence of modeling decisions: an initial Hopf bifurcation approach exhibited systematic frequency collapse to the analysis band boundary, motivating the adoption of a linearized harmonic oscillator model whose second-order damped oscillation physics matches the driven dynamics of synoptic weather systems.

## Project Structure

```
wind_experiment/
├── R/
│   ├── constants.R            ← UKF numerical constants (domain-independent)
│   ├── ukf_engine.R           ← UKF core: propagate_model, UKF_dT, UKF_blend
│   ├── sl_models.R            ← SL models: make_sl_single, make_sl_paired_fixed_aw,
│   │                            hilbert_analytic, omega_from_phase
│   ├── models.R               ← Coupled oscillator models: LHO, gLk, abk
│   ├── optim.R                ← iterative_param_optim, optim_params
│   │                               These 5 files are exact with UKF-MDD.
│   ├── wind_config.R          ← NEW: wind-specific SL_BOUNDS, timing, bounds
│   ├── wind_preprocessing.R   ← NEW: wind data loading, cleaning, Hilbert
│   └── wind_logging.R         ← NEW: per-fit CSV logging, session log
├── data/
│   └── wind.csv               ← USER: place wind data here
├── logs/                      ← Generated: master_log.csv, session_log.txt, etc.
├── plots/                     ← Generated: figures
└── wind_experiment.ipynb      ← Main experiment notebook (R kernel)
```
**Only domain-specific configuration changed:**

| Parameter | BOLD | Wind  | Reason |
|-----------|-------------|-------------|--------|
| `dT` | 2.0 s (TR) | 1.0 day | Data sampling rate |
| `dt` | 0.2 s | 0.1 day | RK4 sub-step (dT/10) |
| `OM_MIN` | 0.126 rad/TR | 0.017 rad/day | Annual cycle |
| `OM_MAX` | 1.257 rad/TR | 2.094 rad/day | 3-day synoptic |
| `T_obs` | ~260 TRs | ~3287 days | Record length |

## Hypotheses

- H1 & Single-station natural frequency & LHO single & $\gamma$ \\
- H2 & Inter-station coupling & LHO two-stage & $k$ \\
- H3 & Seasonal regime transitions & LHO single & $\gamma$ per season \\
- H4 & Yearly stationarity & LHO single & $\gamma$ per year \\
- H5 & Damping sensitivity & LHO (varied $\zeta$) & $\gamma$ and $k$ vs $\zeta$ \\
- H6 & Model comparison (LHO vs pendulum) & Both & $\Delta\chi^2$ \\
- H7 & Cross-domain SL comparison & SL fixed-$\omega$ & $a$ \\
- Exp~A & Gamma ceiling fix ($\gamma_{\max}$: $5 \to 20$) & LHO both stages & $\gamma$, $k$ \\
- Exp~B & Cross-signal robustness & LHO two-stage & $k$ per pair \\
- Exp~C & Full zeta sensitivity (new $\gamma_{\max}$) & LHO both stages & $\gamma$, $k$ vs $\zeta$ \\

## Dependencies

R packages: `pracma`, `MASS`, `Matrix` (installed automatically by the notebook)

---

## Citation

If you use this pipeline or build on this work, please cite:

---

<div align="center">

LHO derived from a [Coupled Oscillator Model](https://doi.org/10.5074/t.2018.002) · [Stuart-Landau](https://en.wikipedia.org/wiki/Stuart%E2%80%93Landau_equation) normal form · [Unscented Kalman Filter](https://github.com/insilico/UKF) · [UKF-MDD](https://github.com/skaraoglu/UKF-MDD)

</div>

