# CHANGES — Digital Patient Simulator (Enhanced Edition)

This file documents all modifications and additions made to the original
Barbiero & Liò (2020) Computational Patient codebase.

## Changes to Original Source (patient_pkg/)

Only one change was made to the original Barbiero source code:

- **`patient_pkg/infection/_equations.py`** — corrected one import path:
  - Before: `from patient.pd._equations import _catalized_AngI`
  - After: `from patient_pkg.pd._equations import _catalized_AngI`
  - Reason: the package folder was renamed from `patient` to `patient_pkg`
    to avoid collision with Python's standard namespace. No mathematical
    content was changed.

## Changes to Solver Settings (main.py)

- **DKD model solver**: changed from `method="RK23"` to `method="LSODA"`
  - Reason: the RAS ODE system is stiff; RK23 (explicit, non-stiff) required
    ~16 seconds per solve on this system. LSODA (adaptive stiff/non-stiff)
    solves the same equations in ~0.08 seconds with identical results.
  - No equations, parameters, or initial conditions were changed.

## Additions (new files)

### backend/main.py
- Wrapped all five Barbiero ODE models as FastAPI POST endpoints
- Added five enhancement endpoints (`/api/enhance/*`)
- Added `.mat` parameter file loader with caching

### backend/extensions.py (new file)
Novel analytical functions, all computed using only the original Barbiero ODE outputs:
- `compute_dprs()` — Digital Patient Risk Score composite index
- `compute_phases()` — Clinical phase timeline decomposition
- `run_dose_response()` — Systematic drug dose sweep
- `run_sensitivity()` — Glucose × Dose parameter grid
- `run_diabetes_indices()` — HOMA-IR, HOMA-B, ISI derivation
- `run_cohort()` — Virtual patient population simulation

### frontend/index.html (new file)
Complete single-file browser UI including:
- Interactive body visualizer with animated organ state
- Gender-selectable patient silhouettes
- ECG strip animation
- All nine result tabs
- Chart.js time-series charts
- Scenario comparison with differential analysis
- Clinical report generation
- JSON export

## Authors

Original work: Pietro Barbiero, Ramon Viñas Torné, Pietro Liò (2020)
Enhancements: Prof. Dr. Utku Köse (2026)
