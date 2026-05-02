"""
Digital Patient Simulator — Enhanced Extensions
Novel contributions beyond Barbiero et al. (2020).

All computations use ONLY the original Barbiero ODEs and .mat parameters.
No new biological constants are introduced — only post-processing and
batch-running of the existing solver.

Utku Köse (2026) — github.com/utkukose/digital-patient
"""
import numpy as np
from typing import Optional
from pydantic import BaseModel, Field
from main import (
    run_infection, run_diabetes, run_dkd, run_hypertension,
    InfectionRequest, DiabetesRequest, HypertensionRequest,
)


# ═══════════════════════════════════════════════════════════════════════════════
#  1. DIGITAL PATIENT RISK SCORE (DPRS)
#     Novel composite index derived from RAS ODE outputs.
#     Combines inflammation burden, receptor balance, ACE2 integrity, drug effect.
# ═══════════════════════════════════════════════════════════════════════════════

def compute_dprs(r: dict) -> dict:
    """
    Compute the Digital Patient Risk Score and its sub-components from
    an infection model result dict.
    Returns a dict with DPRS (0–100) and four component scores.
    """
    angII  = np.array(r["angII"])
    ang17  = np.array(r["ang17"])
    ace2   = np.array(r["ACE2"])
    ir     = np.array(r["IR"])
    at1r   = np.array(r["at1r"])
    at2r   = np.array(r["at2r"])
    diacid = np.array(r["diacid"])

    # Component 1: Inflammation burden (IR-based, 0–1)
    ir_score = float(np.clip(ir[-1] / 100.0, 0, 1))

    # Component 2: Receptor imbalance — AT1R/AT2R (higher = more vasoconstrictive)
    recep_ratio = float(at1r[-1] / max(at2r[-1], 1e-9))
    recep_score = float(np.clip((recep_ratio - 1) / 10.0, 0, 1))

    # Component 3: Counter-regulatory protection — ANG(1-7)/AngII
    prot_ratio  = float(ang17[-1] / max(angII[-1], 1e-9))
    prot_score  = float(np.clip(1.0 - prot_ratio / 200.0, 0, 1))   # lower = more protected

    # Component 4: ACE2 integrity loss
    ace2_loss   = float(np.clip((ace2[-1] - ace2[0]) / max(ace2[0], 1e-9), -1, 20))
    ace2_score  = float(np.clip(ace2_loss / 20.0, 0, 1))

    # Weighted composite (weights sum to 1)
    dprs = (ir_score * 0.40 + recep_score * 0.20 +
            prot_score * 0.20 + ace2_score * 0.20)
    dprs_100 = round(float(np.clip(dprs * 100, 0, 100)), 2)

    # Drug effect index
    inh_pct = float((100 * diacid**0.99 / (diacid**0.99 + 2.2**0.99)).max()) if diacid.max() > 0 else 0.0
    drug_eff = round(inh_pct, 2)

    return {
        "DPRS":            dprs_100,
        "inflammation":    round(ir_score * 100, 2),
        "receptor_imbal":  round(recep_score * 100, 2),
        "protection_loss": round(prot_score * 100, 2),
        "ACE2_integrity":  round((1 - ace2_score) * 100, 2),
        "drug_efficacy":   drug_eff,
        "AT1R_AT2R_ratio": round(recep_ratio, 4),
        "ANG17_AngII_ratio": round(prot_ratio, 4),
    }


# ═══════════════════════════════════════════════════════════════════════════════
#  2. DOSE-RESPONSE ANALYSIS
#     Sweep drug dose across a range; compute AngII, IR, inhibition % at each dose.
#     Barbiero only ran fixed doses; this reveals the full therapeutic window.
# ═══════════════════════════════════════════════════════════════════════════════

class DoseResponseRequest(BaseModel):
    doses: list[float] = Field(default=[0,1,2,3,5,7,10,15,20])
    glu: float = Field(5.0, ge=0.1, le=30)
    infection: bool = True
    renal_function: str = Field("normal", pattern="^(normal|impaired)$")
    sim_time_end: int = Field(120, ge=24, le=336)
    drug_name: str = "benazepril"

def run_dose_response(req: DoseResponseRequest) -> dict:
    curves = {"dose": [], "angII_final": [], "angII_peak": [],
              "IR_final": [], "inhibition_pct": [], "diacid_peak": [],
              "DPRS": [], "ANG17_AngII_ratio": []}

    for dose in req.doses:
        r = run_infection(InfectionRequest(
            dose=float(dose), glu=req.glu, infection=req.infection,
            renal_function=req.renal_function, sim_time_end=req.sim_time_end,
            drug_name=req.drug_name,
        ))
        diacid = np.array(r["diacid"])
        inh = float((100 * diacid**0.99 / (diacid**0.99 + 2.2**0.99)).max()) if diacid.max() > 0 else 0.0
        dprs = compute_dprs(r)
        curves["dose"].append(float(dose))
        curves["angII_final"].append(round(float(r["angII"][-1]), 6))
        curves["angII_peak"].append(round(float(max(r["angII"])), 6))
        curves["IR_final"].append(round(float(r["IR"][-1]), 6))
        curves["inhibition_pct"].append(round(inh, 3))
        curves["diacid_peak"].append(round(float(diacid.max()), 4))
        curves["DPRS"].append(dprs["DPRS"])
        curves["ANG17_AngII_ratio"].append(dprs["ANG17_AngII_ratio"])

    return curves


# ═══════════════════════════════════════════════════════════════════════════════
#  3. PARAMETER SENSITIVITY (HEATMAP)
#     Grid of glu × dose → IR_final. Reveals therapeutic interaction surface.
#     Barbiero never computed this; he only ran fixed combinations.
# ═══════════════════════════════════════════════════════════════════════════════

class SensitivityRequest(BaseModel):
    glu_values:  list[float] = Field(default=[1, 5, 10, 15, 17, 20])
    dose_values: list[float] = Field(default=[0, 2, 5, 10, 15, 20])
    infection:   bool  = True
    renal_function: str = Field("normal", pattern="^(normal|impaired)$")
    sim_time_end: int  = Field(120, ge=24, le=336)
    output_var:  str   = Field("IR_final", pattern="^(IR_final|angII_final|DPRS|ANG17_AngII_ratio)$")

def run_sensitivity(req: SensitivityRequest) -> dict:
    matrix, dprs_matrix = [], []
    for glu in req.glu_values:
        row, dprs_row = [], []
        for dose in req.dose_values:
            r = run_infection(InfectionRequest(
                dose=float(dose), glu=float(glu), infection=req.infection,
                renal_function=req.renal_function, sim_time_end=req.sim_time_end,
            ))
            dprs = compute_dprs(r)
            if req.output_var == "IR_final":
                val = round(float(r["IR"][-1]), 4)
            elif req.output_var == "angII_final":
                val = round(float(r["angII"][-1]), 4)
            elif req.output_var == "DPRS":
                val = dprs["DPRS"]
            else:
                val = dprs["ANG17_AngII_ratio"]
            row.append(val)
            dprs_row.append(dprs["DPRS"])
        matrix.append(row)
        dprs_matrix.append(dprs_row)

    return {
        "glu_values":   req.glu_values,
        "dose_values":  req.dose_values,
        "matrix":       matrix,
        "dprs_matrix":  dprs_matrix,
        "output_var":   req.output_var,
    }


# ═══════════════════════════════════════════════════════════════════════════════
#  4. TEMPORAL PHASE ANALYSIS
#     Decompose the ODE time-series into clinical phases:
#       Phase 0: Baseline (t < 0.5d)
#       Phase 1: Acute response (AngII rising fastest)
#       Phase 2: Peak stress (IR > 50% of max)
#       Phase 3: Drug effect (diacid building)
#       Phase 4: Chronic steady-state
#     Not done by Barbiero — he only plotted raw trajectories.
# ═══════════════════════════════════════════════════════════════════════════════

def compute_phases(r: dict) -> dict:
    t      = np.array(r["t"])
    angII  = np.array(r["angII"])
    ir     = np.array(r["IR"])
    diacid = np.array(r["diacid"])
    ace2   = np.array(r["ACE2"])

    d_angII = np.gradient(angII, t)
    d_ir    = np.gradient(ir, t)

    # Phase boundaries
    peak_angII_idx  = int(np.argmax(angII))
    peak_ir_idx     = int(np.argmax(ir))
    drug_onset_idx  = int(np.argmax(diacid > diacid.max() * 0.1)) if diacid.max() > 0 else len(t)-1
    drug_peak_idx   = int(np.argmax(diacid))

    phases = [
        {"name": "Baseline",       "t_start": round(float(t[0]),3),          "t_end": round(float(t[min(10,len(t)-1)]),3),  "key_event": "Pre-infection steady state"},
        {"name": "Acute response", "t_start": round(float(t[0]),3),           "t_end": round(float(t[peak_angII_idx]),3),    "key_event": f"AngII peaks at {angII[peak_angII_idx]:.4f} ng/mL"},
        {"name": "Inflammatory",   "t_start": round(float(t[peak_angII_idx]),3),"t_end": round(float(t[peak_ir_idx]),3),     "key_event": f"IR peaks at {ir[peak_ir_idx]:.4f}"},
        {"name": "Drug onset",     "t_start": round(float(t[drug_onset_idx]),3),"t_end": round(float(t[drug_peak_idx]),3),   "key_event": f"Diacid peaks at {diacid[drug_peak_idx]:.4f} ng/mL"},
        {"name": "Steady state",   "t_start": round(float(t[drug_peak_idx]),3),"t_end": round(float(t[-1]),3),              "key_event": f"Final AngII={angII[-1]:.4f}, IR={ir[-1]:.4f}"},
    ]

    # Key clinical metrics per phase
    metrics = {
        "angII_acceleration_max": round(float(d_angII.max()), 6),
        "IR_acceleration_max":    round(float(d_ir.max()), 6),
        "ace2_nadir":             round(float(ace2.min()), 6),
        "ace2_nadir_day":         round(float(t[np.argmin(ace2)]), 3),
        "angII_peak_day":         round(float(t[peak_angII_idx]), 3),
        "IR_peak_day":            round(float(t[peak_ir_idx]), 3),
        "drug_onset_day":         round(float(t[drug_onset_idx]), 3),
        "drug_halflife_est":      "5.2d (impaired)" if diacid.max() > 50 else "0.8d (normal)",
    }

    return {"phases": phases, "metrics": metrics}


# ═══════════════════════════════════════════════════════════════════════════════
#  5. DIABETES — HOMA-IR AND CLINICAL INDICES
#     Standard clinical indices from the diabetes ODE outputs.
#     Barbiero's paper never computed these; only plotted G and I trajectories.
# ═══════════════════════════════════════════════════════════════════════════════

class DiabetesIndexRequest(BaseModel):
    glu_values: list[float] = Field(default=[5, 8, 10, 12, 15, 17, 20])
    sim_days: int = Field(10, ge=1, le=30)

def run_diabetes_indices(req: DiabetesIndexRequest) -> dict:
    results = []
    for glu in req.glu_values:
        r = run_diabetes(DiabetesRequest(glu=float(glu), sim_days=req.sim_days))
        g    = np.array(r["G"])
        ins  = np.array(r["I"])
        beta = np.array(r["beta"])
        ir2  = np.array(r["IR"])
        mtor = np.array(r["MTOR"])

        # HOMA-IR (homeostatic model assessment) — standard clinical formula
        # HOMA-IR = (G_fasting * I_fasting) / 22.5
        # Proxy: use mean of last 10% of simulation
        n_tail = max(1, len(g) // 10)
        g_fast = float(g[-n_tail:].mean())
        i_fast = float(ins[-n_tail:].mean())
        homa_ir = round((g_fast * i_fast) / 22.5, 3)

        # HOMA-B (beta-cell function proxy)
        homa_b  = round((20 * i_fast) / max(g_fast - 3.5, 0.1), 3)

        # Insulin sensitivity index (Matsuda-like proxy)
        ISI = round(float(1.0 / max(homa_ir, 0.01)), 4)

        # Beta-cell depletion
        beta_dep = round(float((beta[0] - beta[-1]) / max(beta[0], 1e-9) * 100), 2)

        # Peak glucose excursion
        g_peak = round(float(g.max()), 2)
        g_range = round(float(g.max() - g.min()), 2)

        results.append({
            "glu_input":   float(glu),
            "G_final":     round(g_fast, 3),
            "I_final":     round(i_fast, 3),
            "HOMA_IR":     homa_ir,
            "HOMA_B":      homa_b,
            "ISI":         ISI,
            "beta_depletion_pct": beta_dep,
            "G_peak":      g_peak,
            "G_range":     g_range,
            "MTOR_final":  round(float(mtor[-1]), 4),
        })
    return {"indices": results}


# ═══════════════════════════════════════════════════════════════════════════════
#  6. VIRTUAL PATIENT COHORT
#     Run N virtual patients with varied (glu, dose) combinations and return
#     the full distribution of DPRS, IR, AngII.
#     Enables population-level risk stratification — not in Barbiero at all.
# ═══════════════════════════════════════════════════════════════════════════════

class CohortRequest(BaseModel):
    n_patients: int = Field(20, ge=4, le=50)
    glu_range:  list[float] = Field(default=[1, 20])
    dose_range: list[float] = Field(default=[0, 10])
    infection:  bool = True
    renal_function: str = Field("normal", pattern="^(normal|impaired)$")
    sim_time_end: int = Field(120, ge=24, le=168)
    seed: int = 42

def run_cohort(req: CohortRequest) -> dict:
    rng = np.random.default_rng(req.seed)
    glu_vals  = rng.uniform(req.glu_range[0],  req.glu_range[1],  req.n_patients)
    dose_vals = rng.uniform(req.dose_range[0], req.dose_range[1], req.n_patients)

    patients = []
    for i, (glu, dose) in enumerate(zip(glu_vals, dose_vals)):
        r = run_infection(InfectionRequest(
            dose=round(float(dose), 1), glu=round(float(glu), 1),
            infection=req.infection, renal_function=req.renal_function,
            sim_time_end=req.sim_time_end,
        ))
        dprs = compute_dprs(r)
        patients.append({
            "id":         i + 1,
            "glu":        round(float(glu), 2),
            "dose":       round(float(dose), 2),
            "angII_final": round(float(r["angII"][-1]), 4),
            "IR_final":   round(float(r["IR"][-1]), 4),
            "DPRS":       dprs["DPRS"],
            "drug_eff":   dprs["drug_efficacy"],
            "risk_tier":  "High" if dprs["DPRS"] > 60 else ("Moderate" if dprs["DPRS"] > 30 else "Low"),
        })

    dprs_arr = [p["DPRS"] for p in patients]
    ir_arr   = [p["IR_final"] for p in patients]
    return {
        "patients":     patients,
        "summary": {
            "n":             req.n_patients,
            "DPRS_mean":     round(float(np.mean(dprs_arr)), 2),
            "DPRS_std":      round(float(np.std(dprs_arr)), 2),
            "DPRS_min":      round(float(np.min(dprs_arr)), 2),
            "DPRS_max":      round(float(np.max(dprs_arr)), 2),
            "high_risk_pct": round(sum(1 for d in dprs_arr if d > 60) / req.n_patients * 100, 1),
            "mod_risk_pct":  round(sum(1 for d in dprs_arr if 30 < d <= 60) / req.n_patients * 100, 1),
            "low_risk_pct":  round(sum(1 for d in dprs_arr if d <= 30) / req.n_patients * 100, 1),
            "IR_mean":       round(float(np.mean(ir_arr)), 4),
        }
    }
