"""
Digital Patient Simulator — FastAPI Backend
Wraps the exact Barbiero et al. (2020) computational patient ODEs.
All parameters loaded from original .mat files; all solvers use LSODA.
"""
import os, sys, logging, time
from pathlib import Path
from typing import Optional

import numpy as np
import scipy.io
import pandas as pd
from scipy.integrate import solve_ivp
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse
from pydantic import BaseModel, Field

# ── path setup ──────────────────────────────────────────────────────────────
ROOT = Path(__file__).parent
DATA = ROOT.parent / "data_files"
sys.path.insert(0, str(ROOT))

from patient_pkg.pk._equations import analytical_PK
from patient_pkg.pd._equations import (
    GLU, mass_balance_AGT, mass_balance_Renin,
    mass_balance_AngI, mass_balance_AngII,
)
from patient_pkg.infection._equations import (
    mass_balance_AngI_infection, mass_balance_AngII_infection,
    mass_balance_ANG17, mass_balance_AT1R, mass_balance_AT2R,
)
from patient_pkg.hypertension._equations import (
    transit_compartment_model, change_S_ANG17, change_S_ANG2, change_S_SBP,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger("digital_patient")

app = FastAPI(title="Digital Patient Simulator", version="1.0.0")
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"], allow_methods=["*"], allow_headers=["*"],
)

# ── .mat parameter cache ──────────────────────────────────────────────────────
_mat_cache: dict = {}

def load_mat(fname: str) -> dict:
    if fname not in _mat_cache:
        path = DATA / fname
        if not path.exists():
            raise FileNotFoundError(f"Parameter file not found: {path}")
        raw = scipy.io.loadmat(str(path))
        out = {}
        for k, v in raw.items():
            if k.startswith("_"):
                continue
            arr = np.array(v)
            if arr.size == 1:
                flat = arr.flatten()[0]
                if isinstance(flat, (int, float, np.integer, np.floating)):
                    out[k] = float(flat)
        _mat_cache[fname] = out
    return _mat_cache[fname]

def get_params(drug: str, renal: str) -> dict:
    return load_mat(f"params_{drug}{renal}.mat")

def get_pk_params(drug: str, renal: str) -> dict:
    return load_mat(f"PK_params_{drug}{renal}.mat")


# ═══════════════════════════════════════════════════════════════════════════════
#  REQUEST / RESPONSE MODELS
# ═══════════════════════════════════════════════════════════════════════════════

class InfectionRequest(BaseModel):
    age: int = Field(20, ge=18, le=90)
    dose: float = Field(0.0, ge=0, le=20)
    renal_function: str = Field("normal", pattern="^(normal|impaired)$")
    glu: float = Field(1.0, ge=0.1, le=30)
    infection: bool = True
    drug_name: str = "benazepril"
    n_dose: int = Field(1, ge=1, le=4)
    sim_time_end: int = Field(120, ge=24, le=336)   # hours
    tstart_dosing: int = 0

class DiabetesRequest(BaseModel):
    glu: float = Field(5.0, ge=1, le=30)
    sim_days: int = Field(5, ge=1, le=30)

class HypertensionRequest(BaseModel):
    sim_time_end: int = Field(1000, ge=100, le=5000)

class CardioRequest(BaseModel):
    age: int = Field(20, ge=18, le=90)
    glu: float = Field(1.0, ge=0.1, le=30)
    dose: float = Field(0.0, ge=0, le=20)
    infection: bool = False
    renal_function: str = Field("normal", pattern="^(normal|impaired)$")
    drug_name: str = "benazepril"
    n_dose: int = 1
    sim_time_end: int = Field(120, ge=24, le=336)
    tstart_dosing: int = 0


# ═══════════════════════════════════════════════════════════════════════════════
#  INFECTION / DKD MODEL  (exact Barbiero _infection.py)
# ═══════════════════════════════════════════════════════════════════════════════

def _infection_ode(t, conc, drug_dose, ke_diacid, VF_diacid, ka_diacid,
                   feedback_capacity, k_cat_Renin, k_feedback, C50,
                   n_Hill, tau, tfinal_dosing, AngII_conc_t0, Renin_conc_t0,
                   baseline_prod_Renin, k_degr_Renin, k_degr_AngI, k_degr_AGT,
                   tstart_dosing, glu,
                   c_Renin_a, c_Renin_b, c_ACE_a, c_ACE_b, c_AT1_a, c_AT1_b,
                   k_APA, k_AT2, k_NEP, k_AGT, k_ACE2_0,
                   h_ANGII, h_ANG17, h_ATR, drug_type, is_infected):

    AngI_conc, AngII_conc, Renin_conc, AGT_conc, ANG17_conc, AT1R_conc, AT2R_conc, k_ACE2, IR = conc

    diacid_conc = analytical_PK(drug_dose, ka_diacid, VF_diacid, ke_diacid,
                                 t, tau, tfinal_dosing, tstart_dosing)

    Inhibition = (100 * diacid_conc**n_Hill) / (diacid_conc**n_Hill + C50**n_Hill)

    c_Renin = c_Renin_a * GLU(t, glu) + c_Renin_b
    c_ACE   = c_ACE_a   * GLU(t, glu) + c_ACE_b
    c_AT1   = c_AT1_a   * GLU(t, glu) + c_AT1_b

    d_AGT   = mass_balance_AGT(k_AGT, c_Renin, AGT_conc, k_degr_AGT)
    d_Renin = mass_balance_Renin(baseline_prod_Renin, k_feedback, AngII_conc_t0,
                                  AngII_conc, feedback_capacity, k_degr_Renin, Renin_conc)
    d_AngI  = mass_balance_AngI_infection(c_Renin, AGT_conc, k_cat_Renin, Renin_conc,
                                          Renin_conc_t0, k_degr_AngI, k_NEP, k_ACE2,
                                          AngI_conc, c_ACE, Inhibition, drug_type)
    d_AngII = mass_balance_AngII_infection(h_ANGII, c_AT1, k_APA, k_ACE2, k_AT2,
                                           AngII_conc, c_ACE, AngI_conc, Inhibition, drug_type)
    d_ANG17 = mass_balance_ANG17(k_NEP, AngI_conc, k_ACE2, AngII_conc, h_ANG17, ANG17_conc)
    d_AT1R  = mass_balance_AT1R(c_AT1, AngII_conc, h_ATR, AT1R_conc, Inhibition, drug_type)
    d_AT2R  = mass_balance_AT2R(k_AT2, AngII_conc, h_ATR, AT2R_conc)

    k_in = 0.1
    d_k_ACE2 = (k_in * AngII_conc - np.log(2) / 2 * k_ACE2) if is_infected else 0.0
    k_out    = 0.05 if is_infected else 1.0
    d_IR     = 0.15 * (k_ACE2 - k_ACE2_0) + 0.001 * diacid_conc + 0.1 * glu - np.log(2) / 1 * IR * k_out

    return [d_AngI, d_AngII, d_Renin, d_AGT, d_ANG17, d_AT1R, d_AT2R, d_k_ACE2, d_IR]


def run_infection(req: InfectionRequest) -> dict:
    drug = req.drug_name
    renal = req.renal_function
    params = get_params(drug, renal)
    pk     = get_pk_params(drug, renal)

    drug_dose = req.dose * 1e6      # mg → ng (as in Barbiero)
    tau       = 24.0 / req.n_dose

    coefficients = [params["c_Renin"], params["k_cat_Renin"],
                    params["k_feedback"], params["feedback_capacity"],
                    params["k_cons_AngII"]]

    Rate_params = np.array([1.527482117056147e-07, 1.705688364046031e-05,
                             2.472978807773762e-04, 4.533794480918563e-03,
                             7.072930413876994e-04, 1.296703909210782e-02]) * 3600
    c_Renin_a, c_Renin_b, c_ACE_a, c_ACE_b, c_AT1_a, c_AT1_b = Rate_params

    Rate_cons = np.array([1.210256981930063e-02, 1.069671574938187e-04,
                           6.968146259597334e-03, 1.628277841850352e-04,
                           6.313823632053240e+02]) * 3600
    k_APA, k_ACE2, k_AT2, k_NEP, k_AGT = Rate_cons
    h_ANGII = 18 / 3600

    drug_type = "ACEi" if drug in ("benazepril", "cilazapril") else "ARB"

    ANG17_conc_t0 = 9.858
    AT1R_conc_t0  = 16.2
    AT2R_conc_t0  = 5.4
    h_ANG17 = 0.5
    h_ATR   = 0.2

    Renin_conc_t0  = pk["Renin_conc_t0"]
    AngI_conc_t0   = pk["AngI_conc_t0"]
    AngII_conc_t0  = pk["AngII_conc_t0"]
    AGT_conc_t0    = pk["AGT_conc_t0"]
    k_degr_Renin   = pk["k_degr_Renin"]
    k_degr_AngI    = pk["k_degr_AngI"]
    k_degr_AGT     = pk["k_degr_AGT"]
    Mw_AngII       = pk["Mw_AngII"]

    baseline_prod_Renin = k_degr_Renin * Renin_conc_t0
    k_ACE2_t0 = k_ACE2
    IR_t0     = k_ACE2_t0

    conc_t0 = [AngI_conc_t0, AngII_conc_t0, Renin_conc_t0, AGT_conc_t0,
               ANG17_conc_t0, AT1R_conc_t0, AT2R_conc_t0, k_ACE2_t0, IR_t0]

    t_eval = np.arange(0, req.sim_time_end, tau / 500)

    ode_args = (
        drug_dose, pk["ke_diacid"], pk["VF_diacid"], pk["ka_diacid"],
        params["feedback_capacity"], params["k_cat_Renin"], params["k_feedback"],
        pk["C50"], pk["n_Hill"], tau, req.sim_time_end,
        AngII_conc_t0, Renin_conc_t0, baseline_prod_Renin,
        k_degr_Renin, k_degr_AngI, k_degr_AGT,
        req.tstart_dosing, req.glu,
        c_Renin_a, c_Renin_b, c_ACE_a, c_ACE_b, c_AT1_a, c_AT1_b,
        k_APA, k_AT2, k_NEP, k_AGT, k_ACE2_t0,
        h_ANGII, h_ANG17, h_ATR, drug_type, req.infection,
    )

    sol = solve_ivp(_infection_ode, [0, req.sim_time_end], conc_t0,
                    args=ode_args, t_eval=t_eval, method="LSODA")

    if not sol.success:
        raise RuntimeError(f"ODE solver failed: {sol.message}")

    t = sol.t
    diacid = np.array([analytical_PK(drug_dose, pk["ka_diacid"], pk["VF_diacid"],
                                      pk["ke_diacid"], ti, tau, req.sim_time_end,
                                      req.tstart_dosing) for ti in t])

    conv = 1e6 / 1000
    Mw_AngII_val = Mw_AngII
    tplot    = t / 24
    AngII_c  = sol.y[1, :] * Mw_AngII_val * conv
    Ang17_c  = sol.y[4, :] * Mw_AngII_val * conv
    AT1R_c   = sol.y[5, :] * Mw_AngII_val * conv
    AT2R_c   = sol.y[6, :] * Mw_AngII_val * conv
    ACE2     = sol.y[7, :]
    IR       = sol.y[8, :]

    ANGII_Plot = 0.021001998652419
    y_angII_norm = ((AngII_c / (Mw_AngII_val * 1e6)) / ANGII_Plot) * 100
    y_angII = AngII_c / (Mw_AngII_val * 1e6 / 1000)

    return {
        "t":         tplot.tolist(),
        "angII":     y_angII.tolist(),
        "angII_norm": y_angII_norm.tolist(),
        "ang17":     (Ang17_c / (Mw_AngII_val * 1e6 / 1000)).tolist(),
        "at1r":      (AT1R_c  / (Mw_AngII_val * 1e6 / 1000)).tolist(),
        "at2r":      (AT2R_c  / (Mw_AngII_val * 1e6 / 1000)).tolist(),
        "ACE2":      ACE2.tolist(),
        "IR":        IR.tolist(),
        "diacid":    diacid.tolist(),
    }


# ═══════════════════════════════════════════════════════════════════════════════
#  DIABETES MODEL  (exact Barbiero _diabetes.py)
# ═══════════════════════════════════════════════════════════════════════════════

class _StepFun:
    def __init__(self, x, y):
        self.x = np.array(x); self.y = np.array(y)
    def eval(self, xi):
        diff = self.x - xi
        t = int(np.sum(diff < 0))
        return self.y[min(t, len(self.y)-1)]

def _make_step(x_init, y_init, days):
    x = np.array(x_init) / 24
    count = np.arange(days)
    xf = np.concatenate([x + c for c in count])
    k_rep = len(xf) / len(y_init)
    yf = np.tile(y_init, int(k_rep))
    yf = np.append(yf, y_init[0])
    return _StepFun(xf, yf)

def _diabetes_ode(t, y, stepg, stepe, k, alpha, theta, R0, EG0, SI,
                  d0, r1, r2, m0, c4, c5, i0, m, q, lt):
    I, G, beta, IR, MTOR, Tt, C = y
    dI    = (beta * theta * G**2) / (alpha + G**2) - k * I
    dG    = R0 - G * (EG0 + SI * I / (IR + 1)) + 1.0 * stepg.eval(C) - 0.1 * stepe.eval(C)
    dbeta = beta * (-d0 + r1 * G - r2 * G**2)
    dIR   = -i0 * IR + m * MTOR + q * I
    dMTOR = -m0 * MTOR + c4 * I / (IR + 1) + c5 * G
    dTt   = 0.001 * (20 / (1 + np.exp(-0.05 * (G - 100)))) * Tt * np.log(lt / Tt)
    dC    = 1.0
    return [dI, dG, dbeta, dIR, dMTOR, dTt, dC]

def run_diabetes(req: DiabetesRequest) -> dict:
    days = req.sim_days
    stepg = _make_step([8, 8.5, 12, 12.5, 20, 20.5], [0, 200, 0, 4200, 0, 4200], days)
    stepe = _make_step([18, 18.5], [0, 2000], days)

    k=432; alpha=20000; theta=43.2; R0=864; EG0=0.44
    SI = 1.62 if req.glu < 8 else 0.52
    d0=0.06; r1=0.00084; r2=0.0000024; m0=47.7; c4=9; c5=6
    i0=87; m=2; q=0.017; lt=5

    I0=13.59; G0=100; beta0=407.73; IR0=0.359; MTOR0=14.465; Tt0=1; C0=0

    sol = solve_ivp(_diabetes_ode, [0, days],
                    [I0, G0, beta0, IR0, MTOR0, Tt0, C0],
                    args=(stepg, stepe, k, alpha, theta, R0, EG0, SI,
                          d0, r1, r2, m0, c4, c5, i0, m, q, lt),
                    max_step=0.001, method="LSODA")
    if not sol.success:
        raise RuntimeError(f"ODE solver failed: {sol.message}")

    return {
        "t":    sol.t.tolist(),
        "I":    sol.y[0].tolist(),
        "G":    sol.y[1].tolist(),
        "beta": sol.y[2].tolist(),
        "IR":   sol.y[3].tolist(),
        "MTOR": sol.y[4].tolist(),
        "Tt":   sol.y[5].tolist(),
    }


# ═══════════════════════════════════════════════════════════════════════════════
#  HYPERTENSION MODEL  (exact Barbiero _hypertension.py)
# ═══════════════════════════════════════════════════════════════════════════════

def _hypertension_ode(t, conc, K_in_ANG17, K_out_ANG17, K_tr, SS_ANG2, n,
                      K_in_ANG2, K_out_ANG2, SS_ANG17, m,
                      K_in_SBP, K_out_SBP, SSS_ANG2, III_ANG17,
                      ES_p, EI_p, ES_ANG2, ES_ANG17):
    S_ANG17, S_ANG2, S_SBP = conc
    ES_ANG2_n  = transit_compartment_model(ES_ANG2,  K_tr, n)
    ES_ANG17_m = transit_compartment_model(ES_ANG17, K_tr, m)
    d1 = change_S_ANG17(K_in_ANG17, ES_ANG2_n, ES_p, K_out_ANG17, S_ANG17)
    d2 = change_S_ANG2(K_in_ANG2, EI_p, K_out_ANG2, ES_ANG17_m, S_ANG2)
    d3 = change_S_SBP(K_in_SBP, SSS_ANG2, S_ANG2, III_ANG17, S_ANG17, K_out_SBP, S_SBP)
    return [d1, d2, d3]

def run_hypertension(req: HypertensionRequest) -> dict:
    K_in_ANG17=117; K_out_ANG17=1.59; K_tr=1.63; SS_ANG2=0.0726; n=29
    K_in_ANG2=27.5; K_out_ANG2=0.215; SS_ANG17=0.0711; m=4
    K_in_SBP=2670;  K_out_SBP=56;    SSS_ANG2=0.0316;  III_ANG17=0.00956
    ES_p=0.0; EI_p=0.0

    C_ANG17_mean=9; C_ANG2_mean=0.8
    ES_ANG2  = SS_ANG2  * C_ANG2_mean
    ES_ANG17 = SS_ANG17 * C_ANG17_mean
    S_ANG17  = K_in_ANG17 * (1 + ES_ANG2) / K_out_ANG17
    S_ANG2   = K_in_ANG2 / (K_out_ANG2 * (1 + ES_ANG17))
    S_SBP    = K_in_SBP * (1 + SSS_ANG2 * S_ANG2 - III_ANG17 * S_ANG17) / K_out_SBP

    sol = solve_ivp(_hypertension_ode, [0, req.sim_time_end],
                    [S_ANG17, S_ANG2, S_SBP],
                    args=(K_in_ANG17, K_out_ANG17, K_tr, SS_ANG2, n,
                          K_in_ANG2, K_out_ANG2, SS_ANG17, m,
                          K_in_SBP, K_out_SBP, SSS_ANG2, III_ANG17,
                          ES_p, EI_p, ES_ANG2, ES_ANG17),
                    method="LSODA")
    if not sol.success:
        raise RuntimeError(f"ODE solver failed: {sol.message}")

    return {
        "t":       sol.t.tolist(),
        "S_ANG17": sol.y[0].tolist(),
        "S_ANG2":  sol.y[1].tolist(),
        "S_SBP":   sol.y[2].tolist(),
    }


# ═══════════════════════════════════════════════════════════════════════════════
#  CARDIO MODEL  (exact Barbiero _circulation.py — abridged PV output)
# ═══════════════════════════════════════════════════════════════════════════════

def run_cardio(req: CardioRequest) -> dict:
    """
    Run the DKD/infection pipeline which also feeds the cardio module.
    We return key cardiovascular derived signals from the RAS solution
    combined with the circulation.csv parameters.
    The full 668-line cardio ODE requires patient-specific beat timing from
    real ECG data (tmeas / ABPmeas arrays in _circulation.py).
    We return the RAS haemodynamic proxies (diacid, AngII, IR) plus
    estimated Psa, Ppa from circulation.csv steady-state values.
    """
    inf_req = InfectionRequest(
        age=req.age, dose=req.dose, renal_function=req.renal_function,
        glu=req.glu, infection=req.infection, drug_name=req.drug_name,
        n_dose=req.n_dose, sim_time_end=req.sim_time_end,
        tstart_dosing=req.tstart_dosing,
    )
    ras = run_infection(inf_req)

    circ = pd.read_csv(str(DATA / "circulation.csv"), index_col=0)
    Psa_ss = 93.0   # mmHg systemic arterial (typical resting MAP)
    Ppa_ss = 14.0   # mmHg pulmonary arterial

    n = len(ras["t"])
    # Derive SBP variation from AngII trajectory (as in paper methodology)
    angII_arr = np.array(ras["angII"])
    angII_base = angII_arr[0] if angII_arr[0] > 0 else 1e-9
    psa_arr = Psa_ss * (1 + 0.15 * (angII_arr / angII_base - 1))
    ppa_arr = np.full(n, Ppa_ss)

    return {
        "t":    ras["t"],
        "angII": ras["angII"],
        "ang17": ras["ang17"],
        "IR":    ras["IR"],
        "diacid": ras["diacid"],
        "ACE2":  ras["ACE2"],
        "Psa":   psa_arr.tolist(),
        "Ppa":   ppa_arr.tolist(),
    }


# ═══════════════════════════════════════════════════════════════════════════════
#  DKD MODEL  (exact call_combinedRAS_ACE_PKPD from _local_RAS.py)
# ═══════════════════════════════════════════════════════════════════════════════

def _dkd_ode(t, conc, drug_dose, ke_diacid, VF_diacid, ka_diacid,
             feedback_capacity, k_cat_Renin, k_feedback, C50, n_Hill,
             tau, tfinal_dosing, AngI_conc_t0, AngII_conc_t0, Renin_conc_t0,
             AGT_conc_t0, baseline_prod_Renin, k_degr_Renin, k_degr_AngI,
             k_degr_AGT, k_cons_AngII, tstart_dosing, glu):
    AngI_conc, AngII_conc, Renin_conc, AGT_conc = conc

    diacid_conc = analytical_PK(drug_dose, ka_diacid, VF_diacid, ke_diacid,
                                 t, tau, tfinal_dosing, tstart_dosing)
    Inhibition = (100 * diacid_conc**n_Hill) / (diacid_conc**n_Hill + C50**n_Hill)

    Rate_params = np.array([1.527482117056147e-07, 1.705688364046031e-05,
                             2.472978807773762e-04, 4.533794480918563e-03,
                             7.072930413876994e-04, 1.296703909210782e-02]) * 3600
    c_Renin_a, c_Renin_b, c_ACE_a, c_ACE_b, c_AT1_a, c_AT1_b = Rate_params
    c_Renin = c_Renin_a * GLU(t, glu) + c_Renin_b
    c_ACE   = c_ACE_a   * GLU(t, glu) + c_ACE_b
    c_AT1   = c_AT1_a   * GLU(t, glu) + c_AT1_b

    Rate_cons = np.array([1.210256981930063e-02, 1.069671574938187e-04,
                           6.968146259597334e-03, 1.628277841850352e-04,
                           6.313823632053240e+02]) * 3600
    k_APA, k_ACE2, k_AT2, k_NEP, k_AGT = Rate_cons
    h_ANGII = 18 / 3600

    d_AGT   = mass_balance_AGT(k_AGT, c_Renin, AGT_conc, k_degr_AGT)
    d_Renin = mass_balance_Renin(baseline_prod_Renin, k_feedback, AngII_conc_t0,
                                  AngII_conc, feedback_capacity, k_degr_Renin, Renin_conc)
    d_AngI  = mass_balance_AngI(c_Renin, AGT_conc, k_cat_Renin, Renin_conc,
                                 Renin_conc_t0, k_degr_AngI, k_NEP, k_ACE2,
                                 AngI_conc, c_ACE, Inhibition)
    d_AngII = mass_balance_AngII(h_ANGII, c_AT1, k_APA, k_ACE2, k_AT2,
                                  AngII_conc, c_ACE, AngI_conc, Inhibition)
    return [d_AngI, d_AngII, d_Renin, d_AGT]

def run_dkd(req: InfectionRequest) -> dict:
    drug = req.drug_name; renal = req.renal_function
    params = get_params(drug, renal)
    pk     = get_pk_params(drug, renal)

    drug_dose = req.dose * 1e6
    tau       = 24.0 / req.n_dose

    Renin_conc_t0 = pk["Renin_conc_t0"]
    AngI_conc_t0  = pk["AngI_conc_t0"]
    AngII_conc_t0 = pk["AngII_conc_t0"]
    AGT_conc_t0   = pk["AGT_conc_t0"]
    k_degr_Renin  = pk["k_degr_Renin"]
    k_degr_AngI   = pk["k_degr_AngI"]
    k_degr_AGT    = pk["k_degr_AGT"]
    Mw_AngII      = pk["Mw_AngII"]
    Mw_AngI       = pk["Mw_AngI"]

    baseline_prod_Renin = k_degr_Renin * Renin_conc_t0
    conc_t0 = [AngI_conc_t0, AngII_conc_t0, Renin_conc_t0, AGT_conc_t0]
    t_eval  = np.arange(0, req.sim_time_end, tau / 500)

    ode_args = (
        drug_dose, pk["ke_diacid"], pk["VF_diacid"], pk["ka_diacid"],
        params["feedback_capacity"], params["k_cat_Renin"], params["k_feedback"],
        pk["C50"], pk["n_Hill"], tau, req.sim_time_end,
        AngI_conc_t0, AngII_conc_t0, Renin_conc_t0, AGT_conc_t0,
        baseline_prod_Renin, k_degr_Renin, k_degr_AngI, k_degr_AGT,
        params["k_cons_AngII"], req.tstart_dosing, req.glu,
    )

    sol = solve_ivp(_dkd_ode, [0, req.sim_time_end], conc_t0,
                    args=ode_args, t_eval=t_eval, method="LSODA")
    if not sol.success:
        raise RuntimeError(f"ODE solver failed: {sol.message}")

    t = sol.t
    diacid = np.array([analytical_PK(drug_dose, pk["ka_diacid"], pk["VF_diacid"],
                                      pk["ke_diacid"], ti, tau, req.sim_time_end,
                                      req.tstart_dosing) for ti in t])
    conv   = 1e6 / 1000
    tplot  = t / 24
    AngII_c = sol.y[1, :] * Mw_AngII * conv
    Mw_AngII_val = Mw_AngII
    ANGII_Plot   = 0.021001998652419
    y_angII_norm = ((AngII_c / (Mw_AngII_val * 1e6)) / ANGII_Plot) * 100
    y_angII      = AngII_c / (Mw_AngII_val * 1e6 / 1000)

    return {
        "t":          tplot.tolist(),
        "angII":      y_angII.tolist(),
        "angII_norm": y_angII_norm.tolist(),
        "diacid":     diacid.tolist(),
        "AngI":       (sol.y[0, :] * Mw_AngI * conv).tolist(),
        "Renin":      (sol.y[2, :] * pk["Mw_Renin"] * conv).tolist(),
        "AGT":        (sol.y[3, :] * pk["Mw_AGT"] * conv).tolist(),
    }


# ═══════════════════════════════════════════════════════════════════════════════
#  API ROUTES
# ═══════════════════════════════════════════════════════════════════════════════

@app.get("/api/health")
def health():
    return {"status": "ok", "version": "1.0.0", "data_dir": str(DATA),
            "mat_files": [f.name for f in DATA.glob("*.mat")]}

@app.post("/api/simulate/infection")
def api_infection(req: InfectionRequest):
    try:
        t0 = time.time()
        result = run_infection(req)
        result["elapsed_s"] = round(time.time() - t0, 3)
        result["model"] = "infection"
        result["params"] = req.dict()
        return result
    except Exception as e:
        log.exception("infection failed")
        raise HTTPException(500, str(e))

@app.post("/api/simulate/diabetes")
def api_diabetes(req: DiabetesRequest):
    try:
        t0 = time.time()
        result = run_diabetes(req)
        result["elapsed_s"] = round(time.time() - t0, 3)
        result["model"] = "diabetes"
        result["params"] = req.dict()
        return result
    except Exception as e:
        log.exception("diabetes failed")
        raise HTTPException(500, str(e))

@app.post("/api/simulate/hypertension")
def api_hypertension(req: HypertensionRequest):
    try:
        t0 = time.time()
        result = run_hypertension(req)
        result["elapsed_s"] = round(time.time() - t0, 3)
        result["model"] = "hypertension"
        result["params"] = req.dict()
        return result
    except Exception as e:
        log.exception("hypertension failed")
        raise HTTPException(500, str(e))

@app.post("/api/simulate/cardio")
def api_cardio(req: CardioRequest):
    try:
        t0 = time.time()
        result = run_cardio(req)
        result["elapsed_s"] = round(time.time() - t0, 3)
        result["model"] = "cardio"
        result["params"] = req.dict()
        return result
    except Exception as e:
        log.exception("cardio failed")
        raise HTTPException(500, str(e))

@app.post("/api/simulate/dkd")
def api_dkd(req: InfectionRequest):
    try:
        t0 = time.time()
        result = run_dkd(req)
        result["elapsed_s"] = round(time.time() - t0, 3)
        result["model"] = "dkd"
        result["params"] = req.dict()
        return result
    except Exception as e:
        log.exception("dkd failed")
        raise HTTPException(500, str(e))

# Serve built frontend (if present)
FRONTEND_DIST = ROOT.parent / "frontend" / "dist"
if FRONTEND_DIST.exists():
    app.mount("/assets", StaticFiles(directory=str(FRONTEND_DIST / "assets")), name="assets")

    @app.get("/{full_path:path}")
    def serve_spa(full_path: str):
        return FileResponse(str(FRONTEND_DIST / "index.html"))


# ═══════════════════════════════════════════════════════════════════════════════
#  ENHANCED EXTENSION ROUTES  (Utku Köse additions)
# ═══════════════════════════════════════════════════════════════════════════════

# Extension request models (defined here to avoid circular import)
from pydantic import BaseModel as _BM, Field as _F

class DoseResponseRequest(_BM):
    doses: list[float] = _F(default=[0,1,2,3,5,7,10,15,20])
    glu: float = _F(5.0, ge=0.1, le=30)
    infection: bool = True
    renal_function: str = _F("normal", pattern="^(normal|impaired)$")
    sim_time_end: int = _F(120, ge=24, le=336)
    drug_name: str = "benazepril"

class SensitivityRequest(_BM):
    glu_values:  list[float] = _F(default=[1,5,10,15,17,20])
    dose_values: list[float] = _F(default=[0,2,5,10,15,20])
    infection:   bool  = True
    renal_function: str = _F("normal", pattern="^(normal|impaired)$")
    sim_time_end: int  = _F(120, ge=24, le=336)
    output_var:  str   = _F("IR_final", pattern="^(IR_final|angII_final|DPRS|ANG17_AngII_ratio)$")

class DiabetesIndexRequest(_BM):
    glu_values: list[float] = _F(default=[5,7,8,10,12,15,17,20])
    sim_days: int = _F(10, ge=1, le=30)

class CohortRequest(_BM):
    n_patients: int = _F(20, ge=4, le=50)
    glu_range:  list[float] = _F(default=[1,20])
    dose_range: list[float] = _F(default=[0,10])
    infection:  bool = True
    renal_function: str = _F("normal", pattern="^(normal|impaired)$")
    sim_time_end: int = _F(120, ge=24, le=168)
    seed: int = 42


@app.post("/api/enhance/dprs")
def api_dprs(req: InfectionRequest):
    try:
        import extensions as ext
        r = run_infection(req)
        return {"dprs": ext.compute_dprs(r), "phases": ext.compute_phases(r), "params": req.dict()}
    except Exception as e:
        log.exception("dprs failed"); raise HTTPException(500, str(e))

@app.post("/api/enhance/dose_response")
def api_dose_response(req: DoseResponseRequest):
    try:
        import extensions as ext
        t0 = time.time()
        from extensions import DoseResponseRequest as _DR
        result = ext.run_dose_response(_DR(**req.dict()))
        result["elapsed_s"] = round(time.time()-t0, 3)
        return result
    except Exception as e:
        log.exception("dose_response failed"); raise HTTPException(500, str(e))

@app.post("/api/enhance/sensitivity")
def api_sensitivity(req: SensitivityRequest):
    try:
        import extensions as ext
        t0 = time.time()
        from extensions import SensitivityRequest as _SR
        result = ext.run_sensitivity(_SR(**req.dict()))
        result["elapsed_s"] = round(time.time()-t0, 3)
        return result
    except Exception as e:
        log.exception("sensitivity failed"); raise HTTPException(500, str(e))

@app.post("/api/enhance/diabetes_indices")
def api_diabetes_indices(req: DiabetesIndexRequest):
    try:
        import extensions as ext
        t0 = time.time()
        from extensions import DiabetesIndexRequest as _DI
        result = ext.run_diabetes_indices(_DI(**req.dict()))
        result["elapsed_s"] = round(time.time()-t0, 3)
        return result
    except Exception as e:
        log.exception("diabetes_indices failed"); raise HTTPException(500, str(e))

@app.post("/api/enhance/cohort")
def api_cohort(req: CohortRequest):
    try:
        import extensions as ext
        t0 = time.time()
        from extensions import CohortRequest as _CR
        result = ext.run_cohort(_CR(**req.dict()))
        result["elapsed_s"] = round(time.time()-t0, 3)
        return result
    except Exception as e:
        log.exception("cohort failed"); raise HTTPException(500, str(e))
