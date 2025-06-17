from __future__ import annotations
import json, os, sys
from pathlib import Path
from typing import Tuple, Any

import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

from tqdm import tqdm
import cantera as ct

from mixture_properties import (
    mixture_properties,
    air_mass_flow_for_phi,
    sensible_h_mass,
)


def load_json(path: Path) -> dict:
    """Read a JSON file, accepting an optional UTF-8 BOM."""
    return json.loads(path.read_text(encoding="utf-8-sig"))


# ───────────────────────────────────────────────────────────────────────────────
# USER CONSTANTS
# ───────────────────────────────────────────────────────────────────────────────
T_H2O       = 353.15          # K  (80 °C liquid / vapour injection)
GRID_NODES  = 30              # φ₂ × H₂O resolution per map
X_AIR       = "N2:0.78084, O2:0.20946"
width_lean  = 0.0001          # m – computational domain length
PHI_MIN, PHI_MAX = 0.25, 0.90 # lean combustor equivalence-ratio range
H2O_MIN, H2O_MAX = 0.00, 0.15 # kg s⁻¹ water sweep

# ───────────────────────────────────────────────────────────────────────────────
# 1.  GENERIC LEAN-STAGE SOLVER (with flame‑out handling)
# ───────────────────────────────────────────────────────────────────────────────

def run_case(args: Tuple[float, float, dict,
                         float, float, float,
                         float, float, str, str, float]
            ) -> tuple[Any, Any, float, float, float, float,
                       float, float, float, float, float, float, int]:

    (f, mdot_H2O, data, mdot_exh, h_exh, h_air2, h_H2O,
     TOTAL_AIR, mdot_air1, X_AIR, MECH, width) = args

    # ----- data from rich burn -------------------------------------------------
    P_mix = data["P_inlet_Pa"]
    X_exh = data["X_exh"]

    # Build exhaust composition dict if needed
    g_exh = ct.Solution(MECH)
    g_exh.X = X_exh
    g_exh.HP = h_exh, P_mix

    if isinstance(X_exh, list):
        spec_names = g_exh.species_names
        X_exh_dict = {s: X_exh[i] for i, s in enumerate(spec_names)
                      if X_exh[i] > 0.0}
    else:
        X_exh_dict = X_exh

    # ----- massflow bookkeeping ----------------------------------------------
    AVALIABLE_AIR = TOTAL_AIR * (1.0 - 0.0525)        # bleed fraction
    mdot_cool = f * AVALIABLE_AIR
    mdot_air2 = AVALIABLE_AIR - mdot_air1 - mdot_cool

    # ----- build cantera objects ----------------------------------------------
    g_exh = ct.Solution(MECH)        ; g_exh.X = X_exh ; g_exh.HP = h_exh,  P_mix
    g_air = ct.Solution(MECH)        ; g_air.X = X_AIR ; g_air.HP = h_air2, P_mix
    g_H2O = ct.Water()               ; g_H2O.TP = T_H2O, P_mix

    # sensible enthalpies
    h_exh_s  = sensible_h_mass(g_exh)
    h_air_s  = sensible_h_mass(g_air)
    h_H2O_s  = sensible_h_mass(g_H2O)

    Hdot_in  = (mdot_exh  * h_exh_s +
                mdot_air2 * h_air_s +
                mdot_H2O  * h_H2O_s)

    # premix + water
    X_mix1, mdot_mix1, _ = mixture_properties(
        X1=X_exh_dict, X2=X_AIR, mdot1=mdot_exh, mdot2=mdot_air2)
    X_mix,  mdot_tot, _  = mixture_properties(
        X1=X_mix1, X2="H2O:1", mdot1=mdot_mix1, mdot2=mdot_H2O)

    h_mix = (mdot_exh * g_exh.enthalpy_mass +
             mdot_air2 * g_air.enthalpy_mass +
             mdot_H2O  * g_H2O.enthalpy_mass) / mdot_tot

    gas = ct.Solution(MECH)
    gas.X, gas.HP = X_mix, (h_mix, P_mix)
    gas.transport_model = "mixture-averaged"

    flame = ct.FreeFlame(gas, width=width)
    flame.set_refine_criteria(ratio=5.0, slope=0.05, curve=0.02)
    flame.set_max_grid_points(flame.flame, 10_000)

    # ── flame‑out fallback ⬇︎ ---------------------------------------------------
    flameout = 0  # 0 = normal combustion, 1 = flame‑out (no ignition)
    try:
        flame.solve(loglevel=0, auto=True, refine_grid=True)
    except ct.CanteraError:
        # Flame could not be established – treat as non‑reacting stream
        flameout = 1

    # ----- compute outputs -----------------------------------------------------
    if flameout:
        # No combustion: exit gas equals inlet mixture
        exit_gas = gas
        Tmax           = gas.T
        NOx_ppm        = 0.0
    else:
        exit_gas = ct.Solution(MECH)
        exit_gas.TPX = flame.T[-1], P_mix, flame.X[:, -1]
        Tmax = float(flame.T.max())
        idx = gas.species_index
        X_NO  = flame.X[idx("NO"),  -1] if idx("NO")  >= 0 else 0.0
        X_N2O = flame.X[idx("N2O"), -1] if idx("N2O") >= 0 else 0.0
        NOx_ppm = float((X_NO + X_N2O) * 1e6)

    # Common post‑processing ----------------------------------------------------
    h_mix_s   = sensible_h_mass(gas)
    h_exit_s  = sensible_h_mass(exit_gas)
    Qdot_rel  = mdot_tot * (h_exit_s - h_mix_s)

    # Mix with cooling air to get TIT conditions
    X_mix_out, mdot_TIT, _ = mixture_properties(
        X1=exit_gas.X, X2=X_AIR, mdot1=mdot_tot, mdot2=mdot_cool)
    h_mix_out = (mdot_tot  * exit_gas.enthalpy_mass +
                 mdot_cool * g_air.enthalpy_mass) / mdot_TIT

    gas_TIT = ct.Solution(MECH)
    gas_TIT.X = X_mix_out
    gas_TIT.HP = h_mix_out, P_mix

    gamma = gas_TIT.cp_mass / gas_TIT.cv_mass
    T_total = gas_TIT.T * (1.0 + 0.5 * (gamma - 1.0) * 0.3**2)
    Hdot_out = mdot_TIT * gas_TIT.cp_mass * T_total

    # Species mole fractions we always have (either from flame or gas)
    X_O2  = float(exit_gas["O2"].X)
    X_H2O = float(exit_gas["H2O"].X)

    return (gas.equivalence_ratio(), mdot_H2O,
            Tmax, NOx_ppm, X_O2, X_H2O,
            float(gamma), Hdot_in, Hdot_out, Qdot_rel,
            float(gas_TIT.T), float(mdot_TIT), flameout)


# ───────────────────────────────────────────────────────────────────────────────
# helper: fractional‑air→φ₂ relationship
# ───────────────────────────────────────────────────────────────────────────────

def f_for_phi2(phi2: np.ndarray,
               mdot_air_stoich: float,
               *, TOTAL_AIR: float, mdot_air1: float) -> np.ndarray:
    return (TOTAL_AIR - mdot_air1 - mdot_air_stoich / phi2) / TOTAL_AIR

# ───────────────────────────────────────────────────────────────────────────────
# 2.  MAIN DRIVER
# ───────────────────────────────────────────────────────────────────────────────

def map_single_case(setup: dict, rich: dict, mech: str) -> None:
    """Build φ₂–H₂O map for the given power setting.

    Now records a `flameout` flag (0 = normal, 1 = no combustion) so that
    post‑processing scripts can omit or highlight failed points.
    """
    power = int(setup["power_kW"])

    mdot_H2   = setup["m_H2"]                 # kg s⁻¹ of H₂ from fuel split
    TOTAL_AIR = rich ["TOTAL_AIR"]            # total air from comp. map
    mdot_air1 = rich ["mdot_air1"]            # air already burnt in rich stage
    mdot_exh  = rich ["mdot_exh"]             # kg s⁻¹ exhaust from rich stage
    h_exh     = rich ["h_exh_Jperkg"]         # J kg⁻¹
    T_in      = setup["T_inlet_K"]            # K
    P_mix     = rich ["P_inlet_Pa"]           # Pa

    # lean‑stage inlet air enthalpy
    g_air = ct.Solution(mech); g_air.TPX = T_in, P_mix, X_AIR
    h_air2 = g_air.enthalpy_mass

    # pre‑compute water enthalpies
    steam = ct.Solution("gri30.yaml"); steam.TPX = T_H2O, P_mix, "H2O:1"
    water = ct.Water(); water.TP = T_H2O, P_mix
    h_H2O = steam.enthalpy_mass - (steam.enthalpy_mass - water.enthalpy_mass)

    # ---------- structured sweep meshes ---------------------------------------
    mdot_air_stoich = air_mass_flow_for_phi(mdot_H2, 1.0, X_air=X_AIR)
    phi2_targets = np.linspace(PHI_MIN, PHI_MAX, 50)
    F_SWEEP = f_for_phi2(phi2_targets,
                         mdot_air_stoich,
                         TOTAL_AIR=TOTAL_AIR,
                         mdot_air1=mdot_air1)
    W_SWEEP = np.linspace(H2O_MIN, H2O_MAX, 50)

    # ---------- pack tasks -----------------------------------------------------
    args_common = (rich, mdot_exh, h_exh, h_air2, h_H2O,
                   TOTAL_AIR, mdot_air1, X_AIR, mech, width_lean)
    tasks = [(f, w, *args_common) for f in F_SWEEP for w in W_SWEEP]

    # ---------- run in parallel ------------------------------------------------
    results = []
    with ProcessPoolExecutor() as exe:
        futures = [exe.submit(run_case, t) for t in tasks]
        for fut in tqdm(as_completed(futures), total=len(futures),
                        desc=f"  [{power:4d} kW] lean map"):
            results.append(fut.result())
    results = np.asarray(results, dtype=float)

    # ---------- store ----------------------------------------------------------
    np.savez(f"lean_map_{power}kW.npz",
             phi2      = results[:,0],
             water     = results[:,1],
             Tmax      = results[:,2],
             NOx_wet   = results[:,3],
             X_O2      = results[:,4],
             X_H2O     = results[:,5],
             gamma     = results[:,6],
             H_in      = results[:,7],
             H_out     = results[:,8],
             Qrel      = results[:,9],
             TIT       = results[:,10],
             mdot_TIT  = results[:,11],
             flameout  = results[:,12])

    print(f"→ lean_map_{power}kW.npz  saved.")


# ───────────────────────────────────────────────────────────────────────────────
# 3.  PROGRAM ENTRY‑POINT
# ───────────────────────────────────────────────────────────────────────────────

def main() -> None:

    rich_path  = Path("rich_cases.json")
    setup_path = Path("cases_setup.json")
    if not (rich_path.exists() and setup_path.exists()):
        sys.exit("ERROR: Missing rich_cases.json or cases_setup.json. "
                 "Run scripts 01 and 02 first.")

    rich_json = load_json(rich_path)
    setup_json = load_json(setup_path)

    mech_file  = setup_json["MECH"]

    # build quick index: power_kW → dict
    rich_idx  = {int(c["power_kW"]): c for c in rich_json ["cases"]}
    setup_idx = {int(c["power_kW"]): c for c in setup_json["cases"]}

    produced = []
    for p_kW, setup_case in setup_idx.items():
        if p_kW not in rich_idx:
            print(f"power={p_kW} kW present in setup but missing in rich file.")
            continue
        map_single_case(setup_case, rich_idx[p_kW], mech_file)
        produced.append(p_kW)

    if produced:
        Path("lean_maps_index.txt").write_text(
            "\n".join(f"lean_map_{p}kW.npz" for p in sorted(produced)))
        print("\nDone.  Index written → lean_maps_index.txt")


# ───────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    mp.freeze_support()
    mp.set_start_method("spawn", force=True)
    main()
