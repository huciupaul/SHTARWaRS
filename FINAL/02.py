"""
02_run_rich_stage.py
--------------------
Iterates over the *cases_setup.json* produced by **01_setup_cases.py**, runs a
free premixed rich flame (φ₁ = 5) for each operating point using Cantera, and
serialises the burnt‑gas state plus thermodynamic fluxes required for the
lean‑stage mapper.

Outputs
-------
rich_cases.json   # persistent, power‑indexed dictionary of rich‑stage data

Execution
---------
$ python 02_run_rich_stage.py
"""
from __future__ import annotations

import json
from pathlib import Path
from datetime import datetime

import cantera as ct
ct.suppress_thermo_warnings()

from mixture_properties import (
    mixture_properties,
    kerosene_to_h2,
    air_mass_flow_for_phi,
    MECH,
)

# ---------------------------------------------------------------------------
# 0 – CONSTANTS & INPUTS
# ---------------------------------------------------------------------------
PHI_RICH = 5.0                           # φ₁ maintained across all power points
X_AIR    = "N2:0.78084, O2:0.20946"      # dry air, mole fractions
X_H2     = "H2:1"
WIDTH_RICH = 0.04                       # m – domain length for free flame

# ---------------------------------------------------------------------------
# 1 – LOAD OPERATING POINTS FROM SETUP
# ---------------------------------------------------------------------------
setup_file = Path("cases_setup.json")
if not setup_file.exists():
    raise FileNotFoundError("Run 01_setup_cases.py first; 'cases_setup.json' missing.")

cases_setup = json.loads(setup_file.read_text())

# container for (possibly pre‑existing) results
rich_file = Path("rich_cases.json")
if rich_file.exists():
    rich_data = json.loads(rich_file.read_text())
    # build fast lookup by power so we can update/append
    existing = {c["power_kW"]: c for c in rich_data.get("cases", [])}
else:
    rich_data = {"created": cases_setup["created"], "MECH": MECH, "cases": []}
    existing  = {}

# ---------------------------------------------------------------------------
# 2 – RUN EACH CASE
# ---------------------------------------------------------------------------
for case in cases_setup["cases"]:
    pwr = case["power_kW"]
    if pwr in existing:
        print(f"[rich] {pwr:.1f} kW → already solved; skipping.")
        continue

    print(f"[rich] Solving free flame for {pwr:.1f} kW …")

    T_inlet = case["T_inlet_K"]
    P_inlet = case["P_inlet_Pa"]
    mdot_H2 = case["m_H2"]
    mdot_air1 = case["mdot_air_rich"]
    TOTAL_AIR = case["m_air"]

    # -----------------------------------------------------------------------
    # 2.1 – build premixed reactant state at compressor exit conditions
    # -----------------------------------------------------------------------
    X_mix, mdot_rich_burn, _ = mixture_properties(
        mdot1=mdot_H2, X1=X_H2,
        mdot2=mdot_air1, X2=X_AIR,
        T=T_inlet, P=P_inlet,
    )

    gas = ct.Solution(MECH)
    gas.transport_model = "multicomponent"
    gas.TPX = T_inlet, P_inlet, X_mix

    flame_rich = ct.FreeFlame(gas, width=WIDTH_RICH)
    flame_rich.set_refine_criteria(ratio=6.0, slope=0.01, curve=0.02)
    flame_rich.solve(loglevel=0, auto=True, refine_grid=True)

    # -----------------------------------------------------------------------
    # 2.2 – burnt gas state immediately downstream of the rich stage
    # -----------------------------------------------------------------------
    exit_gas = ct.Solution(MECH)
    exit_gas.TPX = flame_rich.T[-1], P_inlet, flame_rich.X[:, -1]

    mdot_exh = mdot_rich_burn
    h_exh    = exit_gas.enthalpy_mass              # J kg⁻¹

    idx = gas.species_index
    X_NO  = flame_rich.X[idx("NO"),  -1] if idx("NO")  >= 0 else 0.0
    X_N2O = flame_rich.X[idx("N2O"), -1] if idx("N2O") >= 0 else 0.0

    # -----------------------------------------------------------------------
    # 2.3 – assemble case dictionary
    # -----------------------------------------------------------------------
    result = {
        "power_kW":        pwr,
        "T_inlet_K":       T_inlet,
        "P_inlet_Pa":      P_inlet,
        "T_peak_K":        float(flame_rich.T.max()),
        "T_burnt_K":       float(flame_rich.T[-1]),
        "X_exh":           exit_gas.X.tolist(),   # mole fraction array
        "mdot_exh":        mdot_exh,
        "h_exh_Jperkg":    h_exh,
        "TOTAL_AIR":       TOTAL_AIR,
        "mdot_air1":       mdot_air1,
        "grid_points":     len(flame_rich.grid),
        "NOx_ppm":         float((X_NO + X_N2O) * 1e6),
    }

    existing[pwr] = result

# ---------------------------------------------------------------------------
# 3 – WRITE/UPDATE OUTPUT FILE
# ---------------------------------------------------------------------------
rich_data["cases"] = list(existing.values())
rich_data["updated"] = datetime.utcnow().isoformat() + "Z"

rich_file.write_text(json.dumps(rich_data, indent=2))
print(f"[rich] {len(rich_data['cases'])} cases written → {rich_file.resolve()}")
