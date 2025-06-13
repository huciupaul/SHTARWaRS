from __future__ import annotations
import json, os
from pathlib import Path
from typing import Tuple, Any
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import cantera as ct

from mixture_properties import mixture_properties,air_mass_flow_for_phi,sensible_h_mass


# GENERAL DATA AND INPUTS

# INPUTS
T_H2O = 353.15  # K  (80 °C)
GRID_NODES = 38




#DATA
mdot_H2 = 0.029449 #  kg/s fuel @ TOGA ONLY USED TO CALCULATE FUEL WATER RATIO
width_lean = 0.0001  # m
X_AIR = "N2:0.78084, O2:0.20946"

LHV_H2 = 120.0e6   # J kg⁻¹
INPUT_ENTHALPY = 1586918.27373181  # W
TOTAL_CHEMICAL_POWER = mdot_H2*LHV_H2  # 3.53388 MW


# Running one simulation
def run_case(args: Tuple[float, float, dict,
                        float, float, float,
                        float,float, str, str, float]

            ) -> tuple[Any, Any, float, float, float, float, float, float | Any, float, float, float, float] | tuple[
    Any, Any, float, float, float, float, float, float, float, float, float]:

    #colletcting inputs
    (f, mdot_H2O, data, mdot_exh, h_exh, h_air2, h_H2O,
     TOTAL_AIR, mdot_air1, X_AIR, MECH, width) = args

    # Collect the data from the exhaust
    P_mix = data["P_mix"]
    X_exh = data["X_exh"]

    #Split air into Lean and cooling
    AVALIABLE_AIR = TOTAL_AIR*(1-0.0525)
    mdot_cool = f * AVALIABLE_AIR
    mdot_air2 = AVALIABLE_AIR - mdot_air1 - f * AVALIABLE_AIR

    g_exh = ct.Solution(MECH);           g_exh.X = X_exh;  g_exh.HP = h_exh,  P_mix
    g_air = ct.Solution(MECH);           g_air.X = X_AIR;   g_air.HP = h_air2, P_mix
    g_H2O = ct.Water();                  g_H2O.TP = T_H2O, P_mix

    h_exh_s  = sensible_h_mass(g_exh)    # J kg⁻¹  (Δh from 298 K)
    h_air_s  = sensible_h_mass(g_air)
    h_H2O_s  = sensible_h_mass(g_H2O)

    Hdot_in  = (mdot_exh  * h_exh_s +
                mdot_air2 * h_air_s +
                mdot_H2O  * h_H2O_s)      # W

    X_mix1, mdot_mix1, _ = mixture_properties(
        X1=X_exh, X2=X_AIR,
        mdot1=mdot_exh, mdot2=mdot_air2)

    X_mix, mdot_total, _ = mixture_properties(
        X1=X_mix1, X2="H2O:1",
        mdot1=mdot_mix1, mdot2=mdot_H2O)

    h_mix = (mdot_exh * g_exh.enthalpy_mass +
             mdot_air2 * g_air.enthalpy_mass +
             mdot_H2O  * g_H2O.enthalpy_mass) / mdot_total

    gas = ct.Solution(MECH)
    gas.X, gas.HP = X_mix, (h_mix, P_mix)
    gas.transport_model = "mixture-averaged"

    flame = ct.FreeFlame(gas, width=width)
    flame.set_refine_criteria(ratio=5.0, slope=0.05, curve=0.02)
    flame.set_max_grid_points(flame.flame, 10_000)
    flame.solve(loglevel=0, auto=True, refine_grid=True)

    # burnt gas immediately after rich-quench flame
    exit_gas = ct.Solution(MECH)
    exit_gas.TPX = flame.T[-1], P_mix, flame.X[:, -1]

    # sensible enthalpy of premix vs. burnt
    h_mix_sens  = sensible_h_mass(gas)
    h_exit_sens = sensible_h_mass(exit_gas)

    # chemical heat released in this stage
    Qdot_release = mdot_total * (h_exit_sens - h_mix_sens)   # W


    # Dilution with the cooling air  →  TIT mixture
    X_mix_out, mdot_TIT, _ = mixture_properties(
        X1=exit_gas.X, X2=X_AIR,
        mdot1=mdot_total, mdot2=mdot_cool)

    h_mix_out = (mdot_total * exit_gas.enthalpy_mass +
                 mdot_cool  * g_air.enthalpy_mass) / mdot_TIT

    gas_TIT = ct.Solution(MECH)
    gas_TIT.X = X_mix_out
    gas_TIT.HP = h_mix_out, P_mix

    h_TIT_sens = sensible_h_mass(gas_TIT)
    Hdot_out   = mdot_TIT * h_TIT_sens                   # W

    # Fetch data and return
    idx = gas.species_index
    X_NO  = flame.X[idx("NO"),  -1] if idx("NO")  >= 0 else 0.0
    X_N2O = flame.X[idx("N2O"), -1] if idx("N2O") >= 0 else 0.0
    X_O2  = flame.X[idx("O2"),  -1]
    X_H2O = flame.X[idx("H2O"), -1]

    gamma = gas_TIT.cp_mass / gas_TIT.cv_mass

    return (gas.equivalence_ratio(),          # φ₂
            mdot_H2O,                         # water rate
            float(flame.T.max()),             # Tmax
            float((X_NO + X_N2O) * 1e6),      # NOx ppm
            float(X_O2), float(X_H2O),
            float(gamma),
            Hdot_in,                          # sensible IN
            Hdot_out,                         # sensible OUT (after dilution)
            Qdot_release,                     # chem. heat in second stage
            float(gas_TIT.T),                 # TIT
            float(mdot_TIT))                  # total ṁ into turbine






#HELPER FUNCTION TO GET EQUIV. RATIO FROM THE FRACTION OF WATER
def f_for_phi2(phi2: np.ndarray,
                   mdot_air_stoich: float,
                   *,
                   TOTAL_AIR: float,
                   mdot_air1: float) -> np.ndarray:

    return (TOTAL_AIR - mdot_air1 - mdot_air_stoich / phi2) / TOTAL_AIR











# -----------------------------------------------------------
#          MAIN WORKFLOW OF THE CODE; THREADS AND STUFF)
# -----------------------------------------------------------

def main() -> None:

    # gET THE VALUES FROM RICH CC
    rich_file = Path("rich_to_lean.json")
    if not rich_file.exists():
        raise FileNotFoundError("Run the rich-stage script first; "
                                "missing 'rich_to_lean.json'.")
    data = json.loads(rich_file.read_text())

    # unpack the data from Rich CC
    MECH       = data["MECH"]
    mdot_exh   = data["mdot_exh"]
    h_exh      = data["h_exh"]
    h_air2     = data["h_air2"]
    TOTAL_AIR  = data["TOTAL_AIR"]
    mdot_air1  = data["mdot_air1"]




    # WATER VAPOR
    steam = ct.Solution("gri30.yaml")
    steam.TPX = T_H2O, data["P_mix"], "H2O:1"
    h_vap = steam.enthalpy_mass        # J kg⁻¹ (H2O gas, 80 °C)

    # LIQUID WATER
    water = ct.Water()
    water.TP = T_H2O, data["P_mix"]
    h_liq = water.enthalpy_mass

    # explicit latent heat at this pressure
    h_lat = h_vap - h_liq
    h_H2O = h_vap - h_lat              # = h_liq   (but now shown explicitly)



    # TO GET A STRUCTURED MESH BECOUSE I AM STUPID AND BUILT EVERYTHING ON FRACTION OF AIR INSTEAD OF EQUIVALENCE RATIO
    mdot_air_stoich = air_mass_flow_for_phi(mdot_H2, 1.0, X_air=X_AIR)
    PHI_MIN, PHI_MAX = 0.25, 0.90  # range

    phi2_targets = np.linspace(PHI_MIN, PHI_MAX, GRID_NODES)
    F_SWEEP =f_for_phi2(phi2_targets,
                         mdot_air_stoich,
                         TOTAL_AIR=TOTAL_AIR,
                         mdot_air1=mdot_air1)



    #SWEEP THE VALUES OF WATER INPUTTED
    W_SWEEP = np.linspace(0.00,   0.15,  GRID_NODES)    # kg s⁻¹ water


    #PRINT THE CORES
    print(f"\n>>> map with explicit latent subtraction on {os.cpu_count()} cores …")

    #PREPARE AND PACH EACH RUN INTO A TASK
    args_common = (data, mdot_exh, h_exh, h_air2, h_H2O,
                   TOTAL_AIR, mdot_air1, X_AIR, MECH, width_lean)
    tasks = [(f, w, *args_common) for f in F_SWEEP for w in W_SWEEP]

    #INITIALIZE THE VALUES TO STORE
    phi2_vals, w_vals, Tmax_vals, NOx_vals, xo2_vals, xh2o_vals = [], [], [], [], [], []
    cp_vals, gamma_vals,  = [], []
    T_TIT_vals, Mdot_TIT = [],[]
    H_in_vals,H_out_vals,mdot_TIT_vals =[],[],[]
    Qrel_val =[]
    # RUN TASKS IN PARALLEL
    # -------------------------------------------------------------------------
    # 0.  prepare containers
    # -------------------------------------------------------------------------
    phi2_vals, w_vals, Tmax_vals, NOx_vals = [], [], [], []
    xo2_vals, xh2o_vals, gamma_vals = [], [], []
    H_in_vals, H_out_vals, Qrel_vals = [], [], []
    T_TIT_vals, mdot_TIT_vals = [], []
    with ProcessPoolExecutor() as exe:
        futures = [exe.submit(run_case, t) for t in tasks]
        for fut in tqdm(as_completed(futures), total=len(futures)):
            (p, w, tmax, nox, xo2, xh2o, gamma,
             H_in, H_out, Qrel, TIT, mdot_TIT_val) = fut.result()

            phi2_vals.append(p);
            w_vals.append(w)
            Tmax_vals.append(tmax);
            NOx_vals.append(nox)
            xo2_vals.append(xo2);
            xh2o_vals.append(xh2o)
            gamma_vals.append(gamma)
            H_in_vals.append(H_in)
            H_out_vals.append(H_out)
            Qrel_vals.append(Qrel)
            T_TIT_vals.append(TIT);
            mdot_TIT_vals.append(mdot_TIT_val)

    phi2_vals = np.asarray(phi2_vals)
    w_vals = np.asarray(w_vals)
    Tmax_vals = np.asarray(Tmax_vals)
    NOx_vals = np.asarray(NOx_vals)
    xo2_vals = np.asarray(xo2_vals)
    xh2o_vals = np.asarray(xh2o_vals)
    gamma_vals = np.asarray(gamma_vals)
    H_in_vals = np.asarray(H_in_vals)
    H_out_vals = np.asarray(H_out_vals)
    Qrel_vals = np.asarray(Qrel_vals)
    T_TIT_vals = np.asarray(T_TIT_vals)
    mdot_TIT_vals = np.asarray(mdot_TIT_vals)

    NOX_dry = False
    if NOX_dry:
        O2_dry = xo2_vals / (1.0 - xh2o_vals)
        F15 = (20.9 - 15.0) / (20.9 - 100.0 * O2_dry)
        NOx15 = NOx_vals * F15


    np.savez("lean_water_map_80C_good.npz",
             phi2=phi2_vals,
             water=w_vals,
             Tmax=Tmax_vals,
             Tmix=T_TIT_vals,  # TIT after dilution
             NOx_wet=NOx_vals,
             #NOx_15=NOx15,
             X_O2=xo2_vals,
             X_H2O=xh2o_vals,
             gamma_vals=gamma_vals,
             H_in=H_in_vals,
             H_out=H_out_vals,
             Qrel=Qrel_vals,
             mdot_TIT=mdot_TIT_vals)



    print("Map saved → lean_water_map_80C_good.npz")



    # preeliminary plots to check
    fig, ax = plt.subplots()
    sc = ax.scatter(phi2_vals, w_vals, c=NOx_vals, cmap="turbo", s=30)
    ax.set_xlabel("$\\varphi_2$"); ax.set_ylabel("$\\dot m_{H2O}$  [kg/s]")
    plt.colorbar(sc, ax=ax, label="NO$_x$  [ppm]")
    plt.title("NO$_x$ map vs φ₂ and 80 °C **liquid** water injection")
    plt.tight_layout(); plt.show()

    try:
        from scipy.interpolate import griddata
        phi_lin = np.linspace(phi2_vals.min(), phi2_vals.max(), 60)
        w_lin   = np.linspace(w_vals.min(),  w_vals.max(), 60)
        PHI, W = np.meshgrid(phi_lin, w_lin)
        T_grid = griddata((phi2_vals, w_vals), Tmax_vals,
                          (PHI, W), method="cubic")
        plt.figure(); cs = plt.contourf(PHI, W, T_grid,
                                        levels=15, cmap="plasma")
        plt.xlabel("$\\varphi_2$"); plt.ylabel("$\\dot m_{H2O}$  [kg/s]")
        plt.title("$T_{max}$ contour (80 °C liquid water)")
        plt.colorbar(cs, label="K"); plt.tight_layout(); plt.show()
    except Exception:
        pass


if __name__ == "__main__":

    mp.freeze_support() # windows is a little bitch
    mp.set_start_method("spawn", force=True) # same
    main()
