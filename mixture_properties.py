from typing import Tuple, Union, Mapping
import cantera as ct
from RQL import MECH




LHV_KEROSENE = 43.15e6   # J kg⁻¹
LHV_H2       = 120.0e6   # J kg⁻¹


"""
General utilities for mixing two arbitrary gas streams on a *mass‑flow* basis
using Cantera.  Each stream is characterised by:

    • mdot_i –mass‐flow rate [kg s⁻¹]
    • X_i     –mole‑fraction string (Cantera style) or dict

Given (mdot1, X1) and (mdot2, X2), the module returns:

    X_mix_str : str   – compact Cantera‑compatible mole‑fraction string
    mdot_tot  : float – total mass‑flow rate [kg s⁻¹]
    M_mix     : float – mean molar mass of the resulting mixture [kg mol⁻¹]

The algorithm is independent of the thermodynamic mechanism as long as all
species listed in X₁ and X₂ are present.

"""







def _to_mole_fraction_str(x: Union[str, Mapping[str, float]]) -> str:
    """Convert *x* to a Cantera mole‑fraction string."""
    if isinstance(x, str):
        return x
    if isinstance(x, Mapping):
        return ", ".join(f"{k}:{v}" for k, v in x.items())
    raise TypeError("Composition must be str or mapping.")


def mixture_properties(
    mdot1: float,
    X1: Union[str, Mapping[str, float]],
    mdot2: float,
    X2: Union[str, Mapping[str, float]],
    *,
    T: float = 298.15,
    P: float = 101_325.0,
    mech: str = MECH,
) -> Tuple[str, float, float]:
    """Return (X_mix_str, mdot_tot, M_mix) for two arbitrary streams.

    Parameters
    ----------
    mdot1, mdot2
        Mass‑flow rates [kg s⁻¹] (non‑negative).
    X1, X2
        Compositions in Cantera mole‑fraction notation or *dict* mapping
        *species → mole fraction*.  Fractions need not be normalised; Cantera
        normalises internally.
    T, P
        State at which *mean molar mass* is reported [K, Pa].
    mech
        Cantera mechanism providing thermochemical data for all species in
        *X1* ∪ *X2*.

    Notes
    -----
    The mean molar mass of a mixture is

    .. math:: \bar{M} = \sum_i X_i M_i,

    where :math:`X_i` and :math:`M_i` are the mole fraction and molar mass of
    species *i* respectively.
    """
    if mdot1 < 0 or mdot2 < 0:
        raise ValueError("Mass‑flow rates must be non‑negative.")

    gas = ct.Solution(mech)

    # ---------------------------------------------------------------------
    # Stream1
    # ---------------------------------------------------------------------
    gas.TPX = T, P, _to_mole_fraction_str(X1)
    M1 = gas.mean_molecular_weight / 1000.0  # kg mol⁻¹
    n1_total = mdot1 / M1
    n1_vec = n1_total * gas.X

    # ---------------------------------------------------------------------
    # Stream2
    # ---------------------------------------------------------------------
    gas.TPX = T, P, _to_mole_fraction_str(X2)
    M2 = gas.mean_molecular_weight / 1000.0
    n2_total = mdot2 / M2
    n2_vec = n2_total * gas.X

    # ---------------------------------------------------------------------
    # Mix streams
    # ---------------------------------------------------------------------
    n_mix = n1_vec + n2_vec
    n_tot = n_mix.sum()
    if n_tot == 0.0:
        raise ZeroDivisionError("Both mass‑flow rates are zero – undefined mixture.")

    X_mix = n_mix / n_tot
    gas.TPX = T, P, X_mix  # this sets gas.X and recomputes mean molecular weight

    X_mix_str = ", ".join(
        f"{sp}:{xi:.6g}" for sp, xi in zip(gas.species_names, X_mix) if xi > 0.0
    )

    mdot_tot = mdot1 + mdot2
    M_mix = gas.mean_molecular_weight / 1000.0  # kg mol⁻¹

    return X_mix_str, mdot_tot, M_mix

















def kerosene_to_h2(mdot_kerosene: float,
                   *,
                   lhv_kerosene: float = LHV_KEROSENE,
                   lhv_h2: float = LHV_H2) -> float:
    """
        This function inputs the amount of kerosene burnt and outputs the equivalent ammount
        of hydrogen needed for combustion using LHV.
    """
    if mdot_kerosene < 0:
        raise ValueError("Mass-flow must be non-negative.")
    return mdot_kerosene * lhv_kerosene / lhv_h2


def _as_mole_str(x: Union[str, Mapping[str, float]]) -> str:
    """
    Returns a mapping as a format of Cantera readable string library.
    """
    if isinstance(x, str):
        return x
    return ", ".join(f"{k}:{v}" for k, v in x.items())


def air_mass_flow_for_phi(
    mdot_H2: float,
    phi: float,
    *,
    X_air: Union[str, Mapping[str, float]],
    mech: str = MECH,
    T: float = 298.15,
    P: float = 101_325.0,
) -> float:
    """
    Given an equivalence ratio, and a mass flow rate of hydrogen, (and a dictionary of the mole fraction of air desired)
    this function delivers the ammount of mass flow of air necesary.
    """
    if mdot_H2 < 0:
        raise ValueError("mdot_H2 must be non-negative.")
    if phi <= 0:
        raise ValueError("phi must be positive.")

    gas = ct.Solution(mech)
    i_H2 = gas.species_index("H2")
    MW_H2 = gas.molecular_weights[i_H2] / 1000.0  # kg mol⁻¹

    n_H2 = mdot_H2 / MW_H2           # mol s⁻¹ H2
    n_O2_stoich = 0.5 * n_H2         # mol s⁻¹ O2

    MW_O2 = gas.molecular_weights[gas.species_index("O2")] / 1000.0
    mdot_O2_req = n_O2_stoich * MW_O2

    gas.TPX = T, P, _as_mole_str(X_air)
    Y_O2_air = gas.mass_fraction_dict()["O2"]

    mdot_air_stoich = mdot_O2_req / Y_O2_air
    return mdot_air_stoich / phi