"""
pytest unit-tests for `mixture_properties`.

Assumes the implementation lives in `mixture_two_streams.py`; change the
import line if your file is named differently.
"""


import cantera as ct
import pytest

# ──────────────────────────────────────────────────────────────────────────────
# Import the functions under test
# ──────────────────────────────────────────────────────────────────────────────
from mixture_properties import _to_mole_fraction_str

# ------------------------------------------------------------------------------
# 1.  Helper must convert dict → Cantera string exactly
# ------------------------------------------------------------------------------
def test__to_mole_fraction_str_dict():
    comp = {"O2": 0.21, "N2": 0.79}
    s = _to_mole_fraction_str(comp)
    # Order is implementation-defined, but every species must appear once
    assert "O2:0.21" in s and "N2:0.79" in s
    # A second call with a string must be returned verbatim
    s2 = _to_mole_fraction_str("H2:1.0")
    assert s2 == "H2:1.0"


# ------------------------------------------------------------------------------
# 2.  Physically consistent mixing of air + H2
# ------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("mdot_air", "mdot_h2"),
    [(1.0, 0.1), (0.5, 0.5)],  # two mass-flow combinations
)
def test_mixture_air_h2_mean_molar_mass_and_flow(mdot_air: float, mdot_h2: float):
    X_air = "N2:0.78084, O2:0.20946, AR:0.00934, CO2:0.000407"
    X_H2 = "H2:1"

    X_mix_str, mdot_tot, M_mix = mixture_properties(
        mdot_air, X_air, mdot_h2, X_H2, mech="SanDiego.yaml"
    )

    # ---------------- Expected total mass-flow -----------------
    assert mdot_tot == pytest.approx(mdot_air + mdot_h2)

    # ---------------- Expected mean molar mass -----------------
    gas = ct.Solution("SanDiego.yaml")

    # Stream 1 – air
    gas.TPX = 298.15, 101_325.0, X_air
    M_air = gas.mean_molecular_weight / 1000.0  # kg mol⁻¹
    n_air = mdot_air / M_air
    n_air_vec = n_air * gas.X

    # Stream 2 – pure H₂
    i_H2 = gas.species_index("H2")
    M_H2 = gas.molecular_weights[i_H2] / 1000.0
    n_h2 = mdot_h2 / M_H2
    n_h2_vec = [0.0] * gas.n_species
    n_h2_vec[i_H2] = n_h2

    # Mix & compute the reference mean molar mass
    n_mix = n_air_vec + n_h2_vec
    n_tot = n_mix.sum()
    M_expected = (mdot_air + mdot_h2) / n_tot  # kg mol⁻¹
    assert M_mix == pytest.approx(M_expected, rel=1e-8)

    # ---------------- Composition self-consistency --------------
    # Feed the returned string back into Cantera – should reproduce M_mix
    gas.TPX = 298.15, 101_325.0, X_mix_str
    assert gas.mean_molecular_weight / 1000.0 == pytest.approx(M_mix, rel=1e-10)


# ------------------------------------------------------------------------------
# 3.  Negative mass-flow must raise ValueError
# ------------------------------------------------------------------------------
def test_negative_mass_flow_raises():
    with pytest.raises(ValueError):
        mixture_properties(
            mdot1=-0.1,
            X1={"N2": 0.8, "O2": 0.2},
            mdot2=0.1,
            X2="H2:1",
        )
