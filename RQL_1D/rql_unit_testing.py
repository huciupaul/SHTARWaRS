import pytest
import cantera as ct
from RQL_1D.mixture_properties import _to_mole_fraction_str, mixture_properties, kerosene_to_h2, _as_mole_str, air_mass_flow_for_phi, mix_streams_const_HP

# Use a simple Cantera mechanism for testing
TEST_MECH = "gri30.yaml"


#---------_to_mole_fraction_str tests-------------------
def test_to_mole_fraction_str_multiple_species():
    X = {"H2": 0.5, "O2": 0.5}
    result = _to_mole_fraction_str(X)
    assert result in ["H2:0.5, O2:0.5", "O2:0.5, H2:0.5"]

def test_to_mole_fraction_str_empty():
    X = {}
    result = _to_mole_fraction_str(X)
    assert result == ""

def test_to_mole_fraction_str_invalid():
    with pytest.raises(ValueError):
        _to_mole_fraction_str(-1.0)


#---------_as_mole_str tests-------------------
def test_mixture_properties_basic():
    mdot1 = 2.0  # kg/s
    X1 = {"H2": 1.0}
    mdot2 = 8.0  # kg/s
    X2 = {"O2": 1.0}
    X_mix_str, mdot_tot, M_mix = mixture_properties(
        mdot1, X1, mdot2, X2, mech=TEST_MECH
    )
    assert isinstance(X_mix_str, str)
    assert pytest.approx(mdot_tot) == 10.0
    # Check that the mixture contains both H2 and O2
    assert "H2" in X_mix_str and "O2" in X_mix_str
    # Check mean molar mass is between H2 and O2
    gas = ct.Solution(TEST_MECH)
    M_H2 = gas.molecular_weights[gas.species_index("H2")] / 1000.0
    M_O2 = gas.molecular_weights[gas.species_index("O2")] / 1000.0
    assert M_H2 < M_mix < M_O2

def test_mixture_properties_zero_massflow():
    # Both streams zero
    with pytest.raises(ZeroDivisionError):
        mixture_properties(0.0, {"H2": 1.0}, 0.0, {"O2": 1.0}, mech=TEST_MECH)

def test_mixture_properties_negative_massflow():
    # Negative mass flow should raise ValueError
    with pytest.raises(ValueError):
        mixture_properties(-1.0, {"H2": 1.0}, 1.0, {"O2": 1.0}, mech=TEST_MECH)
    with pytest.raises(ValueError):
        mixture_properties(1.0, {"H2": 1.0}, -1.0, {"O2": 1.0}, mech=TEST_MECH)


#---------kerosene_to_h2 tests-------------------
def test_kerosene_to_h2_typical():
    # Typical conversion: 1 kmol kerosene to H2
    n_kerosene = 1.0
    h2_moles = kerosene_to_h2(n_kerosene)
    assert h2_moles > 0
    assert isinstance(h2_moles, float)

def test_kerosene_to_h2_zero():
    # Zero kerosene should yield zero H2
    assert kerosene_to_h2(0.0) == 0.0

def test_kerosene_to_h2_negative():
    # Negative input should raise ValueError
    with pytest.raises(ValueError):
        kerosene_to_h2(-1.0)


#---------_as_mole_str tests-------------------
def test_as_mole_str_multiple_species():
    X = {"H2": 0.5, "O2": 0.5}
    result = _as_mole_str(X)
    # Accept either order
    assert result in ["H2:0.5, O2:0.5", "O2:0.5, H2:0.5"]

def test_as_mole_str_empty():
    X = {}
    result = _as_mole_str(X)
    assert result == ""

def test_as_mole_str_invalid_type():
    with pytest.raises(TypeError):
        _as_mole_str("not a dict")

def test_as_mole_str_invalid_value():
    X = {"H2": -1.0}
    with pytest.raises(ValueError):
        _as_mole_str(X)


#---------air_mass_flow_for_phi tests-------------------
def test_air_mass_flow_for_phi_typical():
    # Typical values for fuel mass flow and equivalence ratio
    mdot_fuel = 1.0  # kg/s
    phi = 0.5
    mech = "gri30.yaml"
    X_fuel = {"H2": 1.0}
    X_air = {"O2": 0.21, "N2": 0.79}
    mdot_air = air_mass_flow_for_phi(mdot_fuel, phi, X_fuel, X_air, mech)
    assert mdot_air > 0
    assert isinstance(mdot_air, float)

def test_air_mass_flow_for_phi_phi_one():
    # phi = 1.0 (stoichiometric)
    mdot_fuel = 1.0
    phi = 1.0
    mech = "gri30.yaml"
    X_fuel = {"H2": 1.0}
    X_air = {"O2": 0.21, "N2": 0.79}
    mdot_air = air_mass_flow_for_phi(mdot_fuel, phi, X_fuel, X_air, mech)
    assert mdot_air > 0
    assert isinstance(mdot_air, float)

def test_air_mass_flow_for_phi_zero_fuel():
    # Zero fuel mass flow should yield zero air mass flow
    mdot_fuel = 0.0
    phi = 0.8
    mech = "gri30.yaml"
    X_fuel = {"H2": 1.0}
    X_air = {"O2": 0.21, "N2": 0.79}
    mdot_air = air_mass_flow_for_phi(mdot_fuel, phi, X_fuel, X_air, mech)
    assert mdot_air == 0.0

def test_air_mass_flow_for_phi_invalid_phi():
    # Negative phi should raise ValueError
    mdot_fuel = 1.0
    phi = -0.5
    mech = "gri30.yaml"
    X_fuel = {"H2": 1.0}
    X_air = {"O2": 0.21, "N2": 0.79}
    with pytest.raises(ValueError):
        air_mass_flow_for_phi(mdot_fuel, phi, X_fuel, X_air, mech)

def test_air_mass_flow_for_phi_invalid_mdot_fuel():
    # Negative fuel mass flow should raise ValueError
    mdot_fuel = -1.0
    phi = 0.8
    mech = "gri30.yaml"
    X_fuel = {"H2": 1.0}
    X_air = {"O2": 0.21, "N2": 0.79}
    with pytest.raises(ValueError):
        air_mass_flow_for_phi(mdot_fuel, phi, X_fuel, X_air, mech)

#---------mix_streams_const_HP tests-------------------
class DummyStream:
    """A minimal dummy stream object for testing mix_streams_const_HP."""
    def __init__(self, gas, mass):
        self.mass = mass
        self.enthalpy = gas.enthalpy_mass * mass
        self.X = gas.X
        self.T = gas.T
        self.phase = gas
        self.phase.source = gas.name  # For ct.Solution(gas.source)
        
@pytest.fixture
def gases():
    # Use Cantera's gri30 mechanism for test
    gas1 = ct.Solution("gri30.yaml")
    gas2 = ct.Solution("gri30.yaml")
    gas3 = ct.Solution("gri30.yaml")
    gas1.TPX = 300, 101325, "H2:1.0"
    gas2.TPX = 400, 101325, "O2:1.0"
    gas3.TPX = 350, 101325, "N2:1.0"
    return gas1, gas2, gas3

def test_mix_streams_const_HP_basic(gases):
    gas1, gas2, gas3 = gases
    stream1 = DummyStream(gas1, 1.0)
    stream2 = DummyStream(gas2, 2.0)
    stream3 = DummyStream(gas3, 3.0)
    P = 101325

    X_mix, T_mix, mdot_mix = mix_streams_const_HP(stream1, stream2, stream3, P)

    # Check mass flow
    assert pytest.approx(mdot_mix, rel=1e-12) == 6.0

    # Check composition: should sum to 1
    assert pytest.approx(sum(X_mix), rel=1e-12) == 1.0

    # Check temperature is within the range of input temperatures
    assert min(stream1.T, stream2.T, stream3.T) <= T_mix <= max(stream1.T, stream2.T, stream3.T)

def test_mix_streams_const_HP_zero_mass(gases):
    gas1, gas2, gas3 = gases
    stream1 = DummyStream(gas1, 0.0)
    stream2 = DummyStream(gas2, 0.0)
    stream3 = DummyStream(gas3, 0.0)
    P = 101325

    with pytest.raises(ZeroDivisionError):
        mix_streams_const_HP(stream1, stream2, stream3, P)

def test_mix_streams_const_HP_single_stream(gases):
    gas1, _, _ = gases
    stream1 = DummyStream(gas1, 1.0)
    # Streams 2 and 3 are zero mass
    class ZeroStream:
        mass = 0.0
        enthalpy = 0.0
        X = gas1.X
        T = gas1.T
        phase = gas1
        phase.source = gas1.name

    stream2 = ZeroStream()
    stream3 = ZeroStream()
    P = 101325

    X_mix, T_mix, mdot_mix = mix_streams_const_HP(stream1, stream2, stream3, P)
    assert pytest.approx(mdot_mix, rel=1e-12) == 1.0
    assert pytest.approx(sum(X_mix), rel=1e-12) == 1.0
    assert pytest.approx(T_mix, rel=1e-6) == stream1.T