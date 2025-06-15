from global_constants import *


def calc_cost(P_fc_tot, P_fc_used, P_eps, m_fc, m_sto, m_h2_nom, m_h2):
    """
    Calculate the cost of the aircraft and its fuel over its lifetime.
    Pilot costs are not included.
    :param P_fc_tot: 100% throttle power of fuel cell [kW]
    :param P_fc_used: max power used by the fuel cell [kW]
    :param P_eps: max power used by the electric propulsion system [kW]
    :param m_fc: mass of the fuel cell [kg]
    :param m_sto: mass of the storage system [kg]
    :param m_h2_nom: nominal mass of hydrogen used [kg]
    :param m_h2: total mass of hydrogen stored [kg]

    :return: total cost of the aircraft and its fuel over its lifetime [EUR]
    """
    fc_cost = P_fc_used * FC_cost_no_bop * 3 + P_fc_tot * (FC_cost - FC_cost_no_bop) + m_fc * FC_disposal_cost + FC_maint_cost * P_fc_used
    sto_cost = m_h2 * Sto_cost + m_sto * Sto_disposal_cost
    eps_cost = P_eps * EPS_cost
    H2_cost = m_h2_nom * flight_lifetime * H2_cost

    AC_cost = AC_dev_cost + AC_purchase_cost + AC_disposal_cost + Beech_maint_cost + Insurance_cost + Crew_cost

    total_cost = fc_cost + sto_cost + eps_cost + H2_cost + AC_cost
    return total_cost

