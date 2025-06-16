from global_constants import *


def calc_cost(P_FC, P_eps, m_sto, m_fc, m_h2):
    """
    Calculate the cost of the aircraft and its fuel over its lifetime.
    Pilot costs are not included.
    :param fc_cost: cost of the fuel cell over lifetime [EUR]
    :param P_eps: max power used by the electric propulsion system [kW]
    :param m_fc: mass of the fuel cell [kg]
    :param m_sto: mass of the storage system [kg]
    :param m_h2_nom: nominal mass of hydrogen used [kg]
    :param m_h2: total mass of hydrogen stored [kg]

    :return: total cost of the aircraft and its fuel over its lifetime [EUR]
    """
    # Initial costs
    fc_init = FC_cost * P_FC
    sto_init = Sto_cost * m_h2
    eps_init = P_eps * EPS_cost

    purchase_cost = fc_init + sto_init + eps_init + AC_purchase_cost + AC_dev_cost

    # Operational costs  per year
    fc_maint = FC_maint_cost * P_FC / years_of_life
    sto_maint = Sto_maint_cost * m_sto / years_of_life
    landing =  7.5 * Landing_tax / years_of_life
    crew = Crew_cost / years_of_life
    fuel_cost = H2_cost * m_h2 * 0.72  # 28% is reserve fuel
    insurance = Insurance_cost / years_of_life

    operational_costs = fc_maint + sto_maint + landing + crew + fuel_cost + insurance + AC_maint_cost

    # Disposal costs
    fc_disposal = FC_disposal_cost * m_fc 
    sto_disposal = Sto_disposal_cost * m_sto
    disposal_cost =  fc_disposal + sto_disposal + AC_disposal_cost
    
    # Total cost
    total_cost = purchase_cost + operational_costs * years_of_life + disposal_cost
    return total_cost

