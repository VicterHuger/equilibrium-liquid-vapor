"""uniquac calculations"""
from typing import List
import math
from utils.uniquac_properties import UniquacProperties
from utils.constants import GAS_CONSTANT_SI, THRESHOLD


def calculate_gibbs_excess_activity_coefficients(number_components: int, temperature: float, liquid_compositions: List[float], uniquac_params: UniquacProperties):
    """Calculate activity coefficients and gibbs excess energy

    Args:
        number_components (int): _description_
        temperature (float): _description_
        liquid_compositions (List[float]): _description_
        uniquac_params (UniquacProperties): _description_

    Returns:
        activity_coefficient (List[float])
    """
# ****************************************************************************
#
#  PROGRAM: UNIQUAC
#
#  PURPOSE:  Cálculo de G de excesso e coeficiente de atividade.
#
# ****************************************************************************

    normalized_uniquac_params: List[List[float]] = [[]]
    for i in range(number_components):
        for j in range(number_components):
            normalized_uniquac_params[i][j] = uniquac_params.uij_0[i][j] + \
                uniquac_params.uij_t[i][j]*(temperature-298.15)

    coordination_number = 10.0

    sum_q: List[float] = []
    surface_fractions: List[float] = []

    for i in range(number_components):
        for j in range(number_components):
            sum_q[i] = sum_q[i] + liquid_compositions[j]*uniquac_params.q[j]
        surface_fractions[i] = liquid_compositions[i] * \
            uniquac_params.q[i]/sum_q[i]

    sum_r: List[float] = []
    volume_fractions: List[float] = []
    for i in range(number_components):
        for j in range(number_components):
            sum_r[i] = sum_r[i] + liquid_compositions[j]*uniquac_params.r[j]
        volume_fractions[i] = liquid_compositions[i] * \
            uniquac_params.r[i]/sum_r[i]

    sum_c1 = 0.0
    sum_c2 = 0.0
    activity_coefficient_combinatorial: List[float] = []
    alngc: List[float] = []

    for i in range(number_components):
        if liquid_compositions[i] < THRESHOLD:
            continue
        sum_c1 = sum_c1 + \
            liquid_compositions[i] * \
            math.log(volume_fractions[i]/liquid_compositions[i])
        sum_c2 = sum_c2 + liquid_compositions[i]*uniquac_params.q[i]*math.log(
            volume_fractions[i]/surface_fractions[i])
        alngc[i] = math.log(volume_fractions[i]/liquid_compositions[i]) + 1.0 - \
            volume_fractions[i]/liquid_compositions[i] - (coordination_number/2.0)*uniquac_params.q[i]*(
                math.log(volume_fractions[i]/surface_fractions[i]) + 1.0 - volume_fractions[i]/surface_fractions[i])
        activity_coefficient_combinatorial[i] = math.exp(alngc[i])

    gibbs_energy_combinatorial_term = GAS_CONSTANT_SI*temperature * \
        (sum_c1 - (coordination_number/2.0)*sum_c2)

    tau: List[List[float]] = [[]]

    for i in range(number_components):
        for j in range(number_components):
            tau[i][j] = math.exp(-normalized_uniquac_params[i][j]/temperature)

    s: List[float] = []

    for i in range(number_components):
        for j in range(number_components):
            s[i] = s[i] + surface_fractions[j]*tau[i][j]

    sum_r1 = 0.0
    sum_r2: List[float] = []
    alngr: List[float] = []
    activity_coefficient_residual: List[float] = []

    for i in range(number_components):
        sum_r1 = sum_r1 + \
            liquid_compositions[i]*uniquac_params.q[i]*math.log(s[i])
        for j in range(number_components):
            sum_r2[i] = sum_r2[i] + tau[i][j]*surface_fractions[j]/s[j]
        alngr[i] = uniquac_params.q[i]*(1.0 - math.log(s[i]) - sum_r2[i])
        activity_coefficient_residual[i] = math.exp(alngr[i])
    gibbs_energy_residual_term = - GAS_CONSTANT_SI*temperature*sum_r1

    excess_gibbs_energy = gibbs_energy_combinatorial_term + gibbs_energy_residual_term
    print('[UNIQUAC - activity coefficient calc] -> Gibbs excess energy:',
          excess_gibbs_energy)

    activity_coefficient: List[float] = []
    for i in range(number_components):
        activity_coefficient[i] = activity_coefficient_combinatorial[i] * \
            activity_coefficient_residual[i]

    return activity_coefficient
