""" Calculate any parameter needed for Peng Robinson Calculations"""
import math
from utils.alpha_model import AlphaModel
from utils.constants import GAS_CONSTANT, THRESHOLD, PI_CONSTANT
from typing import List


def calculate_atrative_repulsive_params(temperature: float, acentric_factor: float, critical_temperature: float, critical_pressure: float, alpha_model: AlphaModel, aat1: float | None, aat2: float | None, aat3: float | None):
    """Function to calculate atrative (a) and repulsive (b) Peng-Robinson params for a component
        PENG E ROBINSON 1978 / ALMEIDA ET AL. 1991

    Args:
        component_index (int): Index of the component in the system
        temperature (float): Temperature of the sytem
        acentric_factor (float): Acentric factor of the component
        critical_temperature (float): Critical temperature of the component
        critical_pressure (float): Critical pressure of the component
        alpha_model (AlphaMOdel): Alpha Model Chosen
        aat1 (float | None): ALMEIDA ET AL. 1991 param
        aat2 (float | None): ALMEIDA ET AL. 1991 param
        aat3 (float | None): ALMEIDA ET AL. 1991 param
    """
    ohm_a = 0.457235528921
    ohm_b = 0.0777960739039

    alpha: float = 0.0

    if alpha_model == AlphaModel.PR87:
        so1: float = 0.0
        w = acentric_factor
        if w <= 0.491:
            so1 = 0.37464 + 1.54226*w - 0.26992*w**2
        else:
            so1 = 0.379642 + 1.487503*w - 0.164423*w**2 + 0.016666*w**3
        alpha = (
            1.0+so1*(1.0-math.sqrt(temperature/critical_temperature)))**2
    elif alpha_model == AlphaModel.AAT:
        alpha = math.exp(aat1*(1-temperature/critical_temperature)*(abs(1-temperature/critical_temperature))**(
            aat2-1)+aat3*(critical_temperature/temperature-1)) if temperature != critical_pressure else 1.0

    attractive_param = ohm_a * alpha * \
        (GAS_CONSTANT**2)*(critical_temperature**2)/critical_pressure
    repulsive_param = ohm_b * GAS_CONSTANT * \
        critical_temperature / critical_pressure

    return attractive_param, repulsive_param


def solve_pr_cubic_equation(a10: float, a11: float, a12: float):
    """solve pr cubic equation

    Args:
        a10 (float): _description_
        a11 (float): _description_
        a12 (float): _description_

    Returns:
        equation_roots: List[float]: _description_
        nr: int
    """
# ****************************************************************************
#
#     SUB-ROTINA: CUBICA
#
#     OBJETIVO: Solução da equação cúbica.
#
# ****************************************************************************
    nr: int = 0
    te = 0.0
    ONE_THIRD = 1.0/3.0
    equation_roots: List[float] = [0.0]*3
    s = 0.0
    z1 = 0.0
    z2 = 0.0
    vacos = 0.0

    q = (3.0 * a11 - (a10 ** 2)) / 9.0
    r = (9.0 * a10 * a11 - 27.0 * a12 - 2.0 * a10 ** 3) / 54.0
    d = q ** 3 + r ** 2
    if abs(d) <= 0.00000001:
        d = 0

    if d == 0:
        nr = 0
        if r >= 0.0:
            s = abs(r) ** ONE_THIRD
        else:
            s = -abs(r) ** ONE_THIRD
        z1 = 2.0 * s - a10 / 3.0
        if z1 <= 0.00000001:
            z1 = 0.0
        z2 = -s - a10 / 3.0
        if z2 <= 0.00000001:
            z2 = 0.0
        equation_roots[0] = max(z1, z2)
        equation_roots[2] = min(z1, z2)
        equation_roots[1] = equation_roots[2]
        return equation_roots, nr

    elif d > 0.0:
        nr = 1

        u = r + math.sqrt(d)
        if u >= 0.0:
            s = abs(u) ** ONE_THIRD
        else:
            s = -abs(u) ** ONE_THIRD
        if abs(s) <= 0.00000001:
            s = 0.0

        u = r - math.sqrt(d)
        if (u >= 0.0):
            te = abs(u) ** ONE_THIRD
        else:
            te = -abs(u) ** ONE_THIRD
        if abs(te) <= 0.00000001:
            te = 0.0

        equation_roots[0] = s + te - a10 / 3.0
        equation_roots[1] = 0.0
        equation_roots[2] = 0.0
        return equation_roots, nr
    elif d < 0.0:
        nr = -1
        u = 2.0 * math.sqrt(-q)
        vacos = r / math.sqrt(-q ** 3)
        if vacos > 1.0:
            vacos = 1.0
        theta = math.acos(vacos) / 3.0
        if abs(theta) <= 0.00000001:
            theta = 0.0
        z1 = u * math.cos(theta) - a10 / 3.0
        if z1 <= 0.00000001:
            z1 = 0.0
        z2 = u * math.cos(theta + 2.0 * PI_CONSTANT / 3.0) - a10 / 3.0
        if z2 <= 0.00000001:
            z2 = 0.0
        z3 = u * math.cos(theta + 4.0 * PI_CONSTANT / 3.0) - a10 / 3.0
        if z3 <= 0.00000001:
            z3 = 0.0
        equation_roots[0] = max(z1, max(z2, z3))
        equation_roots[2] = min(z1, min(z2, z3))
        if equation_roots[0] == z1:
            equation_roots[1] = max(z2, z3)
        if equation_roots[0] == z2:
            equation_roots[1] = max(z1, z3)
        if equation_roots[0] == z3:
            equation_roots[1] = max(z1, z2)
        return equation_roots, nr
    return equation_roots, nr

    # compressability_factor: float = 1.0
    # count: int = 0
    # z0: float = 0.0
    # z2: float = 0.0
    # z3: float = 0.0

    # while abs(compressability_factor - z0) > THRESHOLD and count < 100000:
    #     z0 = compressability_factor
    #     compressability_factor = compressability_factor - (compressability_factor**3 + a12*(
    #         compressability_factor**2) + a11*compressability_factor + a10)/(3*(compressability_factor**2) + 2*a12*compressability_factor + a11)
    #     count = count + 1

    # a22 = a12 + compressability_factor
    # a21 = (compressability_factor**2) + a12*compressability_factor + a11
    # delta = (a22**2)-4*a21

    # if delta >= 0.0:
    #     z2 = (-a22-math.sqrt(delta))/2
    #     z3 = (-a22+math.sqrt(delta))/2
    # else:
    #     z2 = compressability_factor
    #     z3 = compressability_factor
    #     print('solve_pr_cubic_equation: delta lower than zero!')

    # z1 = compressability_factor

    # return z1, z2, z3
