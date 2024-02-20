"""Fugacity calculations by different methods"""
import math
from typing import List
from service.helpers.mixing_rules import calculate_heidemann_kokal_elv, calculate_activity_coefficient_fugchk
from service.helpers.pr_parameters import calculate_atrative_repulsive_params, solve_pr_cubic_equation
from service.helpers.uniquac import calculate_gibbs_excess_activity_coefficients
from utils import phase, component_properties, uniquac_properties, constants, mixing_rule_param


def fugacity_pr_calculate(number_components: int, system_phase: phase.Phase, pressure: float, temperature: float, liquid_compositions: List[float], vapor_compositions: List[float], component_properties_argument: component_properties.ComponentProperties):
    """Calculate fugacity phase using Peng Robinson EoS - PHASES

    Args:
        number_components (int): number of componentes of a system
        system_phase (Phase): phase that fugacity will be calculated
        pressure (float): system pressure
        temperature (float): system temperature
        liquid_compositions (List[float]): liquid composition of system components
        vapor_compositions (List[float]): vapor composition of system components
        component_properties: ComponentProperties
    """
    # **************************************************************************
    #
    #     SUB-ROTINA: PHASES
    #
    #     OBJETIVO: Cálculo da fugacidade das fases via PR EoS
    #
    # **************************************************************************

    atractive_params: List[float] = []
    repulsive_params: List[float] = []
    system_compositions: List[float] = []
    binary_params: List[List[float]] = []
    mixture_atrative_param: float = 0.0
    mixture_repulsive_param: float = 0.0

    for i in range(number_components):
        tc = component_properties_argument.critical_acentric_properties.critical_temperatures[
            i]
        pc = component_properties_argument.critical_acentric_properties.critical_pressures[i]
        w = component_properties_argument.critical_acentric_properties.acentric_factors[i]
        aat1: float | None = None
        aat2: float | None = None
        aat3: float | None = None

        alpha_chosen = component_properties_argument.alpha_model
        if alpha_chosen == component_properties.AlphaModel.AAT:
            aat1 = component_properties_argument.aat_params.aat1[i]
            aat2 = component_properties_argument.aat_params.aat2[i]
            aat3 = component_properties_argument.aat_params.aat3[i]

        atractive_param, repulsive_param = calculate_atrative_repulsive_params(
            temperature, w, tc, pc, alpha_chosen, aat1, aat2, aat3)
        atractive_params.append(atractive_param)
        repulsive_params.append(repulsive_param)

    if system_phase == phase.Phase.LIQUID:
        system_compositions = liquid_compositions
    elif system_phase == phase.Phase.VAPOR:
        system_compositions = vapor_compositions

    for i in range(number_components):
        binary_params.append([])
        for j in range(number_components):
            peng_robinson_binary_param = component_properties_argument.peng_robinson_binary_params[
                i][j]
            binary_params[i].append((
                1.0-peng_robinson_binary_param)*(math.sqrt(atractive_params[i]*atractive_params[j])))
            mixture_atrative_param = mixture_atrative_param + \
                system_compositions[i] * \
                system_compositions[j]*binary_params[i][j]
        mixture_repulsive_param = mixture_repulsive_param + \
            system_compositions[i] * repulsive_params[i]

    a2 = pressure*mixture_atrative_param / \
        ((constants.GAS_CONSTANT*temperature)**2)
    b2 = pressure*mixture_repulsive_param/(constants.GAS_CONSTANT*temperature)
    sig = 1.0 + math.sqrt(2)
    eps = 1.0 - math.sqrt(2)
    r1 = -sig
    r2 = -eps
    a10 = - (mixture_repulsive_param*(r1+r2+1.0) +
             constants.GAS_CONSTANT*temperature/pressure)
    a11 = (mixture_repulsive_param**2)*(r1*r2+r1+r2) + constants.GAS_CONSTANT*temperature * \
        mixture_repulsive_param*(r1+r2)/pressure + \
        mixture_atrative_param/pressure
    a12 = - mixture_repulsive_param*(r1*r2*(mixture_repulsive_param**2) + r1*r2*mixture_repulsive_param *
                                     constants.GAS_CONSTANT*temperature/pressure + mixture_atrative_param/pressure)

    compressability_factor_list, nr = solve_pr_cubic_equation(a10, a11, a12)

    volume = 0.0
    if nr <= 0:
        if system_phase == phase.Phase.LIQUID:
            volume = compressability_factor_list[2]
        if system_phase == phase.Phase.VAPOR:
            volume = compressability_factor_list[0]
    else:
        volume = compressability_factor_list[0]

    compressability_factor = volume / \
        (constants.GAS_CONSTANT*temperature/pressure)
    aas: List[float] = []
    component_fugacities: List[float] = []

    for i in range(number_components):
        aas.append(0.0)
        for j in range(number_components):
            aas[i] = aas[i] + system_compositions[j]*binary_params[i][j]
        component_fugacities.append(system_compositions[i]*pressure*math.exp((repulsive_params[i]/mixture_repulsive_param)*(compressability_factor-1.0)-math.log(compressability_factor-b2)-(
            a2/b2)*(1.0/(sig-eps))*(2.0*aas[i]/mixture_atrative_param - repulsive_params[i]/mixture_repulsive_param)*math.log((compressability_factor+sig*b2)/(compressability_factor+eps*b2))))

    return component_fugacities, compressability_factor, volume


def fugacity_uniquac_calculate(number_molecular_components: int, system_phase: phase.Phase, pressure: float, temperature: float, liquid_compositions: List[float], vapor_compositions: List[float], component_properties_argument: component_properties.ComponentProperties, mixing_rule_param_argument: mixing_rule_param.MixingRuleParam, uniquac_params: uniquac_properties.UniquacProperties):
    """Fugacity calculation with UNIQUAC

    Args:
        number_molecular_components (int): _description_
        system_phase (phase.Phase): _description_
        pressure (float): _description_
        temperature (float): _description_
        liquid_compositions (List[float]): _description_
        vapor_compositions (List[float]): _description_
        component_properties_argument (component_properties.ComponentProperties): _description_
        mixing_rule_param_argument (Mixing Rule Param): _description_
        uniquac_params (uniquac_properties.UniquacProperties): _description_

    Returns:
        compressability_factor, volume, fugacities (tuple(float, float, List[float])): _description_
    """
    # ****************************************************************************
    #
    #  SUBROUTINE: PR0011
    #
    #  PURPOSE:  Cálculo da fugacidade pelo modelo PR-UNIQUAC-HK
    #
    # ****************************************************************************
    # VARIÁVEIS PRINCIPAIS --------------------------------------------------
    # NC    ------ Número de componentes (espécies) no sistema
    # NM    ------ Número de espécies moleculares no sistema
    # GAMA(I) ---- Coeficiente de atividade do componente I
    # T     ------ Temperatura (K)
    # P     ------ Pressão (bar)
    # ZI(I) ------ Fração molar aparente do componente I na fase p
    # ZRI(I) ----- Fração molar real do componente I na fase p
    # Fug_p(I) --- Fugacidade do componente I na fase p

    # VÁRIÁVEIS AUXILIARES --------------------------------------------------
    # PHA

    # INICIALIZAÇÃO ---------------------------------------------------------
    activity_coefficients: List[float] = []
    number_components = number_molecular_components
    apparent_mole_fractions: List[float] = []
    real_mole_fractions: List[float] = []
    atractive_params: List[float] = []
    repulsive_params: List[float] = []

    apparent_mole_fractions = liquid_compositions if system_phase == phase.Phase.LIQUID else vapor_compositions

    real_mole_fractions = apparent_mole_fractions

    # CÁLCULO DOS COEFICIENTES DE ATIVIDADE ---------------------------------
    if mixing_rule_param_argument == mixing_rule_param.MixingRuleParam.HK and system_phase == phase.Phase.LIQUID:
        new_activity_coefficients = calculate_activity_coefficients(
            number_components, number_molecular_components, temperature, real_mole_fractions, uniquac_params)
        if new_activity_coefficients:
            activity_coefficients = new_activity_coefficients

    # CÁLCULO DOS PARÂMETROS ENERGÉTICOS E CO-VOLUMES DAS ESPÉCIES PURAS ----
    atractive_params, repulsive_params = [calculate_atrative_repulsive_params(
        temperature, component_properties_argument.critical_acentric_properties.acentric_factors[
            i],
        component_properties_argument.critical_acentric_properties.critical_temperatures[i],
        component_properties_argument.critical_acentric_properties.critical_pressures[i],
        component_properties_argument.alpha_model, component_properties_argument.aat_params.aat1[
            i],
        component_properties_argument.aat_params.aat2[i], component_properties_argument.aat_params.aat3[i])
        for i in range(number_molecular_components)]

    # for i in range(number_molecular_components):
    #     acentric_factor = component_properties_argument.critical_acentric_properties.acentric_factors[
    #         i]
    #     critical_temperature = component_properties_argument.critical_acentric_properties.critical_temperatures[
    #         i]
    #     critical_pressure = component_properties_argument.critical_acentric_properties.critical_pressures[
    #         i]
    #     alpha_model = component_properties_argument.alpha_model
    #     aat1 = component_properties_argument.aat_params.aat1[i]
    #     aat2 = component_properties_argument.aat_params.aat2[i]
    #     aat3 = component_properties_argument.aat_params.aat3[i]
    #     atractive_params[i], repulsive_params[i] = calculate_atrative_repulsive_params(
    #         temperature, acentric_factor, critical_temperature, critical_pressure, alpha_model, aat1, aat2, aat3)

    # CÁLCULO DA FUGACIDADE DAS ESPÉCIES MOLECULARES -------------------------
    volume = 0.0
    fugacities: List[float] = []
    if mixing_rule_param_argument == mixing_rule_param.MixingRuleParam.HK:
        volume, fugacities = calculate_heidemann_kokal_elv(number_components, number_molecular_components, system_phase,
                                                           temperature, pressure, liquid_compositions, atractive_params, repulsive_params, activity_coefficients)
    elif mixing_rule_param_argument == mixing_rule_param.MixingRuleParam.LCVM:
        volume, fugacities = calculate_activity_coefficient_fugchk(system_phase, number_molecular_components, temperature,
                                                                   pressure, real_mole_fractions, atractive_params, repulsive_params, uniquac_params)
        fugacities = [fugacities[i]*apparent_mole_fractions[i]
                      * pressure for i in range(number_components)]
        # for i in range(number_components):
        #     fugacities[i] = fugacities[i]*apparent_mole_fractions[i]*pressure

    compressability_factor = volume / \
        (constants.GAS_CONSTANT*temperature/pressure)

    return compressability_factor, volume, fugacities


def calculate_activity_coefficients(number_components: int, number_molecular_components: int, temperature: float, real_mole_fractions: List[float], uniquac_params: uniquac_properties.UniquacProperties):
    """Calculate acitivity coefficients or return None

    Args:
        number_components (int): _description_
        number_molecular_components (int): _description_
        temperature (float): _description_
        real_mole_fractions (List[float]): _description_
        uniquac_params (UniquacProperties): _description_

    Returns:
        list | none: updated gama or None
    """

    for i in range(number_molecular_components):
        if abs(real_mole_fractions[i] - 1.0) < constants.THRESHOLD:
            return

    gama = calculate_gibbs_excess_activity_coefficients(
        number_components, temperature, real_mole_fractions, uniquac_params)

    gama = [1.0 if abs(real_mole_fractions[i]) < constants.THRESHOLD else gama[i]
            for i in range(number_components)]

    return gama
