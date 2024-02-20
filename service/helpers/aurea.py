"""File to describe otimization aurea functions"""
from utils.binary_model_param import BinaryModelParam
from utils.component_properties import ComponentProperties
from utils.mixing_rule_param import MixingRuleParam
from utils.uniquac_properties import UniquacProperties
from utils.phase import Phase
from utils.constants import THRESHOLD
from service.helpers.fugacity import fugacity_pr_calculate, fugacity_uniquac_calculate
from typing import List


# *************************************************************************
#
#  PROGRAMA: AUREA
#
#  OBJETIVO: OTIMIZACAO DE MÍNIMO - UNIVARIÁVEL PARA ORVALP
#
# *************************************************************************

def aurea_dew_pressure(liquid_compositions: List[float], vapor_compositions: List[float], temperature: float, pressure: float, newton_raphson_function_value: float, differential_newton_raphson_function_value: float, binary_model_param: BinaryModelParam, number_components: int, mixing_rule_param: MixingRuleParam, component_properties: ComponentProperties, uniquac_properties: UniquacProperties):

    # Inicialização -------------------------------------------------------
    aurea_lower_limit = 0.0
    aurea_upper_limit = 1.0
    aurea_tolerance = 1e-6
    delta = aurea_upper_limit - aurea_lower_limit
    lower_limit_variable = aurea_lower_limit + 0.382*delta
    upper_limit_variable = aurea_upper_limit - 0.382*delta
    lower_function_value = aure_dew_pressure_helper(lower_limit_variable, liquid_compositions, vapor_compositions, temperature, pressure, newton_raphson_function_value,
                                                    differential_newton_raphson_function_value, binary_model_param, number_components, mixing_rule_param, component_properties, uniquac_properties)
    upper_function_value = aure_dew_pressure_helper(upper_limit_variable, liquid_compositions, vapor_compositions, temperature, pressure, newton_raphson_function_value,
                                                    differential_newton_raphson_function_value, binary_model_param, number_components, mixing_rule_param, component_properties, uniquac_properties)

    # Iterações -----------------------------------------------------------

    while delta >= aurea_tolerance:

        if lower_function_value >= upper_function_value:
            aurea_lower_limit = lower_limit_variable
            lower_limit_variable = upper_limit_variable
            lower_function_value = upper_function_value
            delta = aurea_upper_limit - aurea_lower_limit
            upper_limit_variable = aurea_upper_limit - 0.382*delta
            upper_function_value = aure_dew_pressure_helper(upper_limit_variable, liquid_compositions, vapor_compositions, temperature, pressure, newton_raphson_function_value,
                                                            differential_newton_raphson_function_value, binary_model_param, number_components, mixing_rule_param, component_properties, uniquac_properties)
        else:
            aurea_upper_limit = upper_limit_variable
            upper_limit_variable = lower_limit_variable
            upper_function_value = lower_function_value
            delta = aurea_upper_limit - aurea_lower_limit
            lower_limit_variable = aurea_lower_limit + 0.382*delta
            lower_function_value = aure_dew_pressure_helper(lower_limit_variable, liquid_compositions, vapor_compositions, temperature, pressure, newton_raphson_function_value,
                                                            differential_newton_raphson_function_value, binary_model_param, number_components, mixing_rule_param, component_properties, uniquac_properties)

    optimized_variable_value = (lower_limit_variable+upper_limit_variable)/2.0
    return optimized_variable_value

# *************************************************************************
#
#  PROGRAMA: FUNCAO
#
#  OBJETIVO: OTIMIZAÇÃO DE MÍNIMO - FUNÇÃO AUXILIAR PARA ORVALP
#
# *************************************************************************


def aure_dew_pressure_helper(aurea_variable_value: float, liquid_compositions: List[float], vapor_compositions: List[float], temperature: float, pressure_argument: float, newton_raphson_function_value: float, differential_newton_raphson_function_value: float, binary_model_param: BinaryModelParam, number_components: int, mixing_rule_param: MixingRuleParam, component_properties: ComponentProperties, uniquac_properties: UniquacProperties):
    """Helper function to calculate function value of newton raphson for dew pressure

    Args:
        aurea_variable_value (float): _description_
        liquid_compositions (List[float]): _description_
        vapor_compositions (List[float]): _description_
        temperature (float): _description_
        pressure_argument (float): _description_
        newton_raphson_function_value (float): _description_
        differential_newton_raphson_function_value (float): _description_
        binary_model_param (BinaryModelParam): _description_
        number_components (int): _description_
        mixing_rule_param (MixingRuleParam): _description_
        component_properties (ComponentProperties): _description_
        uniquac_properties (UniquacProperties): _description_

    Returns:
        aurea_function_value:float: function value of newton raphson for dew pressure
    """
    pressure = pressure_argument - \
        (newton_raphson_function_value /
         differential_newton_raphson_function_value)*aurea_variable_value

    sum_value = 0.0
    liquid_fugacities: List[float] = [0.0] * number_components
    liquid_fugacity_coefficients: List[float] = [0.0] * number_components
    vapor_fugacities: List[float] = [0.0] * number_components
    vapor_fugacity_coefficients: List[float] = [0.0] * number_components
    fugacity_coefficients_reason: List[float] = [0.0] * number_components

    phase = Phase.LIQUID

    if binary_model_param == BinaryModelParam.PR_UNIQUAC:
        _liquid_compressability, _liquid_volume, liquid_fugacities = fugacity_uniquac_calculate(number_components, phase, pressure, temperature, liquid_compositions,
                                                                                                vapor_compositions, component_properties, mixing_rule_param, uniquac_properties)
    elif binary_model_param == BinaryModelParam.PR:
        liquid_fugacities, _liquid_compressability, _liquid_volume = fugacity_pr_calculate(
            number_components, phase, pressure, temperature, liquid_compositions, vapor_compositions, component_properties)

    # INCLUIDO CALCULO DA FASE VAPOR, # no loop do 5004 o valor de CFV(I) = 0, RAZAO(I) = Inf
    phase = Phase.VAPOR

    if binary_model_param == BinaryModelParam.PR_UNIQUAC:
        _vapor_compressability, _vapor_volume, vapor_fugacities = fugacity_uniquac_calculate(number_components, phase, pressure, temperature, liquid_compositions,
                                                                                             vapor_compositions, component_properties, mixing_rule_param, uniquac_properties)
    elif binary_model_param == BinaryModelParam.PR:
        vapor_fugacities, _vapor_compressability, _vapor_volume = fugacity_pr_calculate(
            number_components, phase, pressure, temperature, liquid_compositions, vapor_compositions, component_properties)

    for i in range(number_components):
        y = vapor_compositions[i]

        if abs(y) > THRESHOLD:
            liquid_fugacity_coefficients[i] = liquid_fugacities[i] / \
                (liquid_compositions[i] * pressure)
            vapor_fugacity_coefficients[i] = vapor_fugacities[i] / \
                (y * pressure)
            fugacity_coefficients_reason[i] = liquid_fugacity_coefficients[i] / \
                vapor_fugacity_coefficients[i]
            sum_value = sum_value + y/fugacity_coefficients_reason[i]

    aurea_function_value = abs(sum_value-1.0)  # TESTAR AO QUADRADO TBM
    # F2 = (SOMA-1.0)  # AQUI o valor da funcao deve ser o valor ORIGINAL
    return aurea_function_value
