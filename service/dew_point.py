"""Service for dew point operations"""

from typing import List
from service.helpers.fugacity import fugacity_pr_calculate, fugacity_uniquac_calculate
from service.helpers.wilson import wilson_correlation
from service.helpers.aurea import aurea_dew_pressure
from utils.phase import Phase
from utils.component_properties import ComponentProperties
from utils.constants import THRESHOLD
from utils.uniquac_properties import UniquacProperties
from utils.exceptions import ConvergenceError
from utils.binary_model_param import BinaryModelParam
from utils.mixing_rule_param import MixingRuleParam


def calc_dew_pressure(number_components: int, vapor_compositions: List[float],
                      temperature: float, binary_model_param: BinaryModelParam, mixing_rule_param: MixingRuleParam, component_properties: ComponentProperties, uniquac_properties: UniquacProperties):
    """Calculate dew pressure for the condition provided

    Args:
        number_components (int): number of components in the system
        vapor_compositions (List[float]):  components vapor composition list
        temperature (float): temperature of the system
        binary_model_param (BinaryModelParam): BinaryModelParam where can be PR or PR_UNIQUAC
        mixing_rule_param (MixingRuleParam): MixingRuleParam where can be HK or LCVM
        component_properties (ComponentProperties): _description_
        uniquac_properties (UniquacProperties): _description_

    Raises:
        ConvergenceError: When count of newton raphson reaches 1000 

    Returns:
        pressure(float): Dew pressure for the scenario provided
    """
    # ****************************************************************************************************************************************************
    #
    #  SUB-ROTINA: ORVALP
    #
    #  OBJETIVO: Cálculo do ponto de orvalho (T e Y conhecidos)
    #
    # ****************************************************************************************************************************************************

    # **************************************************************************
    # Initial Estimations
    sum_list = [0.0]
    sum1_list = [0.0]
    pressure = 300.0  # bar
    count = 1

    liquid_fugacities: List[float] = [0.0] * number_components
    liquid_step_fugacities: List[float] = [0.0] * number_components
    liquid_fugacity_coefficients: List[float] = [0.0] * number_components
    liquid_step_fugacity_coefficients: List[float] = [0.0] * number_components
    vapor_fugacities: List[float] = [0.0] * number_components
    vapor_step_fugacities: List[float] = [0.0] * number_components
    vapor_fugacity_coefficients: List[float] = [0.0] * number_components
    vapor_step_fugacity_coefficients: List[float] = [0.0] * number_components
    fugacity_coefficients_reason: List[float] = [0.0] * number_components
    step_fugacity_coefficient_reason: List[float] = [0.0] * number_components

    liquid_compositions = wilson_correlation(number_components=number_components, pressure=pressure,
                                             temperature=temperature, vapor_compositions=vapor_compositions, component_properties=component_properties)

    newton_raphson_function_value = 1.0

    # Newton Raphson for Buble Pressure
    while abs(newton_raphson_function_value) > 1e-6:
        if count == 10000:
            raise ConvergenceError(
                f"Dew Pressure has not converged - count reached {count} units")
        pressure_step = pressure + 0.0001*pressure
        temperature_step = temperature

        phase = Phase.LIQUID

        if binary_model_param == BinaryModelParam.PR_UNIQUAC:
            _liquid_compressability, _liquid_volume, liquid_fugacities = fugacity_uniquac_calculate(number_components, phase, pressure, temperature, liquid_compositions,
                                                                                                    vapor_compositions, component_properties, mixing_rule_param, uniquac_properties)
        elif binary_model_param == BinaryModelParam.PR:
            liquid_fugacities, _liquid_compressability, _liquid_volume = fugacity_pr_calculate(
                number_components, phase, pressure, temperature, liquid_compositions, vapor_compositions, component_properties)

        phase = Phase.VAPOR

        if binary_model_param == BinaryModelParam.PR_UNIQUAC:
            _vapor_compressability, _vapor_volume, vapor_fugacities = fugacity_uniquac_calculate(number_components, phase, pressure, temperature, liquid_compositions,
                                                                                                 vapor_compositions, component_properties, mixing_rule_param, uniquac_properties)
            _vapor_compressability_step, _vapor_volume_step, vapor_step_fugacities = fugacity_uniquac_calculate(number_components, phase, pressure_step, temperature_step, liquid_compositions,
                                                                                                                vapor_compositions, component_properties, mixing_rule_param, uniquac_properties)
        elif binary_model_param == BinaryModelParam.PR:
            vapor_fugacities, _vapor_compressability, _vapor_volume = fugacity_pr_calculate(
                number_components, phase, pressure, temperature, liquid_compositions, vapor_compositions, component_properties)
            vapor_step_fugacities, _vapor_compressability_step, _vapor_volume_step = fugacity_pr_calculate(
                number_components, phase, pressure_step, temperature_step, liquid_compositions, vapor_compositions, component_properties)

        sum_list.append(0.0)

        for i in range(number_components):
            y = vapor_compositions[i]
            if abs(y) > THRESHOLD:
                liquid_fugacity_coefficients[i] = liquid_fugacities[i] / \
                    (liquid_compositions[i] * pressure)
                vapor_fugacity_coefficients[i] = vapor_fugacities[i] / \
                    (y * pressure)
                vapor_step_fugacity_coefficients[i] = vapor_step_fugacities[i]/(
                    y*pressure_step)
                fugacity_coefficients_reason[i] = liquid_fugacity_coefficients[i] / \
                    vapor_fugacity_coefficients[i]
                sum_list[count-1] = sum_list[count-1] + \
                    y/fugacity_coefficients_reason[i]

        test_1 = 1.0
        while test_1 >= 1.0e-5:
            for i in range(number_components):
                if vapor_compositions[i] > THRESHOLD:
                    liquid_compositions[i] = vapor_compositions[i] / \
                        fugacity_coefficients_reason[i] / sum_list[count-1]

            if count == 10000:
                raise ConvergenceError(
                    f'Dew Pressure has not converged - count === {count}')
            count = count + 1
            sum_list.append(0.0)
            sum1_list.append(0.0)

            phase = Phase.LIQUID
            if binary_model_param == BinaryModelParam.PR_UNIQUAC:
                _liquid_compressability, _liquid_volume, liquid_fugacities = fugacity_uniquac_calculate(number_components, phase, pressure, temperature, liquid_compositions,
                                                                                                        vapor_compositions, component_properties, mixing_rule_param, uniquac_properties)
                _liquid_compressability_step, _liquid_volume_step, liquid_step_fugacities = fugacity_uniquac_calculate(
                    number_components, phase, pressure_step, temperature_step, liquid_compositions, vapor_compositions, component_properties, mixing_rule_param, uniquac_properties)
            elif binary_model_param == BinaryModelParam.PR:
                liquid_fugacities, _liquid_compressability, _liquid_volume = fugacity_pr_calculate(
                    number_components, phase, pressure, temperature, liquid_compositions, vapor_compositions, component_properties)
                liquid_step_fugacities, _liquid_compressability, _liquid_volume = fugacity_pr_calculate(
                    number_components, phase, pressure_step, temperature_step, liquid_compositions, vapor_compositions, component_properties)

            for i in range(number_components):
                y = vapor_compositions[i]
                if abs(y) > THRESHOLD:
                    liquid_fugacity_coefficients[i] = liquid_fugacities[i] / \
                        (liquid_compositions[i]*pressure)
                    fugacity_coefficients_reason[i] = liquid_fugacity_coefficients[i] / \
                        vapor_fugacity_coefficients[i]
                    sum_list[count-1] = sum_list[count-1] + \
                        y/fugacity_coefficients_reason[i]
                    liquid_step_fugacity_coefficients[i] = liquid_step_fugacities[i]/(
                        liquid_compositions[i]*pressure_step)
                    step_fugacity_coefficient_reason[i] = liquid_step_fugacity_coefficients[i] / \
                        vapor_step_fugacity_coefficients[i]
                    sum1_list[count-1] = sum1_list[count-1] + \
                        y / \
                        step_fugacity_coefficient_reason[i]

            test_1 = abs(sum_list[count-1]-sum_list[count-2])

        newton_raphson_function_value = sum_list[count-1] - 1.0
        if abs(newton_raphson_function_value) <= 1.0e-6:
            break

        newton_raphson_function_step_value = sum1_list[count-1] - 1.0
        differential_newton_raphson_function_value = (
            newton_raphson_function_step_value - newton_raphson_function_value)/(pressure_step - pressure)

        for i in range(number_components):
            if abs(vapor_compositions[i]) > THRESHOLD:
                liquid_compositions[i] = vapor_compositions[i] / \
                    fugacity_coefficients_reason[i]/sum_list[count-1]

        newton_raphson_reduce_factor = aurea_dew_pressure(liquid_compositions, vapor_compositions, temperature, pressure, newton_raphson_function_value,
                                                          differential_newton_raphson_function_value, binary_model_param, number_components, mixing_rule_param, component_properties, uniquac_properties)

        pressure = pressure - (newton_raphson_function_value /
                               differential_newton_raphson_function_value)*newton_raphson_reduce_factor
        if pressure < 0.0:
            pressure = 0.05

    return pressure
