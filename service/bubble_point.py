"""Service for bubble point operations"""

from typing import List
from service.helpers.fugacity import fugacity_pr_calculate, fugacity_uniquac_calculate
from utils.phase import Phase
from utils.component_properties import ComponentProperties
from utils.constants import THRESHOLD
from utils.uniquac_properties import UniquacProperties
from utils.exceptions import ConvergenceError
from utils.binary_model_param import BinaryModelParam
from utils.mixing_rule_param import MixingRuleParam


def calc_bubble_pressure(number_components: int, liquid_compositions: List[float],
                         temperature: float, binary_model_param: BinaryModelParam, component_properties: ComponentProperties, uniquac_properties: UniquacProperties, mixing_rule_param=MixingRuleParam.HK):
    """Calculate bubble pressure for the condition provided

    Args:
        number_components (int): number of components in the system
        liquid_compositions (List[float]): _description_
        temperature (float): temperature of the system
        binary_model_param (BinaryModelParam): _description_
        mixing_rule_param (MixingRuleParam): _description_
        component_properties (ComponentProperties): _description_
        uniquac_properties (UniquacProperties): _description_

    Raises:
        ConvergenceError: When count of newton raphson reaches 1000 

    Returns:
        pressure(float): Bubble pressure for the scenario provided
    """
    # ****************************************************************************************************************************************************
    #
    #     SUB-ROTINA: BOLHAP
    #
    #     OBJETIVO: Cálculo do ponto de bolha (T e X conhecidos)
    #
    # ****************************************************************************************************************************************************
    try:

        # Initial Estimations
        sum_list = [0.0]
        sum1_list = [0.0]
        pressure = 100  # bar
        count = 1

        liquid_fugacity_coefficients = [0.0] * number_components
        liquid_step_fugacity_coefficients = [0.0] * number_components
        vapor_fugacities = [0.0] * number_components
        vapor_step_fugacities = [0.0] * number_components
        vapor_fugacity_coefficients = [0.0] * number_components
        vapor_step_fugacity_coefficients = [0.0] * number_components
        fugacity_coefficients_reason = [0.0] * number_components
        step_fugacity_coefficient_reason = [0.0] * number_components

        normalized_number_components = number_components
        for i in range(number_components):
            if abs(liquid_compositions[i]) < THRESHOLD:
                normalized_number_components = normalized_number_components - 1

        vapor_compositions = [0.0 if abs(liquid_compositions[i]) < THRESHOLD else 1.0 /
                              normalized_number_components for i in range(number_components)]

        newton_raphson_function_value = 1

        # Newton Raphson for Buble Pressure
        while abs(newton_raphson_function_value) > 1e-6:
            if count == 1000:
                raise ConvergenceError(
                    "Bubble Pressure has not converged - count reached 1000 units")
            pressure_step = pressure + 0.001*pressure
            temperature_step = temperature

            phase = Phase.LIQUID

            if binary_model_param == BinaryModelParam.PR_UNIQUAC:
                _liquid_compressability, _liquid_volume, liquid_fugacities = fugacity_uniquac_calculate(number_components, phase, pressure, temperature, liquid_compositions,
                                                                                                        vapor_compositions, component_properties, mixing_rule_param, uniquac_properties)
                _liquid_step_compressability, _liquid_step_volume, liquid_step_fugacities = fugacity_uniquac_calculate(
                    number_components, phase, pressure_step, temperature_step, liquid_compositions, vapor_compositions, component_properties, mixing_rule_param, uniquac_properties)
            elif binary_model_param == BinaryModelParam.PR:
                liquid_fugacities, _liquid_compressability, _liquid_volume = fugacity_pr_calculate(
                    number_components, phase, pressure, temperature, liquid_compositions, vapor_compositions, component_properties)
                liquid_step_fugacities, _liquid_step_compressability, _liquid_step_volume = fugacity_pr_calculate(
                    number_components, phase, pressure_step, temperature_step, liquid_compositions, vapor_compositions, component_properties)

            phase = Phase.VAPOR

            if binary_model_param == BinaryModelParam.PR_UNIQUAC:
                _vapor_compressability, _vapor_volume, vapor_fugacities = fugacity_uniquac_calculate(number_components, phase, pressure, temperature, liquid_compositions,
                                                                                                     vapor_compositions, component_properties, mixing_rule_param, uniquac_properties)
            elif binary_model_param == BinaryModelParam.PR:
                vapor_fugacities, _vapor_compressability, _vapor_volume = fugacity_pr_calculate(
                    number_components, phase, pressure, temperature, liquid_compositions, vapor_compositions, component_properties)

            for i in range(number_components):
                x = liquid_compositions[i]
                if abs(x) > THRESHOLD:

                    liquid_fugacity_coefficients[i] = liquid_fugacities[i] / \
                        (x * pressure)
                    liquid_step_fugacity_coefficients[i] = liquid_step_fugacities[i] / (
                        x * pressure_step)
                    vapor_fugacity_coefficients[i] = vapor_fugacities[i] / \
                        (vapor_compositions[i] * pressure)
                    fugacity_coefficients_reason[i] = liquid_fugacity_coefficients[i] / \
                        vapor_fugacity_coefficients[i]
                    sum_list[count-1] = sum_list[count-1] + \
                        fugacity_coefficients_reason[i] * x

            test_1 = 1.0
            while test_1 >= 1.0e-5:
                vapor_composition_total = 0.0
                for i in range(number_components):
                    if liquid_compositions[i] > THRESHOLD:
                        vapor_compositions[i] = fugacity_coefficients_reason[i] * \
                            liquid_compositions[i]/sum_list[count-1]
                    vapor_composition_total = vapor_composition_total + \
                        vapor_compositions[i]
                for i in range(number_components):
                    vapor_compositions[i] = vapor_compositions[i] / \
                        vapor_composition_total

                if count == 1000:
                    raise ConvergenceError(
                        'Bubble Pressure has not converged - count === 1000')
                count = count + 1
                sum_list.append(0.0)
                sum1_list.append(0.0)

                phase = Phase.VAPOR
                if binary_model_param == BinaryModelParam.PR_UNIQUAC:
                    _vapor_compressability, _vapor_volume, vapor_fugacities = fugacity_uniquac_calculate(number_components, phase, pressure, temperature, liquid_compositions,
                                                                                                         vapor_compositions, component_properties, mixing_rule_param, uniquac_properties)
                    _vapor_compressability_step, _vapor_volume_step, vapor_step_fugacities = fugacity_uniquac_calculate(
                        number_components, phase, pressure_step, temperature_step, liquid_compositions, vapor_compositions, component_properties, mixing_rule_param, uniquac_properties)
                elif binary_model_param == BinaryModelParam.PR:
                    vapor_fugacities, _vapor_compressability, _vapor_volume = fugacity_pr_calculate(
                        number_components, phase, pressure, temperature, liquid_compositions, vapor_compositions, component_properties)
                    vapor_step_fugacities, _vapor_compressability, _vapor_volume = fugacity_pr_calculate(
                        number_components, phase, pressure_step, temperature_step, liquid_compositions, vapor_compositions, component_properties)

                for i in range(number_components):
                    if abs(liquid_compositions[i]) > THRESHOLD:
                        vapor_fugacity_coefficients[i] = vapor_fugacities[i] / \
                            (vapor_compositions[i]*pressure)
                        fugacity_coefficients_reason[i] = liquid_fugacity_coefficients[i] / \
                            vapor_fugacity_coefficients[i]
                        sum_list[count-1] = sum_list[count-1] + \
                            fugacity_coefficients_reason[i] * \
                            liquid_compositions[i]
                        vapor_step_fugacity_coefficients[i] = vapor_step_fugacities[i]/(
                            vapor_compositions[i]*pressure_step)
                        step_fugacity_coefficient_reason[i] = liquid_step_fugacity_coefficients[i] / \
                            vapor_step_fugacity_coefficients[i]
                        sum1_list[count-1] = sum1_list[count-1] + \
                            step_fugacity_coefficient_reason[i] * \
                            liquid_compositions[i]

                test_1 = abs(sum_list[count-1]-sum_list[count-2])

            newton_raphson_function_value = sum_list[count-1] - 1.0
            if abs(newton_raphson_function_value) <= 1.0e-6:
                break

            newton_raphson_function_step_value = sum1_list[count-1] - 1.0
            differential_newton_raphson_function_value = (
                newton_raphson_function_step_value - newton_raphson_function_value)/(pressure_step - pressure)

            pressure = pressure - newton_raphson_function_value / \
                differential_newton_raphson_function_value
            if pressure < 0.0:
                pressure = 0.05

        return pressure
    except Exception as err:
        print(f'Erron buble pressure calculation: {err}')
        raise Exception(err)

# Verificar primeira inconsistência de tradução nos services (done)
# Verificar primeira inconsistência de tradução nos helpers
# - PR001 (DONE)
# - PHASES (DONE)
# - UNIQUAC
# - ABPRAL
