"""Utilities classes for storaging component properties"""
from typing import List
from utils.alpha_model import AlphaModel, AATParams


class CriticalAcentricProperties():
    """ Class to storage critical component properties
    Args:
        acentric_factors: (List[float])
        critical_temperatures: (List[float]), 
        critical_pressures: (List[float])
    """

    def __init__(self, acentric_factors: List[float], critical_temperatures: List[float], critical_pressures: List[float]) -> None:
        self.acentric_factors = acentric_factors
        self.critical_temperatures = critical_temperatures
        self.critical_pressures = critical_pressures


class ComponentProperties():
    """Class to storage component properties

    Args:
        peng_robinson_binary_params: (List[List])
        critical_acentric_properties: CriticalAcentricProperties
        alpha_model: AlphaModel,
        aat_params: AATParams
    """

    def __init__(self, peng_robinson_binary_params: List[List[float]], critical_acentric_properties: CriticalAcentricProperties, alpha_model: AlphaModel, aat_params: AATParams | None) -> None:
        self.peng_robinson_binary_params = peng_robinson_binary_params
        self.critical_acentric_properties = critical_acentric_properties
        self.alpha_model = alpha_model
        self.aat_params = aat_params
