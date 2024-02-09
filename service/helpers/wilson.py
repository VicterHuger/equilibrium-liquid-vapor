"""File to calculate initial x(i) - wilson correlation """
from typing import List
from utils.component_properties import ComponentProperties
import math


def wilson_correlation(number_components: int, pressure: float, temperature: float, vapor_compositions: List[float], component_properties: ComponentProperties):
    liquid_compositions = [0.0]*number_components
    for i in range(number_components):
        equilibrium_constant_k = (component_properties.critical_acentric_properties.critical_pressures[i]/pressure)*(math.exp(5.373*(
            1.0+component_properties.critical_acentric_properties.acentric_factors[i])*(1.0-(component_properties.critical_acentric_properties.critical_temperatures[i]/temperature))))
        liquid_compositions[i] = vapor_compositions[i]/equilibrium_constant_k

    liquid_compositions = [liquid_composition/sum(liquid_compositions)
                           for liquid_composition in liquid_compositions]

    return liquid_compositions
