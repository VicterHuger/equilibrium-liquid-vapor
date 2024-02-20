""" Utilities classes for operation with PR Alpha models"""
from enum import Enum
from typing import List


class AlphaModel(Enum):
    """ A Class to set which Peng Robinson alpha model was chosen"""
    PR87 = 1
    AAT = 2


class AATParams():
    """ Class to storage aat params

    Args:
        aat1=(List[float]), 
        aat2=(List[float]), 
        aat3=(List[float])
    """

    def __init__(self, aat1: List[float], aat2: List[float], aat3: List[float]) -> None:
        self.aat1 = aat1
        self.aat2 = aat2
        self.aat3 = aat3
