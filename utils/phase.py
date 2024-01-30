""" Utilities classes for operation with phases"""
from enum import Enum


class Phase(Enum):
    """ A Class to set which phase was chosen. If 1 it is LIQUID, 2 for vapor"""
    LIQUID = 1
    VAPOR = 2
