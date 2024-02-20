""" Utilities classes for operation with mixing rule param"""
from enum import Enum


class MixingRuleParam(Enum):
    """ A Class to set which mixing rule param was chosen. If 1 it is Heidemann Kokal, 2 for LCVM"""
    HK = 1
    LCVM = 2
