""" Utilities classes for operation with binary model param"""
from enum import Enum


class BinaryModelParam(Enum):
    """ A Class to set which binary model param was chosen. If 1 it is Peng Robinson, 2 for Peng Robinson with Uniquac"""
    PR = 1
    PR_UNIQUAC = 2
