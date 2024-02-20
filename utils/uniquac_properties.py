"""Uniquac utils"""
from typing import List


class UniquacProperties():
    """Class to storage Uniquac Properties
    Args:
        uij_0: List[List[float]]
        uij_t: List[List[float]]
        r: List[float]
        q: List[float]
    """

    def __init__(self, uij_0: List[List[float]], uij_t: List[List[float]], r: List[float], q: List[float]) -> None:
        self.uij_0 = uij_0
        self.uij_t = uij_t
        self.r = r
        self.q = q
