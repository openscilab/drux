# -*- coding: utf-8 -*-
"""Drux cubic model implementation."""

from .base_model import DrugReleaseModel
from .messages import ERROR_CUBIC_INITIAL_AMOUNT
from dataclasses import dataclass


@dataclass
class CubicParameters:
    """
    Parameters for the cubic model.

    Attributes:
        M0 (float): initial amount of drug in the solution (most times, M0 = 0) (mg)
        k1 (float): linear release coefficient (mg/s)
        k2 (float): quadratic release coefficient (mg/s^2)
        k3 (float): cubic release coefficient (mg/s^3)
    """

    M0: float
    k1: float
    k2: float
    k3: float


class CubicModel(DrugReleaseModel):
    """Simulator for the cubic drug release model."""

    def __init__(self, k1: float, k2: float, k3: float, M0: float = 0) -> None:
        """
        Initialize the cubic model with the given parameters.

        :param k1: linear release coefficient (mg/s)
        :param k2: quadratic release coefficient (mg/s^2)
        :param k3: cubic release coefficient (mg/s^3)
        :param M0: initial amount of drug in the solution (mg), default is 0
        """
        super().__init__()
        self._parameters = CubicParameters(M0=M0, k1=k1, k2=k2, k3=k3)
        self._plot_parameters["label"] = "Cubic Model"

    def __repr__(self):
        """Return a string representation of the cubic model."""
        return (
            f"drux.CubicModel(k1={self._parameters.k1}, k2={self._parameters.k2}, "
            f"k3={self._parameters.k3}, M0={self._parameters.M0})"
        )

    @staticmethod
    def model_function(t: float, M0: float, k1: float, k2: float, k3: float) -> float:
        """
        Calculate the drug release at time t using the cubic model.

        Formula:
        - M(t) = M0 + k1 * t + k2 * t^2 + k3 * t^3
        :param t: time (s)
        :param M0: initial amount of drug in the solution (mg)
        :param k1: linear release coefficient (mg/s)
        :param k2: quadratic release coefficient (mg/s^2)
        :param k3: cubic release coefficient (mg/s^3)
        """
        Mt = M0 + k1 * t + k2 * t**2 + k3 * t**3

        return Mt

    def _validate_parameters(self) -> None:
        """Validate the parameters of the cubic model."""
        if self._parameters.M0 < 0:
            raise ValueError(ERROR_CUBIC_INITIAL_AMOUNT)
