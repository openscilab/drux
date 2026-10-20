# -*- coding: utf-8 -*-
"""Drux quadratic model implementation."""

from .base_model import DrugReleaseModel
from .messages import ERROR_QUADRATIC_INITIAL_AMOUNT
from dataclasses import dataclass


@dataclass
class QuadraticParameters:
    """
    Parameters for the quadratic model.

    Attributes:
        M0 (float): initial amount of drug in the solution (most times, M0 = 0) (mg)
        k1 (float): linear release coefficient (mg/s)
        k2 (float): quadratic release coefficient (mg/s^2)
    """

    M0: float
    k1: float
    k2: float


class QuadraticModel(DrugReleaseModel):
    """Simulator for the quadratic drug release model."""

    def __init__(self, k1: float, k2: float, M0: float = 0) -> None:
        """
        Initialize the quadratic model with the given parameters.

        :param k1: linear release coefficient (mg/s)
        :param k2: quadratic release coefficient (mg/s^2)
        :param M0: initial amount of drug in the solution (mg), default is 0
        """
        super().__init__()
        self._parameters = QuadraticParameters(M0=M0, k1=k1, k2=k2)
        self._plot_parameters["label"] = "Quadratic Model"

    def __repr__(self):
        """Return a string representation of the quadratic model."""
        return (
            f"drux.QuadraticModel(k1={self._parameters.k1}, k2={self._parameters.k2}, "
            f"M0={self._parameters.M0})"
        )

    @staticmethod
    def model_function(t: float, M0: float, k1: float, k2: float) -> float:
        """
        Calculate the drug release at time t using the quadratic model.

        Formula:
        - M(t) = M0 + k1 * t + k2 * t^2
        :param t: time (s)
        :param M0: initial amount of drug in the solution (mg)
        :param k1: linear release coefficient (mg/s)
        :param k2: quadratic release coefficient (mg/s^2)
        """
        Mt = M0 + k1 * t + k2 * t**2

        return Mt

    def _validate_parameters(self) -> None:
        """Validate the parameters of the quadratic model."""
        if self._parameters.M0 < 0:
            raise ValueError(ERROR_QUADRATIC_INITIAL_AMOUNT)
