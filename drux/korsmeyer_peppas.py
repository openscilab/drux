# -*- coding: utf-8 -*-
"""Drux Korsmeyer-Peppas model implementation."""

from .base_model import DrugReleaseModel
from .messages import (
    ERROR_KORSMEYER_PEPPAS_INITIAL_AMOUNT,
    ERROR_KORSMEYER_PEPPAS_RELEASE_RATE,
    ERROR_KORSMEYER_PEPPAS_RELEASE_EXPONENT,
)
from dataclasses import dataclass


@dataclass
class KorsmeyerPeppasParameters:
    """
    Parameters for the Korsmeyer-Peppas model.

    Attributes:
        M0 (float): initial amount of drug in the solution (most times, M0 = 0) (mg)
        k (float): release rate constant (mg/s^n)
        n (float): release exponent
    """

    M0: float
    k: float
    n: float


class KorsmeyerPeppasModel(DrugReleaseModel):
    """Simulator for the Korsmeyer-Peppas drug release model."""

    def __init__(self, k: float, n: float, M0: float = 0) -> None:
        """
        Initialize the Korsmeyer-Peppas model with the given parameters.

        :param k: release rate constant (mg/s^n)
        :param n: release exponent
        :param M0: initial amount of drug in the solution (mg), default is 0
        """
        super().__init__()
        self._parameters = KorsmeyerPeppasParameters(M0=M0, k=k, n=n)
        self._plot_parameters["label"] = "Korsmeyer-Peppas Model"

    def __repr__(self):
        """Return a string representation of the Korsmeyer-Peppas model."""
        return (
            f"drux.KorsmeyerPeppasModel(k={self._parameters.k}, n={self._parameters.n}, "
            f"M0={self._parameters.M0})"
        )

    @staticmethod
    def model_function(t: float, M0: float, k: float, n: float) -> float:
        """
        Calculate the drug release at time t using the Korsmeyer-Peppas model.

        Formula:
        - M(t) = M0 + k * t^n
        :param t: time (s)
        :param M0: initial amount of drug in the solution (mg)
        :param k: release rate constant (mg/s^n)
        :param n: release exponent
        """
        Mt = M0 + k * t**n

        return Mt

    def _validate_parameters(self) -> None:
        """Validate the parameters of the Korsmeyer-Peppas model."""
        if self._parameters.M0 < 0:
            raise ValueError(ERROR_KORSMEYER_PEPPAS_INITIAL_AMOUNT)
        if self._parameters.k < 0:
            raise ValueError(ERROR_KORSMEYER_PEPPAS_RELEASE_RATE)
        if self._parameters.n <= 0:
            raise ValueError(ERROR_KORSMEYER_PEPPAS_RELEASE_EXPONENT)
