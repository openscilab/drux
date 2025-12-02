# -*- coding: utf-8 -*-
"""Drux first-order model implementation."""
from math import exp
from .base_model import DrugReleaseModel
from .messages import ERROR_FIRST_ORDER_INITIAL_AMOUNT, ERROR_FIRST_ORDER_RELEASE_RATE
from dataclasses import dataclass


@dataclass
class FirstOrderParameters:
    """
    Parameters for the first-order model.

    Attributes:
        M0 (float):  entire releasable amount of drug (normally M0 > 0) (mg)
        k (float):  first-order release rate constant (1/s)
    """

    M0: float
    k: float


class FirstOrderModel(DrugReleaseModel):
    """Simulator for the first-order drug release model."""

    def __init__(self, k: float, M0: float) -> None:
        """
        Initialize the first-order model with the given parameters.

        :param k: first-order release rate constant (1/s)
        :param M0: entire releasable amount of drug (the asymptotic maximum) (mg)
        """
        super().__init__()
        self._parameters = FirstOrderParameters(k=k, M0=M0)
        self._plot_parameters["label"] = "First-Order Model"

    def __repr__(self):
        """Return a string representation of the First-Order model."""
        return f"drux.FirstOrderModel(k={self._parameters.k}, M0={self._parameters.M0})"

    def _model_function(self, t: float) -> float:
        """
        Calculate the drug release at time t using the first-order model.

        Formula:
        - M(t) = M0 * (1 - exp(-k * t))
        :param t: time (s)
        """
        M0 = self._parameters.M0
        k = self._parameters.k

        Mt = M0 * (1 - exp(-k * t))

        return Mt

    def _validate_parameters(self) -> None:
        """Validate the parameters of the first-order model."""
        if self._parameters.M0 < 0:
            raise ValueError(ERROR_FIRST_ORDER_INITIAL_AMOUNT)
        if self._parameters.k < 0:
            raise ValueError(ERROR_FIRST_ORDER_RELEASE_RATE)
