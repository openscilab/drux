# -*- coding: utf-8 -*-
"""Drux zero-order model implementation."""

from .base_model import DrugReleaseModel
from .messages import ERROR_ZERO_ORDER_RELEASE_RATE, ERROR_ZERO_ORDER_INITIAL_AMOUNT
from dataclasses import dataclass


@dataclass
class ZeroOrderParameters:
    """
    Parameters for the Zero-order model.

    Attributes:
        M0 (float):  initial amount of drug in the solution (most times, M0 = 0)
        k0 (float):  zero-order release rate constant (mg/s)
    """

    M0: float
    k0: float


class ZeroOrderModel(DrugReleaseModel):
    """Simulator for the Zero-order drug release model."""

    def __init__(self, k0: float, M0: float = 0) -> None:
        """
        Initialize the zero-order model with the given parameters.

        :param k0: Zero-order release rate constant (mg/s)
        :param M0: Initial amount of drug in the solution (mg), default is 0
        """
        super().__init__()
        self._parameters = ZeroOrderParameters(k0=k0, M0=M0)
        self._plot_parameters["label"] = "Zero-Order Model"

    def __repr__(self):
        """Return a string representation of the Zero-Order model."""
        return f"drux.ZeroOrderModel(k0={self._parameters.k0}, M0={self._parameters.M0})"

    def _model_function(self, t: float) -> float:
        """
        Calculate the drug release at time t using the zero-order model.

        Formula:
        - M(t) = M0 + k0 * t
        :param t: time (s)
        """
        M0 = self._parameters.M0
        k0 = self._parameters.k0

        Mt = M0 + k0 * t

        return Mt

    def _validate_parameters(self) -> None:
        """Validate the parameters of the zero-order model."""
        if self._parameters.M0 < 0:
            raise ValueError(ERROR_ZERO_ORDER_INITIAL_AMOUNT)
        if self._parameters.k0 < 0:
            raise ValueError(ERROR_ZERO_ORDER_RELEASE_RATE)
