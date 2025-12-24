# -*- coding: utf-8 -*-
"""Drux Hopfenberg model implementation."""

from .base_model import DrugReleaseModel
from .messages import (
    ERROR_INVALID_EROSION_CONSTANT,
    ERROR_INVALID_INITIAL_RADIUS,
    ERROR_INVALID_GEOMETRY_FACTOR,
    ERROR_INVALID_CONCENTRATION,
    ERROR_RELEASABLE_AMOUNT,
)
from dataclasses import dataclass


@dataclass
class HopfenbergParameters:
    """
    Parameters for the Hopfenberg model based on surface erosion.

    Attributes:
        k0 (float): Erosion rate constant (mg/(mm^2·s))
        c0 (float): Initial drug concentration in the matrix (mg/mm^3)
        a0 (float): Initial radius or half-thickness of the device (mm)
        n (int): Geometry factor (1=slab, 2=cylinder, 3=sphere)
    """

    M: float
    k0: float
    c0: float
    a0: float
    n: int


class HopfenbergModel(DrugReleaseModel):
    """Simulator for the Hopfenberg drug release model for surface-eroding polymers."""

    def __init__(self, M: float, k0: float, c0: float, a0: float, n: int) -> None:
        """
        Initialize the Hopfenberg model with the given parameters.

        :param M: entire releasable amount of drug (normally M > 0) (mg)
        :param k0: Erosion rate constant (mg/(mm^2·s))
        :param c0: Initial drug concentration in the matrix (mg/mm^3)
        :param a0: Initial radius or half-thickness of the device (mm)
        :param n: Geometry factor (1=slab, 2=cylinder, 3=sphere)
        """
        super().__init__()
        self._parameters = HopfenbergParameters(M=M, k0=k0, c0=c0, a0=a0, n=n)
        self._plot_parameters["label"] = "Hopfenberg Model"

    def __repr__(self):
        """Return a string representation of the Hopfenberg model."""
        return (
            f"drux.HopfenbergModel(M={self._parameters.M}, k0={self._parameters.k0}, "
            f"c0={self._parameters.c0}, a0={self._parameters.a0}, "
            f"n={self._parameters.n})"
        )

    def _model_function(self, t: float) -> float:
        """
        Calculate the fractional drug release at time t using the Hopfenberg model.

        Formula:
        - Mt = M∞(1 - (1 - k0*t / (c0*a0))^n)

        :param t: time (s)
        :return: drug release
        """
        M = self._parameters.M
        k0 = self._parameters.k0
        c0 = self._parameters.c0
        a0 = self._parameters.a0
        n = self._parameters.n

        inner_term = 1 - (k0 * t) / (c0 * a0)

        Mt = M * (1 - (inner_term**n))

        return Mt

    def _validate_parameters(self) -> None:
        """Validate the parameters of the Hopfenberg model."""
        if self._parameters.M < 0:
            raise ValueError(ERROR_RELEASABLE_AMOUNT)
        if self._parameters.k0 < 0:
            raise ValueError(ERROR_INVALID_EROSION_CONSTANT)
        if self._parameters.c0 <= 0:
            raise ValueError(ERROR_INVALID_CONCENTRATION)
        if self._parameters.a0 <= 0:
            raise ValueError(ERROR_INVALID_INITIAL_RADIUS)
        if self._parameters.n not in (1, 2, 3):
            raise ValueError(ERROR_INVALID_GEOMETRY_FACTOR)
