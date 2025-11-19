# -*- coding: utf-8 -*-
"""Drux Higuchi model implementation."""

from .base_model import DrugReleaseModel
from .messages import (
    ERROR_INVALID_DIFFUSION,
    ERROR_INVALID_CONCENTRATION,
    ERROR_INVALID_SOLUBILITY,
    ERROR_SOLUBILITY_HIGHER_THAN_CONCENTRATION,
)
from dataclasses import dataclass
from math import sqrt


@dataclass
class HiguchiParameters:
    """
    Parameters for the Higuchi model based on physical formulation.

    Attributes:
        D (float): Drug diffusivity in the polymer carrier (cm^2/s)
        c0 (float): Initial drug concentration (mg/cm^3)
        cs (float): Drug solubility in the polymer (mg/cm^3)
    """

    D: float
    c0: float
    cs: float


class HiguchiModel(DrugReleaseModel):
    """Simulator for the Higuchi drug release model using analytical expressions based on concentration conditions."""

    def __init__(self, D: float, c0: float, cs: float) -> None:
        """
        Initialize the Higuchi model with the given parameters.

        :param D: Drug diffusivity in the polymer carrier (cm^2/s)
        :param c0: Initial drug concentration (mg/cm^3)
        :param cs: Drug solubility in the polymer (mg/cm^3)
        """
        super().__init__()
        self._parameters = HiguchiParameters(D=D, c0=c0, cs=cs)
        self._plot_parameters["label"] = "Higuchi Model"

    def __repr__(self):
        """Return a string representation of the Higuchi model."""
        return (
            f"drux.HiguchiModel(D={self._parameters.D}, "
            f"c0={self._parameters.c0}, cs={self._parameters.cs})"
        )

    def _model_function(self, t: float) -> float:
        """
        Calculate the drug release at time t using the Higuchi model.

        Formula:
        - General case: Mt = sqrt(D * c0 * (2*c0 - cs) * cs * t)
        :param t: time (s)
        """
        D = self._parameters.D
        c0 = self._parameters.c0
        cs = self._parameters.cs

        Mt = sqrt(D * (2 * c0 - cs) * cs * t)

        return Mt

    def _validate_parameters(self) -> None:
        """Validate the parameters of the Higuchi model."""
        if self._parameters.D <= 0:
            raise ValueError(ERROR_INVALID_DIFFUSION)
        if self._parameters.c0 <= 0:
            raise ValueError(ERROR_INVALID_CONCENTRATION)
        if self._parameters.cs <= 0:
            raise ValueError(ERROR_INVALID_SOLUBILITY)
        if self._parameters.cs > self._parameters.c0:
            raise ValueError(ERROR_SOLUBILITY_HIGHER_THAN_CONCENTRATION)
