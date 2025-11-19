# -*- coding: utf-8 -*-
"""Drux Weibull model implementation."""

from .base_model import DrugReleaseModel
from .messages import (
    ERROR_WEIBULL_SCALE_PARAMETER,
    ERROR_WEIBULL_INITIAL_AMOUNT,
    ERROR_WEIBULL_SHAPE_PARAMETER,
)
from dataclasses import dataclass
from math import exp


@dataclass
class WeibullParameters:
    """
    Parameters for the Weibull model based on physical formulation.

    Attributes:
        M (float): entire releasable amount of drug (normally M > 0) (mg)
        a (float): scale factor
        b (float): shape factor
    """

    M: float
    a: float
    b: float


class WeibullModel(DrugReleaseModel):
    """Simulator for the Weibull drug release model using analytical expressions based on concentration conditions."""

    def __init__(self, M: float, a: float, b: float) -> None:
        """
        Initialize the Weibull model with the given parameters.

        :param M: entire releasable amount of drug (normally M > 0) (mg)
        :param a: scale factor
        :param b: shape factor
        """
        super().__init__()
        self._parameters = WeibullParameters(M=M, a=a, b=b)
        self._plot_parameters["label"] = "Weibull Model"

    def __repr__(self):
        """Return a string representation of the Weibull model."""
        return f"drux.WeibullModel(M={self._parameters.M}, a={self._parameters.a}, b={self._parameters.b})"

    def _model_function(self, t: float) -> float:
        """
        Calculate the drug release at time t using the Weibull model.

        Formula:
        - General case: Mt = M * (1 - exp(-a*t ** b))
        :param t: time (s)
        """
        M = self._parameters.M
        a = self._parameters.a
        b = self._parameters.b

        Mt = M * (1 - exp(-a * t**b))

        return Mt

    def _validate_parameters(self) -> None:
        """Validate the parameters of the Weibull model."""
        if self._parameters.M < 0:
            raise ValueError(ERROR_WEIBULL_INITIAL_AMOUNT)
        if self._parameters.a <= 0:
            raise ValueError(ERROR_WEIBULL_SCALE_PARAMETER)
        if self._parameters.b <= 0:
            raise ValueError(ERROR_WEIBULL_SHAPE_PARAMETER)
