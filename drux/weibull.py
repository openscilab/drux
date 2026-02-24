# -*- coding: utf-8 -*-
"""Drux Weibull model implementation."""

from .base_model import DrugReleaseModel
from .utils import create_parameters_dataclass


class WeibullModel(DrugReleaseModel):
    """Simulator for the Weibull drug release model using analytical expressions based on concentration conditions."""

    def __init__(self, a: float, b: float, M: float = 1) -> None:
        """
        Initialize the Weibull model with the given parameters.

        :param M: entire releasable amount of drug (normally M > 0) (mg)
        :param a: scale factor
        :param b: shape factor
        """
        super().__init__()
        self._model_name = "weibull"
        _params_class = create_parameters_dataclass(self._model_name)
        self._parameters = _params_class(M=M, a=a, b=b)
        self._plot_parameters["label"] = "Weibull Model"

    def __repr__(self):
        """Return a string representation of the Weibull model."""
        return f"drux.WeibullModel(M={self._parameters.M}, a={self._parameters.a}, b={self._parameters.b})"
