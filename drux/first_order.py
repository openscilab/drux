# -*- coding: utf-8 -*-
"""Drux first-order model implementation."""
from .base_model import DrugReleaseModel
from .utils import create_parameters_dataclass


class FirstOrderModel(DrugReleaseModel):
    """Simulator for the first-order drug release model."""

    def __init__(self, k: float, M0: float = 1) -> None:
        """
        Initialize the first-order model with the given parameters.

        :param k: first-order release rate constant (1/s)
        :param M0: entire releasable amount of drug (the asymptotic maximum) (mg)
        """
        super().__init__()
        self._model_name = "first_order"
        _params_class = create_parameters_dataclass(self._model_name)
        self._parameters = _params_class(k=k, M0=M0)
        self._plot_parameters["label"] = "First-Order Model"

    def __repr__(self):
        """Return a string representation of the First-Order model."""
        return f"drux.FirstOrderModel(k={self._parameters.k}, M0={self._parameters.M0})"
