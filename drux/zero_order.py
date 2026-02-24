# -*- coding: utf-8 -*-
"""Drux zero-order model implementation."""

from .base_model import DrugReleaseModel
from .utils import create_parameters_dataclass


class ZeroOrderModel(DrugReleaseModel):
    def __init__(self, k0: float, M0: float = 0):
        super().__init__()
        self._model_name = "zero_order"
        _params_class = create_parameters_dataclass(self._model_name)
        self._parameters = _params_class(k0=k0, M0=M0)
        self._plot_parameters["label"] = "Zero-Order Model"

    def __repr__(self):
        """Return a string representation of the Zero-Order model."""
        return (
            f"drux.ZeroOrderModel(k0={self._parameters.k0}, M0={self._parameters.M0})"
        )
