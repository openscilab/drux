# -*- coding: utf-8 -*-
"""Drux Higuchi model implementation."""

from .base_model import DrugReleaseModel
from .utils import create_parameters_dataclass


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
        self._model_name = "higuchi"
        _params_class = create_parameters_dataclass(self._model_name)
        self._parameters = _params_class(D=D, c0=c0, cs=cs)
        self._plot_parameters["label"] = "Higuchi Model"

    def __repr__(self):
        """Return a string representation of the Higuchi model."""
        return (
            f"drux.HiguchiModel(D={self._parameters.D}, "
            f"c0={self._parameters.c0}, cs={self._parameters.cs})"
        )
