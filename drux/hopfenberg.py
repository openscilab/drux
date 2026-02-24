# -*- coding: utf-8 -*-
"""Drux Hopfenberg model implementation."""

from .base_model import DrugReleaseModel
from .utils import create_parameters_dataclass


class HopfenbergModel(DrugReleaseModel):
    """Simulator for the Hopfenberg drug release model for surface-eroding polymers."""

    def __init__(self, k0: float, c0: float, a0: float, n: int, M: float = 1) -> None:
        """
        Initialize the Hopfenberg model with the given parameters.

        :param M: entire releasable amount of drug (normally M > 0) (mg)
        :param k0: Erosion rate constant (mg/(mm^2·s))
        :param c0: Initial drug concentration in the matrix (mg/mm^3)
        :param a0: Initial radius or half-thickness of the device (mm)
        :param n: Geometry factor (1=slab, 2=cylinder, 3=sphere)
        """
        super().__init__()
        self._model_name = "hopfenberg"
        _params_class = create_parameters_dataclass(self._model_name)
        self._parameters = _params_class(M=M, k0=k0, c0=c0, a0=a0, n=n)
        self._plot_parameters["label"] = "Hopfenberg Model"

    def __repr__(self):
        """Return a string representation of the Hopfenberg model."""
        return (
            f"drux.HopfenbergModel(M={self._parameters.M}, k0={self._parameters.k0}, "
            f"c0={self._parameters.c0}, a0={self._parameters.a0}, "
            f"n={self._parameters.n})"
        )
