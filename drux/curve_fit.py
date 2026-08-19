# -*- coding: utf-8 -*-
"""Drux curve fitting utilities."""

import inspect
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple, Type

import numpy as np
from scipy.optimize import curve_fit

from .base_model import DrugReleaseModel
from .first_order import FirstOrderModel
from .higuchi import HiguchiModel
from .hopfenberg import HopfenbergModel
from .weibull import WeibullModel
from .zero_order import ZeroOrderModel
from .messages import (
    ERROR_UNKNOWN_MODEL,
    ERROR_TIME_RELEASE_LENGTH_MISMATCH,
    ERROR_INSUFFICIENT_DATA_POINTS,
    ERROR_UNKNOWN_FIT_PARAMETER,
    ERROR_NO_FREE_PARAMETERS,
    ERROR_NO_FIT_RESULT,
)

# Single point of registration mapping a model name to its class.
MODEL_CLASSES = {
    "zero_order": ZeroOrderModel,
    "first_order": FirstOrderModel,
    "higuchi": HiguchiModel,
    "weibull": WeibullModel,
    "hopfenberg": HopfenbergModel,
}


def _model_parameter_names(model_class: Type[DrugReleaseModel]) -> List[str]:
    """Return a model's parameter names, in declaration order, excluding time."""
    parameters = inspect.signature(model_class.model_function).parameters
    return [name for name in parameters if name != "t"]


@dataclass
class FitResult:
    """
    Result of fitting a drug release model to experimental data.

    Attributes:
        model_name (str): name of the fitted model
        parameters (Dict[str, float]): known and fitted parameter values, keyed by name
        r_squared (float): coefficient of determination of the fit
        model (DrugReleaseModel): model instance initialized with the fitted parameters
    """

    model_name: str
    parameters: Dict[str, float]
    r_squared: float
    model: DrugReleaseModel


class CurveFit:
    """Fit a registered drug release model to experimental time-series data."""

    def __init__(
        self,
        model_name: str,
        time: Sequence[float],
        release_profile: Sequence[float],
        known_parameters: Optional[Dict[str, float]] = None,
    ) -> None:
        """
        Initialize the curve fitter.

        :param model_name: name of a registered model (see `MODEL_CLASSES`)
        :param time: time points of the experimental data (s)
        :param release_profile: measured drug release at each time point
        :param known_parameters: parameter values to hold fixed instead of fitting
        """
        if model_name not in MODEL_CLASSES:
            raise ValueError(ERROR_UNKNOWN_MODEL.format(model_name, sorted(MODEL_CLASSES)))

        time = np.asarray(time, dtype=float)
        release_profile = np.asarray(release_profile, dtype=float)
        if time.shape != release_profile.shape:
            raise ValueError(ERROR_TIME_RELEASE_LENGTH_MISMATCH)
        if time.size < 2:
            raise ValueError(ERROR_INSUFFICIENT_DATA_POINTS)

        self._model_name = model_name
        self._model_class = MODEL_CLASSES[model_name]
        self._time = time
        self._release_profile = release_profile
        self._known_parameters = known_parameters or {}
        self._parameter_names = _model_parameter_names(self._model_class)

        unknown_keys = sorted(set(self._known_parameters) - set(self._parameter_names))
        if unknown_keys:
            raise ValueError(
                ERROR_UNKNOWN_FIT_PARAMETER.format(unknown_keys, model_name, self._parameter_names)
            )

        self._free_parameters = [p for p in self._parameter_names if p not in self._known_parameters]
        if not self._free_parameters:
            raise ValueError(ERROR_NO_FREE_PARAMETERS)

        self._fit_result: Optional[FitResult] = None

    def _merge_parameters(self, free_values: Sequence[float]) -> Dict[str, float]:
        """Merge fitted values for the free parameters with the known parameters."""
        free_values = iter(free_values)
        return {
            name: self._known_parameters[name] if name in self._known_parameters else next(free_values)
            for name in self._parameter_names
        }

    def _equation(self, t: np.ndarray, *free_values: float) -> np.ndarray:
        """Evaluate the model's own equation for a set of free parameter values."""
        parameters = self._merge_parameters(free_values)
        return np.vectorize(lambda ti: self._model_class.model_function(ti, **parameters))(t)

    def fit(
        self,
        initial_guess: Optional[Sequence[float]] = None,
        bounds: Optional[Tuple[Sequence[float], Sequence[float]]] = None,
    ) -> FitResult:
        """
        Estimate the unknown parameters via non-linear least squares.

        :param initial_guess: initial guess for each free parameter, in the model's
            declared parameter order (default: 1.0 for each)
        :param bounds: (lower, upper) bounds for the free parameters, forwarded to
            `scipy.optimize.curve_fit` (default: (~0, inf))
        """
        n_free = len(self._free_parameters)
        p0 = list(initial_guess) if initial_guess is not None else [1.0] * n_free
        lower_upper = bounds or ([1e-10] * n_free, [np.inf] * n_free)

        fitted_values, _ = curve_fit(
            self._equation, self._time, self._release_profile, p0=p0, bounds=lower_upper, maxfev=10000
        )
        parameters = self._merge_parameters(fitted_values)

        predicted = self._equation(self._time, *fitted_values)
        residual_sum_of_squares = np.sum((self._release_profile - predicted) ** 2)
        total_sum_of_squares = np.sum((self._release_profile - np.mean(self._release_profile)) ** 2)
        r_squared = 1 - residual_sum_of_squares / total_sum_of_squares if total_sum_of_squares > 0 else 0.0

        fit_result = FitResult(
            model_name=self._model_name,
            parameters=parameters,
            r_squared=r_squared,
            model=self._model_class(**parameters),
        )
        self._fit_result = fit_result
        return fit_result

    def get_result(self) -> FitResult:
        """Return the most recent fit result."""
        if self._fit_result is None:
            raise ValueError(ERROR_NO_FIT_RESULT)
        return self._fit_result

    def __repr__(self) -> str:
        """Return a string representation of the curve fitter."""
        if self._fit_result is None:
            return f"drux.CurveFit({self._model_name}, not fitted)"
        params = ", ".join(f"{k}={v:.4f}" for k, v in self._fit_result.parameters.items())
        return f"drux.CurveFit({self._model_name}: {params}, R²={self._fit_result.r_squared:.4f})"
