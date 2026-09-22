# -*- coding: utf-8 -*-
"""Drux curve fitting utilities."""

import inspect
from dataclasses import dataclass
from math import isnan
from numbers import Real
from typing import Any, Dict, List, Optional, Sequence, Tuple, Type

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
    ERROR_UNKNOWN_FIT_ARGUMENT,
    ERROR_KNOWN_FIT_ARGUMENT,
    ERROR_MISSING_FIT_ARGUMENT,
    ERROR_FIT_ARGUMENT_TYPE,
    ERROR_INVALID_FIT_VALUE,
    ERROR_INVALID_FIT_BOUNDS,
    ERROR_FIT_BOUNDS_ORDER,
    ERROR_FIT_GUESS_OUT_OF_BOUNDS,
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

# Defaults applied to every free parameter that the caller does not specify.
DEFAULT_INITIAL_GUESS = 1.0
DEFAULT_BOUNDS = (1e-10, np.inf)


def _fit_value(value: Any, parameter_name: str, argument_name: str) -> float:
    """Return a value of a fit argument as a float.

    :param value: value given for the parameter
    :param parameter_name: name of the parameter, used in the error message
    :param argument_name: name of the fit argument, used in the error message
    """
    if not isinstance(value, Real) or isnan(value):
        raise ValueError(ERROR_INVALID_FIT_VALUE.format(parameter_name=parameter_name, argument_name=argument_name))
    return float(value)


def _model_parameter_names(model_class: Type[DrugReleaseModel]) -> List[str]:
    """Return a model's parameter names, in declaration order, excluding time.

    :param model_class: model class exposing a `model_function` signature
    """
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
        :param release_profile: a sequence of measured drug release at each time point
        :param known_parameters: parameter values to hold fixed instead of fitting
        """
        if model_name not in MODEL_CLASSES:
            raise ValueError(ERROR_UNKNOWN_MODEL.format(model_name=model_name, available_models=sorted(MODEL_CLASSES)))

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
                ERROR_UNKNOWN_FIT_PARAMETER.format(unknown_keys=unknown_keys, model_name=model_name, parameter_names=self._parameter_names)
            )

        self._free_parameters = [p for p in self._parameter_names if p not in self._known_parameters]
        if not self._free_parameters:
            raise ValueError(ERROR_NO_FREE_PARAMETERS)

        self._fit_result: Optional[FitResult] = None

    def _merge_parameters(self, free_values: Sequence[float]) -> Dict[str, float]:
        """Merge fitted values for free parameters with known fixed parameters.

        :param free_values: fitted values in the model's declared free-parameter order
        """
        free_values = iter(free_values)
        return {
            name: self._known_parameters[name] if name in self._known_parameters else next(free_values)
            for name in self._parameter_names
        }

    def _equation(self, t: np.ndarray, *free_values: Any) -> np.ndarray:
        """Evaluate the model equation for a set of free parameter values.

        :param t: time values where the model should be evaluated
        :param free_values: free parameter values passed by the optimizer
        """
        parameters = self._merge_parameters(free_values)
        return np.vectorize(lambda ti: self._model_class.model_function(ti, **parameters))(t)

    def _validate_fit_argument(
            self, values: Optional[Dict[str, Any]], argument_name: str) -> Optional[Dict[str, Any]]:
        """Check that a fit argument gives a value for every free parameter, and for no other parameter.

        :param values: fit argument keyed by parameter name, or None to use the defaults
        :param argument_name: name of the argument, used in the error messages
        """
        if values is None:
            return None
        if not isinstance(values, dict):
            raise ValueError(ERROR_FIT_ARGUMENT_TYPE.format(argument_name=argument_name))

        known_keys = sorted(set(values) & set(self._known_parameters))
        if known_keys:
            raise ValueError(ERROR_KNOWN_FIT_ARGUMENT.format(known_keys=known_keys, argument_name=argument_name))

        unknown_keys = sorted(set(values) - set(self._free_parameters))
        if unknown_keys:
            raise ValueError(
                ERROR_UNKNOWN_FIT_ARGUMENT.format(
                    unknown_keys=unknown_keys,
                    argument_name=argument_name,
                    model_name=self._model_name,
                    free_parameters=self._free_parameters)
            )

        missing_keys = [name for name in self._free_parameters if name not in values]
        if missing_keys:
            raise ValueError(
                ERROR_MISSING_FIT_ARGUMENT.format(
                    missing_keys=missing_keys,
                    argument_name=argument_name,
                    free_parameters=self._free_parameters)
            )
        return values

    def _parameter_bounds(self, bounds: Optional[Dict[str, Any]], name: str) -> Tuple[float, float]:
        """Return the validated `(lower, upper)` bounds of a single free parameter.

        :param bounds: bounds keyed by parameter name, or None to use the defaults
        :param name: name of the free parameter
        """
        if bounds is None:
            return DEFAULT_BOUNDS

        try:
            low, high = bounds[name]
        except (TypeError, ValueError):
            raise ValueError(ERROR_INVALID_FIT_BOUNDS.format(parameter_name=name))

        low = _fit_value(low, name, "bounds")
        high = _fit_value(high, name, "bounds")
        if low >= high:
            raise ValueError(ERROR_FIT_BOUNDS_ORDER.format(parameter_name=name))
        return low, high

    def _parameter_guess(self, initial_guess: Optional[Dict[str, Any]], name: str, low: float, high: float) -> float:
        """Return the validated initial guess of a single free parameter.

        :param initial_guess: initial guess keyed by parameter name, or None to use the default
        :param name: name of the free parameter
        :param low: validated lower bound of the parameter
        :param high: validated upper bound of the parameter
        """
        # A default guess is kept inside the bounds, so that bounds stay usable on their own.
        if initial_guess is None:
            return min(max(DEFAULT_INITIAL_GUESS, low), high)

        guess = _fit_value(initial_guess[name], name, "initial_guess")
        if not low <= guess <= high:
            raise ValueError(ERROR_FIT_GUESS_OUT_OF_BOUNDS.format(parameter_name=name))
        return guess

    def fit(
        self,
        initial_guess: Optional[Dict[str, float]] = None,
        bounds: Optional[Dict[str, Tuple[float, float]]] = None,
    ) -> FitResult:
        """
        Estimate the unknown parameters via non-linear least squares.

        Both arguments are keyed by parameter name, so the caller does not need to know the
        order in which the optimizer receives the free parameters. An argument that is given
        must hold a value for every free parameter; omit the argument to use the defaults.

        :param initial_guess: initial guess per free parameter, keyed by name (default: 1.0 for each)
        :param bounds: `(lower, upper)` bounds per free parameter, keyed by name (default: (~0, inf))

        :raises ValueError: if a parameter name or a value of either argument is invalid
        """
        initial_guess = self._validate_fit_argument(initial_guess, "initial_guess")
        bounds = self._validate_fit_argument(bounds, "bounds")

        lower, upper = zip(*(self._parameter_bounds(bounds, name) for name in self._free_parameters))
        p0 = [
            self._parameter_guess(initial_guess, name, low, high)
            for name, low, high in zip(self._free_parameters, lower, upper)
        ]

        fitted_values, _ = curve_fit(
            self._equation, self._time, self._release_profile, p0=p0, bounds=(list(lower), list(upper)), maxfev=10000
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
        return f"drux.CurveFit({self._model_name}: {params}, R^2={self._fit_result.r_squared:.4f})"
