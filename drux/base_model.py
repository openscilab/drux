"""This file contains the abstract base class for drug release models."""

import numpy as np
import matplotlib.pyplot as plt
from types import SimpleNamespace
from typing import Any, Optional

from .messages import (
    ERROR_TARGET_RELEASE_RANGE,
    ERROR_DURATION_TIME_STEP_POSITIVE,
    ERROR_TIME_STEP_GREATER_THAN_DURATION,
    ERROR_NO_SIMULATION_DATA,
    ERROR_RELEASE_PROFILE_TOO_SHORT,
    ERROR_TARGET_RELEASE_EXCEEDS_MAX,
)

from .registry import get_model_config


class DrugReleaseModel:
    """Base class for drug release models.

    This class provides a common interface and functionality for various
    mathematical models of drug release from delivery systems.

    Model classes are typically generated via :func:`create_model_class`
    from the centralized model registry.  All equation logic,
    parameter validation, and parameter metadata live in the registry
    so that subclasses need no custom overrides.
    """

    def __init__(self):
        """Initialize the drug release model."""
        self._time_points = None
        self._release_profile = None
        self._plot_parameters = {
            "xlabel": "Time (s)",
            "ylabel": "Cumulative Release",
            "title": "Drug Release Profile",
            "label": "Release Profile",
        }
        self._parameters = None
        self._model_name = None

    def _validate_parameters(self) -> None:
        """Validate model parameters using the model registry."""
        if self._model_name is None or self._parameters is None:
            raise ValueError("Model name and parameters must be set")

        config = get_model_config(self._model_name)

        for param_name, rule in config["validation"].items():
            if rule.get("cross_param"):
                # Cross-parameter validation (e.g., cs <= c0)
                if not rule["check"](self._parameters):
                    raise ValueError(rule["error"])
            else:
                # Single parameter validation
                param_value = getattr(self._parameters, param_name)
                if not rule["check"](param_value):
                    raise ValueError(rule["error"])

    def _model_function(self, t: float) -> float:
        """
        Model function that calculates drug release profile over time.

        :param t: time point at which to calculate drug release
        """
        if self._model_name is None or self._parameters is None:
            raise ValueError("Model name and parameters must be set")

        config = get_model_config(self._model_name)
        return config["equation"](self._parameters, t)

    def _get_release_profile(self) -> np.ndarray:
        """Calculate the drug release profile over the specified time points."""
        return np.vectorize(self._model_function)(self._time_points)

    def _validate_plot(self) -> tuple:
        """
        Validate plotting process.

        :raises ValueError: if simulation data is not available
        :raises ValueError: if release profile is too short
        """
        if self._time_points is None or self._release_profile is None:
            raise ValueError(ERROR_NO_SIMULATION_DATA)

        if len(self._release_profile) < 2:
            raise ValueError(ERROR_RELEASE_PROFILE_TOO_SHORT)
        fig, ax = plt.subplots()
        return fig, ax

    def simulate(self, duration: int, time_step: float = 1) -> np.ndarray:
        """
        Simulate drug release over time.

        :param duration: total time for simulation (in seconds)
        :param time_step: time step for simulation (in seconds)
        """
        if duration <= 0 or time_step <= 0:
            raise ValueError(ERROR_DURATION_TIME_STEP_POSITIVE)
        if time_step > duration:
            raise ValueError(ERROR_TIME_STEP_GREATER_THAN_DURATION)
        self._time_points = np.arange(0, duration + time_step, time_step)
        self._validate_parameters()
        self._release_profile = self._get_release_profile()
        return self._release_profile

    def plot(
        self,
        show: bool = True,
        label: Optional[str] = None,
        xlabel: Optional[str] = None,
        ylabel: Optional[str] = None,
        title: Optional[str] = None,
        **kwargs: Any
    ) -> tuple:
        """
        Plot the drug release profile.

        :param show: Whether to display the plot (default: True)
        :param label: The legend label for the release profile curve
        :param xlabel: Label for the x-axis
        :param ylabel: Label for the y-axis
        :param title: Title of the plot
        """
        # Create a new figure and axis if not provided
        fig, ax = self._validate_plot()

        # Plotting the release profile
        ax.plot(
            self._time_points,
            self._release_profile,
            label=label or self._plot_parameters["label"],
            **kwargs
        )
        ax.set_xlabel(xlabel or self._plot_parameters["xlabel"])
        ax.set_ylabel(ylabel or self._plot_parameters["ylabel"])
        ax.set_title(title or self._plot_parameters["title"])
        ax.grid()
        ax.legend()

        # Show the plot if requested
        if show:
            fig.show()

        return fig, ax

    def get_release_rate(self) -> np.ndarray:
        """Calculate the instantaneous release rate (derivative of release profile)."""
        if self._time_points is None or self._release_profile is None:
            raise ValueError(ERROR_NO_SIMULATION_DATA)

        if len(self._release_profile) < 2:
            raise ValueError(ERROR_RELEASE_PROFILE_TOO_SHORT)

        # Calculate the derivative of the release profile
        release_rate = np.gradient(self._release_profile, self._time_points)
        return release_rate

    def time_for_release(self, target_release: float) -> float:
        """
        Estimate time needed to reach a specific release percentage.

        :param target_release: target release fraction (>= 0)

        :raises ValueError: if target_release is negative
        :raises ValueError: if simulation data is not available
        :raises ValueError: if target_release exceeds maximum release
        """
        if self._time_points is None or self._release_profile is None:
            raise ValueError(ERROR_NO_SIMULATION_DATA)

        if target_release < 0:
            raise ValueError(ERROR_TARGET_RELEASE_RANGE)

        if target_release > self._release_profile[-1]:
            raise ValueError(ERROR_TARGET_RELEASE_EXCEEDS_MAX)

        # Find first time point where release >= target
        idx = np.argmax(self._release_profile >= target_release)
        return self._time_points[idx]


def _create_parameters(**kwargs):
    """Create a parameters namespace.

    :param kwargs: Parameter values
    """
    return SimpleNamespace(**kwargs)


def create_model_class(model_name: str, class_name: str, label: str, docstring: str = "") -> type:
    """Generate a :class:`DrugReleaseModel` subclass from the model registry.

    :param model_name: Key in the model registry (e.g. ``"zero_order"``)
    :param class_name: Name of the generated class (e.g. ``"ZeroOrderModel"``)
    :param label: Plot legend label (e.g. ``"Zero-Order Model"``)
    :param docstring: Docstring for the generated class
    """
    config = get_model_config(model_name)
    param_specs = config["params"]

    # Build the ordered list of (name, type, default|REQUIRED) for __init__
    required_params = []
    optional_params = []
    for pname, pinfo in param_specs.items():
        if pinfo["default"] is not None:
            optional_params.append((pname, pinfo["type"], pinfo["default"]))
        else:
            required_params.append((pname, pinfo["type"]))

    def __init__(self, **kwargs):
        # Apply defaults for missing optional params
        for pname, _ptype, pdefault in optional_params:
            kwargs.setdefault(pname, pdefault)
        # Check all required params are present
        for pname, _ptype in required_params:
            if pname not in kwargs:
                raise TypeError(f"Missing required parameter: {pname}")
        DrugReleaseModel.__init__(self)
        self._model_name = model_name
        self._parameters = _create_parameters(**kwargs)
        self._plot_parameters["label"] = label

    def __repr__(self):
        parts = ", ".join(
            f"{pname}={getattr(self._parameters, pname)}"
            for pname in param_specs
        )
        return f"drux.{class_name}({parts})"

    cls = type(class_name, (DrugReleaseModel,), {
        "__init__": __init__,
        "__repr__": __repr__,
        "__doc__": docstring or f"Simulator for the {label.lower()}.",
    })

    return cls
