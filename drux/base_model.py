"""This file contains the abstract base class for drug release models."""

import numpy as np
import matplotlib.pyplot as plt
from abc import ABC
from typing import Any, Optional

from .messages import (
    ERROR_TARGET_RELEASE_RANGE,
    ERROR_DURATION_TIME_STEP_POSITIVE,
    ERROR_TIME_STEP_GREATER_THAN_DURATION,
    ERROR_NO_SIMULATION_DATA,
    ERROR_RELEASE_PROFILE_TOO_SHORT,
    ERROR_TARGET_RELEASE_EXCEEDS_MAX,
)

from .params import MODELS_REGISTRY


class DrugReleaseModel(ABC):
    """
    Abstract base class for drug release models.

    This class provides a common interface and functionality for various
    mathematical models of drug release from delivery systems.

    Subclasses should implement:
    - _model_function(): Core model equation
    - _validate_parameters(): Parameter validation
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
        """
        Validate model parameters using MODELS_REGISTRY.

        Raises ValueError if parameters are invalid.
        """
        if self._model_name is None or self._parameters is None:
            raise ValueError("Model name and parameters must be set")

        config = MODELS_REGISTRY[self._model_name]

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

        config = MODELS_REGISTRY[self._model_name]
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
