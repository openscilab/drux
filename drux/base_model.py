"""This file contains the abstract base class for drug release models."""

import numpy as np
from abc import ABC, abstractmethod
from typing import Any

from .messages import (
    ERROR_TARGET_RELEASE_RANGE,
    ERROR_DURATION_TIME_STEP_POSITIVE,
    ERROR_TIME_STEP_GREATER_THAN_DURATION,
    ERROR_NO_SIMULATION_DATA,
    ERROR_RELEASE_PROFILE_TOO_SHORT,
    ERROR_TARGET_RELEASE_EXCEEDS_MAX,
    MATPLOT_IMPORT_ERROR,
)


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
        self.time_points = None
        self.release_profile = None

    @abstractmethod
    def _validate_parameters(self) -> None:
        """
        Validate model parameters.

        Should raise ValueError if parameters are invalid.
        """
        pass

    @abstractmethod
    def _model_function(self, t: float) -> float:
        """
        Model function that calculates drug release profile over time.

        :param t: time point at which to calculate drug release
        """
        pass

    def _get_release_profile(self) -> np.ndarray:
        """Calculate the drug release profile over the specified time points."""
        return np.vectorize(self._model_function)(self.time_points)

    def _validate_plot(self) -> tuple:
        """
        Validate plotting process.

        :raises ImportError: if matplotlib is not installed
        :raises ValueError: if simulation data is not available
        :raises ValueError: if release profile is too short
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError(MATPLOT_IMPORT_ERROR)

        if self.time_points is None or self.release_profile is None:
            raise ValueError(ERROR_NO_SIMULATION_DATA)

        if len(self.release_profile) < 2:
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
        self.time_points = np.arange(0, duration + time_step, time_step)
        self._validate_parameters()
        self.release_profile = self._get_release_profile()
        return self.release_profile

    def plot(self, show: bool = True, **kwargs: Any) -> tuple:
        """
        Plot the drug release profile.

        :param show: Whether to display the plot (default: True)
        """
        # Create a new figure and axis if not provided
        fig, ax = self._validate_plot()

        # Plotting the release profile
        ax.plot(
            self.time_points, self.release_profile, label="Release Profile", **kwargs
        )
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Cumulative Release")
        ax.set_title("Drug Release Profile")
        ax.grid()
        ax.legend()

        # Show the plot if requested
        if show:
            fig.show()

        return fig, ax

    def get_release_rate(self) -> np.ndarray:
        """Calculate the instantaneous release rate (derivative of release profile)."""
        if self.time_points is None or self.release_profile is None:
            raise ValueError(ERROR_NO_SIMULATION_DATA)

        if len(self.release_profile) < 2:
            raise ValueError(ERROR_RELEASE_PROFILE_TOO_SHORT)

        # Calculate the derivative of the release profile
        release_rate = np.gradient(self.release_profile, self.time_points)
        return release_rate

    def time_for_release(self, target_release: float) -> float:
        """
        Estimate time needed to reach a specific release percentage.

        :param target_release: target release fraction (0-1)

        :raises ValueError: if target_release is not between 0 and 1
        :raises ValueError: if simulation data is not available
        :raises ValueError: if target_release exceeds maximum release
        """
        if self.time_points is None or self.release_profile is None:
            raise ValueError(ERROR_NO_SIMULATION_DATA)

        if target_release < 0 or target_release > 1:
            raise ValueError(ERROR_TARGET_RELEASE_RANGE)

        if target_release > self.release_profile[-1]:
            raise ValueError(ERROR_TARGET_RELEASE_EXCEEDS_MAX)

        # Find first time point where release >= target
        idx = np.argmax(self.release_profile >= target_release)
        return self.time_points[idx]
