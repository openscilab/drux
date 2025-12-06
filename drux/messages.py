# -*- coding: utf-8 -*-
"""Drux messages and error constants."""

# Error messages
ERROR_DURATION_TIME_STEP_POSITIVE = "Duration and time step must be positive values"
ERROR_TIME_STEP_GREATER_THAN_DURATION = "Time step cannot be greater than duration"
ERROR_NO_SIMULATION_DATA = "No simulation data available. Run simulate() first."
ERROR_RELEASE_PROFILE_TOO_SHORT = (
    "Release profile is too short to calculate release rate."
)
ERROR_TARGET_RELEASE_RANGE = "Target release must be non-negative."
ERROR_TARGET_RELEASE_EXCEEDS_MAX = (
    "Target release exceeds maximum release of the simulated duration."
)

# Error messages for Higuchi
ERROR_INVALID_DIFFUSION = "Diffusivity (D) must be positive."
ERROR_INVALID_CONCENTRATION = "Initial drug concentration (c0) must be positive."
ERROR_INVALID_SOLUBILITY = "Solubility (cs) must be positive."
ERROR_SOLUBILITY_HIGHER_THAN_CONCENTRATION = (
    "Solubility (cs) must be lower or equal to initial concentration (c0)."
)

# Error messages for zero-order
ERROR_ZERO_ORDER_RELEASE_RATE = "Release rate (k0) must be non-negative."
ERROR_ZERO_ORDER_INITIAL_AMOUNT = (
    "Initial amount of drug in the solution (M0) must be non-negative."
)

# Error messages for first-order
ERROR_FIRST_ORDER_RELEASE_RATE = "Release rate (k) must be non-negative."
ERROR_FIRST_ORDER_INITIAL_AMOUNT = (
    "Entire releasable amount of drug (M0) must be non-negative."
)

# Error messages for Weibull
ERROR_WEIBULL_SCALE_PARAMETER = "Scale parameter (a) must be positive."
ERROR_WEIBULL_SHAPE_PARAMETER = "Shape parameter (b) must be positive."
ERROR_WEIBULL_INITIAL_AMOUNT = (
    "Entire releasable amount of drug (M) must be non-negative."
)
