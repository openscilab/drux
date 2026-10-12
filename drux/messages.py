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
ERROR_RELEASABLE_AMOUNT = (
    "Entire releasable amount of drug (M) must be non-negative."
)

# Hopfenberg model error messages
ERROR_INVALID_EROSION_CONSTANT = "Erosion rate constant (k0) must be non-negative."
ERROR_INVALID_INITIAL_RADIUS = "Initial radius or half-thickness (a0) must be positive."
ERROR_INVALID_GEOMETRY_FACTOR = "Geometry factor (n) must be 1 (slab), 2 (cylinder), or 3 (sphere)."

# Error messages for curve fitting
ERROR_UNKNOWN_MODEL = "Unknown model '{model_name}'. Available models: {available_models}."
ERROR_TIME_RELEASE_LENGTH_MISMATCH = "Time and release profile must have the same length."
ERROR_INSUFFICIENT_DATA_POINTS = "At least 2 data points are required for fitting."
ERROR_UNKNOWN_FIT_PARAMETER = "Unknown parameter(s) {unknown_keys} for model '{model_name}'. Valid parameters: {parameter_names}."
ERROR_UNKNOWN_FIT_ARGUMENT = (
    "Unknown parameter(s) {unknown_keys} in '{argument_name}' for model '{model_name}'. "
    "Free parameters: {free_parameters}."
)
ERROR_KNOWN_FIT_ARGUMENT = "Parameter(s) {known_keys} in '{argument_name}' are known, so they are not fitted."
ERROR_MISSING_FIT_ARGUMENT = (
    "Missing parameter(s) {missing_keys} in '{argument_name}'. "
    "Give a value for every free parameter: {free_parameters}."
)
ERROR_FIT_ARGUMENT_TYPE = "'{argument_name}' must be a dictionary that maps parameter names to values."
ERROR_INVALID_FIT_VALUE = "Value of parameter '{parameter_name}' in '{argument_name}' must be a real number."
ERROR_INVALID_FIT_BOUNDS = "Bounds of parameter '{parameter_name}' must be a (lower, upper) pair of real numbers."
ERROR_FIT_BOUNDS_ORDER = "Lower bound of parameter '{parameter_name}' must be less than its upper bound."
ERROR_FIT_GUESS_OUT_OF_BOUNDS = "Initial guess of parameter '{parameter_name}' must be between its bounds."
ERROR_NO_FREE_PARAMETERS = "All parameters are known; there is nothing left to fit."
ERROR_NO_FIT_RESULT = "No fit result available. Run fit() first."
