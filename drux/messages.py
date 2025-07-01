# Error messages
ERROR_DURATION_TIME_STEP_POSITIVE = "Duration and time step must be positive values"
ERROR_TIME_STEP_GREATER_THAN_DURATION = "Time step cannot be greater than duration"
ERROR_NO_SIMULATION_DATA = "No simulation data available. Run simulate() first."
ERROR_RELEASE_PROFILE_TOO_SHORT = (
    "Release profile is too short to calculate release rate."
)
ERROR_TARGET_RELEASE_RANGE = "Target release must be between 0 and 1."
ERROR_TARGET_RELEASE_EXCEEDS_MAX = (
    "Target release exceeds maximum release of the simulated duration."
)
MATPLOT_IMPORT_ERROR = "Matplotlib is required for plotting but not installed."

# Error messages for Higuchi
ERROR_INVALID_DIFFUSION = "Diffusivity (D) must be positive."
ERROR_INVALID_CONCENTRATION = "Initial drug concentration (c0) must be positive."
ERROR_INVALID_SOLUBILITY = "Solubility (cs) must be positive."
