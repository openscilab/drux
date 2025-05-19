from drux import DrugReleaseModel, ModelParameters
from drux.messages import (
    ERROR_INVALID_DIFFUSION,
    ERROR_INVALID_CONCENTRATION,
    ERROR_INVALID_THICKNESS,
    ERROR_INVALID_SOLUBILITY,
)
from dataclasses import dataclass
from math import sqrt, pi


@dataclass
class HiguchiParameters(ModelParameters):
    """
    Parameters for the Higuchi model based on physical formulation.

    Attributes:
        D (float): Drug diffusivity in the polymer carrier (cm^2/s)
        c0 (float): Initial drug concentration (mg/cm^3)
        cs (float): Drug solubility in the polymer (mg/cm^3)
        L (float): Film thickness (cm)
    """
    D: float
    c0: float
    cs: float
    L: float

class HiguchiModel(DrugReleaseModel):
    """
    Simulator for the Higuchi drug release model using analytical expressions
    based on concentration conditions:
    """

    def __init__(self, D: float, c0: float, cs: float, L: float) -> None:
        """
        Initialize the Higuchi model with the given parameters.
        :param D: Drug diffusivity in the polymer carrier (cm^2/s)
        :param c0: Initial drug concentration (mg/cm^3)
        :param cs: Drug solubility in the polymer (mg/cm^3)
        :param L: Film thickness (cm)
        """
        super().__init__()
        self.params = HiguchiParameters(D=D, c0=c0, cs=cs, L=L)

    def _model_function(self, t: float) -> float:
        """
        Calculate the drug release at time t using the Higuchi model.
        - Case 1: c0 < cs → Mt = sqrt(D * t / (pi * L^2))
        - Case 2: General case (default): Mt = sqrt(D * c0 * (2*c0 - cs) * cs * t)
        :param t: time (s)
        """
        D = self.params.D
        c0 = self.params.c0
        cs = self.params.cs
        L = self.params.L


        if c0 < cs:
            Mt = sqrt((D * t) / (pi * L**2))
        else:
            Mt = sqrt(D * c0 * (2 * c0 - cs) * cs) * t

        return Mt

    def _validate_parameters(self) -> None:
        """
        Validate the parameters of the Higuchi model.
        """
        if self.params.D <= 0:
            raise ValueError(ERROR_INVALID_DIFFUSION)
        if self.params.c0 <= 0:
            raise ValueError(ERROR_INVALID_CONCENTRATION)
        if self.params.cs <= 0:
            raise ValueError(ERROR_INVALID_SOLUBILITY)
        if self.params.L <= 0:
            raise ValueError(ERROR_INVALID_THICKNESS)