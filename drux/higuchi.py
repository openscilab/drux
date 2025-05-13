from drux import DrugReleaseModel, ModelParameters
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

    def __init__(self, params: HiguchiParameters):
        super().__init__(params)

    def _model_function(self, t: float) -> float:
        """
        Calculate the drug release at time t using the Higuchi model.
        - Case 1: c0 < cs → Mt = sqrt(D * t / (pi * L^2))
        - Case 2: c0 >> cs → Mt = sqrt(D * c0 * cs) * t
        - Case 3: General case (default): Mt = sqrt(D * c0 * (2*c0 - cs) * cs * t)
        :param t: time (s)
        """
        D = self.params.D
        c0 = self.params.c0
        cs = self.params.cs
        L = self.params.L


        if c0 < cs:
            Mt = sqrt((D * t) / (pi * L**2))
        elif c0 > 10 * cs:
            Mt = sqrt(D * c0 * cs) * t
        else:
            Mt = sqrt(D * c0 * (2 * c0 - cs) * cs) * t

        return Mt

    def _validate_parameters(self) -> None:
        """
        Validate the parameters of the Higuchi model.
        """
        if self.params.D <= 0:
            raise ValueError("Diffusivity (D) must be positive.")
        if self.params.c0 <= 0:
            raise ValueError("Initial drug concentration (c0) must be positive.")
        if self.params.cs <= 0:
            raise ValueError("Solubility (cs) must be positive.")
        if self.params.L <= 0:
            raise ValueError("Film thickness (L) must be positive.")