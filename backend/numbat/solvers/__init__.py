"""Mode solvers, SBS gain calculations, and results."""

from . import modecalcs, modes, sbsgain, simresult
from .modecalcs import ACSimulation, EMSimulation, Simulation
from .modes import Mode, ModeAC, ModeEM
from .simresult import ACSimResult, EMSimResult, SimResult

__all__ = [
    "ACSimResult",
    "ACSimulation",
    "EMSimResult",
    "EMSimulation",
    "Mode",
    "ModeAC",
    "ModeEM",
    "SimResult",
    "Simulation",
    "integration",
    "modecalcs",
    "modes",
    "sbsgain",
    "simresult",
]

integration = sbsgain


def __getattr__(name):
    return getattr(sbsgain, name)
