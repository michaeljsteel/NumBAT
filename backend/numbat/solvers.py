"""Mode solvers, SBS gain calculations, and results."""

from numbat_core.solvers import modecalcs, modes, sbsgain, simresult
from numbat_core.solvers.modecalcs import ACSimulation, EMSimulation, Simulation
from numbat_core.solvers.modes import Mode, ModeAC, ModeEM
from numbat_core.solvers.simresult import ACSimResult, EMSimResult, SimResult

__all__ = [
    "ACSimResult", "ACSimulation", "EMSimResult", "EMSimulation", "Mode",
    "ModeAC", "ModeEM", "SimResult", "Simulation", "integration", "sbsgain",
    "modecalcs", "modes", "simresult",
]

# Retain the former public module name for existing scripts and legacy tutorials.
integration = sbsgain


def __getattr__(name):
    """Expose SBS gain symbols for converted ``from integration import ...`` imports."""
    return getattr(sbsgain, name)
