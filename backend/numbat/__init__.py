"""Public NumBAT package."""

from .api import application as _application
from .api.application import *


def __getattr__(name):
    try:
        return getattr(_application, name)
    except AttributeError:
        raise AttributeError(f"module 'numbat' has no attribute {name!r}") from None


def __dir__():
    return sorted(set(globals()) | set(dir(_application)))
