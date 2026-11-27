"""Shared verbosity-controlled output for NumBAT's Python packages."""

from enum import Enum


class Verbosity(Enum):
    """Controls which NumBAT messages are displayed."""

    QUIET = "quiet"
    NORMAL = "normal"
    DETAILED = "detailed"


class OutputLevel(Enum):
    """Classifies a NumBAT message independently of its severity."""

    NORMAL = "normal"
    DETAILED = "detailed"
    ALWAYS = "always"


_verbosity = Verbosity.NORMAL
_verbosity_levels = {Verbosity.QUIET: 0, Verbosity.NORMAL: 1, Verbosity.DETAILED: 2}
_output_levels = {OutputLevel.NORMAL: 1, OutputLevel.DETAILED: 2}


def set_verbosity(verbosity: Verbosity | str) -> None:
    """Set the verbosity of NumBAT's Python output using an enum or string."""
    try:
        selected = Verbosity(verbosity)
    except ValueError:
        raise ValueError("verbosity must be 'quiet', 'normal', or 'detailed'") from None
    global _verbosity
    _verbosity = selected


def get_verbosity() -> str:
    """Return the current NumBAT output verbosity as a string."""
    return _verbosity.value


def nbprint(*args, sep=" ", end="\n", file=None, flush=False,
            level: OutputLevel | str = OutputLevel.NORMAL):
    """Print a NumBAT message if its level is enabled.

    ``level`` is an :class:`OutputLevel` or its string value. ``always`` is
    for warnings, errors, and explicitly requested displays.
    """
    try:
        selected = OutputLevel(level)
    except ValueError:
        raise ValueError("output level must be 'normal', 'detailed', or 'always'") from None
    if selected is OutputLevel.ALWAYS or _verbosity_levels[_verbosity] >= _output_levels[selected]:
        print(*args, sep=sep, end=end, file=file, flush=flush)
