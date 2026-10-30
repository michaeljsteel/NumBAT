"""Application-runtime services for low-level NumBAT modules.

Mesh, structure and plotting code reads output locations, tool paths and plot
preferences through this module instead of importing ``numbat``. The
application object registers itself with :func:`set_application` when it is
created, and path queries are forwarded to it so they always reflect its
current state.
"""

from pathlib import Path

import git

from numbat_core import nbversion, reporting
from numbat_core.plotting.preferences import PlotPrefs

_application = None
_plot_prefs = None
_cached_git_hash = None


def set_application(app):
    """Register the application object that answers path queries."""
    global _application
    _application = app


def _app():
    if _application is None:
        raise RuntimeError(
            'You must create a NumBAT object before calling any other NumBAT functions. '
            'See the tutorials for examples.'
        )
    return _application


def assert_application_created():
    """Exit with an error message if no application object has been created."""
    if _application is None:
        reporting.report_and_exit(
            'You must create a NumBAT object before calling any other NumBAT functions.  '
            'See the tutorials for examples.'
        )


def outdir():
    """Return the current output directory."""
    return _app().outdir()


def outprefix():
    """Return the current output file prefix."""
    return _app().outprefix()


def outpath(prefix=''):
    """Return the full output path for a prefix (or the default prefix)."""
    return _app().outpath(prefix)


def outdir_fields_path(prefix='') -> Path:
    """Return the modal-fields directory for a prefix, creating it if needed."""
    return _app().outdir_fields_path(prefix)


def path_gmsh():
    """Return the path to the gmsh executable."""
    return _app().path_gmsh()


def path_mesh_templates():
    """Return the path to the mesh templates directory."""
    return _app().path_mesh_templates()


SCRATCH_DIRNAME = 'numbat_tmp'


def path_scratch() -> Path:
    """Return the scratch directory for generated mesh files, creating it if needed.

    The directory lives in the current working directory and is shared by all
    runs made from it. Generated files are kept for inspection.
    """
    path = Path.cwd() / SCRATCH_DIRNAME
    path.mkdir(parents=True, exist_ok=True)
    return path


def is_macos():
    """Return True if running on macOS."""
    return _app().is_macos()


def plot_prefs() -> PlotPrefs:
    """Return the global plot preferences, creating defaults if none are set."""
    global _plot_prefs
    if _plot_prefs is None:
        _plot_prefs = PlotPrefs()
    return _plot_prefs


def plotfile_ext():
    """Return the plot file extension, e.g. '.png'."""
    return plot_prefs()._plot_extension


def set_plot_prefs(prefs: PlotPrefs):
    """Set the global plot preferences."""
    global _plot_prefs
    _plot_prefs = prefs


def version():
    """Return the NumBAT version string."""
    return nbversion.NUMBAT_VERSION_STR_MMM


def git_hash():
    """Get the current git commit hash (cached for the session)."""
    global _cached_git_hash
    if _cached_git_hash is None:
        try:
            repo = git.Repo(search_parent_directories=True)
            _cached_git_hash = repo.head.object.hexsha
        except (git.InvalidGitRepositoryError, git.GitCommandError):
            _cached_git_hash = 'unknown'

    return _cached_git_hash if _cached_git_hash != 'unknown' else None
