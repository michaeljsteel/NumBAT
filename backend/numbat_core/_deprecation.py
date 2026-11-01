"""Errors raised by legacy import facades."""


def warn_legacy_import(module_name: str, replacement: str) -> None:
    raise ImportError(
        f"As of version 2.5.0, importing '{module_name}' is no longer supported; use '{replacement}' instead.\n"
        f"See the migration guide in the manual appendix for more information.\n\nTo avoid this error without making updates, you can revert to version 2.3.0 [git commit: 7853637] (no warnings), or version 2.4.1 [git commit: 358a559] (warnings only) by checking out the appropriate tag in the git repository."
    )
