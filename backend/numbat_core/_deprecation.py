"""Helpers for warnings emitted by legacy import facades."""

import warnings


def warn_legacy_import(module_name: str, replacement: str) -> None:
    warnings.warn(
        f"Importing '{module_name}' is deprecated; use '{replacement}' instead.\nYou can avoid these warnings by reverting to version 2.3.0 [git commit 7853637...] or earlier, or by updating your code to use the new import paths.\nThis warning will become an error in a future release.",
        FutureWarning,
        stacklevel=3,
    )
