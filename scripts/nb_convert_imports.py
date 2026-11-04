#!/usr/bin/env python3
"""Convert legacy flat NumBAT imports to the public ``numbat.*`` package.

Usage:
    nb_convert_imports.py [--write] [--no-backup] PATH [PATH ...]

PATH may be a .py or .ipynb file, or a directory (searched recursively).
Without --write the tool only prints a unified diff of what it would change.
With --write, each changed file is rewritten and the original is saved as
``<name>.orig`` unless --no-backup is given.

Only import statements are rewritten, and only simple one-module forms such as
``import materials``, ``import integration as integ`` or ``from nbtypes import
SI_GHz``. Anything else that mentions a legacy module (for example a combined
``import numbat, materials``) is reported for manual attention.
"""

import argparse
import difflib
import json
import re
import sys
from pathlib import Path

# legacy module -> (new package, new module name or None for "from X import names")
IMPORT_MODULES = {
    "materials": ("numbat", "materials"),
    "integration": ("numbat.solvers", "sbsgain"),
    "modecalcs": ("numbat.solvers", "modecalcs"),
    "modes": ("numbat.solvers", "modes"),
    "simresult": ("numbat.solvers", "simresult"),
    "structure": ("numbat.geometry", "structure"),
    "plotgain": ("numbat", "plotting"),
    "nbgmsh": ("numbat", "mesh"),
    "femmesh": ("numbat", "mesh"),
}

# legacy module -> new module for ``from legacy import names``
FROM_MODULES = {
    "nbtypes": "numbat.types",
    "nbthread": "numbat.parallel",
    "numbattools": "numbat.utils",
    "usermesh": "numbat.mesh",
    "plotgain": "numbat.plotting",
    "materials": "numbat.materials",
    "integration": "numbat.solvers",
    "modecalcs": "numbat.solvers.modecalcs",
    "structure": "numbat.geometry",
    "plotting.ipythontools": "numbat.plotting.ipython",
}

LEGACY_NAMES = set(IMPORT_MODULES) | set(FROM_MODULES) | {"numbat"}

RE_IMPORT = re.compile(r"^(?P<ind>\s*)import\s+(?P<mod>[A-Za-z_][\w]*)(?:\s+as\s+(?P<as>\w+))?\s*(?P<tail>#.*)?$")
RE_FROM = re.compile(r"^(?P<ind>\s*)from\s+(?P<mod>[\w.]+)\s+import\s+(?P<rest>.+)$")
RE_ANY_LEGACY = re.compile(
    r"^\s*(?:import|from)\s+(?:[\w.]+\s*,\s*)*(?:" + "|".join(re.escape(m) for m in sorted(LEGACY_NAMES)) + r")\b"
)


def convert_line(line):
    """Return (new_line, needs_review). ``line`` has no trailing newline."""
    m = RE_IMPORT.match(line)
    if m and m["mod"] in IMPORT_MODULES:
        pkg, name = IMPORT_MODULES[m["mod"]]
        alias = m["as"] or (m["mod"] if m["mod"] != name else None)
        tail = f"  {m['tail']}" if m["tail"] else ""
        asclause = f" as {alias}" if alias else ""
        return f"{m['ind']}from {pkg} import {name}{asclause}{tail}", False

    m = RE_FROM.match(line)
    if m and m["mod"] in FROM_MODULES:
        return f"{m['ind']}from {FROM_MODULES[m['mod']]} import {m['rest']}", False

    if RE_ANY_LEGACY.match(line) and "," in line.split("#")[0]:
        return line, True
    return line, False


def convert_source(text):
    out, review = [], []
    for n, line in enumerate(text.split("\n"), 1):
        new, flag = convert_line(line)
        if flag:
            review.append((n, line.strip()))
        out.append(new)
    return "\n".join(out), review


def convert_notebook(text):
    nb = json.loads(text)
    review = []
    for cell in nb.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        src = cell["source"]
        was_list = isinstance(src, list)
        joined = "".join(src) if was_list else src
        new, rev = convert_source(joined)
        review.extend(rev)
        if new != joined:
            cell["source"] = new.splitlines(keepends=True) if was_list else new
    indent = 1 if text.startswith('{\n "') else 2
    trailer = "\n" if text.endswith("\n") else ""
    return json.dumps(nb, indent=indent, ensure_ascii=False) + trailer, review


def process(path, write, backup):
    text = path.read_text(encoding="utf-8")
    new, review = (convert_notebook if path.suffix == ".ipynb" else convert_source)(text)
    for n, line in review:
        print(f"{path}: needs manual review: {line}" + (f" (line {n})" if path.suffix != ".ipynb" else ""))
    if new == text:
        return False
    if write:
        if backup:
            path.with_suffix(path.suffix + ".orig").write_text(text, encoding="utf-8")
        path.write_text(new, encoding="utf-8")
        print(f"converted {path}")
    else:
        sys.stdout.writelines(
            difflib.unified_diff(text.splitlines(True), new.splitlines(True), f"a/{path}", f"b/{path}")
        )
    return True


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+", type=Path)
    ap.add_argument("--write", action="store_true", help="modify files in place (default: show diff only)")
    ap.add_argument("--no-backup", action="store_true", help="do not keep .orig copies when writing")
    args = ap.parse_args(argv)

    files = []
    for p in args.paths:
        if p.is_dir():
            files += sorted(list(p.rglob("*.py")) + list(p.rglob("*.ipynb")))
        else:
            files.append(p)
    nchanged = sum(process(f, args.write, not args.no_backup) for f in files)
    print(f"{nchanged} of {len(files)} file(s) {'converted' if args.write else 'would change'}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
