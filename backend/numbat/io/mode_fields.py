"""Text serialization for sampled mode fields.

The current stable schemas are Cartesian real/imaginary 2D files and
Cartesian-plus-magnitude 1D files.
"""

import numpy as np

from numbat import runtime
from numbat.field_selection import complex_storage_keys, field_symbol, storage_key
from numbat.nbtypes import FieldComponent, FieldSelection, FieldType


def _field_type_label(field_type):
    return "AC" if FieldType.parse(field_type) == FieldType.AC else "EM"


def _write_component(mode_index, filename_prefix, selection, coordinates, fields):
    field_key = storage_key(selection)

    if selection.is_derived:
        filename = f"{filename_prefix}.txt"
        header = f"Mode {mode_index:02d}, {field_key}, {coordinates}"
        np.savetxt(filename, fields[field_key], header=header)
        return

    real_key, imaginary_key = complex_storage_keys(
        selection.field_type, selection.component
    )
    for suffix, key in (("re", real_key), ("im", imaginary_key)):
        filename = f"{filename_prefix}_{suffix}.txt"
        header = f"Mode {mode_index:02d}, {field_key}_{suffix}, {coordinates}"
        np.savetxt(filename, fields[key], header=header)


def write_mode_fields_2d(*, prefix, field_type, mode_index, coordinates, fields):
    """Write Cartesian 2D field components as real and imaginary text files."""
    v_x = coordinates["v_x"]
    v_y = coordinates["v_y"]
    coordinate_summary = (
        f"vx: [{v_x[0]:.6f},{v_x[-1]:.6f}: {len(v_x)}], "
        f"vy: [{v_y[0]:.6f},{v_y[-1]:.6f}: {len(v_y)}]"
    )
    output_directory = runtime.outdir_fields_path(prefix=prefix)
    field_label = _field_type_label(field_type)

    for component in (
        FieldComponent.X,
        FieldComponent.Y,
        FieldComponent.Z,
    ):
        selection = FieldSelection.primary(field_type, component)
        component_label = f"{field_symbol(field_type)}{component.value}"
        filename_prefix = (
            f"{output_directory}/{field_label}_mode_{mode_index:02d}_{component_label}"
        )
        _write_component(
            mode_index, filename_prefix, selection, coordinate_summary, fields
        )


def write_mode_fields_1d(
    *,
    prefix,
    field_type,
    mode_index,
    cut,
    x_axis,
    fields,
):
    """Write a 1D mode cut using the stable Cartesian-plus-magnitude schema.

    The exported derived quantity is ``Fa``. ``Fi`` remains available through
    sampled field dictionaries and plotting but is intentionally omitted here
    until a future export-schema revision.
    """
    output_directory = runtime.outdir_fields_path(prefix=prefix)
    field_label = _field_type_label(field_type)
    filename = f"{output_directory}/{field_label}_mode_{mode_index:02d}_{cut}cut.txt"

    selections = (
        FieldSelection.primary(field_type, FieldComponent.X),
        FieldSelection.secondary(field_type, FieldComponent.X),
        FieldSelection.primary(field_type, FieldComponent.Y),
        FieldSelection.secondary(field_type, FieldComponent.Y),
        FieldSelection.secondary(field_type, FieldComponent.Z),
        FieldSelection.primary(field_type, FieldComponent.Z),
        FieldSelection.magnitude(field_type),
    )
    field_keys = [storage_key(selection) for selection in selections]

    output = np.zeros([len(x_axis), len(field_keys) + 1], dtype=np.float64)
    output[:, 0] = x_axis
    for index, field_key in enumerate(field_keys, start=1):
        output[:, index] = fields[field_key]

    np.savetxt(
        filename,
        output,
        header="x " + " ".join(field_keys),
        fmt="%16.8f",
    )
