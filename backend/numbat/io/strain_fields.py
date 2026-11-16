"""Text serialization for sampled acoustic strain fields."""

import numpy as np

from numbat.runtime import outdir_fields_path
from numbat.strain_selection import StrainComponent, complex_storage_keys


def write_strain_fields_2d(*, prefix, mode_index, coordinates, fields):
    """Write all six strain components as real and imaginary 2D text files.

    Files use the stable name ``AC_strain_mode_{mode:02d}_S{ij}_{re|im}.txt``.
    Their matrices have the same ``(x, y)`` indexing as the sampled plot grid.
    """
    v_x = coordinates["v_x"]
    v_y = coordinates["v_y"]
    coordinate_summary = (
        f"vx: [{v_x[0]:.6f},{v_x[-1]:.6f}: {len(v_x)}], "
        f"vy: [{v_y[0]:.6f},{v_y[-1]:.6f}: {len(v_y)}]"
    )
    output_directory = outdir_fields_path(prefix=prefix)

    for component in StrainComponent:
        real_key, imaginary_key = complex_storage_keys(component)
        component_label = f"S{component.value}"
        filename_prefix = (
            f"{output_directory}/AC_strain_mode_{mode_index:02d}_{component_label}"
        )
        for suffix, key in (("re", real_key), ("im", imaginary_key)):
            np.savetxt(
                f"{filename_prefix}_{suffix}.txt",
                fields[key],
                header=f"Mode {mode_index:02d}, {key}, {coordinate_summary}",
            )


def write_strain_fields_1d(*, prefix, mode_index, cut, x_axis, fields):
    """Write all physical strain parts for one sampled cut in a single file."""
    output_directory = outdir_fields_path(prefix=prefix)
    filename = f"{output_directory}/AC_strain_mode_{mode_index:02d}_{cut}cut.txt"
    field_keys = [
        key
        for component in StrainComponent
        for key in complex_storage_keys(component)
    ]
    output = np.empty((len(x_axis), len(field_keys) + 1), dtype=np.float64)
    output[:, 0] = x_axis
    for index, key in enumerate(field_keys, start=1):
        output[:, index] = fields[key]

    np.savetxt(
        filename,
        output,
        header="x " + " ".join(field_keys),
        fmt="%16.8f",
    )
