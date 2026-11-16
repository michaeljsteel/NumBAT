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
