"""Gmsh mesh-file readers used by NumBAT mesh consumers."""

from pathlib import Path

import numpy as np


def read_msh2_mesh(mesh_path: Path) -> tuple[np.ndarray, list[tuple[tuple[int, int, int], int]], set[tuple[int, int]]]:
    """Read nodes, triangular elements, and line elements from a Gmsh 2.2 mesh."""
    lines = mesh_path.read_text().splitlines()
    node_coordinates = {}
    triangles = []
    line_elements = set()
    line_index = 0

    while line_index < len(lines):
        if lines[line_index] == "$Nodes":
            node_count = int(lines[line_index + 1])
            for node_line in lines[line_index + 2 : line_index + 2 + node_count]:
                node_id, x_coord, y_coord, _ = node_line.split()
                node_coordinates[int(node_id)] = (float(x_coord), float(y_coord))
            line_index += node_count + 3
        elif lines[line_index] == "$Elements":
            element_count = int(lines[line_index + 1])
            for element_line in lines[line_index + 2 : line_index + 2 + element_count]:
                fields = element_line.split()
                element_type = int(fields[1])
                tag_count = int(fields[2])
                node_ids = [int(node_id) for node_id in fields[3 + tag_count :]]
                physical_id = int(fields[3]) if tag_count else 0

                if element_type == 1:
                    line_elements.add(tuple(sorted(node_ids[:2])))
                elif element_type in (2, 9):
                    triangles.append((tuple(node_ids[:3]), physical_id))
            line_index += element_count + 3
        else:
            line_index += 1

    if not node_coordinates or not triangles:
        raise ValueError(f"Mesh file does not contain Gmsh 2.2 nodes and triangular elements: {mesh_path}")

    ordered_node_ids = sorted(node_coordinates)
    node_id_to_index = {node_id: index for index, node_id in enumerate(ordered_node_ids)}
    coordinates = np.array([node_coordinates[node_id] for node_id in ordered_node_ids])
    indexed_triangles = [
        (tuple(node_id_to_index[node_id] for node_id in node_ids), physical_id)
        for node_ids, physical_id in triangles
    ]
    indexed_line_elements = {
        tuple(node_id_to_index[node_id] for node_id in node_ids) for node_ids in line_elements
    }
    return coordinates, indexed_triangles, indexed_line_elements


