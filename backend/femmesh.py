# Copyright (C) 2017-2025  Michael Steel, Bjorn Sturmberg, Kokou Dossou.

# NumBAT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import itertools

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import copy
from numbattools import np_min_max
from pathlib import Path


from math import sqrt

import nbgmsh
from nbtypes import SI_to_gmpercc, SI_um
import numbat
from plottools import save_and_close_figure
import reporting
import plotting






# def omake_interper_f_2d(tri_triang6p, finder, vx_out, vy_out, nx, ny):
#    return lambda femsol: matplotlib.tri.LinearTriInterpolator(
#    tri_triang6p, femsol, trifinder=finder)(vx_out, vy_out).reshape(
#        nx, ny)


# TODO: an attempt to allow matplotlib triangulation through multiprocessing
# ProcessPoolExecutor, because neither Triangulations or nested functions are
# pickleable
class TriInterpWrapper:
    def __init__(self, tri_triang6p, finder, vx_out, vy_out, nx=0, ny=0):
        self.tri_triang6p = tri_triang6p
        self.finder = finder
        self.vx_out = vx_out
        self.vy_out = vy_out
        self.nx = nx  # Nonzero if 2D context
        self.ny = ny  # Nonzero if 2D context

    def __call__(self, femsol):
        #print('\n\nCalling TriInterpWrapper:', type(self.tri_triang6p))
        #print(dir(self.tri_triang6p))
        #print(self.tri_triang6p._cpp_triangulation)

        #del self.tri_triang6p._cpp_triangulation
        #self.tri_triang6p._cpp_triangulation = None

        #mptri = matplotlib.tri.LinearTriInterpolator(
        #                                            self.tri_triang6p, femsol, trifinder=self.finder)(
        #                                                    self.vx_out, self.vy_out).reshape(self.nx, self.ny)

        mptri = matplotlib.tri.LinearTriInterpolator(
                                                    self.tri_triang6p, femsol, trifinder=self.finder)(
                                                            self.vx_out, self.vy_out)

        if self.nx:
            mptri = mptri.reshape(self.nx, self.ny)


        return mptri

def make_interper_f_2d(tri_triang6p, finder, vx_out, vy_out, nx, ny):
    """vx_out and vy_out are flattened 1D lists of x and y coords from a 2D grid of dimension nx x ny."""

#    def mif2d(femsol):
#        return matplotlib.tri.LinearTriInterpolator(
#            tri_triang6p, femsol, trifinder=finder
#        )(vx_out, vy_out).reshape(nx, ny)
#
#    return mif2d

    return TriInterpWrapper(tri_triang6p, finder, vx_out, vy_out, nx, ny)



def make_interper_f_1d(tri_triang6p, finder, vx_out, vy_out):
    """vx_out and vy_out are 1D lists of x and y coords from a 1D sampling line."""

    # def mif1d(femsol):
    #     return matplotlib.tri.LinearTriInterpolator(
    #         tri_triang6p, femsol, trifinder=finder
    #     )(vx_out, vy_out)

    # return mif1d
    return TriInterpWrapper(tri_triang6p, finder, vx_out, vy_out)




class FemMesh:
    """Mirrors the FEM mesh constructed in Fortran.

    Used by FEMScalarFieldPlotter to understand FEM grid to sample from
    Simulation (EM and AC) to communicate with Fortran side.

    The properties can be constructed eitehr by reading an existing .mail file or by
    taking the properties from the Fortran side after a simulation has been run.
    """


    def __init__(self):
        ## filename of original Gmsh-derived .mail file
        self.mesh_mail_fname = ""

        # Number of mesh nodes (points) in the FEM calculation
        # (Different to original underlying mesh because of elt sub-triangulation)
        self.n_msh_pts = 0

        # Number of elements in .msh mesh file
        self.n_msh_el = 0

        #
        # Quantities made by python
        #

        # Array[0:n_msh_el] of index of elt's material into list of active em or materials.
        #    Values in range:  em simulation - [1, optical_props.n_mats_em]
        #                      ac simulation - [1, elastic_props.n_mats_ac]
        # Rename to .elt_2_material_index
        self.v_el_2_mat_idx = None

        # An elastic el_conv_tbl_n  (a dictionary!), completely diff to self.type_el
        # Rename to .active_material_index
        # Belongs in ElasticProperties, not FemMesh
        self.typ_el_AC = None

        #
        # Quantities made by fortran
        #
        # Map of each element to its 6 nodes by node index (1..6). shape = (6, n_msh_el])
        self.elnd_to_mshpt = None

        # Line or surface index of a node [1..num_gmsh_types],     shape = (n_msh_pts,1)
        self.node_physindex = None

        # physical scaled x-y, locations of every node             shape= (n_msh_pts,2)
        self.v_nd_xy = None

        # Nodes per each element (is always 6)
        self.n_nodes = 6

        self.ac_mesh_from_em = True  # Always True

        # Dict map of elastically active mesh elements into original em mesh elements.
        #    Length: n_msh_el(ac).   Values in [1..n_msh_el(em)]
        self.el_convert_tbl = None

        # measured from nbgmsh to have them early before we know v_nd_xy
        self.extents = ( None )

        # Seems never used
        # self.el_convert_tbl_inv = None   # Inversion of el_convert_tbl
        # self.node_convert_tbl = None

    def build_from_gmsh_mail(self, wguide):
        # Takes list of all material refractive indices
        # Discards any that are zero

        # (MJS: Not sure about the counting from 1, possibly needed for fortran?)

        self.mesh_mail_fname = wguide.mesh_mail_fname
        mesh = nbgmsh.MailData(wguide.mesh_mail_fname)

        # mesh = wguide.get_mail_mesh_data()

        self.v_el_2_mat_idx = mesh.v_elts_mat

        # Electromagnetic so the mesh properties come straight from the Mail file.
        self.n_msh_pts = mesh.n_msh_pts  # TODO: get rid of these
        self.n_msh_el = mesh.n_msh_elts

        self.elnd_to_mshpt = mesh.v_elts[:, :6].T
        self.v_nd_xy = np.zeros([2, self.n_msh_pts])

        # Messy: Mail file does not include the domain scaling which is in nm, not microns
        # Everything including yvalues is scaled by domain_x, not domain_y
        self.v_nd_xy[0, :] = mesh.v_x * wguide.domain_x * 0.001
        self.v_nd_xy[1, :] = (mesh.v_y * wguide.domain_x * 0.001)

        self.extents = [
            np.min(self.v_nd_xy[0, :]),
            np.max(self.v_nd_xy[0, :]),
            np.min(self.v_nd_xy[1, :]),
            np.max(self.v_nd_xy[1, :]),
        ]

        self._report_em_mesh_properties(wguide)

    def _report_em_mesh_properties(self, wguide):
        opt_props = wguide.optical_props
        print(
            f"\n The final EM sim mesh has {self.n_msh_pts} nodes, {self.n_msh_el} elements and {opt_props.n_mats_em} element types (materials)."
        )

        #matvals = list(wguide.d_materials.values())[: opt_props.n_mats_em]
        matvals = wguide.get_optical_materials()

        print(" The material index table is:", opt_props.el_conv_table_n, "\n")
        print(f" There are {opt_props.n_mats_em} active materials:")
        for im, m in enumerate(matvals):
            # +1 because materials are reported by their Fortran index
            print(f'  {m.material_name+",":20} n = {opt_props.v_refindexn[im]:.5f}, mat. index = {im+1}.')


    def store_fortran_em_mesh_properties(self, type_el, node_physindex, elnd_to_mshpt, v_nd_xy):
        """Called after Fortran calc_modes_em() to store its mesh properties."""

        print("storing em")
        self.v_el_2_mat_idx = type_el
        self.elnd_to_mshpt = elnd_to_mshpt
        self.v_nd_xy = v_nd_xy
        self.node_physindex = node_physindex

        # print("EM mesh properties:")
        # print("  type_el", list(self.v_el_2_mat_idx), self.v_el_2_mat_idx.shape)
        # print("  elt2nodes index map", self.elnd_to_mshpt, self.elnd_to_mshpt.shape)
        # print( #    "  node_physindex index map", self.node_physindex, self.node_physindex.shape)

    def store_fortran_ac_mesh_properties(self, type_el, elnd_to_mshpt, v_nd_xy):
        """Called after Fortran calc_modes_ac() to store its mesh properties."""

        self.v_el_2_mat_idx = type_el
        self.elnd_to_mshpt = elnd_to_mshpt
        self.v_nd_xy = v_nd_xy
        # print("AC after sim mesh properties:")
        # print("  type_el", list(self.v_el_2_mat_idx), self.v_el_2_mat_idx.shape)
        # print("  elt2nodes index map", self.elnd_to_mshpt)

    def ac_build_from_em(self, wguide, em_fem):
        """Builds the FemMesh for an acoustic simulation based on an existing EM mesh
        and the known elastic properties.

        Called from ACSimulation.__init__()
        """

        # Build a table of materials only containing elastic properties referencing the original full list in Struct
        # Needed for restricting meshes to elastic only for example
        # el_conv_table maps the number of the acoustically active material in original material
        # list to the position in the list of acoustically active materials
        # eg [vacuum, silicon, glass, vacuum, chalc] ->  {2:1,3:2,5:3}

        # rename  el_conv_table -> mat_conv_table

        opt_props = wguide.optical_props
        el_props = wguide.elastic_props


        # Take existing msh from EM FEM and manipulate mesh to exclude vacuum areas.
        # simres_EM = self.simres_EM
        # if simres_EM:  # Invariably the case

        n_msh_el = em_fem.n_msh_el
        # type_el = em_fem.v_el_2_mat_idx       # material index of each element into list self.v_refindexn (unit-based)
        elnd_to_mshpt = em_fem.elnd_to_mshpt
        v_nd_xy = em_fem.v_nd_xy

        type_el_AC = []  # material index for each element (length = n_msh_el)

        #  fortran ordered table
        elnd_to_mshpt_AC_tmp = np.zeros(np.shape(elnd_to_mshpt), dtype=np.int64)

        el_convert_tbl = {}

        # Find all elements that have elastic properties (basically drop vacuum borders)
        n_msh_el_AC = 0
        for el in range(n_msh_el):
            mat_idx = em_fem.v_el_2_mat_idx[el]

            if el_props.is_elastic_material_index(mat_idx):

                type_el_AC.append(
                    el_props.active_material_index(mat_idx)
                )  # could do this and if test with try/catch

                el_convert_tbl[n_msh_el_AC] = el
                elnd_to_mshpt_AC_tmp[:, n_msh_el_AC] = elnd_to_mshpt[:, el]
                n_msh_el_AC += 1

        # inverse mapping using map and reversed.  Relies on map being a bijection
        # el_convert_tbl_inv = dict(map(reversed, el_convert_tbl.items()))

        # Now we have the elastic elements, need to get all the nodes encompassed by these elements
        # Each node appearing once. (To essentially make the top half of the mail file.)

        nodes_AC_mul = []  # Find every elastic node allowing multiple entries
        for el in range(n_msh_el_AC):  # for each elastic element
            # for i in range(6):        # add the absolute indices of its 6 nodes
            #    node_lst_tmp.append(elnd_to_mshpt_AC_tmp[i][el])
            nodes_AC_mul.extend(elnd_to_mshpt_AC_tmp[:6, el])

        # Now remove the multiple nodes
        nodes_AC = list(set(nodes_AC_mul))
        n_msh_pts_AC = len(nodes_AC)  # number of mesh points in acoustic grid

        # Now map the unique nodes to start from zero
        # d_nodes_2_acnodes: map from unique node index to ints 0:len(unique_nodes)
        d_nodes_2_acnodes = {}  # rename to  d_nodes_2_acnodes
        for i in range(n_msh_pts_AC):
            d_nodes_2_acnodes[nodes_AC[i]] = i

        # Creating finalised elnd_to_mshpt.
        # TODO: Would be nice to do this with slicing, but d_nodes_2_acnodes is a dict not a numpy array so tricky
        # but can be done by making another array long enough to hold all the unique_nodes counting from zero?
        elnd_to_mshpt_AC = []
        for i in range(6):
            el_tbl = []
            for el in range(n_msh_el_AC):
                el_tbl.append(d_nodes_2_acnodes[elnd_to_mshpt_AC_tmp[i][el]])
            elnd_to_mshpt_AC.append(el_tbl)

        elnd_to_mshpt_AC = (
            np.array(elnd_to_mshpt_AC) + 1
        )  # list to np array and adjust to fortran indexing

        # Find the physical x-y coordinates of the chosen AC nodes.
        v_nd_xy_AC = np.zeros((2, n_msh_pts_AC))
        for node in nodes_AC:
            v_nd_xy_AC[:, d_nodes_2_acnodes[node]] = v_nd_xy[:, node - 1]

        # AC FEM uses Neumann B.C.s so node_physindex is totally irrelevant!
        # # Find nodes on boundaries of materials
        # node_array = -1*np.ones(n_msh_pts)
        # interface_nodes = []
        # for el in range(n_msh_el):
        #     for i in range(6):
        #         node = elnd_to_mshpt[i][el]
        #         # Check if first time seen this node
        #         if node_array[node - 1] == -1: # adjust to python indexing
        #             node_array[node - 1] = type_el[el]
        #         else:
        #             if node_array[node - 1] != type_el[el]:
        #                 interface_nodes.append(node)
        # interface_nodes = list(set(interface_nodes))

        node_physindex_AC = np.zeros(n_msh_pts_AC)

        self.el_convert_tbl = el_convert_tbl
        # self.el_convert_tbl_inv = el_convert_tbl_inv
        # self.node_convert_tbl = node_convert_tbl

        self.n_msh_pts = n_msh_pts_AC
        self.n_msh_el = n_msh_el_AC
        self.elnd_to_mshpt = elnd_to_mshpt_AC
        self.v_el_2_mat_idx = type_el_AC
        self.v_nd_xy = v_nd_xy_AC
        self.node_physindex = node_physindex_AC  # TODO: Does this ever get filled?

        self._report_ac_mesh_properties(el_props)

    def _report_ac_mesh_properties(self, el_props):

        print(
            f"\n The final elastic sim mesh has {self.n_msh_pts} nodes, {self.n_msh_el} mesh elements and {len(el_props.typ_el_AC)} element types (materials)."
        )

        print(" The material index table is:", el_props.typ_el_AC, "\n")
        print(f" There are {el_props.n_mats_ac} active materials:")
        for mat in el_props.v_active_mats:
            vphases = mat.Vac_phase()
            print(
                f'   {mat.material_name+",":20} rho = {mat.rho*SI_to_gmpercc:.3f} g/cc, v_i = [{vphases[0]:.4f}, {vphases[1]:.4f}, {vphases[2]:.4f}] km/s.'
            )

        # print("AC mesh properties:")
        # print("  type_el", self.v_el_2_mat_idx)
        # print("  typ_el_AC", el_props.typ_el_AC)
        # print("  el_convert_tbl", self.el_convert_tbl)
        # print("  elt2nodes index map", self.elnd_to_mshpt)


    def make_sub_triangulation(self):
        # In elnd_to_mshpt
        # Nodes around a triangle element are numbered as corners: 0 1 2,  midpts: 3,4,5
        # This induces a 4-triangle sub-triangulation of each element, with clockwise vertices
        # (0 3 5), (1, 4, 3), (2, 5, 4),  (3, 4, 5)

        # v_triang6p contains an ideal list of these for a set of elements labelled 1.n_msh_el
        # The nodes for each element have their own unique labels, even when they coincide
        # with those of other neighbouring elements
        #
        # [
        #  0 + [(0 3 5), (1, 4, 3), (2, 5, 4),  (3, 4, 5)],
        #  6 + [(0 3 5), (1, 4, 3), (2, 5, 4),  (3, 4, 5)],
        #  12 + [(0 3 5), (1, 4, 3), (2, 5, 4),  (3, 4, 5)],
        #  ...
        #  6*n_msh_el + [(0 3 5), (1, 4, 3), (2, 5, 4),  (3, 4, 5)] ]

        v_triang6p = []
        for idx in range(0, 6 * self.n_msh_el, 6):
            triangles = [
                [idx + 0, idx + 3, idx + 5],
                [idx + 1, idx + 4, idx + 3],
                [idx + 2, idx + 5, idx + 4],
                [idx + 3, idx + 4, idx + 5],
            ]
            v_triang6p.extend(triangles)

        # v_triang1p maps the same triangles to the actual nodes defined by Gmsh, of which there are only n__pt

        tabnod_py = self.elnd_to_mshpt.T - 1  # shift fortran to python indexing

        v_triang1p = []
        # elnd_to_mshpt = self.elnd_to_mshpt
        for i_el in np.arange(self.n_msh_el):
            triangles = [
                [tabnod_py[i_el, 0], tabnod_py[i_el, 3], tabnod_py[i_el, 5]],
                [tabnod_py[i_el, 1], tabnod_py[i_el, 4], tabnod_py[i_el, 3]],
                [tabnod_py[i_el, 2], tabnod_py[i_el, 5], tabnod_py[i_el, 4]],
                [tabnod_py[i_el, 3], tabnod_py[i_el, 4], tabnod_py[i_el, 5]],
            ]
            v_triang1p.extend(triangles)

        return v_triang6p, v_triang1p

    def _save_triangulation_plots(self, triang1p, triang6p, v_nd_xy):
        fig, axs = plt.subplots(1, 1)
        axs[0].triplot(triang1p, linewidth=0.5)
        # axs[1].triplot(triang6p, linewidth=.5)
        # for ax in axs:
        axs[0].set_aspect(1.0)
        axs[0].scatter(v_nd_xy[0, :], v_nd_xy[1, :], s=2, c="red")

        axs[0].set_xlabel(r"$x$ [μm]")
        axs[0].set_ylabel(r"$y$ [μm]")

        pref = numbat.NumBATApp().outprefix()
        fname = pref + f"-{'ac' if self.sim_result.is_AC else 'em'}_triplots.png"
        save_and_close_figure(fig, fname)

    def make_interpolator_for_grid(self, vx_out, vy_out, nx, ny):
        """Constructs interpolator to map triangular fem grid scalar functions onto requested rectangular grid.

        Args:
            vx_out (array(float)): Flattened x-values of rectangular grid
            vy_out (array(float)): Flattened x-values of rectangular grid
            nx (int): x-dimension of rectangular grid
            ny (int): y-dimension of rectangular grid

        Returns:
            make_interpolator_for_grid: The interpolator ready for use
        """

        v_triang6p, v_triang1p = self.make_sub_triangulation()

        # This is for testing only. Normally turn off
        check_tris = False
        if check_tris:
            check_triangulation(self.v_nd_xy[0, :], self.v_nd_xy[1, :], v_triang1p)

        v_x6p, v_y6p = self.get_fullmesh_nodes_xy()

        # triangulations:  x and y coords of all points, list of triangles defined by triples of indices of the points

        # Plots show that these are equivalent meshes with different mesh point orderings
        # triang6p: tabnod_py[i_el, i_node] ordering: sequence of numbers reading out the elnd_to_mshpt
        # triang1p: tabnod_py[i_el, i_node] ordering: straight node ordering 0, 1, 2, ..5, 6+(0, 1, 2, ..5), 12+ 0, 1, 2, ..5
        tri_triang6p = matplotlib.tri.Triangulation(v_x6p, v_y6p, v_triang6p)
        tri_triang1p = matplotlib.tri.Triangulation(
            self.v_nd_xy[0, :], self.v_nd_xy[1, :], v_triang1p
        )

        draw_triangulation = False
        if draw_triangulation:
            pl_tri_triang6p = matplotlib.tri.Triangulation(
                v_x6p / SI_um, v_y6p / SI_um, v_triang6p
            )
            pl_tri_triang1p = matplotlib.tri.Triangulation(
                self.v_nd_xy[0, :] / SI_um, self.v_nd_xy[1, :] / SI_um, v_triang1p
            )

            self._save_triangulation_plots(
                pl_tri_triang1p, pl_tri_triang6p, self.v_nd_xy / SI_um
            )

        # The trifinder only works with triang1p, not triang6p.
        # The latter is apparently an 'invalid triangulation'.
        # Why?!  Perhaps it's clockwise, when anticlock is required?
        # finder = matplotlib.tri.TrapezoidMapTriFinder(tri_triang1p)
        finder = tri_triang1p.get_trifinder()
        # finder = tri_triang6p.get_trifinder()

        # The solutions we plug in are in 6p ordering so using tri_triang6p for the interperloator makes sense
        # But why is the trifinder based on 1p?  Does it make a difference?

        if ny > 1:  # a 2D sampling, not a line cut
            interper_f = make_interper_f_2d(tri_triang6p, finder, vx_out, vy_out, nx, ny)
        else:
            interper_f = make_interper_f_1d(tri_triang6p, finder, vx_out, vy_out)

        return interper_f

    def rect_grid(self, n_points):

        x_min, x_max = self.extents[0:2]
        y_min, y_max = self.extents[2:]

        area = abs((x_max - x_min) * (y_max - y_min))
        n_pts_x = int(n_points * abs(x_max - x_min) / np.sqrt(area))
        n_pts_y = int(n_points * abs(y_max - y_min) / np.sqrt(area))

        v_regx = np.linspace(x_min, x_max, n_pts_x)
        v_regy = np.linspace(y_min, y_max, n_pts_y)
        return v_regx, v_regy

    def element_to_material_index(self, i_el):
        # -1 because FEM material indices are fortran unit-indexed
        return (self.v_el_2_mat_idx[i_el] - 1)



    def get_fullmesh_nodes_xy(self):
        """Returns vectors of x and y physical positions from the 6 nodes of each element
        in the standard order.

        Many nodes will appear multiple times because they are edges (twice) or vertex
        nodes (>=2).

        Used only by the raw mode plotter in plotmoderaw.py.
        """

        v_x6p = np.zeros(6 * self.n_msh_el)
        v_y6p = np.zeros(6 * self.n_msh_el)

        tabnod_py = self.elnd_to_mshpt.T - 1  # shift fortran to python indexing

        i = 0
        for i_el in range(self.n_msh_el):
            for i_node in range(6):

                i_ex = tabnod_py[i_el, i_node]
                # v_x6p[i] = self.v_nd_xy[0, i_ex]
                # v_y6p[i] = self.v_nd_xy[1, i_ex]
                v_x6p[i], v_y6p[i] = self.v_nd_xy[:, i_ex]
                i += 1

        return v_x6p, v_y6p




class FEMScalarFieldPlotter:
    """Class to make 1D and 2D plots of a scalar quantity represented on a FEM mesh."""


    def __init__(self, wguide, n_points=500):
        """Build objects to hold a scalar field represented on a FEM mesh.

        The FEM mesh is constructed from the mesh file mesh_mail_fname.
        n_points is the nominal linear resolution of the rectangular grid output,
        partitioned so that x-y pixels are roughly square.
        That is, in the final grid, nx * ny ~= n_points^2, and dx ~= dy.

        """

        self.fem_mesh = FemMesh()
        self.fem_mesh.build_from_gmsh_mail(wguide)

        self.quantity_name = ""
        self.file_suffix = ""
        self.n_points = n_points

        self.scalar_fields = []

        # get approx square pixel rectangular grid
        self.v_x, self.v_y = self.fem_mesh.rect_grid(n_points)

        # m_x, m_y = np.meshgrid(self.v_x, self.v_y)
        # v_x_flat = m_x.flatten('F')
        # v_y_flat = m_y.flatten('F')
        # self.interper = self.fem_mesh.make_interpolator_for_grid(v_x_flat, v_y_flat, len(self.v_x), len(self.v_y))

        self.x_lab = r"$x$ [μm]"
        self.y_lab = r"$y$ [μm]"
        self.d_lab = r"$d$ [μm]"

    def setup_scalar_properties(self, nm_eng, unit, nm_math, fname_suffix):
        self.dim = 1
        self.nm_eng = nm_eng
        self.nm_math = nm_math
        self.unit = unit
        self.fname_suffix = fname_suffix

    def setup_vector_properties(
        self,
        dim,
        nm_eng,
        unit,
        nm_math,
        nm_math_comps,
        fname_suffix,
        fname_suffix_comps,
    ):
        self.dim = dim
        self.nm_eng = nm_eng
        self.nm_math = nm_math
        self.unit = unit
        self.nm_math_comps = nm_math_comps
        self.fname_suffix = fname_suffix
        self.fname_suffix_comps = fname_suffix_comps

    def set_quantity_name(self, nm, suf):
        self.quantity_name = nm
        self.file_suffix = suf

    def n_elts(self):  # NAME PROBLEMATIC
        return self.fem_mesh.n_msh_el

    def fill_quantity_by_material_index(self, mati_to_scalars):
        """Takes an array mati_to_scalar that maps the index of a material to the desired scalar property, eg refractive index.

        mati_to_scalar is indexed from 0 with 0 being the vacuum material
        """

        # make both scalar and vector inputs usable
        if self.dim > 1:
            mati_to_vals = mati_to_scalars
        else:
            mati_to_vals = np.expand_dims(copy.deepcopy(mati_to_scalars), axis=1)

        n_elts = self.fem_mesh.n_msh_el
        self.scalar_fields = []
        for i in range(self.dim):
            t_scalar_field = np.zeros([6, n_elts], dtype=np.float64)

            for i_el in range(n_elts):
                matel = self.element_to_material_index(i_el)
                t_scalar_field[:, i_el] = np.real(mati_to_vals[matel, i])

            t_scalar_field = t_scalar_field.flatten("F")
            self.scalar_fields.append(t_scalar_field)

    def element_to_material_index(self, i_el):
        return self.fem_mesh.element_to_material_index(i_el)

    def make_plot_1D_cut(self, pref, s_cut, val1, val2=None):
        """Make a 1D plot of the scalar along a line.

        s_cut is a string determining the cut direction, with allowed values: 'x', 'y', or 'line'.

        For s_cut='x', the function is plotted along y at fixed x=val1, where val1 is a float.

        For s_cut='y', the function is plotted along x at fixed x=val1, where val1 is a float.

        For s_cut='line', the function is plotted along a straight line from val1 to val2, where the latter are float tuples (x0,y0) and (x1,y1).


        """

        # full out-directory parth
        prefix = str(Path(numbat.NumBATApp().outdir(), pref))

        if s_cut.lower() == "x":
            self._make_plot_xcut(prefix, val1)
        elif s_cut.lower() == "y":
            self._make_plot_ycut(prefix, val1)
        else:
            self._make_plot_1D(prefix, val1, val2)

    def _make_plot_ycut(self, prefix, y0):
        v_x_flat = self.v_x
        v_y_flat = y0 + np.zeros(len(self.v_x))
        self.interper = self.fem_mesh.make_interpolator_for_grid(
            v_x_flat, v_y_flat, len(self.v_x), 1
        )

        ## TODO: explain the need for the transpose in the imshow call below
        # epslo, epshi = np_min_max(m_scalar)

        fig, ax = plt.subplots()
        ax.set_xlabel(self.x_lab)
        ylab = self.nm_eng + " " + self.nm_math + " " + self.unit
        ax.set_ylabel(ylab)

        for i in range(self.dim):
            v_scalar = self.interper(self.scalar_fields[i])
            ax.plot(self.v_x, v_scalar)

        plotting.save_and_close_figure(
            fig, prefix + "-" + self.fname_suffix + "_xcut.png"
        )

    def _make_plot_xcut(self, prefix, x0):
        v_y_flat = self.v_y
        v_x_flat = x0 + np.zeros(len(self.v_y))

        self.interper = self.fem_mesh.make_interpolator_for_grid(
            v_x_flat, v_y_flat, len(self.v_y), 1
        )

        fig, ax = plt.subplots()
        ax.set_xlabel(self.y_lab)

        ylab = self.nm_eng + " " + self.nm_math + " " + self.unit
        ax.set_ylabel(ylab)

        for i in range(self.dim):
            v_scalar = self.interper(self.scalar_fields[i])
            ax.plot(self.v_y, v_scalar)

        plotting.save_and_close_figure(
            fig, prefix + "-" + self.fname_suffix + "_ycut.png"
        )

    def _make_plot_1D(self, prefix, pt0, pt1):

        x0, y0 = pt0
        x1, y1 = pt1
        v_x_flat = np.linspace(x0, x1, self.n_points)
        v_y_flat = np.linspace(y0, y1, self.n_points)
        v_d = np.sqrt((v_x_flat - v_x_flat[0]) ** 2 + (v_y_flat - v_y_flat[0]) ** 2)

        self.interper = self.fem_mesh.make_interpolator_for_grid(
            v_x_flat, v_y_flat, len(v_x_flat), 1
        )

        fig, ax = plt.subplots()
        ax.set_xlabel(self.d_lab)
        ylab = self.nm_eng + " " + self.nm_math + " " + self.unit
        ax.set_ylabel(ylab)

        for i in range(self.dim):
            v_scalar = self.interper(self.scalar_fields[i])
            ax.plot(v_d, v_scalar)

        plotting.save_and_close_figure(
            fig, prefix + "-" + self.fname_suffix + "_linecut.png"
        )

    def make_plot_2D(self, pref, aspect=1.0, with_cb=True):

        # full out-directory path
        prefix = str(Path(numbat.NumBATApp().outdir(), pref))

        m_x, m_y = np.meshgrid(self.v_x, self.v_y)
        v_x_flat = m_x.flatten("F")
        v_y_flat = m_y.flatten("F")
        self.interper = self.fem_mesh.make_interpolator_for_grid(
            v_x_flat, v_y_flat, len(self.v_x), len(self.v_y)
        )

        for i in range(self.dim):
            m_scalar = self.interper(self.scalar_fields[i])

            # TODO: explain the need for the transpose in the imshow call below
            epslo, epshi = np_min_max(m_scalar)

            fig, ax = plt.subplots()
            cmap = "cool"
            im = ax.imshow(
                m_scalar.T,
                cmap=cmap,
                vmin=1.0,
                vmax=epshi,
                origin="lower",
                extent=[
                    np.min(self.v_x),
                    np.max(self.v_x),
                    np.min(self.v_y),
                    np.max(self.v_y),
                ],
            )
            # cf=ax.contourf(m_regx, m_regy, v_regindex, cmap=cmap, vmin=1.0, vmax=np.nanmax(v_regindex))
            ax.set_xlabel(self.x_lab)
            ax.set_ylabel(self.y_lab)
            ax.set_aspect(aspect)
            im.set_clim(1, np.nanmax(m_scalar))
            if with_cb:
                ticks = np.linspace(epslo, epshi, 5)
                fmts = (
                    (f"{epslo:.4f}",),
                    map(lambda x: f"{x:.1f}", ticks[1:-1]),
                    (f"{epshi:.4f}",),
                )
                fmts = list(itertools.chain.from_iterable(fmts))

                fmt = mticker.FixedFormatter(fmts)
                cax = ax.inset_axes([1.04, 0.1, 0.03, 0.8])
                cb = fig.colorbar(im, cax=cax, ticks=ticks, format=fmt)

                if self.dim == 1:
                    fname = prefix + "-" + self.fname_suffix + ".png"
                    qlab = self.nm_eng + " " + self.nm_math + " " + self.unit
                else:
                    fname = (
                        prefix
                        + "-"
                        + self.fname_suffix
                        + f"_{self.fname_suffix_comps[i]}.png"
                    )
                    qlab = self.nm_eng + " " + self.nm_math_comps[i] + " " + self.unit

                fig.suptitle(qlab)
                cb.ax.tick_params(labelsize=12)
                cb.outline.set_linewidth(0.5)
                cb.outline.set_color("gray")

            plotting.save_and_close_figure(fig, fname)












# Checks of mesh and triangles satisfy conditions for triangulation
# Quadratic algorithm. Use on the smallest grid possible
def check_triangulation(vx, vy, triangs):
    # are points unique
    print("\n\nChecking triangulation goodness")
    npts = len(vx)
    dsepmin = 1e6
    dsi = 0
    dsj = 0
    for i in range(npts):
        for j in range(i + 1, npts):
            dsep = sqrt((vx[i] - vx[j]) ** 2 + (vy[i] - vy[j]) ** 2)
            if dsep < dsepmin:
                dsepmin = dsep
                dsi = i
                dsj = j

    print("  Closest space of triangle points was", dsepmin)
    if dsepmin < 1e-11:
        msg = f"Point collision at {dsi}, {dsj}: ({vx[dsi]},{vy[dsi]}) =  ({vx[dsj]},{vy[dsj]})."
        msg += "\nIt seems the mesh grid reordering has failed."
        reporting.report_and_exit(msg)

    # is list of triangles unique
    s_vtri = set()
    clean = True
    for tri in triangs:
        stri = str(tri)
        if stri in s_vtri:
            print("        Double triangle at", stri)
            clean = False
        else:
            s_vtri.add(stri)
    if clean:
        print("  No doubled triangles found")
    else:
        print("  Found doubled triangles")