
import numpy as np

from nbtypes import SI_to_gmpercc


class FemMesh:
    def __init__(self):
        self.mesh_mail_fname = ""  # filename of original Gmsh-derived .mail file

        self.n_msh_pts = 0  # Number of mesh nodes (points) in the FEM calculation
                            # (Different to origianl underlying mesh because of elt sub-triangulation)
        self.n_msh_el = 0   # Number of elements in .msh mesh file

        # made by python
        self.v_el_2_mat_idx = None  # Array[0:n_msh_el] of index of elt's material into list of active em or materials.
        #    Values in range:  em simulation - [1, optical_props.n_mats_em]
        #                      ac simulation - [1, elastic_props.n_mats_ac]
        # Rename to .elt_2_material_index

        self.typ_el_AC = None  # An elastic el_conv_tbl_n  (a dictionary!), completely diff to self.type_el
        # Rename to .active_material_index
        # Belongs in ElasticProperties, not FemMesh

        # made by fortran
        self.table_nod = None  # Map of each element to its 6 nodes by node index (1..6). shape = (6, n_msh_el])
        self.node_physindex = None  # Line or surface index of a node [1..num_gmsh_types],     shape = (n_msh_pts,1)
        self.xy_nodes = None  # physical scaled x-y, locations of every node             shape= (n_msh_pts,2)

        self.n_nodes = 6  # Nodes per each element (is always 6)
        self.ac_mesh_from_em = True  # Always True

        self.el_convert_tbl = None  # Dict map of elastically active mesh elements into original em mesh elements.
        #    Length: n_msh_el(ac).   Values in [1..n_msh_el(em)]

        # Seems never used
        # self.el_convert_tbl_inv = None   # Inversion of el_convert_tbl
        # self.node_convert_tbl = None

    def build_from_gmsh_mail(self, struc):
        # Takes list of all material refractive indices
        # Discards any that are zero

        # (MJS: Not sure about the counting from 1, possibly needed for fortran?)

        self.mesh_mail_fname = struc.mesh_mail_fname
        opt_props = struc.optical_props

        mesh = struc.get_mail_mesh_data()

        # Electromagnetic so the mesh properties come straight from the Mail file.
        self.n_msh_pts = mesh.n_msh_pts  # TODO: get rid of these
        self.n_msh_el = mesh.n_msh_elts

        print(
            f"\n The final EM sim mesh has {self.n_msh_pts} nodes, {self.n_msh_el} elements and {opt_props.n_mats_em} element types (materials)."
        )

        matvals = list(struc.d_materials.values())[: opt_props.n_mats_em]

        print(" The material index table is:", opt_props.el_conv_table_n, "\n")
        print(f" There are {opt_props.n_mats_em} active materials:")
        for im, m in enumerate(matvals):
            print(
                f'  {m.material_name+",":20} n = {opt_props.v_refindexn[im]:.5f}, mat. index = {im+1}.'
            )  # +1 because materials are reported by their Fortran index

    def store_em_mode_outputs(self, type_el, node_physindex, table_nod, xy_nodes):
        self.v_el_2_mat_idx = type_el
        self.table_nod = table_nod
        self.xy_nodes = xy_nodes
        self.node_physindex = node_physindex

        #print("EM mesh properties:")
        #print("  type_el", list(self.v_el_2_mat_idx), self.v_el_2_mat_idx.shape)
        #print("  elt2nodes index map", self.table_nod, self.table_nod.shape)
        #print( #    "  node_physindex index map", self.node_physindex, self.node_physindex.shape)

    def store_ac_mode_outputs(self, type_el, table_nod, xy_nodes):
        self.v_el_2_mat_idx = type_el
        self.table_nod = table_nod
        self.xy_nodes = xy_nodes
        #print("AC after sim mesh properties:")
        #print("  type_el", list(self.v_el_2_mat_idx), self.v_el_2_mat_idx.shape)
        #print("  elt2nodes index map", self.table_nod)

    def ac_build_from_em(self, structure, em_fem):

        # Build a table of materials only containing elastic properties referencing the original full list in Struct
        # Needed for restricting meshes to elastic only for example
        # el_conv_table maps the number of the acoustically active material in original material
        # list to the position in the list of acoustically active materials
        # eg [vacuum, silicon, glass, vacuum, chalc] ->  {2:1,3:2,5:3}

        # rename  el_conv_table -> mat_conv_table

        opt_props = structure.optical_props
        el_props = structure.elastic_props

        el_props.extract_elastic_mats(structure, opt_props)

        # d_mats_AC = {}

        # Take existing msh from EM FEM and manipulate mesh to exclude vacuum areas.
        # simres_EM = self.simres_EM
        # if simres_EM:  # Invariably the case

        n_msh_el = em_fem.n_msh_el
        # type_el = em_fem.v_el_2_mat_idx       # material index of each element into list self.v_refindexn (unit-based)
        table_nod = em_fem.table_nod
        xy_nodes = em_fem.xy_nodes

        type_el_AC = []  # material index for each element (length = n_msh_el)
        table_nod_AC_tmp = np.zeros(
            np.shape(table_nod), dtype=np.int64
        )  #  fortran ordered table
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
                table_nod_AC_tmp[:, n_msh_el_AC] = table_nod[:, el]
                n_msh_el_AC += 1

        # inverse mapping using map and reversed.  Relies on map being a bijection
        # el_convert_tbl_inv = dict(map(reversed, el_convert_tbl.items()))

        # Now we have the elastic elements, need to get all the nodes encompassed by these elements
        # Each node appearing once. (To essentially make the top half of the mail file.)

        nodes_AC_mul = []  # Find every elastic node allowing multiple entries
        for el in range(n_msh_el_AC):  # for each elastic element
            # for i in range(6):        # add the absolute indices of its 6 nodes
            #    node_lst_tmp.append(table_nod_AC_tmp[i][el])
            nodes_AC_mul.extend(table_nod_AC_tmp[:6, el])

        # Now remove the multiple nodes
        nodes_AC = list(set(nodes_AC_mul))
        n_msh_pts_AC = len(nodes_AC)  # number of mesh points in acoustic grid

        # Now map the unique nodes to start from zero
        # d_nodes_2_acnodes: map from unique node index to ints 0:len(unique_nodes)
        d_nodes_2_acnodes = {}  # rename to  d_nodes_2_acnodes
        for i in range(n_msh_pts_AC):
            d_nodes_2_acnodes[nodes_AC[i]] = i

        # Creating finalised table_nod.
        # TODO: Would be nice to do this with slicing, but d_nodes_2_acnodes is a dict not a numpy array so tricky
        # but can be done by making another array long enough to hold all the unique_nodes counting from zero?
        table_nod_AC = []
        for i in range(6):
            el_tbl = []
            for el in range(n_msh_el_AC):
                el_tbl.append(d_nodes_2_acnodes[table_nod_AC_tmp[i][el]])
            table_nod_AC.append(el_tbl)

        table_nod_AC = (
            np.array(table_nod_AC) + 1
        )  # list to np array and adjust to fortran indexing

        # Find the physical x-y coordinates of the chosen AC nodes.
        xy_nodes_AC = np.zeros((2, n_msh_pts_AC))
        for node in nodes_AC:
            xy_nodes_AC[:, d_nodes_2_acnodes[node]] = xy_nodes[:, node - 1]

        # AC FEM uses Neumann B.C.s so node_physindex is totally irrelevant!
        # # Find nodes on boundaries of materials
        # node_array = -1*np.ones(n_msh_pts)
        # interface_nodes = []
        # for el in range(n_msh_el):
        #     for i in range(6):
        #         node = table_nod[i][el]
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
        self.table_nod = table_nod_AC
        self.v_el_2_mat_idx = type_el_AC
        self.xy_nodes = xy_nodes_AC
        self.node_physindex = node_physindex_AC  # TODO: Does this ever get filled?

        print(
            f"\n The final elastic sim mesh has {self.n_msh_pts} nodes, {self.n_msh_el} mesh elements and {len(el_props.typ_el_AC)} element types (materials)."
        )

        print(" The material index table is:", el_props.typ_el_AC, "\n")
        print(f" There are {el_props.n_mats_ac} active materials:")
        for mat in el_props.v_active_mats:
            print(
                f'   {mat.material_name+",":20} rho = {mat.rho*SI_to_gmpercc:.3f} g/cc.'
            )

        #print("AC mesh properties:")
        #print("  type_el", self.v_el_2_mat_idx)
        #print("  typ_el_AC", el_props.typ_el_AC)
        #print("  el_convert_tbl", self.el_convert_tbl)
        #print("  elt2nodes index map", self.table_nod)

    def get_fullmesh_nodes_xy(self):
        '''Returns vectors of x and y physical positions from the 6 nodes of each element
        in the standard order.
        Many nodes will appear multiple times because they are edges (twice) or vertex
        nodes (>=2).
        '''

        v_x6p = np.zeros(6*self.n_msh_el)
        v_y6p = np.zeros(6*self.n_msh_el)

        tabnod_py = self.table_nod.T - 1  # shift fortran to python indexing


        i = 0
        for i_el in range(self.n_msh_el):
            for i_node in range(6):

                i_ex = tabnod_py[i_el, i_node]
                #v_x6p[i] = self.xy_nodes[0, i_ex]
                #v_y6p[i] = self.xy_nodes[1, i_ex]
                v_x6p[i], v_y6p[i] = self.xy_nodes[:, i_ex]
                i += 1

        return v_x6p, v_y6p

    def make_sub_triangulation(self):
        # In table_nod
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
        for idx in range(0, 6*self.n_msh_el, 6):
            triangles = [[idx+0, idx+3, idx+5],
                         [idx+1, idx+4, idx+3],
                         [idx+2, idx+5, idx+4],
                         [idx+3, idx+4, idx+5]]
            v_triang6p.extend(triangles)

        # v_triang1p maps the same triangles to the actual nodes defined by Gmsh, of which there are only n__pt

        tabnod_py = self.table_nod.T - 1  # shift fortran to python indexing

        v_triang1p = []
        # table_nod = self.table_nod
        for i_el in np.arange(self.n_msh_el):
            triangles = [[tabnod_py[i_el, 0], tabnod_py[i_el, 3], tabnod_py[i_el, 5]],
                         [tabnod_py[i_el, 1], tabnod_py[i_el, 4], tabnod_py[i_el, 3]],
                         [tabnod_py[i_el, 2], tabnod_py[i_el, 5], tabnod_py[i_el, 4]],
                         [tabnod_py[i_el, 3], tabnod_py[i_el, 4], tabnod_py[i_el, 5]]]
            v_triang1p.extend(triangles)

        return v_triang6p, v_triang1p