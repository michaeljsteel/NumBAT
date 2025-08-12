import numpy as np

class OpticalProps:
    """EM properties in unit-indexed forms suitable for fortran."""

    def __init__(self, d_mats, n_mats_em, loss):
        """Create OpticalProps object.

        Parameters:
            v_mats_em : list of Material structure
            n_mats_em : int - number of active materials
            loss : bool - whether to include loss in the refractive index

        """

        v_mats_em = list(d_mats.values())
        self.n_mats_em = n_mats_em         #

        matvals = v_mats_em[:n_mats_em]    # Material structure of those

        # Set up mapping tables for refractive indices
        # (Why we need el_conv_table_n is mystery)
        # el_conv_table_n maps the number of the material to the position in the nonzero v_refindexn
        # el_conv_table_n[ith material] = index into v_refindexn  of non-zero refractive indices
        # Except for zero index materials,
        #  it will always be {1:1, 2:2, 3:3, .., num_mats:num_mats}

        # trivial identity map of index to active material
        self.el_conv_table_n = {i:i for i in range(1, self.n_mats_em+1)}   #_n for refractive index n


        self.v_refindexn =np.array([m.refindex_n for m in matvals])

        if not loss:  # turn of the loss but keep the values complex for fortran passing
            self.v_refindexn = self.v_refindexn.real * complex(1.0, 0.0)




# TODO: could replace wguide with d_materials and OpticalProps
class ElasticProps:
    '''Elastic tensors in unit-indexed forms suitable for fortran'''

    def __init__(self, wguide, symmetry_flag):

        # construct list of materials with acoustic properties
        # Any material not given v_acoustic_mats assumed to be vacuum.
        self.v_acoustic_mats = [m for m in wguide.d_materials.values() if m.has_elastic_properties()]
        self.n_mats_ac = len(self.v_acoustic_mats)

        self.v_acoustic_refindexn = [m.refindex_n for m in self.v_acoustic_mats]
        self.v_acoustic_eps_eff = [m.refindex_n**2 for m in self.v_acoustic_mats]

        # density  shape = [n_mats_ac]
        self.rho = np.zeros(self.n_mats_ac)

        # stiffness tensor in 6x6 Voigt notation [3 x 3  x  n_mats_ac]
        self.c_IJ = np.zeros((6, 6, self.n_mats_ac))

        # stiffness tensor as rank 4 ijkz tensor [3x3x3  x  n_mats_ac]
        self.c_ijkz = np.zeros((3, 3, 3, self.n_mats_ac))

        # stiffness tensor as rank 4 zjkl tensor [3x3x3  x  n_mats_ac]
        self.c_zjkl = np.zeros((3, 3, 3, self.n_mats_ac))

        # photelastic tensor as rank 4 ijkl tensor  # [3x3x3x3  x  n_mats_ac]
        self.p_ijkl = np.zeros((3, 3, 3, 3, self.n_mats_ac))

        # eta tensor as rank 4 ijkl tensor [3x3x3x3  x  n_mats_ac]
        self.eta_ijkl = np.zeros((3, 3, 3, 3, self.n_mats_ac))

        self.fill_tensors(self.v_acoustic_mats, symmetry_flag)

        self._extract_elastic_mats(wguide)


    def _extract_elastic_mats(self, wguide):

        opt_props = wguide.optical_props

        el_conv_table = {}
        oldloc = 1
        newloc = 1

        #No need to examine any materials beyond the max in the EM simulation (they are all vacuum anyway)
        for mat in list(wguide.d_materials.values())[:opt_props.n_mats_em]:
            if mat.has_elastic_properties():
                el_conv_table[oldloc] = newloc
                newloc += 1
            oldloc += 1

        # Table mapping {opt_mat_index -> acoustic_mat_index}
        # ie for k, v in self.typ_al_AC:  v_acoustic_mat[v] = v_opt_mat[k], but indexed from 1
        #self.typ_el_AC = {opt_props.el_conv_table_n[k] : v for k, v in el_conv_table.items()}
        self.typ_el_AC = el_conv_table


    def is_elastic_material_index(self, idx):
        return idx in self.typ_el_AC

    def active_material_index(self, idx):
        return self.typ_el_AC[idx]

    def fill_tensors(self, v_acoustic_mats, symmetry_flag):

        # map a zero-indexed 3x3 elt to unit indexed 6x1 form.  eg x,x == 0,0 == 1
        # TODO: use a zero-indexed form of toVoigt map
        voigt_map = {(0, 0): 1, (1, 1): 2, (2, 2): 3, (2, 1): 4,
                     (2, 0): 5, (0, 1): 6, (1, 2): 4, (0, 2): 5, (1, 0): 6}


        # Build zero-based material tensors from unit-based
        for k_typ in range(self.n_mats_ac):
            if v_acoustic_mats[k_typ]:
                t_ac = v_acoustic_mats[k_typ]
                t_ac_c_IJ = t_ac.stiffness_c_IJ
                t_ac_p_IJ = t_ac.photoel_p_IJ
                t_ac_eta_IJ = t_ac.viscosity_eta_IJ

                self.rho[k_typ] = t_ac.rho

                if False and symmetry_flag:  # is it actually worth making this saving?
                    print('Surprise: using symmetry_flag tensor buildings.')
                    self.c_IJ[0, 0, k_typ] = t_ac_c_IJ[1, 1]
                    self.c_IJ[1, 1, k_typ] = t_ac_c_IJ[1, 1]
                    self.c_IJ[2, 2, k_typ] = t_ac_c_IJ[1, 1]
                    self.c_IJ[0, 1, k_typ] = t_ac_c_IJ[1, 2]
                    self.c_IJ[0, 2, k_typ] = t_ac_c_IJ[1, 2]
                    self.c_IJ[1, 0, k_typ] = t_ac_c_IJ[1, 2]
                    self.c_IJ[1, 2, k_typ] = t_ac_c_IJ[1, 2]
                    self.c_IJ[2, 0, k_typ] = t_ac_c_IJ[1, 2]
                    self.c_IJ[2, 1, k_typ] = t_ac_c_IJ[1, 2]
                    self.c_IJ[3, 3, k_typ] = t_ac_c_IJ[4, 4]
                    self.c_IJ[4, 4, k_typ] = t_ac_c_IJ[4, 4]
                    self.c_IJ[5, 5, k_typ] = t_ac_c_IJ[4, 4]

                    self.c_ijkz[2, 2, 2, k_typ] = t_ac_c_IJ[1, 1]
                    self.c_ijkz[2, 0, 0, k_typ] = t_ac_c_IJ[1, 2]
                    self.c_ijkz[2, 1, 1, k_typ] = t_ac_c_IJ[1, 2]
                    self.c_ijkz[1, 1, 2, k_typ] = t_ac_c_IJ[4, 4]
                    self.c_ijkz[1, 2, 1, k_typ] = t_ac_c_IJ[4, 4]
                    self.c_ijkz[0, 0, 2, k_typ] = t_ac_c_IJ[4, 4]
                    self.c_ijkz[0, 2, 0, k_typ] = t_ac_c_IJ[4, 4]

                    self.c_zjkl[2, 2, 2, k_typ] = t_ac_c_IJ[1, 1]
                    self.c_zjkl[2, 0, 0, k_typ] = t_ac_c_IJ[2, 1]
                    self.c_zjkl[2, 1, 1, k_typ] = t_ac_c_IJ[2, 1]
                    self.c_zjkl[1, 1, 2, k_typ] = t_ac_c_IJ[4, 4]
                    self.c_zjkl[1, 2, 1, k_typ] = t_ac_c_IJ[4, 4]
                    self.c_zjkl[0, 0, 2, k_typ] = t_ac_c_IJ[4, 4]
                    self.c_zjkl[0, 2, 0, k_typ] = t_ac_c_IJ[4, 4]


                    self.p_ijkl[0, 0, 0, 0, k_typ] = t_ac_p_IJ[1, 1]
                    self.p_ijkl[1, 1, 1, 1, k_typ] = t_ac_p_IJ[1, 1]
                    self.p_ijkl[2, 2, 2, 2, k_typ] = t_ac_p_IJ[1, 1]
                    self.p_ijkl[0, 0, 1, 1, k_typ] = t_ac_p_IJ[1, 2]
                    self.p_ijkl[0, 0, 2, 2, k_typ] = t_ac_p_IJ[1, 2]
                    self.p_ijkl[1, 1, 0, 0, k_typ] = t_ac_p_IJ[1, 2]
                    self.p_ijkl[1, 1, 2, 2, k_typ] = t_ac_p_IJ[1, 2]
                    self.p_ijkl[2, 2, 0, 0, k_typ] = t_ac_p_IJ[1, 2]
                    self.p_ijkl[2, 2, 1, 1, k_typ] = t_ac_p_IJ[1, 2]
                    self.p_ijkl[1, 2, 1, 2, k_typ] = t_ac_p_IJ[4, 4]
                    self.p_ijkl[1, 2, 2, 1, k_typ] = t_ac_p_IJ[4, 4]
                    self.p_ijkl[2, 1, 1, 2, k_typ] = t_ac_p_IJ[4, 4]
                    self.p_ijkl[2, 1, 2, 1, k_typ] = t_ac_p_IJ[4, 4]
                    self.p_ijkl[0, 2, 0, 2, k_typ] = t_ac_p_IJ[4, 4]
                    self.p_ijkl[0, 2, 2, 0, k_typ] = t_ac_p_IJ[4, 4]
                    self.p_ijkl[2, 0, 0, 2, k_typ] = t_ac_p_IJ[4, 4]
                    self.p_ijkl[2, 0, 2, 0, k_typ] = t_ac_p_IJ[4, 4]
                    self.p_ijkl[0, 1, 0, 1, k_typ] = t_ac_p_IJ[4, 4]
                    self.p_ijkl[0, 1, 1, 0, k_typ] = t_ac_p_IJ[4, 4]
                    self.p_ijkl[1, 0, 0, 1, k_typ] = t_ac_p_IJ[4, 4]
                    self.p_ijkl[1, 0, 1, 0, k_typ] = t_ac_p_IJ[4, 4]

                    self.eta_ijkl[0, 0, 0, 0, k_typ] = t_ac_eta_IJ[1, 1]
                    self.eta_ijkl[1, 1, 1, 1, k_typ] = t_ac_eta_IJ[1, 1]
                    self.eta_ijkl[2, 2, 2, 2, k_typ] = t_ac_eta_IJ[1, 1]
                    self.eta_ijkl[0, 0, 1, 1, k_typ] = t_ac_eta_IJ[1, 2]
                    self.eta_ijkl[0, 0, 2, 2, k_typ] = t_ac_eta_IJ[1, 2]
                    self.eta_ijkl[1, 1, 0, 0, k_typ] = t_ac_eta_IJ[1, 2]
                    self.eta_ijkl[1, 1, 2, 2, k_typ] = t_ac_eta_IJ[1, 2]
                    self.eta_ijkl[2, 2, 0, 0, k_typ] = t_ac_eta_IJ[1, 2]
                    self.eta_ijkl[2, 2, 1, 1, k_typ] = t_ac_eta_IJ[1, 2]
                    self.eta_ijkl[1, 2, 1, 2, k_typ] = t_ac_eta_IJ[4, 4]
                    self.eta_ijkl[1, 2, 2, 1, k_typ] = t_ac_eta_IJ[4, 4]
                    self.eta_ijkl[2, 1, 1, 2, k_typ] = t_ac_eta_IJ[4, 4]
                    self.eta_ijkl[2, 1, 2, 1, k_typ] = t_ac_eta_IJ[4, 4]
                    self.eta_ijkl[0, 2, 0, 2, k_typ] = t_ac_eta_IJ[4, 4]
                    self.eta_ijkl[0, 2, 2, 0, k_typ] = t_ac_eta_IJ[4, 4]
                    self.eta_ijkl[2, 0, 0, 2, k_typ] = t_ac_eta_IJ[4, 4]
                    self.eta_ijkl[2, 0, 2, 0, k_typ] = t_ac_eta_IJ[4, 4]
                    self.eta_ijkl[0, 1, 0, 1, k_typ] = t_ac_eta_IJ[4, 4]
                    self.eta_ijkl[0, 1, 1, 0, k_typ] = t_ac_eta_IJ[4, 4]
                    self.eta_ijkl[1, 0, 0, 1, k_typ] = t_ac_eta_IJ[4, 4]
                    self.eta_ijkl[1, 0, 1, 0, k_typ] = t_ac_eta_IJ[4, 4]

                else:
                    for i in range(6):
                        for j in range(6):
                            self.c_IJ[i, j, k_typ] = t_ac_c_IJ[i+1, j+1]  # TODO: replace with Voigt.value() ?

                    for i in [0, 1, 2]:
                        for j in [0, 1, 2]:
                            I = voigt_map[(i, j)]
                            for k in [0, 1, 2]:
                                Jz = voigt_map[(k, 2)]
                                self.c_ijkz[i, j, k, k_typ] = t_ac_c_IJ[I, Jz]
                                self.c_zjkl[k, i, j, k_typ] = t_ac_c_IJ[Jz, I]

                                for l in [0, 1, 2]:
                                    J = voigt_map[(k, l)]
                                    self.p_ijkl[i, j, k, l, k_typ] = t_ac_p_IJ[I, J]
                                    self.eta_ijkl[i, j, k, l, k_typ] = t_ac_eta_IJ[I, J]

