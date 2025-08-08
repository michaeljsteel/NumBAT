import sys
sys.path.append('../../backend')
import numpy as np
import scipy.linalg as spla

import materials

#mat_LiNbO3 = materials.make_material("LiNbO3_1973_Auld")
mat_LiNbO3 = materials.make_material("LiNbO3_2023_Rodrigues")

mat_LiNbO3_y = mat_LiNbO3.copy()
#mat_LiNbO3_y.set_orientation('y-cut')
mat_LiNbO3_x = mat_LiNbO3.copy()
mat_LiNbO3_x.set_orientation('negx-cut')
mat_LiNbO3_z = mat_LiNbO3.copy()
mat_LiNbO3_z.set_orientation('z-cut')

print(mat_LiNbO3_y.full_str())
mat_LiNbO3_x.enable_piezoelectric_effects()
mat_LiNbO3_y.enable_piezoelectric_effects()
mat_LiNbO3_z.enable_piezoelectric_effects()

sys.exit(0)
with np.printoptions(precision=3, suppress=True):

    stiffc_IJy1=mat_LiNbO3_y.stiffness_c_IJ
    stiffc_IJz1=mat_LiNbO3_z.stiffness_c_IJ

    cE_IJy = mat_LiNbO3_y._piezo._tens_cE_IJ
    cE_IJz = mat_LiNbO3_z._piezo._tens_cE_IJ

    diJy=mat_LiNbO3_y._piezo._tens_d_iJ
    diJz=mat_LiNbO3_z._piezo._tens_d_iJ

    eiJy=mat_LiNbO3_y._piezo._tens_e_iJ
    eiJz=mat_LiNbO3_z._piezo._tens_e_iJ

    epiJy=mat_LiNbO3_y._piezo._tens_relepsS_ij
    epiJz=mat_LiNbO3_z._piezo._tens_relepsS_ij



    print('Raw stiffnesses')
    print(stiffc_IJy1, stiffc_IJz1)

    print('\n\nPiezo stiffnesses')
    print(cE_IJy, cE_IJz)

    print('\n\npiezo d')
    print(diJy, diJz)

    print('\n\npiezo e')
    print(eiJy, eiJz)

    print('\n\nepsS')
    print(epiJy, epiJz)

    #kap1 = (1,0,0)
    #kap2 = (1,0,0)
    kap1 = (0,1,0)
    kap2 = (0,0,-1)

    stiffc_IJy2=mat_LiNbO3_y.get_stiffness_for_kappa(kap1)
    stiffc_IJz2=mat_LiNbO3_z.get_stiffness_for_kappa(kap2)



    print('\n\nstiffs y')
    print(stiffc_IJy1, stiffc_IJy2)
    print('eigs1\n',spla.eig(stiffc_IJy1.mat)[0])
    print('eigs2\n',spla.eig(stiffc_IJy2.mat)[0])

    print('\n\nstiffs z')
    print(stiffc_IJz1, stiffc_IJz2)
    print('eigs1\n',spla.eig(stiffc_IJz1.mat)[0])
    print('eigs2\n',spla.eig(stiffc_IJz2.mat)[0])

