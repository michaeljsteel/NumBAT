import sys
from PIL import Image
import time



sys.path.append("../../../backend/")

import objects
import materials
#import mode_calcs
#import integration
import plotting
from numbattools import join_figs
from fortran import NumBAT

mat_bkg = materials.make_material("Vacuum")
mat_air = materials.make_material("Air")
mat_a = materials.make_material("Si_2016_Smith")
mat_b = materials.make_material("SiO2_2016_Smith")
mat_c = materials.make_material("SiO2_2016_Smith")
mat_d = materials.make_material("SiO2_2016_Smith")
mat_e = materials.make_material("SiO2_2016_Smith")


#def join_figs(fn1, fn2, fnout):
# 
#    images = [Image.open(x) for x in [fn2, fn1]]
#
#    nsz0 = [images[1].size[0], 
#            int(images[0].size[1]*
#            images[1].size[0]/
#            images[0].size[0])]
#    images[0]=images[0].resize(nsz0)
#
#    widths, heights = zip(*(i.size for i in images))
#
#    total_width = sum(widths)
#    max_height = max(heights)
#
#    new_im = Image.new('RGB', (total_width, max_height))
#
#    x_offset = 0
#    yoffs  = [0,0]
#    if images[0].size[1]>images[1].size[1]:
#        yoffs=[0, int((images[0].size[1]-images[1].size[1])/2)]
#    else:
#        yoffs=[int((images[1].size[1]-images[0].size[1])/2), 0]
#
#    for i, im in enumerate(images):
#      new_im.paste(im, (x_offset,yoffs[i]))
#      x_offset += im.size[0]
#
#    new_im.save(fnout)
#

def do_oneincl():
#
# oneincl_mesh.geo
# 
    unitcell_x = 4000
    unitcell_y = 3000
    inc_a_x = 1500
    inc_a_y = 750
    # plain rect waveguide
    wguide1 = objects.Structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,'rectangular',
                            material_bkg=mat_bkg, material_a=mat_a,
                            lc_bkg=.1, lc_refine_1=10, lc_refine_2=10)
    
    
    wguide2 = objects.Structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,'circular',
                            material_bkg=mat_bkg, material_a=mat_a,
                            lc_bkg=.1, lc_refine_1=10, lc_refine_2=10)
    wguide1.plot_mesh('rect_wg')
    wguide2.plot_mesh('circ_wg')
    
    #time.sleep(3)  # to let gmsh finish writing plots
    #join_figs('rect_wg_mesh_nodes.png', '../../msh_type_lib/1.png', 'rect_wg_mesh_both.png')
    #join_figs('circ_wg_mesh_nodes.png', '../../msh_type_lib/1_circular.png', 'circ_wg_mesh_both.png')
  

def do_twoincl():
    #
    # twoincl_mesh.geo
    # 
    
    # rect inclusions
    unitcell_x = 4000
    unitcell_y = 3000
    inc_a_x = 1500
    inc_a_y = 750
    inc_b_x = int(inc_a_x/2)
    inc_b_y = inc_a_y*2
    incs_y_offset = -400
    two_inc_sep=500
    wguide1 = objects.Structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,'rectangular',
                                inc_b_x=inc_b_x, inc_b_y=inc_b_y, two_inc_sep=two_inc_sep,
                               incs_y_offset=incs_y_offset,
                            material_bkg=mat_air, material_a=mat_a,
                            lc_bkg=.1, lc_refine_1=10, lc_refine_2=10)
    # circ inclusions
    wguide2 = objects.Structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,'circular',
                                inc_b_x=inc_b_x, inc_b_y=inc_b_y, two_inc_sep=two_inc_sep,
                            material_bkg=mat_a, material_a=mat_a,
                            lc_bkg=.1, lc_refine_1=10, lc_refine_2=10)
    wguide2.plot_mesh('twoincl_circ_wg')
    
    #time.sleep(5)
    #join_figs('twoincl_rect_wg_mesh_nodes.png', '../../msh_type_lib/2.png', 'twoincl_rect_wg_mesh_both.png')
    #join_figs('twoincl_circ_wg_mesh_nodes.png', '../../msh_type_lib/twoincl_circ_wg_mesh_geom.png', 'twoincl_circ_wg_mesh_both.png')
   
def do_rib():
    # rib structures inclusions
    unitcell_x = 4000
    unitcell_y = 3000
    inc_a_x = 1000
    inc_a_y = 400

    slab_a_x = 3000
    slab_a_y = 400
    slab_b_y = 100
    coat_x = 50
    coat_y = 50
    coat2_x = 70
    coat2_y = 70

    wguide1 = objects.Structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,'rib',
                                slab_a_x=slab_a_x, slab_a_y=slab_a_y, 
                            material_bkg=mat_bkg, material_a=mat_a, 
                                material_b=mat_b,
                            lc_bkg=.1, lc_refine_1=10, lc_refine_2=10, lc_refine_3=10)
    wguide1.plot_mesh('rib_wg')

    wguide2 = objects.Structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,'rib_coated',
                                slab_a_x=slab_a_x, slab_a_y=slab_a_y, 
                                coat_x=coat_x, coat_y=coat_y, 
                            material_bkg=mat_bkg, material_a=mat_a, 
                                material_b=mat_b,
                            lc_bkg=.1, lc_refine_1=10, lc_refine_2=10, lc_refine_3=10, lc_refine_4=20)
    wguide2.plot_mesh('rib_coated_wg')


    wguide3 = objects.Structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,'rib_double_coated',
                                slab_a_x=slab_a_x, slab_a_y=slab_a_y, 
                                coat_x=coat_x, coat_y=coat_y, 
                                coat2_x=coat2_x, coat2_y=coat2_y, 
                            material_bkg=mat_bkg, material_a=mat_a, material_b=mat_b,
                                material_c=mat_c, material_d=mat_d, material_e=mat_e,
                            lc_bkg=.1, lc_refine_1=10, lc_refine_2=10, lc_refine_3=10, 
                                lc_refine_4=20, lc_refine_5=20)
    wguide3.check_mesh()
    wguide3.plot_mesh('rib_double_coated_wg')

    #join_figs('rib_wg_mesh_nodes.png', '../../msh_type_lib/rib.png', 'rib_wg_mesh_both.png')
    #join_figs('rib_coated_wg_mesh_nodes.png', '../../msh_type_lib/rib_coated.png', 'rib_coated_wg_mesh_both.png')
    #join_figs('rib_double_coated_wg_mesh_nodes.png', '../../msh_type_lib/rib_double_coated.png', 'rib_double_coated_wg_mesh_both.png')
    
def do_slot():
    # rib structures inclusions
    unitcell_x = 4000
    unitcell_y = 3000
    slab_a_x = 3000
    slab_a_y = 400

    inc_a_x = 100  # slot width/pillar separation
    inc_a_y = 400  # slot depth/pillar height

    inc_b_x = 300  # pillar width
    coat_y = 50  # pillar coat thickness

    wguide1 = objects.Structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,'slot',
                                slab_a_x=slab_a_x, slab_a_y=slab_a_y, 
                                inc_b_x=inc_b_x,
                            material_bkg=mat_bkg, material_a=mat_bkg, 
                                material_b=mat_b,
                                material_c=mat_c,
                            lc_bkg=.1, lc_refine_1=10, lc_refine_2=10, lc_refine_3=10)
    wguide1.plot_mesh('slot_wg')

    unitcell_x = 4000
    unitcell_y = 1000

    wguide2 = objects.Structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,'slot_coated',
                                slab_a_x=slab_a_x, slab_a_y=slab_a_y, 
                                inc_b_x=inc_b_x,
                                coat_y=coat_y,
                            material_bkg=mat_bkg, material_a=mat_bkg, 
                                material_b=mat_b,
                                material_c=mat_c,
                                material_d=mat_d,
                            lc_bkg=.1, lc_refine_1=10, lc_refine_2=10, 
                                lc_refine_3=10, lc_refine_4=20)
    wguide2.plot_mesh('slot_coated_wg')
    wguide1.check_mesh()
    wguide2.check_mesh()

    #join_figs('slot_wg_mesh_nodes.png', '../../msh_type_lib/slot.png', 'slot_wg_mesh_both.png')
    #join_figs('slot_coated_wg_mesh_nodes.png', '../../msh_type_lib/slot_coated.png', 
    #          'slot_coated_wg_mesh_both.png')


def do_onion():

    # A concentric Bragg fibre
    unitcell_x = 1500
    unitcell_y = 1500

    # layer thicknesses
    d1=50
    d2=75

    wguide1 = objects.Structure(unitcell_x,d1,unitcell_y,d1,'onion',
                                inc_b_x=d2, inc_c_x=d1, inc_d_x=d2, inc_e_x=d1, inc_f_x=d2,
                                inc_g_x=d1, inc_h_x=d2, inc_i_x=d1, inc_j_x=d2, inc_k_x=d1,
                                inc_l_x=d2, inc_m_x=d1, inc_n_x=d2, inc_o_x=d1,
                            material_bkg=mat_bkg, material_a=mat_a, material_b=mat_b,
                                material_c=mat_a, material_d=mat_b, material_e=mat_a, material_f=mat_b,
                                material_g=mat_a, material_h=mat_b, material_i=mat_a, material_j=mat_b,
                                material_k=mat_a, material_l=mat_b, material_m=mat_a, material_n=mat_b,
                                material_o=mat_a,
                            lc_bkg=.1, lc_refine_1=5, lc_refine_2=5, lc_refine_3=5)
    wguide1.plot_mesh('onion_wg')

    # Single mode fiber 
    unitcell_x = 300
    unitcell_y = 300
    # layer thicknesses
    d1=20
    d2=127/2.

    wguide2 = objects.Structure(unitcell_x,d1,unitcell_y,d1,'onion2',
                                inc_b_x=d2, 
                            material_bkg=mat_bkg, material_a=mat_a, material_b=mat_b,
                            lc_bkg=.05, lc_refine_1=5, lc_refine_2=5, lc_refine_3=5)
    wguide2.plot_mesh('onion2_wg')


    # Single mode fiber with explicit cladding
    unitcell_x = 300
    unitcell_y = 300
    # layer thicknesses
    d1=20
    d2=127/2.
    d3=20

    wguide3 = objects.Structure(unitcell_x,d1,unitcell_y,d1,'onion3',
                                inc_b_x=d2, 
                                inc_c_x=d3, 
                            material_bkg=mat_bkg, material_a=mat_a, 
                                material_b=mat_b, material_c=mat_c,
                            lc_bkg=.05, lc_refine_1=5, lc_refine_2=5, lc_refine_3=5)
    wguide3.plot_mesh('onion3_wg')


    join_figs('onion_wg_mesh_nodes.png', '../../msh_type_lib/onion.png', 'onion_wg_mesh_both.png')
    join_figs('onion2_wg_mesh_nodes.png', '../../msh_type_lib/onion2_wg_mesh_geom.png', 'onion2_wg_mesh_both.png')
    join_figs('onion3_wg_mesh_nodes.png', '../../msh_type_lib/onion3_wg_mesh_geom.png', 'onion3_wg_mesh_both.png')

    unitcell_x = 1500
    unitcell_y = 1500
    wguide1 = objects.Structure(unitcell_x,d1,unitcell_y,d1,'circ_onion',
                                inc_b_x=d2, inc_c_x=d1, inc_d_x=d2, inc_e_x=d1, inc_f_x=d2,
                                inc_g_x=d1, inc_h_x=d2, inc_i_x=d1, inc_j_x=d2, inc_k_x=d1,
                                inc_l_x=d2, inc_m_x=d1, inc_n_x=d2, inc_o_x=d1,
                            material_bkg=mat_bkg, material_a=mat_a, material_b=mat_b,
                                material_c=mat_a, material_d=mat_b, material_e=mat_a, material_f=mat_b,
                                material_g=mat_a, material_h=mat_b, material_i=mat_a, material_j=mat_b,
                                material_k=mat_a, material_l=mat_b, material_m=mat_a, material_n=mat_b,
                                material_o=mat_a,
                            lc_bkg=.05, lc_refine_1=5, lc_refine_2=5, lc_refine_3=5)
    wguide1.plot_mesh('circ_onion_wg')

    # Single mode fiber 
    unitcell_x = 300
    unitcell_y = 300
    # layer thicknesses
    d1=20
    d2=127/2.

    wguide2 = objects.Structure(unitcell_x,d1,unitcell_y,d1,'circ_onion2',
                                inc_b_x=d2, 
                            material_bkg=mat_bkg, material_a=mat_a, material_b=mat_b,
                            lc_bkg=.05, lc_refine_1=5, lc_refine_2=5, lc_refine_3=5)
    wguide2.plot_mesh('circ_onion2_wg')

    wguide3 = objects.Structure(unitcell_x,d1,unitcell_y,d1,'circ_onion3',
                                inc_b_x=d2, 
                                inc_c_x=d3, 
                            material_bkg=mat_bkg, material_a=mat_a, 
                                material_b=mat_b, material_c=mat_c,
                            lc_bkg=.05, lc_refine_1=5, lc_refine_2=5, lc_refine_3=5)
    wguide3.plot_mesh('circ_onion3_wg')

def do_trapezoid():
    unitcell_x = 2000
    unitcell_y = 2000
    # layer thicknesses
    inc_a_x = 500
    inc_a_y = 300
    slab_a_x = 1000
    slab_a_y = 300

    wguide3 = objects.Structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,'trapezoidal_rib',
                                slab_a_x=slab_a_x, slab_a_y=slab_a_y,
                            material_bkg=mat_bkg, material_a=mat_a, 
                                material_b=mat_b,
                            lc_bkg=.05, lc_refine_1=10, lc_refine_2=10)
    wguide3.plot_mesh('trapezoidal_rib_wg')
    #wguide3.check_mesh()
    join_figs('trapezoidal_rib_wg_mesh_nodes.png', '../../msh_type_lib/trapezoidal_rib_wg_mesh_geom.png', 'trapezoidal_rib_wg_mesh_both.png')


def do_pedestal():
    unitcell_x = 3000
    unitcell_y = 2000

    # layer thicknesses
    inc_a_x = 500
    inc_a_y = 100
    inc_b_x = 350

    pillar_x = 20
    pillar_y = 200

    slab_a_x = 2000
    slab_a_y = 300


    wguide1 = objects.Structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,'pedestal',
                                inc_b_x = inc_b_x, slab_a_x=slab_a_x, slab_a_y=slab_a_y,
                            material_bkg=mat_bkg, material_a=mat_a, material_b=mat_b,
                            pillar_x=pillar_x, pillar_y=pillar_y,
                            lc_bkg=.1, lc_refine_1=10, lc_refine_2=10)
    wguide1.plot_mesh('pedestal_wg')
    #wguide3.check_mesh()
    join_figs('pedestal_wg_mesh_nodes.png', '../../msh_type_lib/pedestal_wg_mesh_geom.png', 
              'pedestal_wg_mesh_both.png')

def do_main():
    #do_oneincl()
    #do_twoincl()
    #do_rib()
    #do_slot()
    do_onion()
    #do_trapezoid()
    #do_pedestal()


if __name__=='__main__':
    do_main()
