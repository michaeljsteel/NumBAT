import sys
from PIL import Image, ImageDraw, ImageFont
import time
import math



sys.path.append("../../../backend/")

import numbat
import materials

mat_bkg = materials.make_material("Vacuum")
mat_air = materials.make_material("Air")
mat_a = materials.make_material("Si_2016_Smith")
mat_b = materials.make_material("SiO2_2016_Smith")
mat_c = materials.make_material("SiO2_2016_Smith")
mat_d = materials.make_material("SiO2_2016_Smith")
mat_e = materials.make_material("SiO2_2016_Smith")

fnt_mat = ImageFont.truetype("Pillow/Tests/fonts/FreeMono.ttf", 18)

col_mat = 'black'
col_lc = 'darkgreen'
col_dim = 'brown'

def add_mat_lab(im, bw, bh, lab, x, y):
    d = ImageDraw.Draw(im)
    d.text((x,y), lab, font=fnt_mat, fill=col_mat, anchor='mm')

def add_lc(im, bw, bh, lab, x, y, compass):

    arlen = .05*bw
    sqlen = .05*bw/math.sqrt(2)
    ax = x
    ay = y

    match compass:
        case 'N':
            ay -= alen
            anc = 'mb'
        case 'S':
            ay += alen
            anc = 'mt'
        case 'E':
            ax += alen
            anc = 'lm'
        case 'W':
            ax -= alen
            anc = 'rm'
        case 'NW':
            ax -= sqlen
            ay -= sqlen
            anc = 'rb'
        case 'NE':
            ax += sqlen
            ay -= sqlen
            anc = 'lb'
        case 'SW':
            ax -= sqlen
            ay += sqlen
            anc = 'rt'
        case 'SE':
            ax += sqlen
            ay += sqlen
            anc = 'lt'

    d = ImageDraw.Draw(im)
    d.line((x,y,ax,ay), fill=col_lc, width=2)
    d.text((ax,ay), lab, font=fnt_mat, fill=col_lc, anchor=anc)

def add_dim(im, bw, bh, lab, x1, y1, x2, y2, compass):

    alen = bh*.02
    tx = (x1 + x2)/2
    ty = (y1+y2)/2

    match compass:
        case 'N':
            ty -= alen
            anc = 'mb'
        case 'S':
            ty += alen
            anc = 'mt'
        case 'E':
            tx += alen
            anc = 'lm'
        case 'W':
            tx -= alen
            anc = 'rm'

    d = ImageDraw.Draw(im)
    d.line((x1,y1,x2,y2), fill=col_dim)
    d.text((tx,ty), lab, font=fnt_mat, fill=col_dim, anchor=anc)

def do_oneincl(nbapp):
#
# oneincl_mesh.geo
# 
    unitcell_x = 4000
    unitcell_y = 3000
    inc_a_x = 1500
    inc_a_y = 750
    # plain rect waveguide
    wguide1 = nbapp.make_structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,'rectangular',
                            material_bkg=mat_bkg, material_a=mat_a,
                            lc_bkg=.1, lc_refine_1=10, lc_refine_2=10)
    
    
    wguide2 = nbapp.make_structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,'circular',
                            material_bkg=mat_bkg, material_a=mat_a,
                            lc_bkg=.1, lc_refine_1=10, lc_refine_2=10)
    #wguide1.plot_mesh('rect_wg')
    #wguide2.plot_mesh('circ_wg')

    
    unitcell_x = 4000
    unitcell_y = 3000
    inc_a_x = 2000
    inc_b_x = 600
    inc_b_y = 1100

    wguide3 = nbapp.make_structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,'triangular',
                            material_bkg=mat_bkg, material_a=mat_a,
                            inc_b_x = inc_b_x, inc_b_y = inc_b_y,
                            lc_bkg=.1, lc_refine_1=10, lc_refine_2=10)


    frt = 'triangular'
    wguide3.plot_mesh(frt+'_wg')


    with Image.open(frt+'_wg-mesh.png').convert('RGBA') as im:
        sz = im.size
        bl = im.size[0]*.025
        br = im.size[0]*.475
        bt = im.size[1]*.06
        bb = im.size[1]*.93
        bw = br-bl
        bh = bb-bt

        scalx = 1/unitcell_x
        scaly = 1/unitcell_y

        bumpx = -.0*bw
        bumpy = +.01*bh
        x0 = (bl+br)/2  + bumpx
        y0 = (bt+bb)/2 + bumpy
        ymid = bt+bh*(.5+.5*inc_b_y*scaly) + bumpy

        #add_mat_lab(im, bw, bh, 'x', x0, ymid)
        #add_mat_lab(im, bw, bh, 'x', bl, bt)
        #add_mat_lab(im, bw, bh, 'x', br, bt)
        #add_mat_lab(im, bw, bh, 'x', bl, bb)
        #add_mat_lab(im, bw, bh, 'x', br, bb)

        add_mat_lab(im, bw, bh, 'mat_b1', x0-bw/15, ymid-bh/8)
        add_mat_lab(im, bw, bh, 'mat_bg', x0-bw/3, y0-bh/6)
        add_mat_lab(im, bw, bh, 'mat_bg', x0+bw/3, y0-bh/6)
        add_mat_lab(im, bw, bh, 'mat_bg', x0, ymid+bh/6)

        add_lc(im, bw, bh, 'lc_bg', x0-bw/6, bt, 'SW')
        add_lc(im, bw, bh, 'lc_1',  x0-bw/6, ymid, 'SW')

        x1 = x0-.5*bw*inc_a_x*scalx; x2 = x1+bw*inc_a_x*scalx
        y1 = ymid+.03*bh; y2 = ymid+.03*bh
        add_dim(im, bw, bh, 'inc_a_x', x1, y1, x2, y2, 'S')

        x1 = x0-.5*bw*inc_a_x*scalx; x2 = x1+bw*inc_b_x*scalx
        y1 = ymid-.03*bh; y2 = ymid-.03*bh
        add_dim(im, bw, bh, 'inc_b_x', x1, y1, x2, y2, 'N')

        x1 = x0+.5*bw*inc_a_x*scalx; x2 = x1
        y1 = ymid; y2 = y1-bh*inc_b_y*scaly
        add_dim(im, bw, bh, 'inc_b_y', x1, y1, x2, y2, 'E')


        im.show()
        im.save(frt+'_wg-mesh-annotated.png')


def do_twoincl(nbapp):
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
    wguide1 = nbapp.make_structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,'rectangular',
                                inc_b_x=inc_b_x, inc_b_y=inc_b_y, two_inc_sep=two_inc_sep,
                               incs_y_offset=incs_y_offset,
                            material_bkg=mat_air, material_a=mat_a,
                            lc_bkg=.1, lc_refine_1=10, lc_refine_2=10)
    # circ inclusions
    wguide2 = nbapp.make_structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,'circular',
                                inc_b_x=inc_b_x, inc_b_y=inc_b_y, two_inc_sep=two_inc_sep,
                            material_bkg=mat_a, material_a=mat_a,
                            lc_bkg=.1, lc_refine_1=10, lc_refine_2=10)
    wguide2.plot_mesh('twoincl_circ_wg')
    
   
def do_rib(nbapp):
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

    wguide1 = nbapp.make_structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,'rib',
                                slab_a_x=slab_a_x, slab_a_y=slab_a_y, 
                            material_bkg=mat_bkg, material_a=mat_a, 
                                material_b=mat_b,
                            lc_bkg=.1, lc_refine_1=10, lc_refine_2=10, lc_refine_3=10)
    wguide1.plot_mesh('rib_wg')

    wguide2 = nbapp.make_structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,'rib_coated',
                                slab_a_x=slab_a_x, slab_a_y=slab_a_y, 
                                coat_x=coat_x, coat_y=coat_y, 
                            material_bkg=mat_bkg, material_a=mat_a, 
                                material_b=mat_b,
                            lc_bkg=.1, lc_refine_1=10, lc_refine_2=10, lc_refine_3=10, lc_refine_4=20)
    wguide2.plot_mesh('rib_coated_wg')


    wguide3 = nbapp.make_structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,'rib_double_coated',
                                slab_a_x=slab_a_x, slab_a_y=slab_a_y, 
                                coat_x=coat_x, coat_y=coat_y, 
                                coat2_x=coat2_x, coat2_y=coat2_y, 
                            material_bkg=mat_bkg, material_a=mat_a, material_b=mat_b,
                                material_c=mat_c, material_d=mat_d, material_e=mat_e,
                            lc_bkg=.1, lc_refine_1=10, lc_refine_2=10, lc_refine_3=10, 
                                lc_refine_4=20, lc_refine_5=20)
    #wguide3.check_mesh()
    wguide3.plot_mesh('rib_double_coated_wg')

    
def do_slot(nbapp):
    # rib structures inclusions
    unitcell_x = 4000
    unitcell_y = 3000
    slab_a_x = 3000
    slab_a_y = 400

    inc_a_x = 100  # slot width/pillar separation
    inc_a_y = 400  # slot depth/pillar height

    inc_b_x = 300  # pillar width
    coat_y = 50  # pillar coat thickness

    wguide1 = nbapp.make_structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,'slot',
                                slab_a_x=slab_a_x, slab_a_y=slab_a_y, 
                                inc_b_x=inc_b_x,
                            material_bkg=mat_bkg, material_a=mat_bkg, 
                                material_b=mat_b,
                                material_c=mat_c,
                            lc_bkg=.1, lc_refine_1=10, lc_refine_2=10, lc_refine_3=10)
    wguide1.plot_mesh('slot_wg')

    unitcell_x = 4000
    unitcell_y = 1000

    wguide2 = nbapp.make_structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,'slot_coated',
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
    #wguide1.check_mesh()
    #wguide2.check_mesh()



def do_onion(nbapp):

    # A concentric Bragg fibre
    unitcell_x = 1500
    unitcell_y = 1500

    # layer thicknesses
    d1=50
    d2=75

    wguide1 = nbapp.make_structure(unitcell_x,d1,unitcell_y,d1,'onion',
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

    wguide2 = nbapp.make_structure(unitcell_x,d1,unitcell_y,d1,'onion2',
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

    wguide3 = nbapp.make_structure(unitcell_x,d1,unitcell_y,d1,'onion3',
                                inc_b_x=d2, 
                                inc_c_x=d3, 
                            material_bkg=mat_bkg, material_a=mat_a, 
                                material_b=mat_b, material_c=mat_c,
                            lc_bkg=.05, lc_refine_1=5, lc_refine_2=5, lc_refine_3=5)
    wguide3.plot_mesh('onion3_wg')



    unitcell_x = 1500
    unitcell_y = 1500
    wguide1 = nbapp.make_structure(unitcell_x,d1,unitcell_y,d1,'circ_onion',
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

    wguide2 = nbapp.make_structure(unitcell_x,d1,unitcell_y,d1,'circ_onion2',
                                inc_b_x=d2, 
                            material_bkg=mat_bkg, material_a=mat_a, material_b=mat_b,
                            lc_bkg=.05, lc_refine_1=5, lc_refine_2=5, lc_refine_3=5)
    wguide2.plot_mesh('circ_onion2_wg')

    wguide3 = nbapp.make_structure(unitcell_x,d1,unitcell_y,d1,'circ_onion3',
                                inc_b_x=d2, 
                                inc_c_x=d3, 
                            material_bkg=mat_bkg, material_a=mat_a, 
                                material_b=mat_b, material_c=mat_c,
                            lc_bkg=.05, lc_refine_1=5, lc_refine_2=5, lc_refine_3=5)
    wguide3.plot_mesh('circ_onion3_wg')

def do_trapezoid(nbapp):
    unitcell_x = 2000
    unitcell_y = 2000
    # layer thicknesses
    inc_a_x = 500
    inc_a_y = 300
    slab_a_x = 1000
    slab_a_y = 300

    wguide3 = nbapp.make_structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,'trapezoidal_rib',
                                slab_a_x=slab_a_x, slab_a_y=slab_a_y,
                            material_bkg=mat_bkg, material_a=mat_a, 
                                material_b=mat_b,
                            lc_bkg=.05, lc_refine_1=10, lc_refine_2=10)
    wguide3.plot_mesh('trapezoidal_rib_wg')
    #wguide3.check_mesh()

def do_pedestal(nbapp):
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


    wguide1 = nbapp.make_structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,'pedestal',
                                inc_b_x = inc_b_x, slab_a_x=slab_a_x, slab_a_y=slab_a_y,
                            material_bkg=mat_bkg, material_a=mat_a, material_b=mat_b,
                            pillar_x=pillar_x, pillar_y=pillar_y,
                            lc_bkg=.1, lc_refine_1=10, lc_refine_2=10)
    wguide1.plot_mesh('pedestal_wg')

def do_main():

    nbapp = numbat.NumBATApp()

    do_oneincl(nbapp)
    do_twoincl(nbapp)
    do_rib(nbapp)
    do_slot(nbapp)
    do_onion(nbapp)
    do_trapezoid(nbapp)
    do_pedestal(nbapp)


if __name__=='__main__':
    do_main()
