import sys
from PIL import Image, ImageDraw, ImageFont
import math
from pathlib import Path

sys.path.append("../../../backend/")

import numbat
import materials
from plottools import join_figs

mat_bkg = materials.make_material("Vacuum")
mat_air = materials.make_material("Air")
mat_a = materials.make_material("Si_2016_Smith")
mat_b = materials.make_material("As2S3_2023_Steel")
mat_c = materials.make_material("SiO2_2016_Smith")
mat_d = materials.make_material("Si3N4_2021_Steel")
mat_e = materials.make_material("GaAs_2021_Poulton")

#fnt_mat = ImageFont.truetype("Pillow/Tests/fonts/FreeMono.ttf", 18)
fnt_mat = ImageFont.truetype(r'C:\Windows\Fonts\Times.ttf', 18)


col_mat = 'black'
col_lc = 'darkgreen'
col_dim = 'brown'

fn_suff_raw = '_wg'
fn_ext_raw = '_wg-mesh.png'
fn_ext_anno = '_wg-mesh-annotated.png'

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
            ay -= arlen
            anc = 'mb'
        case 'S':
            ay += arlen
            anc = 'mt'
        case 'E':
            ax += arlen
            anc = 'lm'
        case 'W':
            ax -= arlen
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

def get_sizes(im,x0,x1,y0,y1,un_x, un_y):
    '''Controls the view into the image in pixels.
    Open the combined grid and mesh image and find the coordinates
    of the grid plot corners. Express x0, x1, y0, y1 as fractions of the
    whole window.'''

    sz = im.size
    bl = im.size[0]*x0
    br = im.size[0]*x1
    bt = im.size[1]*y0
    bb = im.size[1]*y1
    bw = br-bl
    bh = bb-bt

    scalx = bw/un_x
    scaly = bh/un_y

    #print(sz, bl, br, bt, bb, bw, bh)
    return (bl, br, bt, bb, bw, bh, scalx, scaly)

def get_sizesb(im, un_x, un_y, border=20):
    '''Controls the view into the image in pixels.
    Open the combined grid and mesh image and find the coordinates
    of the grid plot corners. Express x0, x1, y0, y1 as fractions of the
    whole window.'''

    sz = im.size
    bl = border                   # new left border
    br = im.size[0] - border      # new right border
    bt = border
    bb = im.size[1] - border      # new top border
    bw = br-bl                    # new bottom border
    bh = bb-bt

    scalx = bw/un_x
    scaly = bh/un_y

    return (bl, br, bt, bb, bw, bh, scalx, scaly)

def do_oneincl_rect(nbapp):
    un_x = 4000; un_y = 3000; inc_a_x = 1500; inc_a_y = 750

    # plain rect waveguide
    wguide = nbapp.make_structure( 'rectangular', un_x, un_y, inc_a_x, inc_a_y,
                            material_bkg=mat_bkg, material_a=mat_a,
                            lc_bkg=.1, lc_refine_1=10, lc_refine_2=10)
    frt = 'rect'

    # 'rect_wg-mesh.png'
    fn_wire, fn_mesh, fn_join = wguide.plot_mesh(frt+fn_suff_raw,
                                                 combo_plot=False)
    # 'rect_wg-ref_index.png'
    fn_refindex = wguide.plot_refractive_index_profile(frt+fn_suff_raw)

    fn_anno = frt + fn_ext_anno
    do_oneincl_rect_annotate(frt, un_x, un_y, inc_a_x, inc_a_y,
                             fn_wire, fn_mesh, fn_anno)

def do_oneincl_rect_annotate(frt, un_x, un_y, inc_a_x, inc_a_y,
                             fin_wire, fin_mesh, fout_anno):

    with Image.open(fin_wire).convert('RGBA') as im:
        (bl, br, bt, bb, bw, bh, scalx, scaly) = get_sizesb(im, un_x, un_y)

        x0 = (bl+br)/2        # horz. middle
        y0 = (bt+bb)/2        # vert. middle
        dx = (br-bl)/20       # useful horz. delta
        dy = (bt-bb)/20       # useful vert. delta

        add_mat_lab(im, bw, bh, 'mat_bkg', x0-2*dx, y0+dy)
        add_mat_lab(im, bw, bh, 'mat_a', x0-dx*6, y0+dy*5)

        add_lc(im, bw, bh, 'lc_bkg', x0-dx*3, bt-dy, 'SW')
        add_lc(im, bw, bh, 'lc_1',  x0-inc_a_x*scalx*.5, y0-dy, 'SW')

        xs = x0-.53*inc_a_x*scalx; xf=xs
        ys = y0+.5*inc_a_y*scaly; yf = y0-.5*inc_a_y*scaly
        add_dim(im, bw, bh, 'inc_a_y', xs, ys, xf, yf, 'W')

        xs = x0-.5*inc_a_x*scalx; xf = x0+.5*inc_a_x*scalx;
        ys = y0-.53*inc_a_y*scaly; yf = ys
        add_dim(im, bw, bh, 'inc_a_x', xs, ys, xf, yf, 'N')

        im.show()
        ftmp = fin_wire[:-4] + '_tmp.png'
        im.save(frt+fn_ext_anno)
        im.save(ftmp)
        join_figs([ftmp, fin_mesh], fout_anno, clip=None, trimwhite=False, border=20, delete_inputs=True)
        Path(fin_wire).unlink()




def do_oneincl_circ(nbapp):
    un_x = 4000; un_y = 3000; inc_a_x = 1500; inc_a_y = 750

    wguide = nbapp.make_structure('circular', un_x, un_y, inc_a_x, inc_a_y,
                            material_bkg=mat_bkg, material_a=mat_a,
                            lc_bkg=.1, lc_refine_1=10, lc_refine_2=10)
    frt = 'circ'
    wguide.plot_mesh(frt+fn_suff_raw)
    wguide.plot_refractive_index_profile(frt+fn_suff_raw)

    do_oneincl_circ_annotate(frt, un_x, un_y, inc_a_x, inc_a_y)

def do_oneincl_circ_annotate(frt, un_x, un_y, inc_a_x, inc_a_y):

    with Image.open(frt+fn_ext_raw).convert('RGBA') as im:
        #(bl, br, bt, bb, bw, bh, scalx, scaly) = get_sizes(im, 0.025, 0.475,  0.06, 0.93,  un_x, un_y)
        (bl, br, bt, bb, bw, bh, scalx, scaly) = get_sizesb(im, un_x, un_y)

        x0 = (bl+br)/2
        y0 = (bt+bb)/2
        dx = (br-bl)/20
        dy = (bt-bb)/20

        add_mat_lab(im, bw, bh, 'mat_bkg', x0-dx,   y0+dy)
        add_mat_lab(im, bw, bh, 'mat_a',   x0-5*dx, y0+dy*5)


        add_lc(im, bw, bh, 'lc_bkg', x0-2*dx, bt-dy, 'SW')
        add_lc(im, bw, bh, 'lc_2',   x0-dx, y0+bh/8, 'SW')
        add_lc(im, bw, bh, 'lc_3',   x0,      y0, 'SE')

        xs = x0-.53*inc_a_x*scalx; xf=xs
        ys = y0+.5*inc_a_y*scaly; yf = y0-.5*inc_a_y*scaly
        add_dim(im, bw, bh, 'inc_a_y', xs, ys, xf, yf, 'W')

        xs = x0-.5*inc_a_x*scalx; xf = x0+.5*inc_a_x*scalx;
        ys = y0-.53*inc_a_y*scaly; yf = ys
        add_dim(im, bw, bh, 'inc_a_x', xs, ys, xf, yf, 'N')



        im.show()
        im.save(frt+fn_ext_anno)

def do_oneincl_triang(nbapp):
    un_x = 4000; un_y = 3000; bw = 2000; pxo = 750; ph = 600;

    wguide = nbapp.make_structure('triangular', un_x, un_y,
                                   base_width=bw, peak_xoff = pxo, peak_height=ph,
                                   material_bkg=mat_bkg, material_a=mat_a,
                                   lc_bkg=.1, lc_refine_1=10, lc_refine_2=10)


    frt = 'triangular'
    wguide.plot_mesh(frt+fn_suff_raw)
    wguide.plot_refractive_index_profile(frt+fn_suff_raw)

    do_oneincl_triang_annotate(frt, un_x, un_y, bw, pxo, ph)

def do_oneincl_triang_annotate(frt, un_x, un_y, basewid, pxo, ph):
    with Image.open(frt+fn_ext_raw).convert('RGBA') as im:
        #(bl, br, bt, bb, bw, bh, scalx, scaly) = get_sizes(im, 0.025, 0.475,  0.06, 0.93,  un_x, un_y)
        (bl, br, bt, bb, bw, bh, scalx, scaly) = get_sizesb(im, un_x, un_y)


        dx = (br-bl)/20
        dy = (bt-bb)/20   # negative value so adding dy moves curso up screen

        bumpx = -.0*bw
        bumpy = +.01*bh

        x0 = (bl+br)/2 + bumpx
        y0 = (bt+bb)/2 + bumpy
        ymid = bt+bh*.5+.5*ph*scaly + bumpy  # bh is middle of triangle, ymid is baseline

        add_mat_lab(im, bw, bh, 'mat_a',   x0-bw/15, ymid-bh/8)
        add_mat_lab(im, bw, bh, 'mat_bkg', x0-bw/3,  y0-bh/6)
        add_mat_lab(im, bw, bh, 'mat_bkg', x0+bw/3,  y0-bh/6)
        add_mat_lab(im, bw, bh, 'mat_bkg', x0,       ymid+bh/6)

        add_lc(im, bw, bh, 'lc_bkg', x0-bw/6, bt-dy, 'SW')
        add_lc(im, bw, bh, 'lc_1',   x0-bw/3, ymid,  'SW')

        x1 = x0-.5*basewid*scalx; x2 = x1+basewid*scalx
        y1 = ymid+.03*bh; y2 = ymid+.03*bh
        add_dim(im, bw, bh, 'base_width', x1, y1, x2, y2, 'S')

        x1 = x0-.5*basewid*scalx; x2 = x1+pxo*scalx
        y1 = ymid-.03*bh; y2 = ymid-.03*bh
        add_dim(im, bw, bh, 'peak_xoff', x1, y1, x2, y2, 'N')

        x1 = x0+.5*ph*scalx; x2 = x1
        y1 = ymid; y2 = y1-ph*scaly
        add_dim(im, bw, bh, 'peak_height', x1, y1, x2, y2, 'E')


        im.show()
        im.save(frt+fn_ext_anno)




def do_oneincl(nbapp):
    do_oneincl_rect(nbapp)
    #do_oneincl_circ(nbapp)
    #do_oneincl_triang(nbapp)


def do_twoincl(nbapp):
    #
    # twoincl_mesh.geo
    #

    # rect inclusions
    un_x = 4000
    un_y = 3000
    inc_a_x = 1500
    inc_a_y = 750
    inc_b_x = int(inc_a_x/2)
    inc_b_y = inc_a_y*2
    incs_y_offset = -400
    two_inc_sep=500
    wguide1 = nbapp.make_structure('rectangular', un_x,un_y,inc_a_x,inc_a_y,
                                inc_b_x=inc_b_x, inc_b_y=inc_b_y, two_inc_sep=two_inc_sep,
                               incs_y_offset=incs_y_offset,
                            material_bkg=mat_air, material_a=mat_a, material_b=mat_b,
                            lc_bkg=.1, lc_refine_1=10, lc_refine_2=10)

    wguide1.plot_mesh('twoincl_rect_wg')
    wguide1.plot_refractive_index_profile('twoincl_rect_wg')
    # circ inclusions
    wguide2 = nbapp.make_structure('circular', un_x,un_y,inc_a_x,inc_a_y,
                                inc_b_x=inc_b_x, inc_b_y=inc_b_y, two_inc_sep=two_inc_sep,
                            material_bkg=mat_air, material_a=mat_a, material_b=mat_b,
                            lc_bkg=.1, lc_refine_1=10, lc_refine_2=10)
    wguide2.plot_mesh('twoincl_circ_wg')
    wguide2.plot_refractive_index_profile('twoincl_circ_wg')


def do_rib(nbapp):
    # rib structures inclusions
    un_x = 4000
    un_y = 3000
    inc_a_x = 1000
    inc_a_y = 400

    slab_a_x = 3000
    slab_a_y = 400
    slab_b_h = 100
    coat_x = 70
    coat_y = 70
    coat2_x = 140
    coat2_y = 140

    wguide1 = nbapp.make_structure('rib',un_x,un_y, rib_w=inc_a_x, rib_h=inc_a_y,
                                slab_w=slab_a_x, slab_h=slab_a_y,
                            material_bkg=mat_bkg, material_a=mat_a,
                                material_b=mat_b,
                            lc_bkg=.1, lc_refine_1=10, lc_refine_2=10)
    wguide1.plot_mesh('rib_wg')
    wguide1.plot_refractive_index_profile('rib_wg')

    wguide2 = nbapp.make_structure('rib_coated', un_x, un_y, rib_w=inc_a_x, rib_h=inc_a_y,
                                   slab_w=slab_a_x, slab_h=slab_a_y, coat_w=coat_x, coat_h=coat_y,
                                   material_bkg=mat_bkg, material_a=mat_a, material_b=mat_b, material_c=mat_c,
                                   lc_bkg=.1, lc_refine_1=10, lc_refine_2=10)
    wguide2.plot_mesh('rib_coated_wg')
    wguide2.plot_refractive_index_profile('rib_coated_wg')


    wguide3 = nbapp.make_structure('rib_double_coated',un_x,un_y, rib_w=inc_a_x, rib_h=inc_a_y,
                                   slab_w=slab_a_x, slab_h=slab_a_y, slab_b_h=slab_b_h,
                                   coat_w=coat_x, coat_h=coat_y, coat2_w=coat2_x, coat2_h=coat2_y,
                                   material_bkg=mat_bkg, material_a=mat_a, material_b=mat_b, material_c=mat_c,
                                   material_d=mat_d, material_e=mat_e,
                                   lc_bkg=.1, lc_refine_1=10, lc_refine_2=10, lc_refine_3=10, lc_refine_4=20, lc_refine_5=20)
    #wguide3.check_mesh()
    wguide3.plot_mesh('rib_double_coated_wg')
    wguide3.plot_refractive_index_profile('rib_double_coated_wg')


def do_slot(nbapp):
    # rib structures inclusions
    un_x = 4000
    un_y = 3000
    slab_a_x = 3000
    slab_a_y = 400

    inc_a_x = 100  # slot width/pillar separation
    inc_a_y = 400  # slot depth/pillar height

    inc_b_x = 300  # pillar width
    coat_y = 50  # pillar coat thickness

    wguide1 = nbapp.make_structure('slot', un_x,un_y,
                                slab_w=slab_a_x, slab_h=slab_a_y,
                                   rib_w=inc_b_x,
                                   rib_h=inc_a_y,
                                slot_w=inc_b_x,
                            material_bkg=mat_bkg, material_a=mat_bkg,
                                material_b=mat_b,
                                material_c=mat_c,
                            lc_bkg=.1, lc_refine_1=10, lc_refine_2=10)
    wguide1.plot_mesh('slot_wg')

    un_x = 4000
    un_y = 1000

    wguide2 = nbapp.make_structure( 'slot_coated', un_x, un_y,
                                slab_w=slab_a_x, slab_h=slab_a_y,
                                   rib_w=inc_b_x,
                                   rib_h=inc_a_y,
                                slot_w=inc_b_x,
                                coat_t=coat_y,
                            material_bkg=mat_bkg, material_a=mat_bkg,
                                material_b=mat_b,
                                material_c=mat_c,
                                material_d=mat_d,
                            lc_bkg=.1, lc_refine_1=10, lc_refine_2=10,
                                lc_refine_3=10)
    wguide2.plot_mesh('slot_coated_wg')
    #wguide1.check_mesh()
    #wguide2.check_mesh()



def do_onion(nbapp):

    # A concentric Bragg fibre
    un_x = 1000
    un_y = 1000

    # layer thicknesses
    d1=50
    d2=127/2
    d3=50

    wguide = nbapp.make_structure('onion',un_x, un_y, d1,
                                inc_b_x=d2, inc_c_x=d1, inc_d_x=d2, inc_e_x=d1, inc_f_x=d2,
                                inc_g_x=d1, inc_h_x=d2, inc_i_x=d1, inc_j_x=d2, inc_k_x=d1,
                                inc_l_x=d2, inc_m_x=d1, inc_n_x=d2, inc_o_x=d1,
                            material_bkg=mat_bkg, material_a=mat_a, material_b=mat_b,
                                material_c=mat_a, material_d=mat_b, material_e=mat_a, material_f=mat_b,
                                material_g=mat_a, material_h=mat_b, material_i=mat_a, material_j=mat_b,
                                material_k=mat_a, material_l=mat_b, material_m=mat_a, material_n=mat_b,
                                material_o=mat_a,
                            lc_bkg=.1, lc_refine_1=5, lc_refine_2=5)
    frt = 'onionN'
    wguide.plot_mesh(frt+fn_suff_raw)
    wguide.plot_refractive_index_profile(frt+fn_suff_raw)

    #do_onion23n_annotate(frt, un_x, un_y, d1, d2, 4)



    # Single mode fiber
    un_x = 300
    un_y = 300
    # layer thicknesses
    d1=20
    d2=127/2.

    wguide = nbapp.make_structure('onion2',un_x, un_y,
                                   inc_a_x=d1, inc_b_x=d2,
                            material_bkg=mat_bkg, material_a=mat_a, material_b=mat_b,
                            lc_bkg=.05, lc_refine_2=5)

    frt = 'onion2'
    wguide.plot_mesh(frt+fn_suff_raw)
    wguide.plot_refractive_index_profile(frt+fn_suff_raw)
    do_onion23n_annotate(frt, un_x, un_y, d1, d2, 2)

    # Single mode fiber with explicit cladding
    un_x = 300
    un_y = 300
    # layer thicknesses
    d1=20
    d2=127/2.
    d3=d1

    wguide = nbapp.make_structure('onion3',un_x,un_y,
                                   inc_a_x=d1, inc_b_x=d2, inc_c_x=d3,
                                    material_bkg=mat_bkg, material_a=mat_a,
                                material_b=mat_b, material_c=mat_c,
                            lc_bkg=.05, lc_refine_2=5)
    frt = 'onion3'
    wguide.plot_mesh(frt+fn_suff_raw)
    wguide.plot_refractive_index_profile(frt+fn_suff_raw)
    do_onion23n_annotate(frt, un_x, un_y, d1, d2, 3)

    un_x = 1500
    un_y = 1500
    wguide = nbapp.make_structure('circ_onion',un_x,un_y,
                                   inc_a_x=d1,
                                inc_b_x=d2, inc_c_x=d1, inc_d_x=d2, inc_e_x=d1, inc_f_x=d2,
                                inc_g_x=d1, inc_h_x=d2, inc_i_x=d1, inc_j_x=d2, inc_k_x=d1,
                                inc_l_x=d2, inc_m_x=d1, inc_n_x=d2, inc_o_x=d1,
                            material_bkg=mat_bkg, material_a=mat_a, material_b=mat_b,
                                material_c=mat_a, material_d=mat_b, material_e=mat_a, material_f=mat_b,
                                material_g=mat_a, material_h=mat_b, material_i=mat_a, material_j=mat_b,
                                material_k=mat_a, material_l=mat_b, material_m=mat_a, material_n=mat_b,
                                material_o=mat_a,
                            lc_bkg=.05, lc_refine_2=5)

    frt = 'circ_onionN'
    wguide.plot_mesh(frt+fn_suff_raw)
    wguide.plot_refractive_index_profile(frt+fn_suff_raw)
    #o_circ_onionN_annotate(frt, un_x, un_y, d1, d2, d3)


    # Single mode fiber
    un_x = 300
    un_y = 300
    # layer thicknesses
    d1=20
    d2=127/2.
    d3=d1

    wguide = nbapp.make_structure('circ_onion2',un_x,un_y,
                                   inc_a_x=d1, inc_b_x=d2,
                            material_bkg=mat_bkg, material_a=mat_a, material_b=mat_b,
                            lc_bkg=.05, lc_refine_2=5)

    frt = 'circ_onion2'
    wguide.plot_mesh(frt+fn_suff_raw)
    wguide.plot_refractive_index_profile(frt+fn_suff_raw)
    do_circ_onion23_annotate(frt, un_x, un_y, d1, d2,2)

    wguide = nbapp.make_structure('circ_onion3',un_x,un_y,
                                   inc_a_x=d1, inc_b_x=d2, inc_c_x=d3,
                            material_bkg=mat_bkg, material_a=mat_a,
                                material_b=mat_b, material_c=mat_c,
                            lc_bkg=.05, lc_refine_2=5)

    frt = 'circ_onion3'
    wguide.plot_mesh(frt+fn_suff_raw)
    wguide.plot_refractive_index_profile(frt+fn_suff_raw)
    do_circ_onion23_annotate(frt, un_x, un_y, d1, d2, 3)



def do_onion23n_annotate(frt, un_x, un_y, d1, d2, lev):
    with Image.open(frt+fn_ext_raw).convert('RGBA') as im:
        (bl, br, bt, bb, bw, bh, scalx, scaly) = get_sizes(im, 0.09116, 0.408,  0.0753, 0.9246,  un_x, un_y)

    # get base coords
    bumpx = -.0*bw
    bumpy = +.01*bh

    d1 *= scalx
    d2 *= scalx

    x0 = (bl+br)/2 + bumpx
    y0 = (bt+bb)/2 + bumpy
    dx = (br-bl)/20
    dy = (bt-bb)/20   # negative value so adding dy moves curso up screen


    add_mat_lab(im, bw, bh, 'mat_bkg', x0-4*dx, y0+7*dy)
    add_mat_lab(im, bw, bh, 'mat_a',   x0-1.25*dx, y0+.75*dy)
    add_mat_lab(im, bw, bh, 'mat_b',   x0+1.5*dx, y0+3*dy)


    add_lc(im, bw, bh, 'lc_bkg', x0-bw/6, bt, 'SW')
    add_lc(im, bw, bh, 'lc_2',   x0-d1/2, y0,  'SW')
    add_lc(im, bw, bh, 'lc_2',   x0-d1/2-d2, y0,  'SW')
    add_lc(im, bw, bh, 'lc_2',   x0-d1/2-d2-d1, y0,  'NW')


    c1x1 = x0-d1/2; c1x2 = x0+d1/2
    y1 = y0-dy; y2=y0-dy
    add_dim(im, bw, bh, 'inc_a_x', c1x1, y1, c1x2, y2, 'S')

    c2x1 = c1x2; c2x2 = c2x1+d2
    y1 = y0+dy; y2=y0+dy
    add_dim(im, bw, bh, 'inc_b_x', c2x1, y1, c2x2, y2, 'N')

    if lev>=3:
        c3x1 = c2x2; c3x2 = c3x1+d1
        y1 = y0-dy/2; y2=y0-dy/2
        add_dim(im, bw, bh, 'inc_c_x', c3x1, y1, c3x2, y2, 'S')
        add_mat_lab(im, bw, bh, 'mat_c',   x0+1.5*dx, y0-d1-d2)

    if lev>=4:
        c4x1 = c3x2; c4x2 = c4x1+d2
        y1 = y0+dy/2; y2=y0+dy/2
        add_dim(im, bw, bh, 'inc_d_x', c4x1, y1, c4x2, y2, 'N')
        add_mat_lab(im, bw, bh, 'mat_d',   x0+1.5*dx, y0-d1-d2-d1)

    im.show()
    im.save(frt+fn_ext_anno)




def do_circ_onion23_annotate(frt, un_x, un_y, d1, d2, lev):
    with Image.open(frt+fn_ext_raw).convert('RGBA') as im:
        (bl, br, bt, bb, bw, bh, scalx, scaly) = get_sizes(im, 0.09116, 0.408,  0.0753, 0.9246,  un_x, un_y)

    # get base coords
    bumpx = -.0*bw
    bumpy = +.01*bh

    d1 *= scalx
    d2 *= scalx

    x0 = (bl+br)/2 + bumpx
    y0 = (bt+bb)/2 + bumpy
    dx = (br-bl)/20
    dy = (bt-bb)/20   # negative value so adding dy moves curso up screen


    add_mat_lab(im, bw, bh, 'mat_bkg', x0-4*dx, y0+7*dy)
    add_mat_lab(im, bw, bh, 'mat_a',   x0-1.25*dx, y0+.75*dy)
    add_mat_lab(im, bw, bh, 'mat_b',   x0+1.5*dx, y0+3*dy)


    add_lc(im, bw, bh, 'lc_bkg', x0-bw/6, bt, 'SW')
    add_lc(im, bw, bh, 'lc_2',   x0-d1/2, y0,  'SW')
    add_lc(im, bw, bh, 'lc_2',   x0-d1/2-d2, y0,  'SW')
    add_lc(im, bw, bh, 'lc_2',   x0-d1/2-d2-d1, y0,  'NW')


    c1x1 = x0-d1/2; c1x2 = x0+d1/2
    y1 = y0-dy; y2=y0-dy
    add_dim(im, bw, bh, 'inc_a_x', c1x1, y1, c1x2, y2, 'S')

    c2x1 = c1x2; c2x2 = c2x1+d2
    y1 = y0+dy; y2=y0+dy
    add_dim(im, bw, bh, 'inc_b_x', c2x1, y1, c2x2, y2, 'N')

    if lev>=3:
        c3x1 = c2x2; c3x2 = c3x1+d1
        y1 = y0-dy/2; y2=y0-dy/2
        add_dim(im, bw, bh, 'inc_c_x', c3x1, y1, c3x2, y2, 'S')
        add_mat_lab(im, bw, bh, 'mat_c',   x0+1.5*dx, y0-d1-d2)

    if lev>=4:
        c4x1 = c3x2; c4x2 = c4x1+d2
        y1 = y0+dy/2; y2=y0+dy/2
        add_dim(im, bw, bh, 'inc_d_x', c4x1, y1, c4x2, y2, 'N')
        add_mat_lab(im, bw, bh, 'mat_d',   x0+1.5*dx, y0-d1-d2-d1)

    im.show()
    im.save(frt+fn_ext_anno)






def do_trapezoid(nbapp):
    un_x = 2000
    un_y = 2000
    # layer thicknesses
    inc_a_x = 500
    inc_a_y = 300
    slab_a_x = 1000
    slab_a_y = 300

    wguide3 = nbapp.make_structure( 'trapezoidal_rib', un_x, un_y,
                                slab_width=slab_a_x, slab_thickness=slab_a_y,
                                   rib_top_width=inc_a_x/2, rib_base_width=inc_a_x,
                                   rib_height=inc_a_y,
                            material_bkg=mat_bkg, material_a=mat_a,
                                material_b=mat_b, material_c=mat_c,
                            lc_bkg=.05, lc_refine_1=10, lc_refine_2=10)
    #wguide3.plot_mesh('trapezoidal_rib_wg')
    wguide3.plot_mesh('trapezoidal_rib-mesh-annotated')
    #wguide3.check_mesh()

def do_pedestal(nbapp):
    un_x = 3000
    un_y = 2000

    # layer thicknesses
    ped_top_w = 350
    ped_base_w = 500
    ped_h = 100

    pillar_x = 20
    pillar_y = 200

    slab_a_x = 2000
    slab_a_y = 300


    wguide1 = nbapp.make_structure('pedestal', un_x, un_y,
                                   ped_top_w = ped_top_w,
                                   ped_base_w = ped_base_w,
                                   ped_h = ped_h,
                                   slab_a_x=slab_a_x, slab_a_y=slab_a_y,
                            material_bkg=mat_bkg, material_a=mat_a,
                                   material_b=mat_b,
                                   material_c=mat_c,
                            pillar_x=pillar_x, pillar_y=pillar_y,
                            lc_bkg=.1, lc_refine_1=10, lc_refine_2=10)

    wguide1.plot_mesh('pedestal_wg')


def do_main():

    nbapp = numbat.NumBATApp()

    do_oneincl(nbapp)
    #do_twoincl(nbapp)
    #do_rib(nbapp)
    #do_slot(nbapp)
    #do_onion(nbapp)
    #do_trapezoid(nbapp)
    #do_pedestal(nbapp)


if __name__=='__main__':
    do_main()
