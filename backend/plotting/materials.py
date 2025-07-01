import tempfile
import subprocess

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as mplcm
import matplotlib.ticker as ticker

import numpy as np
import copy
import reporting

from bulkprops import solve_christoffel

def make_axes_square(ext0, ax, flip_x=False, flip_y=False):
    ext = 1.1 * ext0
    xlims = (ext, -ext) if flip_x else (-ext, ext)
    ylims = (ext, -ext) if flip_y else (-ext, ext)
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.set_aspect("equal")

def plot_bulk_dispersion_ivp(material, pref, label=None, show_poln=True,
                            flip_x=False, flip_y=False):

    ax_vp, ax_sl, ax_vg, ax_ivp_3d = None, None, None, None

    fig, ax_sl = plt.subplots(1, 1, dpi=300)
    setup_bulk_dispersion_2D_plot(ax_vp, ax_sl, ax_vg, ax_ivp_3d)

    cm = "cool"  # Color map for polarisation coding
    add_bulk_slowness_curves_to_axes(material, pref, fig, ax_vp, ax_sl, ax_vg, cm,
                                    show_poln, flip_x, flip_y)

    if label is not None:
        ax_sl.text( -0.1, 1.1, label, fontsize=14, style="italic", transform=ax_sl.transAxes)

    fname = pref + "-bulkdisp-ivp.png"
    fig.savefig(fname)
    plt.close(fig)
    return fname

def plot_bulk_dispersion_2D_all(material, pref, label=None, show_poln=True,
                            flip_x=False, flip_y=False):
    """Draw slowness surface 1/v_p(kappa) and ray surface contours in the horizontal
    (x-z) plane for the crystal axes current orientation.

    Solving the Christoffel equation: D C D^T u = -\rho v_p^2 u, for eigenvalue v_p
    and eigengector u.

    C is the Voigt form stiffness.
    D = [
    [kapx  0   0   0  kapz  kapy  ]
    [0   kapy  0   kapz 0   kapx  ]
    [0   0   kapz  kapy kapx  0]] where kap=(cos phi, 0, sin phi).

    """

    #fig, axs = setup_bulk_dispersion_2D_plot()

    #ax_sl, ax_vp, ax_vg, ax_ivp_3d = axs
    fig, axs = plt.subplots(2, 2, figsize=(7, 6), dpi=300)
    fig.subplots_adjust(hspace=0.35, wspace=0)

    ax_vp, ax_sl, ax_vg = axs[0, 0], axs[0, 1], axs[1, 0]

    axs[1, 1].set_axis_off()  # Hide axis 2,2

    axs[1, 1].remove()
    ax_ivp_3d = fig.add_subplot(2, 2, 4, projection="3d")

    setup_bulk_dispersion_2D_plot(ax_vp, ax_sl, ax_vg, ax_ivp_3d)


    cm = "cool"  # Color map for polarisation coding
    add_bulk_slowness_curves_to_axes(material, pref, fig, ax_vp, ax_sl, ax_vg, cm,
                                     show_poln, flip_x, flip_y)

    # markers to reproduce Auld Fig 7.2
    #c11=material.stiffness_c_IJ[1,1]
    #c12=material.stiffness_c_IJ[1,2]
    #c44=material.stiffness_c_IJ[4,4]
    #rho=material.rho
    #anA=2*c44/(c11-c12)
    #l1=np.sqrt(rho/(c11+c44*(1-1/anA)))
    #l2=np.sqrt(rho/c11)
    #l3=np.sqrt(rho/c44)
    #l4=np.sqrt(rho*anA/c44)
    #ax_sl.plot([l1/np.sqrt(2)*1e3], [l1/np.sqrt(2)*1e3], 'x')
    #ax_sl.plot([0], [l2*1e3], 'x')
    #ax_sl.plot([l3/np.sqrt(2)*1e3], [l3/np.sqrt(2)*1e3], 'x')
    #ax_sl.plot([l4/np.sqrt(2)*1e3], [l4/np.sqrt(2)*1e3], 'x')

    # markers to reproduce Auld Fig 7.3
    #c11=material._stiffness_c_IJ_orig[1,1]
    #c12=material._stiffness_c_IJ_orig[1,2]
    #c44=material._stiffness_c_IJ_orig[4,4]
    #rho=material.rho
    #anA=2*c44/(c11-c12)
#
#    l1=np.sqrt(rho/c44)*1e3
#    l3=np.sqrt(rho*anA/c44)*1e3
#    l4=np.sqrt(rho/c11)*1e3
#    l5=np.sqrt(3*rho*anA/((2+anA)*c44))*1e3
#    l6=np.sqrt(rho/(c11+c44*((anA-1)/anA)))*1e3
#    l7=np.sqrt(3*rho/(3*c11+4*c44*((anA-1)/anA)))*1e3
#    th=np.atan(1/np.sqrt(2))
#    ax_sl.plot([0], [l1], 'x')
#    ax_sl.plot([l1], [0], 'x')
#    ax_sl.plot([l3], [0], 'o')
#    ax_sl.plot([0], [-l4], 'v')
#    ax_sl.plot([-l6], [0], 'o')
#    ax_sl.plot([l5*np.cos(th)], [l5*np.sin(th)], 'o')
#    ax_sl.plot([l7*np.cos(th)], [l7*np.sin(th)], 'o')


    #label = self.material_name
    if label is not None:
        ax_vp.text( -0.1, 1.1, label, fontsize=14, style="italic", transform=ax_sl.transAxes)

    add_3d_dispersion_curves_to_axes(material, ax_ivp_3d)

    fname = pref + "-bulkdisp-all.png"
    fig.savefig(fname)
    plt.close(fig)
    return fname



def plot_bulk_dispersion_3D(material, pref):
    """Generate isocontour surfaces of the bulk dispersion in 3D k-space."""

    fig, axs = plt.subplots(1, 2, subplot_kw={"projection": "3d"})
    ax_vp, ax_vg = axs

    add_3d_dispersion_curves_to_axes(material, ax_vp, ax_vg)

    fname = pref + "-bulkdisp3D.png"
    fig.savefig(fname)
    plt.close(fig)
    return fname


def add_3d_dispersion_curves_to_axes(material, ax_ivp=None, ax_vg=None):
    ''' Draw phase and group velocity surfaces on 3D axes.

    :param ax_ivp: 3D axes for the slowness (reciprocal phase velocity).
    :param ax_vg: 3D axes for the group velocity.
    '''
    axs = []

    if ax_ivp is not None: axs.append(ax_ivp)
    if ax_vg is not None: axs.append(ax_vg)

    # Make data
    tpts = 50
    ppts = 100
    vphi = np.linspace(0, 2 * np.pi, ppts)
    vtheta = np.linspace(0, np.pi, tpts)

    ivx = np.zeros([tpts, ppts, 3])
    ivy = np.zeros([tpts, ppts, 3])
    ivz = np.zeros([tpts, ppts, 3])

    ivgx = np.zeros([tpts, ppts, 3])
    ivgy = np.zeros([tpts, ppts, 3])
    ivgz = np.zeros([tpts, ppts, 3])

    for ip, phi in enumerate(vphi):
        for itheta, theta in enumerate(vtheta):
            vkap = np.array(
                [
                    np.sin(theta) * np.cos(phi),
                    np.sin(theta) * np.sin(phi),
                    np.cos(theta),
                ]
            )
            v_vphase, vecs, v_vgroup = solve_christoffel(vkap, material.stiffness_c_IJ, material.rho)

            # slowness curve  eta(vkap) = 1/v_phase(vkap)
            ivx[itheta, ip, :] = vkap[0] / v_vphase
            ivy[itheta, ip, :] = vkap[1] / v_vphase
            ivz[itheta, ip, :] = vkap[2] / v_vphase

            ivgx[itheta, ip, :] = v_vgroup[:, 0]
            ivgy[itheta, ip, :] = v_vgroup[:, 1]
            ivgz[itheta, ip, :] = v_vgroup[:, 2]

    for i in range(3):
        if ax_ivp:
            ax_ivp.plot_surface( ivx[:, :, i], ivy[:, :, i], ivz[:, :, i], alpha=0.25)

        if ax_vg:
            ax_vg.plot_surface( ivgx[:, :, i], ivgy[:, :, i], ivgz[:, :, i], alpha=0.25)

    if ax_ivp:
        ax_ivp.set_xlabel(r"$1/v_x^{(p)}$ [s/km]", fontsize=8, labelpad=1)
        ax_ivp.set_ylabel(r"$1/v_y^{(p)}$ [s/km]", fontsize=8, labelpad=1)
        ax_ivp.set_zlabel(r"$1/v_z^{(p)}$ [s/km]", fontsize=8, labelpad=1)

    if ax_vg:
        ax_vg.set_xlabel(r"$v_x^{(g)}$ [km/s]", fontsize=8, labelpad=1)
        ax_vg.set_ylabel(r"$v_y^{(g)}$ [km/s]", fontsize=8, labelpad=1)
        ax_vg.set_zlabel(r"$v_z^{(g)}$ [km/s]", fontsize=8, labelpad=1)

    for ax in axs:
        for a in ("x", "y", "z"):
            ax.tick_params(axis=a, labelsize=8, pad=0)
        for t_ax in [ax.xaxis, ax.yaxis, ax.zaxis]:
            t_ax.line.set_linewidth(0.5)

        # ax.set_aspect('equal')





def plot_material_photoelastic_IJ(prefix, v_comps, mat):
    npts = 200

    fig, ax = plt.subplots(dpi=200, subplot_kw={'projection':'polar'}, figsize=(3,3))
    # fig.subplots_adjust(hspace=.35, wspace=0)

    d_p_vecs = {}  # map of "IJ" strings -> (el_IJ, np.zeros(npts))

    for s_elt in v_comps:
        if len(s_elt) !=2:
            reporting.report_and_exit(f'Bad length for photoelastic tensor index: {s_elt}.')

        el_IJ=(int(s_elt[0]), int(s_elt[1]))

        if not set(el_IJ) <= set(range(1,7)): # all elements are in range 1-6
        #if el_IJ[0] not in range(1,7) or el_IJ[1] not in range(1,7): #TODO express as a set operation
            reporting.report_and_exit(f'Bad range for photoelastic tensor index: {s_elt}.')

        d_p_vecs[s_elt] = (el_IJ, np.zeros(npts))

    mat0 = copy.deepcopy(mat)
    haty=np.array([0,1,0])

    v_phi = np.linspace(0.0, np.pi * 2, npts)
    for iphi, phi in enumerate(v_phi): # for all angles around the circle
        t_mat = copy.deepcopy(mat0)
        t_mat.rotate(haty, phi)

        for (k,v) in d_p_vecs.items(): # extract the desired p_IJ components
            (I,J) = v[0]
            v[1][iphi] = t_mat.photoel_p_IJ[I,J]

    for (k,v) in d_p_vecs.items():
        (I,J) = v[0]
        v_pIJ = v[1]

        lab = '$p_{' + f'{I},{J}' + '}$'
        plt.polar(v_phi, v_pIJ,  label=lab, lw=.5)

    ax.set_rmax(0.4)
    ax.set_rmin(-0.2)
    ax.set_rticks(np.arange(-0.2,0.4+.0001,.1))
    ax.tick_params(labelsize=8)

    ax.set_rlabel_position(90)
    ax.set_thetagrids(np.arange(0,360-1e-5,30))
    ax.grid(True)

    #ax.set_ylim(-.2,.4)
    #ax.legend(fontsize=6)

     # Move the legend to the right outside the plot
    # bbox_to_anchor: (x, y) coordinates in axes fraction.
    # (1.05, 1) means just outside the top right corner of the plot area.
    # loc='upper left' specifies where the anchor point of the legend box is.
    # In this case, the upper left corner of the legend box is anchored at (1.05, 1).
    ax.legend(fontsize=6, bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0.)

    #plt.tight_layout(rect=[0, 0, 1, 1]) # Adjust tight_layout to make space for legend
                                       # [left, bottom, right, top] in figure coordinates.
                                       # Adjust 'right' to be less than 1 to give space.
                                       # Or better, let tight_layout figure it out, but
                                       # sometimes manual adjustment is needed.



def make_crystal_axes_plot(pref, crystal_axes):
    """Build crystal coordinates diagram using call to external asymptote application."""

    fn = tempfile.NamedTemporaryFile(suffix=".asy", mode="w+t", delete=False)

    asy_cmds = asy_draw_crystal_axes(crystal_axes)
    fn.write(asy_cmds)
    fn.close()

    # run .asy
    subprocess.run(["asy", fn.name, "-o", f"{pref}-crystal"], check=False)
    fname_out = f"{pref}-crystal.png"
    return fname_out


def asy_draw_crystal_axes(crystal_axes):
    (va, vb, vc) = crystal_axes
    s_avec = "(" + ",".join(map(str, va)) + ")"
    s_bvec = "(" + ",".join(map(str, vb)) + ")"
    s_cvec = "(" + ",".join(map(str, vc)) + ")"

    s1 = r"""
settings.outformat='png';
settings.render=8;
import three;
import graph3;

size(2cm,0);
defaultpen(fontsize(7pt));
defaultpen(.2);

real axlen=1.25;
int arrsize=3;
real blen=.5;

//currentprojection=orthographic(1,1,1);
currentprojection=oblique;


draw(O--2X, black, Arrow3(arrsize), L=Label("$\hat{x}$", position=EndPoint));
draw(O--2Y, black, Arrow3(arrsize), L=Label("$\hat{y}$", position=EndPoint));
draw(O--3Z, black, Arrow3(arrsize), L=Label("$\hat{z}$", position=EndPoint));

draw(O-- -2X, gray);
draw(O-- -2Y, gray);
draw(O-- -2Z, gray);


//label("$\hat{x}$", 3X*1.1);
//label("$\hat{y}$", 3Y*1.1);
//label("$\hat{z}$", 3Z*1.1);

draw(box((-1,-.5,-2)*blen,(1,.5,2)*blen),blue);
"""

    s2 = f"""triple avec={s_avec};
triple bvec={s_bvec};
triple cvec={s_cvec};
"""

    s3 = """triple corig=(0,.0,0)*blen;
draw(corig--avec+corig, red, Arrow3(arrsize), L=Label("$c_x$", position=EndPoint));
draw(corig--bvec+corig, red, Arrow3(arrsize), L=Label("$c_y$", position=EndPoint));
draw(corig--cvec+corig, red, Arrow3(arrsize), L=Label("$c_z$", position=EndPoint));

triple k0=(1,-1,-1);
triple k1=k0+(0,0,2);

draw(k0--k1,green, Arrow3(arrsize), L=Label("$k$"));
"""

# used to be triple corig=(0,.5,.2)*blen;
    return s1 + s2 + s3


def setup_bulk_dispersion_2D_plot(ax_vp, ax_sl, ax_vg, ax_ivp3d):
    """Plots both slowness and ray normal contours."""

    if ax_vp is not None:
        ax_vp.set_xlabel(r"$v^{(p)}_{x}$ [km/s]")
        ax_vp.set_ylabel(r"$v^{(p)}_{z}$ [km/s]")

    if ax_sl is not None:
        ax_sl.set_xlabel(r"$1/v^{(p)}_{x}$ [s/km]")
        ax_sl.set_ylabel(r"$1/v^{(p)}_{z}$ [s/km]")

    if ax_vg is not None:
        ax_vg.set_xlabel(r"$v^{(g)}_{x}$ [km/s]")
        ax_vg.set_ylabel(r"$v^{(g)}_{z}$ [km/s]")

    axs = [ax for ax in (ax_vp, ax_sl, ax_vg) if ax is not None]

    for ax in axs:  # Don't write to axis 2,2
        ax.axhline(0, c="gray", lw=0.5)
        ax.axvline(0, c="gray", lw=0.5)
        ax.tick_params(width=0.5)
        for item in ( [ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(10)
        for t_ax in ["top", "bottom", "left", "right"]:
            ax.spines[t_ax].set_linewidth(0.5)
    #axs = ax_sl, ax_vp, ax_vg, ax_ivp3d
    #return fig, axs


def setup_bulk_dispersion_2D_plot_2x1():
    """Plots both slowness and ray normal contours."""

    fig, axs = plt.subplots(1, 1, figsize=(6, 4))
    # fig.subplots_adjust(hspace=.35, wspace=0)
    axs = (axs,)
    # ax_sl, ax_vg = axs
    ax_sl = axs[0]

    # ax_sl, ax_vp, ax_vg = axs[0,0], axs[0,1], axs[1,0]

    # axs[1,1].set_axis_off()  # Hide axis 2,2

    # axs[1,1].remove()
    # ax_ivp3d = fig.add_subplot(2,2,4, projection='3d')
    ax_sl.set_xlabel(r"$1/v^{(p)}_{x}$ [s/km]")
    ax_sl.set_ylabel(r"$1/v^{(p)}_{z}$ [s/km]")
    # ax_vp.set_xlabel(r'$v^{(p)}_{x}$ [s/km]')
    # ax_vp.set_ylabel(r'$v^{(p)}_{z}$ [s/km]')
    # ax_vg.set_xlabel(r'$v^{(g)}_{x}$ [km/s]')
    # ax_vg.set_ylabel(r'$v^{(g)}_{z}$ [km/s]')

    for ax in axs:  # Don't write to axis 2,2
        ax.axhline(0, c="gray", lw=0.5)
        ax.axvline(0, c="gray", lw=0.5)
        ax.tick_params(width=0.5)
        for item in (
            [ax.title, ax.xaxis.label, ax.yaxis.label]
            + ax.get_xticklabels()
            + ax.get_yticklabels()
        ):
            item.set_fontsize(12)
        for t_ax in ["top", "bottom", "left", "right"]:
            ax.spines[t_ax].set_linewidth(0.5)
    # axs = ax_sl, ax_vp, ax_vg, ax_ivp3d
    return fig, axs


def add_bulk_slowness_curves_to_axes(material, pref, fig, ax_vp, ax_sl, ax_vg, cm,
                                     show_poln=True, flip_x=False, flip_y=False):
    npolpts = 28
    npolskip = 10  # make bigger
    npts = npolpts * npolskip  # about 1000
    v_kphi = np.linspace(0.0, np.pi * 2, npts)
    v_vel = np.zeros([npts, 3])
    v_velc = np.zeros([npts, 3])
    v_vgx = np.zeros([npts, 3])
    v_vgz = np.zeros([npts, 3])


    cmm = mpl.colormaps[cm]
    with open(pref + "-bulkdisp.dat", "w") as fout:
        fout.write(
            "#phi    kapx     kapz      vl          vs1         vs2         vlx      vly       vlz     vs1x     vs1y      vs1z       vs2x     vs2y     vs2z      k.v1    k.v2   k.v3\n"
        )

        kapcomp = np.zeros(3)
        ycomp = np.zeros(3)
        for ik, kphi in enumerate(v_kphi):
            # kapx = np.cos(kphi)
            # kapz = np.sin(kphi)
            # kapy = 0.0
            vkap = np.array([np.cos(kphi), 0.0, np.sin(kphi)])

            fout.write(f"{kphi:.4f}  {vkap[0]:+.4f}  {vkap[2]:+.4f}  ")

            # solve_christoffel returns:
            # eigvecs are sorted by phase velocity
            # v_vphase[m]:   |vphase| of modes m=1 to 3
            # vecs[:,m]:     evecs of modes m=1 to 3
            # v_vgroup[m,:]  vgroup of mode m, second index is x,y,z
            t_stiffness = material.get_stiffness_for_kappa(vkap)

            #if ik==0:
            #    print('got stiff', t_stiffness)

            v_vphase, vecs, v_vgroup = solve_christoffel(vkap, t_stiffness, material.rho)

            v_vel[ik, :] = v_vphase  # phase velocity
            v_vgx[ik, :] = v_vgroup[:, 0]  # group velocity components
            v_vgz[ik, :] = v_vgroup[:, 2]

            ycomp = np.abs(vecs[1, :])  # $\unity \cdot u_i$
            kapcomp = np.abs( np.matmul(vkap, vecs))  # component of vkap along each evec
            v_velc[ik, :] = kapcomp  # phase velocity color by polarisation

            for iv in range(3):
                fout.write(f"{v_vphase[iv]*1000:10.4f}  ")
            for iv in range(3):
                fout.write( f"{vecs[0,iv]:7.4f}  {vecs[1,iv]:7.4f}   {vecs[2,iv]:7.4f}  ")
            fout.write(f"{kapcomp[0]:6.4f}  {kapcomp[1]:6.4f} {kapcomp[2]:6.4f}")

            fout.write("\n")

            if show_poln:
                # Draw polarisation ball and stick notations
                irad = 0.07 / v_vel[0, 0]  # length of polarisation sticks
                rad = 0.07 * v_vel[0, 0]  # length of polarisation sticks
                lwstick = 0.9
                srad = 3  # diameter of polarisation dots
                if ik % npolskip == 0:
                    for i in range(3):
                        radsl = 1 / v_vel[ik, i]
                        radvp = v_vel[ik, i]
                        polc = cmm(kapcomp[i])
                        polc = "gray"  # all black for now

                        ptm = radsl * np.array([np.cos(kphi), np.sin(kphi)])
                        pt0 = np.real(ptm - vecs[0:3:2, i] * irad)
                        pt1 = np.real(ptm + vecs[0:3:2, i] * irad)
                        if ax_sl:
                            ax_sl.plot( (pt0[0], pt1[0]), (pt0[1], pt1[1]), c=polc, lw=lwstick)
                            ax_sl.plot( ptm[0], ptm[1], "o", c=polc, markersize=srad * ycomp[i])

                        ptm = radvp * np.array([np.cos(kphi), np.sin(kphi)])
                        pt0 = np.real(ptm - vecs[0:3:2, i] * rad)
                        pt1 = np.real(ptm + vecs[0:3:2, i] * rad)

                        if ax_vp:
                            ax_vp.plot( (pt0[0], pt1[0]), (pt0[1], pt1[1]), c=polc, lw=lwstick)
                            ax_vp.plot( ptm[0], ptm[1], "o", c=polc, markersize=srad * ycomp[i])

    # the main curves for v_p, 1/v_p and v_g
    lw=.5
    for i in range(3):
        if ax_vp:
            ax_vp.scatter( np.cos(v_kphi) * v_vel[:, i], np.sin(v_kphi) * v_vel[:, i],
                      c=v_velc[:, i], vmin=0, vmax=1, s=0.5, cmap=cm, lw=lw)

        if ax_sl:
            ax_sl.scatter( np.cos(v_kphi) / v_vel[:, i], np.sin(v_kphi) / v_vel[:, i],
                      c=v_velc[:, i], vmin=0, vmax=1, s=0.5, cmap=cm, lw=lw)

        if ax_vg:
            ax_vg.scatter( v_vgx[:, i], v_vgz[:, i], c=v_velc[:, i], vmin=0, vmax=1, s=0.5, cmap=cm, lw=lw)

    # Tick location seems to need help here
    active_axes = []
    for ax in (ax_vp, ax_vg):
        if ax is not None:
            active_axes.extend([ax.xaxis, ax.yaxis])

    for tax in active_axes:
        tax.set_major_locator( ticker.MultipleLocator( 2.0 ))# , offset=0

    if ax_sl is not None:
        make_axes_square(np.abs(1 / v_vel).max(), ax_sl, flip_x, flip_y)
    if ax_vp is not None:
        make_axes_square(np.abs(v_vel).max(), ax_vp, flip_x, flip_y)
    if ax_vg is not None:
        make_axes_square(max(np.abs(v_vgx).max(), np.abs(v_vgz).max()), ax_vg, flip_x, flip_y)

    # Add radial speed grid
    v_theta = np.linspace(0, 2*np.pi, 300)
    for vr in range(1,8):
        if ax_vp is not None:
            ax_vp.plot( np.cos(v_theta) * vr , np.sin(v_theta) * vr, ':', c='gray', lw=.25)
        if ax_vg is not None:
            ax_vg.plot( np.cos(v_theta) * vr , np.sin(v_theta) * vr, ':', c='gray', lw=.25)

    for ivr in np.arange(0,.4,.05):
        if ax_sl is not None:
            ax_sl.plot( np.cos(v_theta) * ivr , np.sin(v_theta) * ivr, ':', c='gray', lw=.25)

    # fig.colorbar(mplcm.ScalarMappable(cmap=cm), ax=ax_vp, shrink=.5,
    #             pad=.025, location='top', label='$\hat{e} \cdot \hat{\kappa}$')


def add_bulk_slowness_curves_to_axes_2x1(
        material, pref, fig, ax_sl, ax_vp, cm, mat1or2
    ):
        npolpts = 28
        npolskip = 10  # make bigger
        npts = npolpts * npolskip  # about 1000
        v_kphi = np.linspace(0.0, np.pi * 2, npts)
        v_vel = np.zeros([npts, 3])
        v_velc = np.zeros([npts, 3])
        # v_vgx = np.zeros([npts, 3])
        # v_vgz = np.zeros([npts, 3])

        cmm = mpl.colormaps[cm]

        kapcomp = np.zeros(3)
        ycomp = np.zeros(3)
        for ik, kphi in enumerate(v_kphi):
            vkap = np.array([np.cos(kphi), 0.0, np.sin(kphi)])

            # sol   ve_christoffel returns:
            # eigvecs are sorted by phase velocity
            # v_vphase[m]:   |vphase| of modes m=1 to 3
            # vecs[:,m]:     evecs of modes m=1 to 3
            # v_vgroup[m,:]  vgroup of mode m, second index is x,y,z
            v_vphase, vecs, v_vgroup = solve_christoffel(vkap, material.stiffness_c_IJ, material.rho)

            v_vel[ik, :] = v_vphase  # phase velocity
            # v_vgx[ik, :] = v_vgroup[:,0]  # group velocity components
            # v_vgz[ik, :] = v_vgroup[:,2]

            ycomp = np.abs(vecs[1, :])  # $\unity \cdot u_i$
            kapcomp = np.abs(np.matmul(vkap, vecs))  # component of vkap along each evec
            # This causes both shear waves to have the same colour if there are two pure
            v_velc[ik, :] = kapcomp  # phase velocity color by polarisation.

            # Draw polarisation ball and stick notations
            irad = 0.07 / v_vel[0, 0]  # length of polarisation sticks
            rad = 0.07 * v_vel[0, 0]  # length of polarisation sticks
            lwstick = 0.9
            srad = 5  # diameter of polarisation dots
            if ik % npolskip == 0:
                for i in range(3):
                    radsl = 1 / v_vel[ik, i]
                    radvp = v_vel[ik, i]
                    polc = cmm(kapcomp[i])
                    polc = "k"  # all black for now

                    ptm = radsl * np.array([np.cos(kphi), np.sin(kphi)])
                    pt0 = np.real(ptm - vecs[0:3:2, i] * irad)
                    pt1 = np.real(ptm + vecs[0:3:2, i] * irad)
                    ax_sl.plot((pt0[0], pt1[0]), (pt0[1], pt1[1]), c=polc, lw=lwstick)
                    ax_sl.plot(ptm[0], ptm[1], "o", c=polc, markersize=srad * ycomp[i])

                    ptm = radvp * np.array([np.cos(kphi), np.sin(kphi)])
                    pt0 = np.real(ptm - vecs[0:3:2, i] * rad)
                    pt1 = np.real(ptm + vecs[0:3:2, i] * rad)

                    # ax_vp.plot((pt0[0], pt1[0]), (pt0[1], pt1[1]), c=polc, lw=lwstick)
                    # ax_vp.plot(ptm[0], ptm[1], 'o', c=polc, markersize=srad*ycomp[i])

        # the main curves for 1/v_p and v_g
        for i in range(3):
            ax_sl.scatter(
                np.cos(v_kphi) / v_vel[:, i],
                np.sin(v_kphi) / v_vel[:, i],
                c=v_velc[:, i],
                vmin=0,
                vmax=1,
                s=0.5,
                cmap=cm,
            )

            # ax_vp.scatter(np.cos(v_kphi)*v_vel[:, i], np.sin(v_kphi) *
            #            v_vel[:, i], c=v_velc[:, i], vmin=0, vmax=1, s=0.5, cmap=cm)

            # ax_vg.scatter(v_vgx[:,i], v_vgz[:,i],  c=v_velc[:, i], vmin=0, vmax=1, s=0.5, cmap=cm)

        # Tick location seems to need help here
        # for tax in [ax_vp.xaxis, ax_vp.yaxis, ax_vg.xaxis, ax_vg.yaxis]:
        #   tax.set_major_locator(ticker.MultipleLocator(2.0, offset=0))

        make_axes_square(np.abs(1 / v_vel).max(), ax_sl)
        # make_axes_square(np.abs(v_vel).max(), ax_vp)
        # make_axes_square(max(np.abs(v_vgx).max(), np.abs(v_vgz).max()), ax_vg)

        cbar = fig.colorbar(
            mplcm.ScalarMappable(cmap=cm),
            ax=ax_sl,
            shrink=0.5,
            pad=0.025,
            location="right",
        )
        cbar.ax.tick_params(labelsize=6, width=0.25)

        cbar.outline.set_linewidth(1)
        cbar.set_label(
            label=f"Mat {mat1or2} " + r"$\hat{e} \cdot \hat{\kappa}$", fontsize=10
        )


def compare_bulk_dispersion(mat1, mat2, pref):

    fig, axs = setup_bulk_dispersion_2D_plot_2x1()

    # ax_sl, ax_vg = axs
    ax_sl = axs[0]
    ax_vg = None

    cm1 = "cool"  # Color map for polarisation coding
    cm2 = "autumn"  # Color map for polarisation coding

    add_bulk_slowness_curves_to_axes_2x1(mat1, pref + "_mat1", fig, ax_sl, ax_vg, cm1, 1)
    add_bulk_slowness_curves_to_axes_2x1(mat2, pref + "_mat2", fig, ax_sl, ax_vg, cm2, 2)

    ax_sl.text(
        0.05,
        1.15,
        f"Mat 1: {mat1.material_name}",
        fontsize=14,
        style="italic",
        transform=ax_sl.transAxes,
    )
    ax_sl.text(
        0.05,
        1.05,
        f"Mat 2: {mat2.material_name}",
        fontsize=14,
        style="italic",
        transform=ax_sl.transAxes,
    )

    plt.savefig(pref + "-compare-bulkdisp.png")


