import numpy as np

import matplotlib as plt
import scipy


import nbgmsh
import reporting
import femmesh
import plotting.plottools

def _make_dielectric_plotter(wguide, optical_props, as_epsilon, n_points):


    v_neffeps = optical_props.v_refindexn  # mapping from material index to refractive index

    if as_epsilon:
        v_neffeps = v_neffeps**2
        nm_eng = 'Dielectric constant'
        nm_math=r'$\epsilon(\vec x)$'
        fname_suffix='dielectric_constant'
    else:
        nm_eng = 'Refractive index'
        nm_math=r'$n(\vec x)$'
        fname_suffix='refractive_index'

    fsfp = femmesh.FEMScalarFieldPlotter(wguide, n_points)

    unit=''

    fsfp.setup_scalar_properties(nm_eng, unit, nm_math, fname_suffix)
    fsfp.fill_quantity_by_material_index(v_neffeps)

    return fsfp

def get_structure_plotter_stiffness(wguide, c_I, c_J, n_points=500):
    """Make plotter for arbitrary 1D and 2D slices of the elastic stiffness."""

    if c_I not in range(1,7) or c_J not in range(1,7):
        reporting.report_and_exit('Stiffness tensor indices c_I, c_J must be in the range 1..6.')

    v_stiff = np.zeros(5) # fill me

    fsfp = femmesh.FEMScalarFieldPlotter(wguide, n_points)
    qname = 'Stiffness $c_{'+f'{c_I},{c_J}' +'}$'
    suffname = f'stiffness_c_{c_I}{c_J}'
    fsfp.set_quantity_name(qname, suffname)
    fsfp.fill_scalar_by_material_index(v_stiff)
    return fsfp


def get_structure_plotter_acoustic_velocity(wguide, d_materials, vel_index=0, n_points=500):
    """Make plotter for arbitrary 1D and 2D slices of the elastic acoustic phase speed.

    Args:
        vel_index (0,1,2): Index of the elastic mode phase speed to plot.

    Currently only works for isotropic materials."""

    v_mats = list(d_materials.values())
    v_acvel = np.zeros([len(v_mats),3])

    for i in range(len(v_mats)): # get all 3 phase velocities for each material
        if v_mats[i].has_elastic_properties():
            v_acvel[i,:] = v_mats[i].Vac_phase()

    fsfp = femmesh.FEMScalarFieldPlotter(wguide, n_points)
    fsfp.setup_vector_properties(3, 'Elastic velocity', '[km/s]', r'$v_i$',
                                    [r'$v_0$', r'$v_1$', r'$v_2$'],
                                    'elastic_velocity', ['v0', 'v1', 'v2'])

    fsfp.fill_quantity_by_material_index(v_acvel)
    return fsfp



def plot_refractive_index_profile_rough(mesh_mail_fname, d_materials,
                                        prefix, n_points = 200, as_epsilon=False):
        print('\n\nPlotting ref index')

        mail = nbgmsh.MailData(mesh_mail_fname)
        v_x, v_y = mail.v_centx, mail.v_centy
        v_elt_indices = mail.v_elts[:,-1]  # the elt number column
        v_refindex = 0*v_x

        print('Mesh props', mail.n_msh_pts, mail.n_msh_elts)

        uniq_elts = set(list(v_elt_indices))


        # find ref index at each centroid
        v_elt_refindex = np.zeros(len(uniq_elts))

        for i in range(len(v_elt_refindex)):
            v_elt_refindex[i] = np.real(list(d_materials.values())[i].refindex_n)


        for i,elt in enumerate(v_elt_indices):
            v_refindex[i] = v_elt_refindex[elt-1]  # the type of element is labelled by gmsh from 1.

        # Now we have an irregular x,y,n array to interpolate onto.


        # Construct a regular rect array with n_pts_x * n_pts_y ~ n_points**2
        # and with approximately square pixels
        x_min = np.min(v_x)
        x_max = np.max(v_x)
        y_min = np.min(v_y)
        y_max = np.max(v_y)

        area = abs((x_max-x_min)*(y_max-y_min))
        n_pts_x = int(n_points*abs(x_max-x_min)/np.sqrt(area))
        n_pts_y = int(n_points*abs(y_max-y_min)/np.sqrt(area))

        v_regx = np.linspace(x_min, x_max, n_pts_x)
        v_regy = np.linspace(y_min, y_max, n_pts_y)
        m_regx, m_regy = np.meshgrid(v_regx, v_regy)

        xy_in = np.array([v_x, v_y]).T
        xy_out = np.vstack([m_regx.ravel(), m_regy.ravel()]).T


        v_regindex = scipy.interpolate.griddata(xy_in, v_refindex, xy_out).reshape([n_pts_y, n_pts_x])
        fig, ax = plt.subplots()

        #v_regindex = np.where(v_regindex==0, 1, v_regindex)
        v_regindex = np.nan_to_num(v_regindex, nan=1.0)

        if as_epsilon:
            v_regindex = v_regindex**2
            fig.suptitle('Dielectric constant')
        else:
            fig.suptitle('Refractive index')

        cmap='cool'
        cf=ax.imshow(v_regindex, cmap=cmap, vmin=1.0, vmax=np.nanmax(v_regindex), origin='lower',
                     extent = [x_min, x_max, y_min, y_max])
        #cf=ax.contourf(m_regx, m_regy, v_regindex, cmap=cmap, vmin=1.0, vmax=np.nanmax(v_regindex))
        ax.set_xlabel(r'$x$ [μm]')
        ax.set_ylabel(r'$y$ [μm]')
        cb = fig.colorbar(cf)
        cf.set_clim(1,np.nanmax(v_regindex))
        cb.outline.set_linewidth(.5)
        cb.outline.set_color('gray')


        plottools.save_and_close_figure(fig, prefix+'refn.png')




