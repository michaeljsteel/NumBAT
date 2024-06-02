

# contour plots a single component specified by cc
def delme_plot_(ax, m_X, m_Y, l_fields, plps, cc):
    decorator = plps['decorator']
    cmap_signed = 'seismic'
    cmap_unsigned = 'OrRd'

    # Adjustments to the visible plot domain
    xlmi = plps.get('xlim_min', 0)
    xlma = plps.get('xlim_max', 0)
    ylmi = plps.get('ylim_min', 0)
    ylma = plps.get('ylim_max', 0)

    comp_label = cc.get_label()  # for plotting
    signed = cc.is_signed_field()

    cmap = cmap_unsigned
    if signed:
        cmap = cmap_signed

    # Extract and tidy the field
    field = (l_fields[cc._f_code]).copy().T  # transpose because arrays

    if cc.is_abs():
        field = np.abs(field)**2  # TODO: cleanup: plot |u| as |u|^2

    # if the data is all noise, just plot zeros
    plot_threshold = 1e-8
    if np.max(np.abs(field[~np.isnan(field)])) < plot_threshold:
        field = np.zeros(np.shape(field))

    # set imshow plot (and tick) range to match the input x and y domain
    if True or plps['ticks']:
        x0 = m_X[0, 0]
        x1 = m_X[0, -1]
        y0 = m_Y[0, 0]
        y1 = m_Y[-1, 0]

        xm = x0+x1
        ym = y0+y1
        # Convert to length units in microns
        extents = np.array([x0-xm/2, x1-xm/2, y0-ym/2, y1-ym/2])*1e6
    else:
        extents = None

    interp = None
    # interp='bilinear';

    # User requested specific limits for each component x, y or z
    if decorator.get_cmap_limits(cc._xyz) != None:
        (act_zlo, act_zhi) = decorator.get_cmap_limits(cc._xyz)
        tsnorm = mplcolors.TwoSlopeNorm(
            vmin=act_zlo, vmax=act_zhi, vcenter=(act_zlo+act_zhi)/2)

        im = ax.imshow(field, origin='lower', extent=extents,
                       interpolation=interp, cmap=cmap, norm=tsnorm)
    else:
        act_zlo = np.nanmin(field)
        act_zhi = np.nanmax(field)
        vma = max(abs(act_zlo), abs(act_zhi))

        if signed:
            vmi = -vma
        else:
            vmi = 0
        im = ax.imshow(field, origin='lower', extent=extents,
                       interpolation=interp, cmap=cmap, vmin=vmi, vmax=vma)

    if xlmi > 0 or xlma > 0:
        xmin, xmax = ax.get_xlim()
        width_x = xmax-xmin
        ax.set_xlim(xmin+xlmi*width_x, xmax-xlma*width_x)

    if ylmi > 0 or ylma > 0:
        ymin, ymax = ax.get_ylim()
        width_y = ymax-ymin
        ax.set_ylim(ymin+ylmi*width_y, ymax-ylma*width_y)

    plot_set_ticks(ax, plps, decorator)

    plot_set_title(ax, comp_label, plps, decorator)

    plot_set_axes_style(ax, plps, decorator)

    # plot_set_contours(ax, plops, decorator)

    # contours and colorbars
    # colorbar
    do_cbar = plps['colorbar']
    do_contours = plps['contours']

    if do_cbar or do_contours:
        if plps['contour_lst']:
            cbarticks = plps['contour_lst']
        else:
            nt = plps.get('num_ticks', 7)
            cbarticks = np.linspace(act_zlo, act_zhi, nt)

    if do_cbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size=decorator.get_axes_property('cbar_size'),
                                  pad=decorator.get_axes_property('cbar_pad'))
        cbar = plt.colorbar(im, cax=cax)
        cbar.set_ticks(cbarticks)
        cbarlabels = [f'{t:.2f}' for t in cbarticks]
        cbar.set_ticklabels(cbarlabels)
        cbar.ax.tick_params(labelsize=decorator.get_font_size('cbar_tick'))

    if do_contours:
        if np.max(np.abs(field[~np.isnan(field)])) > plot_threshold:
            CS2 = ax.contour(m_X, m_Y, field.T, levels=cbarticks,
                             #colors=mycolors[::-1],   # FIXME
                             linewidths=(1.5,))
            if do_cbar:
                cbar.add_lines(CS2)

    if decorator is not None:
        decorator.extra_axes_commands(ax)


def del_plot_one_component_quiver(ax, m_X, m_Y, v_fields, plps, cc):
    decorator = plps['decorator']

    quiver_points = plps.get('quiver_points', 20)

    # Adjustments to the visible plot domain
    xlmi = plps.get('xlim_min', 0)
    xlma = plps.get('xlim_max', 0)
    ylmi = plps.get('ylim_min', 0)
    ylma = plps.get('ylim_max', 0)

    n_pts_y, n_pts_x = m_X.shape

    # this could probably be chosen nicer
    quiver_skip_x = int(round(min(n_pts_x, n_pts_y) /
                        quiver_points * (1-xlmi-xlma)))
    # this could probably be chosen nicer
    quiver_skip_y = int(round(min(n_pts_x, n_pts_y) /
                        quiver_points * (1-ylmi-ylma)))

    # TODO: ensure the points are symmetric about the centre of the domain
    # is quiver_skip_x and _y around the right way given the .T operation?

    v_x_q = m_X.T[0::quiver_skip_x, 0::quiver_skip_y]
    v_y_q = m_Y.T[0::quiver_skip_x, 0::quiver_skip_y]

    # TODO: why no transpose on these fields?
    m_ReEx_q = v_fields['Fxr'][0::quiver_skip_x, 0::quiver_skip_y]
    m_ReEy_q = v_fields['Fyr'][0::quiver_skip_x, 0::quiver_skip_y]
    m_ImEx_q = v_fields['Fxi'][0::quiver_skip_x, 0::quiver_skip_y]
    m_ImEy_q = v_fields['Fyi'][0::quiver_skip_x, 0::quiver_skip_y]

    # convert to microns
    v_x_q_um = v_x_q*1e6
    v_y_q_um = v_y_q*1e6

    # centre at zero
    xm = v_x_q_um[-1, 0]+v_x_q_um[0, 0]
    ym = v_y_q_um[0, -1]+v_y_q_um[0, 0]
    v_x_q_um -= xm/2
    v_y_q_um -= ym/2


# Ignore all imaginary values. If there are significant imag values,
# then instaneous vector plots don't make much sense anyway
    m_arrcolour = np.sqrt(m_ReEx_q*m_ReEx_q + m_ReEy_q*m_ReEy_q)
    ax.quiver(v_x_q_um, v_y_q_um, m_ReEx_q, m_ReEy_q, m_arrcolour,
              linewidths=(0.2,), edgecolors=('k'), pivot='mid', headlength=5)  # length of the arrows

    ax.set_aspect('equal')
    # this step is needed because quiver doesn't seem
    ax.set_xlim(v_x_q_um[0, 0], v_x_q_um[-1, 0])
    # to use its input x and y vectors to set range limits
    ax.set_ylim(v_y_q_um[0, 0], v_y_q_um[0, -1])
    # clean this up so as to avoid seemingly circular calls following

    plot_set_ticks(ax, plps, decorator)

    comp_label = cc.get_label()
    plot_set_title(ax, comp_label, plps, decorator)

    plot_set_axes_style(ax, plps, decorator)

    if True or xlmi > 0 or xlma > 0:
        xmin, xmax = ax.get_xlim()
        width_x = xmax-xmin
        ax.set_xlim(xmin+xlmi*width_x, xmax-xlma*width_x)

    if True or ylmi > 0 or ylma > 0:
        ymin, ymax = ax.get_ylim()
        width_y = ymax-ymin
        ax.set_ylim(ymin+ylmi*width_y, ymax-ylma*width_y)

    if decorator is not None:
        decorator.extra_axes_commands(ax)


# def plot_strain_mode_i(self, ival):
    #     # TODO: this interpolation looks very old. Can we get strain directly from fortran?
    #     # Interpolate onto rectangular Cartesian grid

    #     m_Fx = self.m_ReFx + 1j*self.m_ImFx
    #     m_Fy = self.m_ReFy + 1j*self.m_ImFy
    #     m_Fz = self.m_ReFz + 1j*self.m_ImFz
    #     dx = self.v_x[1]-self.v_x[0]
    #     dy = self.v_y[1]-self.v_y[0]

    #     print('finding gradients')
    #     # TODO: Check that the axis choice checks out
    #     del_x_Fx = np.gradient(m_Fx, dx, axis=0)
    #     del_y_Fx = np.gradient(m_Fx, dy, axis=1)
    #     del_x_Fy = np.gradient(m_Fy, dx, axis=0)
    #     del_y_Fy = np.gradient(m_Fy, dy, axis=1)
    #     del_x_Fz = np.gradient(m_Fz, dx, axis=0)
    #     del_y_Fz = np.gradient(m_Fz, dy, axis=1)
    #     del_z_Fx = 1j*self.sim.q_AC*m_Fx
    #     del_z_Fy = 1j*self.sim.q_AC*m_Fy
    #     del_z_Fz = 1j*self.sim.q_AC*m_Fz

    #     return

    #     self.v_x = np.linspace(x_min, x_max, self.n_pts_x)
    #     # For now, get these avlues from v_x already figured out earlier.
    #     x_min = self.v_x[0]
    #     x_max = self.v_x[-1]
    #     n_pts_x = len(self.v_x)
    #     y_min = self.v_y[0]
    #     y_max = self.v_y[-1]
    #     n_pts_y = len(self.v_y)

    #     xy = list(zip(self.v_x6p, self.v_y6p))

    #     # This seems to be equivalent  to taking grid_x = self.m_Y, grid_y = self.m_X
    #     # CONFIRM!
    #     # grid_x, grid_y = np.mgrid[x_min:x_max:n_pts_x*1j, y_min:y_max:n_pts_y*1j]  #OLD CODE
    #     grid_x, grid_y = self.m_Y, self.m_X  # NEW CODE

    #     m_ReFx = interpolate.griddata(
    #         xy, v_Fx6p.real, (grid_x, grid_y), method='linear')
    #     m_ReFy = interpolate.griddata(
    #         xy, v_Fy6p.real, (grid_x, grid_y), method='linear')
    #     m_ReFz = interpolate.griddata(
    #         xy, v_Fz6p.real, (grid_x, grid_y), method='linear')
    #     m_ImFx = interpolate.griddata(
    #         xy, v_Fx6p.imag, (grid_x, grid_y), method='linear')
    #     m_ImFy = interpolate.griddata(
    #         xy, v_Fy6p.imag, (grid_x, grid_y), method='linear')
    #     m_ImFz = interpolate.griddata(
    #         xy, v_Fz6p.imag, (grid_x, grid_y), method='linear')
    #     m_AbsF = interpolate.griddata(
    #         xy, v_F6p.real, (grid_x, grid_y), method='linear')

    #     dx = grid_x[-1, 0] - grid_x[-2, 0]
    #     dy = grid_y[0, -1] - grid_y[0, -2]

    #     m_Fx = m_ReFx + 1j*m_ImFx
    #     m_Fy = m_ReFy + 1j*m_ImFy
    #     m_Fz = m_ReFz + 1j*m_ImFz
    #     m_Fx = m_Fx.reshape(n_pts_x, n_pts_y)
    #     m_Fy = m_Fy.reshape(n_pts_x, n_pts_y)
    #     m_Fz = m_Fz.reshape(n_pts_x, n_pts_y)
    #     m_AbsF = m_AbsF.reshape(n_pts_x, n_pts_y)

    #     m_ReFx = np.real(m_Fx)
    #     m_ReFy = np.real(m_Fy)
    #     m_ReFz = np.real(m_Fz)
    #     m_ImFx = np.imag(m_Fx)
    #     m_ImFy = np.imag(m_Fy)
    #     m_ImFz = np.imag(m_Fz)

    #     del_x_Fx = np.gradient(m_Fx, dx, axis=0)
    #     del_y_Fx = np.gradient(m_Fx, dy, axis=1)
    #     del_x_Fy = np.gradient(m_Fy, dx, axis=0)
    #     del_y_Fy = np.gradient(m_Fy, dy, axis=1)
    #     del_x_Fz = np.gradient(m_Fz, dx, axis=0)
    #     del_y_Fz = np.gradient(m_Fz, dy, axis=1)
    #     del_z_Fx = 1j*sim_wguide.q_AC*m_Fx
    #     del_z_Fy = 1j*sim_wguide.q_AC*m_Fy
    #     del_z_Fz = 1j*sim_wguide.q_AC*m_Fz

    #     # Flip y order as imshow has origin at top left
    #     del_mat = np.array([del_x_Ex[:, ::-1].real, del_x_Ey[:, ::-1].real, del_x_Ez[:, ::-1].real, del_x_Ex[:, ::-1].imag, del_x_Ey[:, ::-1].imag, del_x_Ez[:, ::-1].imag, del_y_Ex[:, ::-1].real, del_y_Ey[:, ::-1].real, del_y_Ez[:, ::-1].real,
    #                        del_y_Ex[:, ::-1].imag, del_y_Ey[:, ::-1].imag, del_y_Ez[:, ::-1].imag, del_z_Ex[:, ::-1].real, del_z_Ey[:, ::-1].real, del_z_Ez[:, ::-1].real, del_z_Ex[:, ::-1].imag, del_z_Ey[:, ::-1].imag, del_z_Ez[:, ::-1].imag])
    #     v_labels = ["Re($S_{xx}$)", "Re($S_{xy}$)", "Re($S_{xz}$)", "Im($S_{xx}$)", "Im($S_{xy}$)", "Im($S_{xz}$)", "Re($S_{yx}$)", "Re($S_{yy}$)", "Re($S_{yz}$)",
    #                 "Im($S_{yx}$)", "Im($S_{yy}$)", "Im($S_{yz}$)", "Re($S_{zx}$)", "Re($S_{zy}$)", "Re($S_{zz}$)", "Im($S_{zx}$)", "Im($S_{zy}$)", "Im($S_{zz}$)"]

    #     # stress field plots
    #     plt.clf()
    #     fig = plt.figure(figsize=(15, 30))
    #     for i_p, plot in enumerate(del_mat):
    #         ax = plt.subplot(6, 3, i_p+1)
    #         im = plt.imshow(plot.T)
    #         # no ticks
    #         plt.xticks([])
    #         plt.yticks([])
    #         # limits
    #         if xlim_min > 0:
    #             ax.set_xlim(xlim_min*n_points, (1-xlim_max)*n_points)
    #         if ylim_min > 0:
    #             ax.set_ylim((1-ylim_min)*n_points, ylim_max*n_points)
    #         # titles
    #         plt.title(v_labels[i_p], fontsize=decorator.get_font_size(
    #             'subplot_title'))
    #         # colorbar
    #         divider = make_axes_locatable(ax)
    #         cax = divider.append_axes("right", size="5%", pad=0.1)
    #         cbar = plt.colorbar(im, cax=cax, format='%.2e')
    #         if num_ticks:
    #             cbarticks = np.linspace(
    #                 np.min(plot), np.max(plot), num=num_ticks)
    #         elif ylim_min != 0:
    #             if xlim_min/ylim_min > 3:
    #                 cbarlabels = np.linspace(np.min(plot), np.max(plot), num=3)
    #             if xlim_min/ylim_min > 1.5:
    #                 cbarlabels = np.linspace(np.min(plot), np.max(plot), num=5)
    #             else:
    #                 cbarlabels = np.linspace(np.min(plot), np.max(plot), num=7)
    #         else:
    #             cbarlabels = np.linspace(np.min(plot), np.max(plot), num=7)
    #         cbar.set_ticks(cbarlabels)
    #         cbarlabels = ['%.2f' % t for t in cbarlabels]
    #         cbar.set_ticklabels(cbarlabels)
    #         if contours:
    #             if contour_lst:
    #                 cbarticks = contour_lst
    #             if np.max(np.abs(plot[~np.isnan(plot)])) > plot_threshold:
    #                 CS2 = ax.contour(
    #                     m_X, m_Y, plot.T, levels=cbarticks, colors=colors[::-1], linewidths=(1.5,))
    #             cbar.add_lines(CS2)
    #         cbar.ax.tick_params(labelsize=decorator.get_font_size('cbar_tick'))
    #     fig.set_tight_layout(True)
    #     n_str = ''
    #     if np.imag(sim_wguide.Eig_values[ival]) < 0:
    #         k_str = r'$\Omega/2\pi = %(re_k)f %(im_k)f i$ GHz' % \
    #             {'re_k': np.real(sim_wguide.Eig_values[ival]*1e-9),
    #              'im_k': np.imag(sim_wguide.Eig_values[ival]*1e-9)}
    #     else:
    #         k_str = r'$\Omega/2\pi = %(re_k)f + %(im_k)f i$ GHz' % \
    #             {'re_k': np.real(sim_wguide.Eig_values[ival]*1e-9),
    #              'im_k': np.imag(sim_wguide.Eig_values[ival]*1e-9)}
    #     plt.suptitle('Mode #' + str(ival) + '   ' + k_str + '   ' +
    #                  n_str, fontsize=decorator.get_font_size('title'))

    #     if pdf_png == 'png':
    #         plt.savefig('%(pre)sfields/%(s)s_S_field_%(i)i%(add)s.png' %
    #                     {'pre': prefix, 's': EM_AC, 'i': ival, 'add': suffix})
    #     elif pdf_png == 'pdf':
    #         plt.savefig('%(pre)sfields/%(s)s_S_field_%(i)i%(add)s.pdf' %
    #                     {'pre': prefix, 's': EM_AC, 'i': ival, 'add': suffix}, bbox_inches='tight')
    #     if not keep_plots_open:
    #         plt.close()


# for i_el in range(sim.n_msh_el):
        #     # triangles
        #     idx = np.arange(6*i_el, 6*(i_el+1))
        #     triangles = [[idx[0], idx[3], idx[5]],
        #                  [idx[1], idx[4], idx[3]],
        #                  [idx[2], idx[5], idx[4]],
        #                  [idx[3], idx[4], idx[5]]]
        #     self.v_triang6p.extend(triangles)
