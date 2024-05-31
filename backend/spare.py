

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
