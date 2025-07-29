import numpy as np
import csv
from scipy import interpolate
from scipy.constants import epsilon_0 as SI_permittivity_eps0


def grad_u(dx, dy, u_mat, q_AC):
    """Take the gradient of field as well as of conjugate of field."""

    m_ux = u_mat[0]
    m_uy = u_mat[1]
    m_uz = u_mat[2]
    del_x_ux = np.gradient(m_ux, dx, axis=0)
    del_y_ux = np.gradient(m_ux, dy, axis=1)
    del_x_uy = np.gradient(m_uy, dx, axis=0)
    del_y_uy = np.gradient(m_uy, dy, axis=1)
    del_x_uz = np.gradient(m_uz, dx, axis=0)
    del_y_uz = np.gradient(m_uz, dy, axis=1)
    del_z_ux = 1j * q_AC * m_ux
    del_z_uy = 1j * q_AC * m_uy
    del_z_uz = 1j * q_AC * m_uz
    del_x_ux_star = np.gradient(np.conj(m_ux), dx, axis=0)
    del_y_ux_star = np.gradient(np.conj(m_ux), dy, axis=1)
    del_x_uy_star = np.gradient(np.conj(m_uy), dx, axis=0)
    del_y_uy_star = np.gradient(np.conj(m_uy), dy, axis=1)
    del_x_uz_star = np.gradient(np.conj(m_uz), dx, axis=0)
    del_y_uz_star = np.gradient(np.conj(m_uz), dy, axis=1)
    del_z_ux_star = -1j * q_AC * np.conj(m_ux)
    del_z_uy_star = -1j * q_AC * np.conj(m_uy)
    del_z_uz_star = -1j * q_AC * np.conj(m_uz)

    del_u_mat = np.array(
        [
            [del_x_ux, del_x_uy, del_x_uz],
            [del_y_ux, del_y_uy, del_y_uz],
            [del_z_ux, del_z_uy, del_z_uz],
        ]
    )
    del_u_mat_star = np.array(
        [
            [del_x_ux_star, del_x_uy_star, del_x_uz_star],
            [del_y_ux_star, del_y_uy_star, del_y_uz_star],
            [del_z_ux_star, del_z_uy_star, del_z_uz_star],
        ]
    )

    return del_u_mat, del_u_mat_star


def comsol_fields(data_file, n_points, mode_index=0):
    """Load Comsol field data on (assumed) grid mesh."""

    with open(data_file, "rt", encoding="ascii") as csvfile:
        spamreader = csv.reader(csvfile, delimiter=" ")  # , quotechar='|')
        for _ in range(9):  # skip header
            next(spamreader)
        x_coord = []
        y_coord = []
        u_x = []
        u_y = []
        u_z = []
        for row in spamreader:
            row = [_f for _f in row if _f]
            row = [float(x) for x in row]
            x_coord.append(row[0])
            y_coord.append(row[1])
            u_x.append(row[(mode_index * 6) + 2] + 1j * row[(mode_index * 6) + 3])
            u_y.append(row[(mode_index * 6) + 4] + 1j * row[(mode_index * 6) + 5])
            u_z.append(row[(mode_index * 6) + 6] + 1j * row[(mode_index * 6) + 7])

    x_coord = np.array(x_coord).reshape(n_points, n_points)
    y_coord = np.array(y_coord).reshape(n_points, n_points)
    x_coord = np.swapaxes(x_coord, 0, 1)
    y_coord = np.swapaxes(y_coord, 0, 1)
    u_x = np.array(u_x).reshape(n_points, n_points)
    u_y = np.array(u_y).reshape(n_points, n_points)
    u_z = np.array(u_z).reshape(n_points, n_points)
    u_x = np.swapaxes(u_x, 0, 1)
    u_y = np.swapaxes(u_y, 0, 1)
    u_z = np.swapaxes(u_z, 0, 1)
    field_mat = np.array([u_x, u_y, u_z])

    return x_coord, y_coord, field_mat


def interp_py_fields(
    sim_EM_pump,
    sim_EM_Stokes,
    sim_AC,
    q_AC,
    n_points,
    EM_mode_index_pump=0,
    EM_mode_index_Stokes=0,
    AC_mode_index=0,
):
    """Interpolate fields from FEM mesh to square grid."""

    # Trim EM fields to non-vacuum area where AC modes are defined
    # n_modes_EM = sim_EM_pump.n_modes
    # n_modes_AC = sim_AC.n_modes
    n_msh_el_AC = sim_AC.n_msh_el
    ncomps = 3
    nnodes = 6
    trimmed_EM_field_p = np.zeros((ncomps, nnodes, n_msh_el_AC), dtype=complex)
    trimmed_EM_field_S = np.zeros((ncomps, nnodes, n_msh_el_AC), dtype=complex)
    trimmed_EM_n = np.zeros((1, nnodes, n_msh_el_AC), dtype=complex)
    for el in range(n_msh_el_AC):
        new_el = sim_AC.el_convert_tbl[el]
        for n in range(nnodes):
            for x in range(ncomps):
                trimmed_EM_field_p[x, n, el] = sim_EM_pump.fem_evecs[
                    x, n, EM_mode_index_pump, new_el
                ]
                trimmed_EM_field_S[x, n, el] = sim_EM_Stokes.fem_evecs[
                    x, n, EM_mode_index_Stokes, new_el
                ]
            trimmed_EM_n[0, n, el] = sim_EM_pump.ls_material[0, n, new_el]

    # field mapping
    x_tmp = []
    y_tmp = []
    for i in np.arange(sim_AC.n_msh_pts):
        x_tmp.append(sim_AC.v_mshpt_xy[0, i])
        y_tmp.append(sim_AC.v_mshpt_xy[1, i])
    x_min = np.min(x_tmp)
    x_max = np.max(x_tmp)
    y_min = np.min(y_tmp)
    y_max = np.max(y_tmp)
    # area = abs((x_max-x_min)*(y_max-y_min))
    n_pts_x = n_points
    n_pts_y = n_points
    v_x = np.zeros(n_pts_x * n_pts_y)
    v_y = np.zeros(n_pts_x * n_pts_y)
    i = 0
    for x in np.linspace(x_min, x_max, n_pts_x):
        for y in np.linspace(y_max, y_min, n_pts_y):
            v_x[i] = x
            v_y[i] = y
            i += 1
    v_x = np.array(v_x)
    v_y = np.array(v_y)

    # unrolling data for the interpolators
    m_elnd_to_mshpt = sim_AC.m_elnd_to_mshpt.T
    v_mshpt_xy = sim_AC.v_mshpt_xy.T

    # dense triangulation with multiple points
    v_x6p = np.zeros(6 * sim_AC.n_msh_el)
    v_y6p = np.zeros(6 * sim_AC.n_msh_el)
    v_ux6p = np.zeros(6 * sim_AC.n_msh_el, dtype=np.complex128)
    v_uy6p = np.zeros(6 * sim_AC.n_msh_el, dtype=np.complex128)
    v_uz6p = np.zeros(6 * sim_AC.n_msh_el, dtype=np.complex128)
    v_Ex6p_E_p = np.zeros(6 * sim_AC.n_msh_el, dtype=np.complex128)
    v_Ey6p_E_p = np.zeros(6 * sim_AC.n_msh_el, dtype=np.complex128)
    v_Ez6p_E_p = np.zeros(6 * sim_AC.n_msh_el, dtype=np.complex128)
    v_Ex6p_E_S = np.zeros(6 * sim_AC.n_msh_el, dtype=np.complex128)
    v_Ey6p_E_S = np.zeros(6 * sim_AC.n_msh_el, dtype=np.complex128)
    v_Ez6p_E_S = np.zeros(6 * sim_AC.n_msh_el, dtype=np.complex128)
    v_n = np.zeros(6 * sim_AC.n_msh_el, dtype=np.complex128)
    # v_triang6p = []

    i = 0
    for i_el in np.arange(sim_AC.n_msh_el):
        for i_node in np.arange(6):
            # index for the coordinates
            i_ex = m_elnd_to_mshpt[i_el, i_node] - 1
            # values
            v_x6p[i] = v_mshpt_xy[i_ex, 0]
            v_y6p[i] = v_mshpt_xy[i_ex, 1]
            v_ux6p[i] = sim_AC.fem_evecs[0, i_node, AC_mode_index, i_el]
            v_uy6p[i] = sim_AC.fem_evecs[1, i_node, AC_mode_index, i_el]
            v_uz6p[i] = sim_AC.fem_evecs[2, i_node, AC_mode_index, i_el]
            v_Ex6p_E_p[i] = trimmed_EM_field_p[0, i_node, i_el]
            v_Ey6p_E_p[i] = trimmed_EM_field_p[1, i_node, i_el]
            v_Ez6p_E_p[i] = trimmed_EM_field_p[2, i_node, i_el]
            v_Ex6p_E_S[i] = trimmed_EM_field_S[0, i_node, i_el]
            v_Ey6p_E_S[i] = trimmed_EM_field_S[1, i_node, i_el]
            v_Ez6p_E_S[i] = trimmed_EM_field_S[2, i_node, i_el]
            v_n[i] = trimmed_EM_n[0, i_node, i_el]
            i += 1

    xy = list(zip(v_x6p, v_y6p))
    grid_x, grid_y = np.mgrid[
        x_min : x_max : n_pts_x * 1j, y_min : y_max : n_pts_y * 1j
    ]
    # pump mode
    m_ReEx_E = interpolate.griddata(
        xy, v_Ex6p_E_p.real, (grid_x, grid_y), method="cubic"
    )
    m_ReEy_E = interpolate.griddata(
        xy, v_Ey6p_E_p.real, (grid_x, grid_y), method="cubic"
    )
    m_ReEz_E = interpolate.griddata(
        xy, v_Ez6p_E_p.real, (grid_x, grid_y), method="cubic"
    )
    m_ImEx_E = interpolate.griddata(
        xy, v_Ex6p_E_p.imag, (grid_x, grid_y), method="cubic"
    )
    m_ImEy_E = interpolate.griddata(
        xy, v_Ey6p_E_p.imag, (grid_x, grid_y), method="cubic"
    )
    m_ImEz_E = interpolate.griddata(
        xy, v_Ez6p_E_p.imag, (grid_x, grid_y), method="cubic"
    )
    m_Ex_E = m_ReEx_E + 1j * m_ImEx_E
    m_Ey_E = m_ReEy_E + 1j * m_ImEy_E
    m_Ez_E = m_ReEz_E + 1j * m_ImEz_E
    m_Ex_E = m_Ex_E.reshape(n_pts_x, n_pts_y)
    m_Ey_E = m_Ey_E.reshape(n_pts_x, n_pts_y)
    m_Ez_E = m_Ez_E.reshape(n_pts_x, n_pts_y)
    E_mat_p = np.array([m_Ex_E, m_Ey_E, m_Ez_E])
    # Stokes mode
    m_ReEx_E = interpolate.griddata(
        xy, v_Ex6p_E_S.real, (grid_x, grid_y), method="cubic"
    )
    m_ReEy_E = interpolate.griddata(
        xy, v_Ey6p_E_S.real, (grid_x, grid_y), method="cubic"
    )
    m_ReEz_E = interpolate.griddata(
        xy, v_Ez6p_E_S.real, (grid_x, grid_y), method="cubic"
    )
    m_ImEx_E = interpolate.griddata(
        xy, v_Ex6p_E_S.imag, (grid_x, grid_y), method="cubic"
    )
    m_ImEy_E = interpolate.griddata(
        xy, v_Ey6p_E_S.imag, (grid_x, grid_y), method="cubic"
    )
    m_ImEz_E = interpolate.griddata(
        xy, v_Ez6p_E_S.imag, (grid_x, grid_y), method="cubic"
    )
    m_Ex_E = m_ReEx_E + 1j * m_ImEx_E
    m_Ey_E = m_ReEy_E + 1j * m_ImEy_E
    m_Ez_E = m_ReEz_E + 1j * m_ImEz_E
    m_Ex_E = m_Ex_E.reshape(n_pts_x, n_pts_y)
    m_Ey_E = m_Ey_E.reshape(n_pts_x, n_pts_y)
    m_Ez_E = m_Ez_E.reshape(n_pts_x, n_pts_y)
    E_mat_S = np.array([m_Ex_E, m_Ey_E, m_Ez_E])
    # AC mode
    m_Reux = interpolate.griddata(xy, v_ux6p.real, (grid_x, grid_y), method="cubic")
    m_Reuy = interpolate.griddata(xy, v_uy6p.real, (grid_x, grid_y), method="cubic")
    m_Reuz = interpolate.griddata(xy, v_uz6p.real, (grid_x, grid_y), method="cubic")
    m_Imux = interpolate.griddata(xy, v_ux6p.imag, (grid_x, grid_y), method="cubic")
    m_Imuy = interpolate.griddata(xy, v_uy6p.imag, (grid_x, grid_y), method="cubic")
    m_Imuz = interpolate.griddata(xy, v_uz6p.imag, (grid_x, grid_y), method="cubic")
    m_ux = m_Reux + 1j * m_Imux
    m_uy = m_Reuy + 1j * m_Imuy
    m_uz = m_Reuz + 1j * m_Imuz
    m_ux = m_ux.reshape(n_pts_x, n_pts_y)
    m_uy = m_uy.reshape(n_pts_x, n_pts_y)
    m_uz = m_uz.reshape(n_pts_x, n_pts_y)
    u_mat = np.array([m_ux, m_uy, m_uz])

    dx = grid_x[-1, 0] - grid_x[-2, 0]
    dy = grid_y[0, -1] - grid_y[0, -2]
    del_u_mat, del_u_mat_star = grad_u(dx, dy, u_mat, q_AC)

    m_Ren = interpolate.griddata(xy, v_n.real, (grid_x, grid_y), method="cubic")
    m_Imn = interpolate.griddata(xy, v_n.imag, (grid_x, grid_y), method="cubic")
    m_n = m_Ren + 1j * m_Imn
    m_n = m_n.reshape(n_pts_x, n_pts_y)

    return (
        n_pts_x,
        n_pts_y,
        dx,
        dy,
        E_mat_p,
        E_mat_S,
        u_mat,
        del_u_mat,
        del_u_mat_star,
        m_n,
    )


def grid_integral(
    m_n,
    sim_AC_structure,
    sim_AC_Omega_AC,
    n_pts_x,
    n_pts_y,
    dx,
    dy,
    E_mat_p,
    E_mat_S,
    u_mat,
    del_u_mat,
    del_u_mat_star,
    AC_mode_index,
):
    """Quadrature integration of AC energy density, AC loss (alpha), and PE gain."""

    # AC energy density integral
    F_AC_energy = 0
    for i in range(3):
        integrand_AC = np.conj(u_mat[i]) * u_mat[i] * sim_AC_structure.rho
        # do a 1-D integral over every row
        I_en = np.zeros(n_pts_x)
        for r in range(n_pts_x):
            I_en[r] = np.trapz(np.real(integrand_AC[r, :]), dx=dy)
        # then an integral over the result
        F_AC_energy += np.trapz(I_en, dx=dx)
        # Adding imag comp
        I_en = np.zeros(n_pts_x)
        for r in range(n_pts_x):
            I_en[r] = np.trapz(np.imag(integrand_AC[r, :]), dx=dy)
        F_AC_energy += 1j * np.trapz(I_en, dx=dx)
    energy_py = 2 * F_AC_energy * sim_AC_Omega_AC[AC_mode_index] ** 2

    # AC loss (alpha) integral
    F_alpha = 0
    for i in range(3):
        for k in range(3):
            for l in range(3):
                for j in range(3):
                    integrand = (
                        del_u_mat[i, j]
                        * del_u_mat_star[k, l]
                        * sim_AC_structure.elastic_props.eta_ijkl[i, j, k, l]
                    )
                    I_en = np.zeros(n_pts_x)
                    for r in range(n_pts_x):
                        I_en[r] = np.trapz(np.real(integrand[r, :]), dx=dy)
                    F_alpha += np.trapz(I_en, dx=dx)
                    I_en = np.zeros(n_pts_x)
                    for r in range(n_pts_x):
                        I_en[r] = np.trapz(np.imag(integrand[r, :]), dx=dy)
                    F_alpha += 1j * np.trapz(I_en, dx=dx)
    alpha_py = np.real(F_alpha * sim_AC_Omega_AC[AC_mode_index] ** 2 / energy_py)

    # PE gain integral

    F_PE = 0
    for i in range(3):
        for k in range(3):
            for l in range(3):
                for j in range(3):
                    # integrand_PE = acoustic_eps_effs[0]**2 * E_mat_p[j]*np.conj(E_mat_S[i])*sim_AC_structure.elastic_props.p_ijkl[i,j,k,l]*del_u_mat_star[k,l]
                    integrand_PE = (
                        m_n**4
                        * E_mat_p[j]
                        * np.conj(E_mat_S[i])
                        * sim_AC_structure.elastic_props.p_ijkl[i, j, k, l]
                        * del_u_mat_star[k, l]
                    )
                    I_en = np.zeros(n_pts_x)
                    for r in range(n_pts_x):
                        I_en[r] = np.trapz(np.real(integrand_PE[r, :]), dx=dy)
                    F_PE += np.trapz(I_en, dx=dx)
                    I_en = np.zeros(n_pts_x)
                    for r in range(n_pts_x):
                        I_en[r] = np.trapz(np.imag(integrand_PE[r, :]), dx=dy)
                    F_PE += 1j * np.trapz(I_en, dx=dx)
    Q_PE_py = F_PE * SI_permittivity_eps0

    return energy_py, alpha_py, Q_PE_py


def gain_python(
    sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC, comsol_data_file, comsol_mode_indices=1
):
    """Calculate interaction integrals and SBS gain in python.
    Load in acoustic mode displacement and calculate gain from this also.
    """

    print("gain python is out of action")
    return
    n_modes_EM = sim_EM_pump.n_modes
    # n_modes_AC = sim_AC.n_modes
    EM_mode_index_pump = 0
    EM_mode_index_Stokes = 0

    n_points = 100
    n_points_comsol_data = 100

    # acoustic_eps_effs =[]
    # for el_typ in range(sim_EM_pump.structure.n_mats_em):
    #     if el_typ+1 in sim_AC.typ_el_AC:
    #         acoustic_eps_effs.append(sim_EM_pump.v_refindexn[el_typ]**2)

    energy_py = np.zeros(comsol_mode_indices, dtype=np.complex128)
    alpha_py = np.zeros(comsol_mode_indices)
    Q_PE_py = np.zeros(
        (len(sim_EM_pump.eigs_kz), len(sim_EM_Stokes.eigs_kz), comsol_mode_indices),
        dtype=np.complex128,
    )
    energy_comsol = np.zeros(comsol_mode_indices, dtype=np.complex128)
    alpha_comsol = np.zeros(comsol_mode_indices)
    Q_PE_comsol = np.zeros(
        (len(sim_EM_pump.eigs_kz), len(sim_EM_Stokes.eigs_kz), comsol_mode_indices),
        dtype=np.complex128,
    )

    for AC_mode_index in range(comsol_mode_indices):  # Comsol data only contains some AC modes
        # Interpolate NumBAT FEM fields onto grid
        (
            n_pts_x,
            n_pts_y,
            dx,
            dy,
            E_mat_p,
            E_mat_S,
            u_mat,
            del_u_mat,
            del_u_mat_star,
            m_n,
        ) = interp_py_fields(
            sim_EM_pump,
            sim_EM_Stokes,
            sim_AC,
            q_AC,
            n_points,
            EM_mode_index_pump=EM_mode_index_pump,
            EM_mode_index_Stokes=EM_mode_index_Stokes,
            AC_mode_index=AC_mode_index,
        )

        # Carry out integration
        (
            energy_py[AC_mode_index],
            alpha_py[AC_mode_index],
            Q_PE_py[EM_mode_index_pump, EM_mode_index_Stokes, AC_mode_index],
        ) = grid_integral(
            m_n,
            sim_AC.structure,
            sim_AC.Omega_AC,
            n_pts_x,
            n_pts_y,
            dx,
            dy,
            E_mat_p,
            E_mat_S,
            u_mat,
            del_u_mat,
            del_u_mat_star,
            AC_mode_index,
        )

        # Load Comsol FEM fields onto grid - acoustic displacement fields
        x_coord, y_coord, u_mat_comsol = comsol_fields(
            comsol_data_file, n_points_comsol_data, mode_index=AC_mode_index
        )
        dx_comsol = x_coord[-1, 0] - x_coord[-2, 0]
        dy_comsol = y_coord[0, -1] - y_coord[0, -2]
        del_u_mat_comsol, del_u_mat_star_comsol = grad_u(dx, dy, u_mat_comsol, q_AC)

        # Carry out integration
        n_pts_x_comsol = n_points_comsol_data
        n_pts_y_comsol = n_points_comsol_data
        (
            energy_comsol[AC_mode_index],
            alpha_comsol[AC_mode_index],
            Q_PE_comsol[EM_mode_index_pump, EM_mode_index_Stokes, AC_mode_index],
        ) = grid_integral(
            m_n,
            sim_AC.structure,
            sim_AC.Omega_AC,
            n_pts_x_comsol,
            n_pts_y_comsol,
            dx_comsol,
            dy_comsol,
            E_mat_p,
            E_mat_S,
            u_mat_comsol,
            del_u_mat_comsol,
            del_u_mat_star_comsol,
            AC_mode_index,
        )

    # Note this is only the PE contribution to gain.
    gain_PE_py = (
        2
        * sim_EM_pump.omega_EM
        * sim_AC.Omega_AC[:comsol_mode_indices]
        * np.real(Q_PE_py * np.conj(Q_PE_py))
    )
    normal_fact_py = np.zeros((n_modes_EM, n_modes_EM, comsol_mode_indices), dtype=complex)
    gain_PE_comsol = (
        2
        * sim_EM_pump.omega_EM
        * sim_AC.Omega_AC[:comsol_mode_indices]
        * np.real(Q_PE_comsol * np.conj(Q_PE_comsol))
    )
    normal_fact_comsol = np.zeros((n_modes_EM, n_modes_EM, comsol_mode_indices), dtype=complex)
    for i in range(n_modes_EM):
        P1 = sim_EM_pump.EM_mode_power[i]
        for j in range(n_modes_EM):
            P2 = sim_EM_Stokes.EM_mode_power[j]
            for k in range(comsol_mode_indices):
                P3_py = energy_py[k]
                normal_fact_py[i, j, k] = P1 * P2 * P3_py * alpha_py[k]
                P3_comsol = energy_comsol[k]
                normal_fact_comsol[i, j, k] = P1 * P2 * P3_comsol * alpha_comsol[k]
    SBS_gain_PE_py = np.real(gain_PE_py / normal_fact_py)
    SBS_gain_PE_comsol = np.real(gain_PE_comsol / normal_fact_comsol)

    return SBS_gain_PE_py, alpha_py, SBS_gain_PE_comsol, alpha_comsol
