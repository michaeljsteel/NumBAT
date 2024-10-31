#include "numbat_decl.h"
! Calculate the elastic energy v_energy integral of an AC mode with itself using
! analytic expressions for basis function v_energys on linear elements.
!

! E = 2 \Omega^2 \int dxdy \rho |u|^2


subroutine ac_mode_energy_analytic_old (n_modes, n_msh_el, n_msh_pts, &
   elnd_to_mshpt, v_nd_xy, n_elt_mats, v_el_material,  &
   rho, Omega_AC, soln_ac_u, v_energy, errco, emsg)

   use numbatmod

   use class_TriangleIntegrators


   integer(8) n_modes, md_i
   integer(8) n_msh_el, n_msh_pts, n_elt_mats
   integer(8) v_el_material(n_msh_el)
   integer(8) elnd_to_mshpt(P2_NODES_PER_EL,n_msh_el)
   double precision v_nd_xy(2,n_msh_pts)
   complex(8) rho(n_elt_mats)

   complex(8) Omega_AC(n_modes)
   complex(8) soln_ac_u(3,P2_NODES_PER_EL,n_modes,n_msh_el)
   complex(8), dimension(n_modes) :: v_energy
   integer(8), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg

   ! Locals
   integer(8) nod_el_p(P2_NODES_PER_EL)
   double precision nds_xy(2,P2_NODES_PER_EL)
   complex(8) bas_ovrlp(P2_NODES_PER_EL,P2_NODES_PER_EL)
   complex(8) ac_U, ac_Ustar
   integer(8) j, j1, typ_e
   integer(8) i_el, xyz_i
   integer(8) bf_i, bf_j
   complex(8) rho_el

   !type(PyFrontEnd) frontend
   type(NBError) nberr
   !type(QuadIntegrator) quadint


!f2py intent(in) n_modes, n_msh_el, n_msh_pts, P2_NODES_PER_EL, elnd_to_mshpt
!f2py intent(in) v_el_material, x, n_elt_mats, rho
!f2py intent(in) soln_ac_u, Omega_AC
!
!f2py depend(elnd_to_mshpt) P2_NODES_PER_EL, n_msh_el
!f2py depend(v_el_material) n_msh_pts
!f2py depend(v_nd_xy) n_msh_pts
!f2py depend(soln_ac_u) P2_NODES_PER_EL, n_modes, n_msh_el
!f2py depend(Omega_AC) n_modes
!f2py depend(rho) n_elt_mats
!
!f2py intent(out) v_energy
!

   errco = 0
   call nberr%reset()

   !call quadint%setup_reference_quadratures()

   !call frontend%init_from_py(n_msh_el, n_msh_pts, elnd_to_mshpt, v_nd_xy, nberr)
   RET_ON_NBERR_UNFOLD(nberr)


   v_energy = D_ZERO


   do i_el=1,n_msh_el
      typ_e = v_el_material(i_el)

      do j=1,P2_NODES_PER_EL
         j1 = elnd_to_mshpt(j,i_el)
         nod_el_p(j) = j1
         nds_xy(:,j) = v_nd_xy(:,j1)

      enddo
      rho_el = rho(typ_e)

      !call frontend%nodes_at_el(i_el, nds_xy)

      call mat_el_energy_rho (nds_xy, rho_el, bas_ovrlp)

      ! Having calculated v_energy of basis functions on element
      ! now multiply by specific field values for modes of interest.
      do md_i=1,n_modes

         do bf_i=1,P2_NODES_PER_EL
            do bf_j=1,P2_NODES_PER_EL
               do xyz_i=1,3

                  ac_Ustar = conjg(soln_ac_u(xyz_i,bf_i,md_i,i_el))
                  ac_U = soln_ac_u(xyz_i,bf_j,md_i,i_el)

                  v_energy(md_i) = v_energy(md_i) + ac_Ustar*ac_U*bas_ovrlp(bf_i,bf_j)

               enddo
            enddo
         enddo

      enddo
   enddo

   v_energy = 2.0 * Omega_AC**2 * v_energy

end subroutine ac_mode_energy_analytic_old


! Calculate the elastic energy v_energy integral of an AC mode with itself using
! analytic expressions for basis function v_energys on linear elements.
!

! E = 2 \Omega^2 \int dxdy \rho |u|^2


subroutine ac_mode_energy_analytic (n_modes, n_msh_el, n_msh_pts, &
   elnd_to_mshpt, v_nd_xy, n_elt_mats, v_el_material,  &
   rho, Omega_AC, soln_ac_u, v_energy_r, errco, emsg)

   use numbatmod
   use class_TriangleIntegrators
   use class_BasisFunctions

   integer(8) n_modes, md_i
   integer(8) n_msh_el, n_msh_pts, n_elt_mats
   integer(8) v_el_material(n_msh_el)
   integer(8) elnd_to_mshpt(P2_NODES_PER_EL,n_msh_el)
   double precision v_nd_xy(2,n_msh_pts)
   complex(8) rho(n_elt_mats)

   complex(8) Omega_AC(n_modes)
   complex(8) soln_ac_u(3,P2_NODES_PER_EL,n_modes,n_msh_el)

   double precision, dimension(n_modes), intent(out) :: v_energy_r

   integer(8), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg

   ! Locals
   complex(8), dimension(n_modes)  :: v_energy
   double precision nds_xy(2,P2_NODES_PER_EL)
   double precision m_int_p2p2(P2_NODES_PER_EL,P2_NODES_PER_EL)
   complex(8) ac_U, ac_Ustar
   integer(8) typ_e
   integer(8) i_el, xyz_i, bf_i, bf_j
   complex(8) rho_el

   type(NBError) nberr
   type(PyFrontEnd) frontend
   type(BasisFunctions) basfuncs


!f2py intent(in) n_modes, n_msh_el, n_msh_pts, P2_NODES_PER_EL, elnd_to_mshpt
!f2py intent(in) v_el_material, x, n_elt_mats, rho
!f2py intent(in) soln_ac_u, Omega_AC

!f2py depend(elnd_to_mshpt) P2_NODES_PER_EL, n_msh_el
!f2py depend(v_el_material) n_msh_pts
!f2py depend(v_nd_xy) n_msh_pts
!f2py depend(soln_ac_u) P2_NODES_PER_EL, n_modes, n_msh_el
!f2py depend(Omega_AC) n_modes
!f2py depend(rho) n_elt_mats


   errco = 0
   call nberr%reset()

   call frontend%init_from_py(n_msh_el, n_msh_pts, elnd_to_mshpt, v_nd_xy, nberr)
   RET_ON_NBERR_UNFOLD(nberr)


   v_energy = D_ZERO


   do i_el=1,n_msh_el
      typ_e = v_el_material(i_el)
      rho_el = rho(typ_e)

      call frontend%nodes_at_el(i_el, nds_xy)

      call basfuncs%set_affine_for_elt(nds_xy, nberr)
      RET_ON_NBERR(nberr)

      call basfuncs%get_triint_p2_p2(m_int_p2p2)

      ! Having calculated v_energy of basis functions on element
      ! now multiply by specific field values for modes of interest.
      do md_i=1,n_modes

         do bf_i=1,P2_NODES_PER_EL
            do bf_j=1,P2_NODES_PER_EL
               do xyz_i=1,3

                  ac_Ustar = conjg(soln_ac_u(xyz_i, bf_i, md_i, i_el))
                  ac_U = soln_ac_u(xyz_i, bf_j, md_i, i_el)

                  v_energy(md_i) = v_energy(md_i) &
                     + rho_el * ac_Ustar * ac_U * m_int_p2p2(bf_i, bf_j)

               enddo
            enddo
         enddo

      enddo
   enddo

   !v_energy = 2.0 * Omega_AC**2 * v_energy
   v_energy_r = real(2.0 * Omega_AC**2 * v_energy)

end subroutine ac_mode_energy_analytic
