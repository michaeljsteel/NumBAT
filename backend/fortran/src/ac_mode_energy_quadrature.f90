#include "numbat_decl.h"

! Calculate the elastic energy v_energy integral of an AC mode with itself using
! numerical quadrature.

! E = 2 \Omega^2 \int dxdy \rho |u|^2

subroutine AC_mode_energy_quadrature_impl (n_modes, n_msh_elts, n_msh_pts, v_mshpt_xy, &
   m_elnd_to_mshpt, n_elt_mats, v_elt_material,  &
   rho, Omega_AC, soln_ac_u, v_energy_r, nberr)

   use numbatmod
   use class_TriangleIntegrators
   use class_BasisFunctions

   integer(8) n_modes, md_i
   integer(8) n_msh_elts, n_msh_pts, n_elt_mats
   integer(8) v_elt_material(n_msh_elts)
   integer(8) m_elnd_to_mshpt(P2_NODES_PER_EL,n_msh_elts)
   double precision v_mshpt_xy(2,n_msh_pts)
   complex(8) soln_ac_u(3,P2_NODES_PER_EL,n_modes,n_msh_elts)
   complex(8) Omega_AC(n_modes)
   double precision, dimension(n_modes), intent(out) :: v_energy_r

   complex(8) rho(n_elt_mats)
   type(NBError) nberr

   ! Locals

   double precision el_nds_xy(2,P2_NODES_PER_EL), t_xy(2)

   complex(8) bas_ovrlp(P2_NODES_PER_EL,P2_NODES_PER_EL)
   complex(8), dimension(n_modes) :: v_energy
   complex(8) U, Ustar
   integer(8) typ_e, iq
   integer(8) i_el, k_eq
   integer(8) bf_i, bf_j

   logical  is_curved
   integer(8) n_curved

   type(QuadIntegrator) quadint
   type(PyFrontEnd) frontend
   type(BasisFunctions) basfuncs

   double precision t_quadwt, t_qwt

   call quadint%setup_reference_quadratures()

   call frontend%init_from_py(n_msh_elts, n_msh_pts, m_elnd_to_mshpt, v_mshpt_xy, nberr)
   RET_ON_NBERR(nberr)

   v_energy = C_ZERO

   do i_el=1,n_msh_elts
      typ_e = v_elt_material(i_el)

      call frontend%nodes_at_el(i_el, el_nds_xy)

      call basfuncs%set_affine_for_elt(el_nds_xy, nberr)
      RET_ON_NBERR(nberr)

      is_curved = frontend%elt_is_curved()
      if (is_curved) n_curved = n_curved + 1


      bas_ovrlp = D_ZERO

! For each quadrature point evaluate v_energy of Lagrange polynomials
! or derivative of Lagrange polynomials
      do iq=1,quadint%n_quad

         call quadint%get_quad_point(iq, t_xy, t_quadwt)
         RET_ON_NBERR(nberr)

         call basfuncs%evaluate_at_position(i_el, t_xy, is_curved, el_nds_xy, nberr)
         RET_ON_NBERR(nberr)

         t_qwt = t_quadwt * abs(basfuncs%det)


         ! Calculate v_energy of basis functions at quadrature point,
         ! which is a superposition of P2 polynomials for each function (field).
         do bf_i=1,P2_NODES_PER_EL
            do bf_j=1,P2_NODES_PER_EL
               bas_ovrlp(bf_i,bf_j) = bas_ovrlp(bf_i,bf_j) + t_qwt * basfuncs%phi_P2_ref(bf_i) * basfuncs%phi_P2_ref(bf_j)
            enddo
         enddo
      enddo

      ! Having calculated v_energy of basis functions on element
      ! now multiply by specific field values for modes of interest.
      do md_i=1,n_modes
         do bf_i=1,P2_NODES_PER_EL
            do bf_j=1,P2_NODES_PER_EL
               do k_eq=1,3
                  Ustar = conjg(soln_ac_u(k_eq,bf_i,md_i,i_el))
                  U = soln_ac_u(k_eq,bf_j,md_i,i_el)
                  v_energy(md_i) = v_energy(md_i) + rho(typ_e) * Ustar * U * bas_ovrlp(bf_i,bf_j)
               enddo
            enddo
         enddo
      enddo

   enddo

   v_energy_r = real(2.0 * Omega_AC**2 * v_energy)

end subroutine AC_mode_energy_quadrature_impl
