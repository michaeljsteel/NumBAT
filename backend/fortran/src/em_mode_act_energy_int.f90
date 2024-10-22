#include "numbat_decl.h"

! Calculate the energy of all EM modes using numerical quadrature.
!  This is the energy of both the electric and magnetic fields.
!
!  E_em = \int eps/2 |E|^2 dV + \int mu/2 |H|^2 dV
!       = 2 eps_0 \int dV eps_r |E|^2

subroutine em_mode_act_energy_int (n_modes, n_msh_el, n_msh_pts, &
   elnd_to_mesh, v_nd_xy, n_elt_mats, el_material, &
   v_refindex, soln_em_e, m_energy, errco, emsg)

   use numbatmod
   use class_TriangleIntegrators

   integer(8) n_modes, n_msh_el, n_msh_pts
   integer(8) elnd_to_mesh(P2_NODES_PER_EL,n_msh_el)
   double precision v_nd_xy(2,n_msh_pts)
   integer(8) n_elt_mats
   integer(8) el_material(n_msh_el)
   complex(8) v_refindex(n_elt_mats)
   complex(8) soln_em_e(3,P2_NODES_PER_EL+7,n_modes,n_msh_el)
   complex(8), dimension(n_modes), intent(out) :: m_energy
   integer(8), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg

   ! Locals

   type(NBError) nberr

   complex(8) t_eps
   complex(8) bas_ovrlap(P2_NODES_PER_EL)
   integer(8)  typ_e
   integer(8) i_el, ival
   integer(8) nd_i, i_eq, iq
   logical  is_curved, do_P3
   integer(8) n_curved
   double precision nds_xy(2,P2_NODES_PER_EL)

   complex(8) em_E, em_Estar

   type(QuadIntegrator) quadint
   type(PyFrontEnd) frontend


   double precision t_quadwt


!f2py depend(elnd_to_mesh) P2_NODES_PER_EL, n_msh_el
!f2py depend(v_nd_xy) n_msh_pts
!f2py depend(soln_em_e) P2_NODES_PER_EL, n_modes, n_msh_el
!f2py depend(v_refindex) n_elt_mats
!f2py depend(el_material) n_msh_el


   do_P3 = .false.
   m_energy = D_ZERO
   n_curved = 0
   call nberr%reset()

   call quadint%setup_reference_quadratures()

   call frontend%init_from_py(n_msh_el, n_msh_pts, elnd_to_mesh, v_nd_xy, nberr)
   RET_ON_NBERR_UNFOLD(nberr)

   do i_el=1,n_msh_el
      typ_e = el_material(i_el)
      t_eps = v_refindex(typ_e)**2  ! Calculate permittivity

      call frontend%nodes_at_el(i_el, nds_xy)

      is_curved = frontend%elt_is_curved()
      if (is_curved) n_curved = n_curved + 1

      bas_ovrlap = 0.0d0

      do iq=1,quadint%n_quad

         call quadint%build_transforms_at(iq, nds_xy, is_curved, do_P3, errco, emsg)
         RETONERROR(errco)

         t_quadwt = quadint%get_current_quadweight()

         ! Calculate m_energy of basis functions at quadrature point,
         ! which is a superposition of P2 polynomials for each function (field).
         ! TODO: why are there no cross terms here?
         do nd_i=1,P2_NODES_PER_EL
            bas_ovrlap(nd_i) = bas_ovrlap(nd_i) + &
               t_quadwt * quadint%phi_P2_ref(nd_i) * quadint%phi_P2_ref(nd_i)
         enddo
      enddo

      ! Having calculated m_energy of basis functions on element
      ! now multiply by specific field values for each mode.
      do ival=1,n_modes
         do nd_i=1,P2_NODES_PER_EL
            do i_eq=1,3
               em_Estar = conjg(soln_em_e(i_eq,nd_i,ival,i_el))
               em_E = soln_em_e(i_eq,nd_i,ival,i_el)
               m_energy(ival) = m_energy(ival) + t_eps * em_Estar * em_E * bas_ovrlap(nd_i)
            enddo
         enddo
      enddo
   enddo

   m_energy = 2.0 * SI_EPS_0 * m_energy

end subroutine em_mode_act_energy_int
