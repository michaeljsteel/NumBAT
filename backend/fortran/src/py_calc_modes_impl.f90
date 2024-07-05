#include "numbat_decl.h"


subroutine set_boundary_conditions(bdy_cdn, n_msh_pts, n_msh_el, mesh_xy, d_nodes_per_el, &
   type_nod, table_nod, n_ddl, neq, & !ip_type_N_E_F, ip_eq, &
   d_dwork, type_N_E_F, m_eqs, debug, &
   iperiod_N, iperiod_N_E_F, inperiod_N, inperiod_N_E_F )

   use numbatmod

   integer(8) :: bdy_cdn, neq, n_msh_pts, n_msh_el, n_ddl, d_nodes_per_el
   !integer(8) :: ip_type_N_E_F, ip_eq, jp_x_n_e_f
   !integer(8) :: ip_period_N, ip_nperiod_N, ip_period_N_E_F, ip_nperiod_N_E_F
   integer(8) :: debug
   double precision mesh_xy(2,n_msh_pts)
   integer(8) type_nod(n_msh_pts)
   integer(8) table_nod(d_nodes_per_el, n_msh_el)
   integer(8) type_N_E_F(2, n_ddl)
   integer(8) m_eqs(3, n_ddl)

   integer(8) :: iperiod_N(n_msh_pts), iperiod_N_E_F(n_ddl)
   integer(8) :: inperiod_N(n_msh_pts), inperiod_N_E_F(n_ddl)



   double precision, dimension(2,n_ddl) :: d_dwork

   double precision, dimension(2,2) :: lat_vecs

   if ( bdy_cdn .eq. BCS_DIRICHLET .or.  bdy_cdn .eq. BCS_NEUMANN) then

      call bound_cond ( bdy_cdn, n_ddl, neq, type_N_E_F, m_eqs)

   elseif( bdy_cdn .eq. BCS_PERIODIC) then  !  Periodic  conditions (never in NumBAT)
      if (debug .eq. 1) then
         write(*,*) "###### periodic_node"
      endif

      call lattice_vec (n_msh_pts, mesh_xy, lat_vecs, debug)

      call periodic_node(n_msh_el, n_msh_pts, d_nodes_per_el, type_nod, mesh_xy, iperiod_N, &
         inperiod_N, table_nod, lat_vecs)

      if (debug .eq. 1) then
         write(*,*) "set_boundary_conditions: ###### periodic_N_E_F"
      endif

      call periodic_N_E_F (n_ddl, type_N_E_F, d_dwork, iperiod_N_E_F, &
         inperiod_N_E_F, lat_vecs)

      call periodic_cond ( bdy_cdn, n_ddl, neq, type_N_E_F, &
         iperiod_N_E_F, m_eqs, debug)

   endif

end subroutine set_boundary_conditions

!end module nbinterfaces
