#include "numbat_decl.h"


subroutine set_boundary_conditions(bdy_cdn, n_msh_pts, n_msh_el,  d_nodes_per_el, n_ddl, &
   mesh_raw,  entities, neq, m_eqs, debug, &
   iperiod_N, iperiod_N_E_F, inperiod_N, inperiod_N_E_F )

   use numbatmod
   use class_MeshRaw

   integer(8) :: bdy_cdn, neq, n_msh_pts, n_msh_el, n_ddl, d_nodes_per_el
   integer(8) :: debug

   type(MeshRaw) :: mesh_raw
   type(MeshEntities) :: entities


   integer(8) m_eqs(3, n_ddl)

   integer(8) :: iperiod_N(n_msh_pts), iperiod_N_E_F(n_ddl)
   integer(8) :: inperiod_N(n_msh_pts), inperiod_N_E_F(n_ddl)

   double precision, dimension(2,2) :: lat_vecs

   if ( bdy_cdn .eq. BCS_DIRICHLET .or.  bdy_cdn .eq. BCS_NEUMANN) then

      call bound_cond ( bdy_cdn, n_ddl, entities%v_phys_i, neq, m_eqs)

   elseif( bdy_cdn .eq. BCS_PERIODIC) then  !  Periodic  conditions (never in NumBAT)

      call lattice_vec (n_msh_pts, mesh_raw%v_nd_xy, lat_vecs, debug)

      call periodic_node(n_msh_el, n_msh_pts, d_nodes_per_el, mesh_raw%node_phys_i, mesh_raw%v_nd_xy, iperiod_N, &
         inperiod_N, mesh_raw%elnd_to_mesh, lat_vecs)

      if (debug .eq. 1) then
         write(*,*) "set_boundary_conditions: ###### periodic_N_E_F"
      endif

      call periodic_N_E_F (n_ddl, entities%v_phys_i, entities%v_xy, iperiod_N_E_F, &
         inperiod_N_E_F, lat_vecs)

      call periodic_cond ( bdy_cdn, n_ddl, neq, entities%v_phys_i, &
         iperiod_N_E_F, m_eqs, debug)

   endif

end subroutine set_boundary_conditions

!end module nbinterfaces
