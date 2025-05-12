#include "numbat_decl.h"


subroutine set_boundary_conditions(bdy_cdn,  &
   mesh,  entities, n_dof, m_eqs, debug, &
   iperiod_N, iperiod_N_E_F, inperiod_N, inperiod_N_E_F, errco, emsg )

   use numbatmod
   use class_Mesh

   integer(8) :: bdy_cdn, n_dof
   integer(8) :: debug

   type(MeshEM) :: mesh
   type(MeshEntities) :: entities


   integer(8) m_eqs(3, entities%n_entities)

   integer(8) :: iperiod_N(mesh%n_msh_pts), iperiod_N_E_F(entities%n_entities)
   integer(8) :: inperiod_N(mesh%n_msh_pts), inperiod_N_E_F(entities%n_entities)
   integer(8),  intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg


   !locals
   double precision, dimension(2,2) :: lat_vecs

   if ( bdy_cdn .eq. BCS_DIRICHLET .or.  bdy_cdn .eq. BCS_NEUMANN) then

      call bound_cond_em (bdy_cdn, entities%n_entities, entities%v_ety_props, n_dof, m_eqs, errco, emsg)

   elseif( bdy_cdn .eq. BCS_PERIODIC) then  !  Periodic  conditions (never in NumBAT)

      call periodic_lattice_vec (mesh%n_msh_pts, mesh%v_mshpt_xy, lat_vecs, debug)

      call periodic_node(mesh%n_msh_elts, mesh%n_msh_pts, &
      P2_NODES_PER_EL, mesh%node_phys_i, mesh%v_mshpt_xy, iperiod_N, &
         inperiod_N, mesh%m_elnd_to_mshpt, lat_vecs)

      call periodic_N_E_F (entities%n_entities, entities%v_ety_props, entities%v_xy, iperiod_N_E_F, &
         inperiod_N_E_F, lat_vecs)

      call periodic_cond ( bdy_cdn, entities%n_entities, n_dof, entities%v_ety_props, &
         iperiod_N_E_F, m_eqs, debug)

   endif

end subroutine
