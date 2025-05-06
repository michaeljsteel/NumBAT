#include "numbat_decl.h"

! Calculates unnormalised moving boundary coupling Q_MB
!
! \nhat is the surface normal vector from material a into material b
! Q_MB = \int_C  dr
!         (\vec u^* \dot \hatn)
!           \cross
!  [   eps0  (eps_a - eps_b)     (\hatn \cross \vecEs)^* \dot (\hatn \cross \vecEp)
!   - 1/eps0 (1/eps_a - 1/epsb)  (\hatn \dot \vecDs)^*   \dot (\hatn \dot \vecDp)  ]
!


subroutine moving_boundary (nval_em_p, nval_em_s, nval_ac, ival_p, ival_s, ival_ac, &
   n_msh_elts, n_msh_pts, elnd_to_mshpt, v_mshpt_xy, &
   n_elt_mats, v_elt_material, typ_select_in, typ_select_out, &
   soln_em_p, soln_em_s, soln_ac_u, v_eps_rel, Q_MB, errco, emsg)

   use numbatmod
   integer(8) n_msh_elts, n_msh_pts,n_elt_mats
   integer(8) nval_em_p, nval_em_s, nval_ac, ival_p, ival_s, ival_ac

   integer(8) v_elt_material(n_msh_elts)
   integer(8) elnd_to_mshpt(P2_NODES_PER_EL, n_msh_elts)
   double precision v_mshpt_xy(2, n_msh_pts)

   integer(8) typ_select_in, typ_select_out

   complex(8) soln_em_p(3, P2_NODES_PER_EL, nval_em_p, n_msh_elts)
   complex(8) soln_em_s(3, P2_NODES_PER_EL, nval_em_s, n_msh_elts)
   complex(8) soln_ac_u(3, P2_NODES_PER_EL, nval_ac, n_msh_elts)

   complex(8) v_eps_rel(n_elt_mats)
   complex(8), intent(out) :: Q_MB(nval_em_p, nval_em_s, nval_ac)
   integer(8), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) ::  emsg

   !     Local variables
   integer(8) nb_visited(n_msh_pts)
   integer(8) ls_edge_endpoint(2,n_msh_pts)
   integer(8) edge_direction(n_msh_pts)
   integer(8) i_el, inod, typ_e
   integer(8) inod_1, inod_2, inod_3, ls_inod(3)
   integer(8) j, j_1, j_2, j_3
   integer(8) nb_edges, nb_interface_edges
   integer(8) edge_endpoints(2,3), opposite_node(3)
   double precision xy_1(2), xy_2(2), xy_3(2), ls_xy(2,3)
   double precision edge_vec(2), edge_perp(2), vec_0(2)
   double precision edge_length, r_tmp
   !complex(8) ls_n_dot(3)
   !complex(8)  ls_n_cross(3,3)
   !complex(8) vec(3,3)
   !integer(8) ivals_ac, ivals_s, ivals_p

   complex(8) evec_p(3), evec_sc(3), uvec_ac(3)

   complex(8) n_dot_d(2)
   complex(8) eps_a, eps_b, tmp1, tmp2
   double precision p2_p2_p2_1d(3,3,3)

   integer(8) v_ival_p(nval_em_p), v_ival_s(nval_em_s), v_ival_ac(nval_ac)
   integer(8) ivs, ivp, ivac, t_ival_s, t_ival_p, t_ival_ac

   complex(8) n_cross_ev_p(3), n_cross_ev_sc(3)
   complex(8) n_dot_ev_p, n_dot_ev_sc, n_dot_uv_ac


   !f2py intent(in) nval_em_p, nval_em_s, nval_ac
   !f2py intent(in) ival_p, ival_s, ival_ac, n_elt_mats
   !f2py intent(in) n_msh_elts, n_msh_pts, P2_NODES_PER_EL, elnd_to_mshpt, debug
   !f2py intent(in) v_elt_material, x, soln_em_p, soln_em_s, soln_ac_u
   !f2py intent(in) typ_select_in, typ_select_out, v_eps_rel, debug

   !f2py depend(elnd_to_mshpt) P2_NODES_PER_EL, n_msh_elts
   !f2py depend(v_elt_material) n_msh_pts
   !f2py depend(x) n_msh_pts
   !f2py depend(soln_em_p) P2_NODES_PER_EL, nval_em_p, n_msh_elts
   !f2py depend(soln_em_s) P2_NODES_PER_EL, nval_em_s, n_msh_elts
   !f2py depend(soln_ac_u) P2_NODES_PER_EL, nval_ac, n_msh_elts
   !f2py depend(v_eps_rel) n_elt_mats
   !
   !fo2py intent(out) Q_MB


   ! typ_select_in: Only the elements i_el with v_elt_material(i_el)=typ_select_in will be analysed
   ! When nb_visited(j) is not zero: nb_visited(j) indicates the number of element the edge j belongs



   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Rules for finding edge ends and opposite nodes on the standard
   ! six node anticlockwise triangle
   !             3
   !           6    5
   !        1    4     2

   edge_endpoints(1,1) = 1   ! Edge 1 has endpoints at nodes 1 and 2
   edge_endpoints(2,1) = 2
   opposite_node(1) = 3      !        and is opposite node 3

   edge_endpoints(1,2) = 2   ! Edge 2 has endpoints at nodes 2 and 3
   edge_endpoints(2,2) = 3
   opposite_node(2) = 1

   edge_endpoints(1,3) = 3   ! Edge 3 has endpoints at nodes 3 and 1
   edge_endpoints(2,3) = 1
   opposite_node(3) = 2


   nb_visited = 0
   ls_edge_endpoint = 0
   ls_edge_endpoint =0
   edge_direction = 0

   Q_MB = D_ZERO



   ! build arrays holding which modes to be calculated
   call fill_ival_arrays(v_ival_p, nval_em_p, ival_p)
   call fill_ival_arrays(v_ival_s, nval_em_s, ival_s)
   call fill_ival_arrays(v_ival_ac, nval_ac, ival_ac)

   do i_el=1,n_msh_elts
      typ_e = v_elt_material(i_el)
      if(typ_e == typ_select_in) then
         !   Scan the edges
         do inod=4,6
            j = elnd_to_mshpt(inod,i_el)
            !   Will indicate the number of
            nb_visited(j) = nb_visited(j) + 1
         enddo
      endif
   enddo

   nb_edges = 0
   nb_interface_edges = 0
   do inod=1,n_msh_pts
      if (nb_visited(inod) >= 1) then
         nb_edges = nb_edges + 1
      endif
      if (nb_visited(inod) == 1) then
         nb_interface_edges = nb_interface_edges + 1
      endif
   enddo


   ! Outward pointing normal vector to the interface edges
   do i_el=1,n_msh_elts
      typ_e = v_elt_material(i_el)
      if (typ_e .ne. typ_select_in) cycle  ! not the inside material

      !   Scan the edges
      do inod=4,6
         j = elnd_to_mshpt(inod,i_el)
         if (nb_visited(j) .ne. 1) cycle ! not an active edge

         inod_1 = edge_endpoints(1,inod-3)
         inod_2 = edge_endpoints(2,inod-3)
         ls_edge_endpoint(1,j) = elnd_to_mshpt(inod_1,i_el)
         ls_edge_endpoint(2,j) = elnd_to_mshpt(inod_2,i_el)

         xy_1 = v_mshpt_xy(:, elnd_to_mshpt(inod_1,i_el))
         xy_2 = v_mshpt_xy(:, elnd_to_mshpt(inod_2,i_el))

         ! edge_vec: vector parallel to the edge
         edge_vec = xy_2 - xy_1

         ! Normalisation of edge_vec
         r_tmp = sqrt(edge_vec(1)**2+edge_vec(2)**2)
         edge_vec = edge_vec / r_tmp

         ! edge_vec: vector perpendicular to the edge (rotation of edge_vec by -pi/2)
         edge_perp(1) = edge_vec(2)
         edge_perp(2) = -edge_vec(1)

         ! Node opposite to the edge inod
         inod_3 = opposite_node(inod-3)
         xy_3(:) = v_mshpt_xy(:,elnd_to_mshpt(inod_3,i_el))
         vec_0 = xy_3 - xy_1

         ! Scalar product of edge_perp and vec_0:
         r_tmp = edge_perp(1)*vec_0(1)+edge_perp(2)*vec_0(2)
         ! if r_tmp < 0: then edge_perp is oriented in the outward direction

         if( r_tmp < 0) then
            edge_direction(j) = 1
         elseif( r_tmp > 0) then
            edge_direction(j) = -1
         else
            errco = NBERR_BAD_MB_EDGES
            write(emsg,*) "edge_orientation: illegal, edge_perp is perpendicular to vec_0"
            return
         endif
      enddo
   enddo




   ! Numerical integration
   do i_el=1,n_msh_elts

      typ_e = v_elt_material(i_el)
      if(typ_e .ne. typ_select_in) then
         cycle
      endif

      eps_a = v_eps_rel(typ_e)
      if (typ_select_out .eq. -1) then
         eps_b = 1.0d0
      else
         eps_b = v_eps_rel(typ_select_out)
      endif

      !   Scan the edges
      do inod=4,6
         j = elnd_to_mshpt(inod,i_el)

         if (ls_edge_endpoint(1,j) .eq. 0) then ! Not an edge
            cycle
         endif

         inod_1 = ls_edge_endpoint(1,j)
         inod_2 = ls_edge_endpoint(2,j)
         xy_1(:) = v_mshpt_xy(:, inod_1)
         xy_2(:) = v_mshpt_xy(:, inod_2)
         xy_3(:) = v_mshpt_xy(:, j)

         ! List of the nodes coordinates
         ls_xy(:,1) = xy_1  ! x-coord. of node 1
         ls_xy(:,2) = xy_2  !             node 2
         ls_xy(:,3) = xy_3  !             node 3

         edge_vec(:) = ls_xy(:,2) - ls_xy(:,1)

         ! Normalisation of edge_vec
         r_tmp = sqrt(edge_vec(1)**2+edge_vec(2)**2)
         edge_vec = -1*edge_direction(j)*edge_vec / r_tmp

         ! edge_vec: vector perpendicular to the edge (rotation of edge_vec by -pi/2)
         edge_perp(1) = -1*edge_vec(2)
         edge_perp(2) = edge_vec(1)


         r_tmp = (ls_xy(1,2) - ls_xy(1,1))**2 + (ls_xy(2,2) - ls_xy(2,1))**2
         edge_length = sqrt(r_tmp)

         ! Some integration factor
         call find_overlaps_p2_p2_p2_1d (p2_p2_p2_1d, edge_length)

         ! Identification number of the two end-points and mid-edge point
         ls_inod(1) = edge_endpoints(1,inod-3)
         ls_inod(2) = edge_endpoints(2,inod-3)
         ls_inod(3) = inod



         !             !
         ! ! \nhat is the surface normal vector from material a into material b
         ! ! Q_MB = \int_C  dr
         ! !         (\vec u^* \dot \hatn)
         ! !           \cross
         ! !  [   eps0  (eps_a - eps_b)     (\hatn \cross \vecEs)^* \dot (\hatn \cross \vecEp)
         ! !   - 1/eps0 (1/eps_a - 1/epsb)  (\hatn \dot \vecDs)^*   \dot (\hatn \dot \vecdp)  ]

         do ivs = 1,nval_em_s
            t_ival_s = v_ival_s(ivs)
            if (t_ival_s .eq. 0) exit


            ! Nodes of the edge
            do j_1=1,3
               ! (x,y,z)-components of the electric field
               evec_sc = conjg(soln_em_s(:, ls_inod(j_1), t_ival_s, i_el))
               ! ls_n_dot(1): Normal component of vec(:,1)
               ! ls_n_dot(1) = evec_sc(1) * edge_perp(1) + evec_sc(2) * edge_perp(2)
               ! ls_n_dot(1) = v2_dot_v3(edge_perp, evec_sc)
               ! ls_n_cross(1,1) = evec_sc(3) * edge_perp(2)
               ! ls_n_cross(2,1) = -1*evec_sc(3) * edge_perp(1)
               ! ls_n_cross(3,1) = evec_sc(2) * edge_perp(1) - evec_sc(1) * edge_perp(2)

               n_dot_ev_sc = v2_dot_v3(edge_perp, evec_sc)
               call v2_cross_v3(edge_perp, evec_sc, n_cross_ev_sc)


               do ivp = 1,nval_em_p
                  t_ival_p = v_ival_p(ivp)
                  if (t_ival_p .eq. 0) exit

                  do j_2=1,3
                     ! (x,y,z)-components of the electric field
                     evec_p = soln_em_p(:, ls_inod(j_2), t_ival_p, i_el)
                     ! ls_n_dot(2): Normal component of vec(:,2)
                     ! ls_n_dot(2) = v2_dot_v3(edge_perp, evec_p)
                     ! ls_n_cross(1,2) = evec_p(3) * edge_perp(2)
                     ! ls_n_cross(2,2) = -1*evec_p(3) * edge_perp(1)
                     ! ls_n_cross(3,2) = evec_p(2) * edge_perp(1) - evec_p(1) * edge_perp(2)

                     n_dot_ev_p = v2_dot_v3(edge_perp, evec_p)
                     call v2_cross_v3(edge_perp, evec_p, n_cross_ev_p)


                     do ivac = 1,nval_ac
                        t_ival_ac = v_ival_ac(ivac)
                        if (t_ival_ac .eq. 0) exit

                        do j_3=1,3
                           ! (x,y,z)-components of the acoustic field
                           uvec_ac = soln_ac_u(:,ls_inod(j_3),t_ival_ac,i_el)

                           ! ls_n_dot(3): scalar product of vec(:,3) and normal vector edge_perp
                           !ls_n_dot(3) = uvec_ac(1) * edge_perp(1) + uvec_ac(2) * edge_perp(2)
                           n_dot_uv_ac = v2_dot_v3(edge_perp, uvec_ac)


                           !tmp1 = (eps_a - eps_b)*SI_EPS_0
                           !tmp1 = tmp1*((ls_n_cross(1,1))*ls_n_cross(1,2) + (ls_n_cross(2,1))*ls_n_cross(2,2) + (ls_n_cross(3,1))*ls_n_cross(3,2))
                           tmp1 = (eps_a - eps_b)*SI_EPS_0 * cv3_dot_cv3(n_cross_ev_sc, n_cross_ev_p)

                           !n_dot_d(1) = SI_EPS_0*eps_a * ls_n_dot(1)
                           !n_dot_d(2) = SI_EPS_0*eps_a * ls_n_dot(2)
                           n_dot_d(1) = SI_EPS_0*eps_a * n_dot_ev_sc
                           n_dot_d(2) = SI_EPS_0*eps_a * n_dot_ev_p

                           tmp2 = (1.0d0/eps_b-1.0d0/eps_a)/SI_EPS_0 * n_dot_d(1) * n_dot_d(2)
                           r_tmp = p2_p2_p2_1d(j_1, j_2, j_3)


                           Q_MB(t_ival_p, t_ival_s, t_ival_ac) = Q_MB(t_ival_p, t_ival_s, t_ival_ac) &
                              + r_tmp*conjg(n_dot_uv_ac)*(tmp1 + tmp2)

                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo

end subroutine moving_boundary



subroutine fill_ival_arrays(v_ival, n_ivals, ival)

   integer(8) n_ivals, ival
   integer(8) v_ival(n_ivals)
   integer(8) i

   v_ival = 0
   if (ival .gt. 0) then
      v_ival(1) = ival
   else
      do i=1,n_ivals
         v_ival(i) = i
      enddo
   endif

end subroutine
