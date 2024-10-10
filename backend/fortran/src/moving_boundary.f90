! Calculates moving boundary coupling Q_MB
!
! \nhat is the surface normal vector from material a into material b
! Q_MB = \int_C  dr
!         (\vec u^* \dot \hatn)
!           \cross
!  [   eps0  (eps_a - eps_b)     (\hatn \cross \vecEs)^* \dot (\hatn \cross \vecEp)
!   - 1/eps0 (1/eps_a - 1/epsb)  (\hatn \dot \vecDs)^*   \dot (\hatn \dot \vecdp)  ]
!


subroutine moving_boundary (nval_EM_p, nval_EM_S, nval_AC, ival_p, &
   ival_s, ival_ac, n_msh_el, n_msh_pts, nodes_per_el, table_nod, type_el, x, &
   nb_typ_el, typ_select_in, typ_select_out, &
   soln_EM_p, soln_EM_S, soln_AC, eps_lst, debug, Q_MB)

   use numbatmod
   integer(8) n_msh_el, n_msh_pts, nodes_per_el, nb_typ_el
   integer(8) type_el(n_msh_el)
   integer(8) table_nod(6,n_msh_el)
   double precision x(2,n_msh_pts)
   integer(8) nval_EM_p, nval_EM_S, nval_AC
   integer(8) ival_p, ival_s, ival_ac

   integer(8) typ_select_in, typ_select_out
   complex(8) soln_EM_p(3,nodes_per_el,nval_EM_p,n_msh_el)
   complex(8) soln_EM_S(3,nodes_per_el,nval_EM_S,n_msh_el)
   complex(8) soln_AC(3,nodes_per_el,nval_AC,n_msh_el)
   complex(8) eps_lst(nb_typ_el)
   complex(8) Q_MB(nval_EM_S, nval_EM_p, nval_AC)

   !     Local variables
   integer(8) debug
   integer(8) nb_visited(n_msh_pts)
   integer(8) ls_edge_endpoint(2,n_msh_pts)
   integer(8) edge_direction(n_msh_pts)
   integer(8) iel, inod, typ_e
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

   integer(8) v_ival_p(nval_EM_p), v_ival_s(nval_EM_s), v_ival_ac(nval_AC)
   integer(8) ivs, ivp, ivac, t_ival_s, t_ival_p, t_ival_ac
   complex(8) n_cross_ev_p(3), n_cross_ev_sc(3)
   complex(8) n_dot_ev_p, n_dot_ev_sc, n_dot_uv_ac


   !f2py intent(in) nval_EM_p, nval_EM_S, nval_AC
   !f2py intent(in) ival_p, ival_s, ival_ac, nb_typ_el
   !f2py intent(in) n_msh_el, n_msh_pts, nodes_per_el, table_nod, debug
   !f2py intent(in) type_el, x, soln_EM_p, soln_EM_S, soln_AC
   !f2py intent(in) typ_select_in, typ_select_out, eps_lst, debug

   !f2py depend(table_nod) nodes_per_el, n_msh_el
   !f2py depend(type_el) n_msh_pts
   !f2py depend(x) n_msh_pts
   !f2py depend(soln_EM_p) nodes_per_el, nval_EM_p, n_msh_el
   !f2py depend(soln_EM_S) nodes_per_el, nval_EM_S, n_msh_el
   !f2py depend(soln_AC) nodes_per_el, nval_AC, n_msh_el
   !f2py depend(eps_lst) nb_typ_el
   !
   !f2py intent(out) Q_MB


   ! typ_select_in: Only the elements iel with type_el(iel)=typ_select_in will be analysed
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
   call fill_ival_arrays(v_ival_p, nval_EM_p, ival_p)
   call fill_ival_arrays(v_ival_s, nval_EM_s, ival_s)
   call fill_ival_arrays(v_ival_ac, nval_ac, ival_ac)


   do iel=1,n_msh_el
      typ_e = type_el(iel)
      if(typ_e == typ_select_in) then
         !   Scan the edges
         do inod=4,6
            j = table_nod(inod,iel)
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

   if (debug .eq. 1) then
      write(*,*)
      write(*,*) "edge_orientation: n_msh_pts, n_msh_el = ", n_msh_pts, n_msh_el
      write(*,*) "edge_orientation: nb_edges = ", nb_edges
      write(*,*) "nb_interface_edges = ", nb_interface_edges
   endif



   ! Outward pointing normal vector to the interface edges
   do iel=1,n_msh_el
      typ_e = type_el(iel)
      if(typ_e .ne. typ_select_in) then
         cycle
      endif

      !   Scan the edges
      do inod=4,6
         j = table_nod(inod,iel)
         if (nb_visited(j) .ne. 1) then ! not an active edge
            cycle
         endif

         inod_1 = edge_endpoints(1,inod-3)
         inod_2 = edge_endpoints(2,inod-3)
         ls_edge_endpoint(1,j) = table_nod(inod_1,iel)
         ls_edge_endpoint(2,j) = table_nod(inod_2,iel)
         xy_1(1) = x(1,table_nod(inod_1,iel))
         xy_1(2) = x(2,table_nod(inod_1,iel))
         xy_2(1) = x(1,table_nod(inod_2,iel))
         xy_2(2) = x(2,table_nod(inod_2,iel))
         ! edge_vec: vector parallel to the edge
         edge_vec(1) = xy_2(1) - xy_1(1)
         edge_vec(2) = xy_2(2) - xy_1(2)
         ! Normalisation of edge_vec
         r_tmp = sqrt(edge_vec(1)**2+edge_vec(2)**2)
         edge_vec(1) = edge_vec(1) / r_tmp
         edge_vec(2) = edge_vec(2) / r_tmp
         ! edge_vec: vector perpendicular to the edge (rotation of edge_vec by -pi/2)
         edge_perp(1) = edge_vec(2)
         edge_perp(2) = -edge_vec(1)
         ! Node opposite to the edge inod
         inod_3 = opposite_node(inod-3)
         xy_3(1) = x(1,table_nod(inod_3,iel))
         xy_3(2) = x(2,table_nod(inod_3,iel))
         vec_0(1) = xy_3(1) - xy_1(1)
         vec_0(2) = xy_3(2) - xy_1(2)
         ! Scalar product of edge_perp and vec_0:
         r_tmp = edge_perp(1)*vec_0(1)+edge_perp(2)*vec_0(2)
         ! if r_tmp < 0: then edge_perp is oriented in the outward direction

         if( r_tmp < 0) then
            edge_direction(j) = 1
         elseif( r_tmp > 0) then
            edge_direction(j) = -1
         else
            write(*,*) "edge_orientation: illegal:"
            write(*,*) "edge_perp is perpendicular to vec_0"
            write(*,*) "edge_orientation: Aborting..."
            stop
         endif
      enddo
   enddo




   ! Numerical integration
   do iel=1,n_msh_el

      typ_e = type_el(iel)
      if(typ_e .ne. typ_select_in) then
         cycle
      endif

      eps_a = eps_lst(typ_e)
      if (typ_select_out .eq. -1) then
         eps_b = 1.0d0
      else
         eps_b = eps_lst(typ_select_out)
      endif

      !   Scan the edges
       do inod=4,6
          j = table_nod(inod,iel)

          if (ls_edge_endpoint(1,j) .eq. 0) then ! Not an edge
             cycle
          endif

         inod_1 = ls_edge_endpoint(1,j)
         inod_2 = ls_edge_endpoint(2,j)
         xy_1(1) = x(1,inod_1)
         xy_1(2) = x(2,inod_1)
         xy_2(1) = x(1,inod_2)
         xy_2(2) = x(2,inod_2)
         xy_3(1) = x(1,j)
         xy_3(2) = x(2,j)

         ! List of the nodes coordinates
         ls_xy(1,1) = xy_1(1)  ! x-coord. of node 1
         ls_xy(2,1) = xy_1(2)  ! y-coord. of node 1
         ls_xy(1,2) = xy_2(1)  !             node 2
         ls_xy(2,2) = xy_2(2)
         ls_xy(1,3) = xy_3(1)  !             node 3
         ls_xy(2,3) = xy_3(2)

         edge_vec(1) = ls_xy(1,2) - ls_xy(1,1)
         edge_vec(2) = ls_xy(2,2) - ls_xy(2,1)

         ! Normalisation of edge_vec
         r_tmp = sqrt(edge_vec(1)**2+edge_vec(2)**2)
         edge_vec(1) = -1*edge_direction(j)*edge_vec(1) / r_tmp
         edge_vec(2) = -1*edge_direction(j)*edge_vec(2) / r_tmp

         ! edge_vec: vector perpendicular to the edge (rotation of edge_vec by -pi/2)
         edge_perp(1) = -1*edge_vec(2)
         edge_perp(2) = edge_vec(1)


         r_tmp = (ls_xy(1,2) - ls_xy(1,1))**2 + (ls_xy(2,2) - ls_xy(2,1))**2
         edge_length = sqrt(r_tmp)

         ! Some integration factor
         call mat_p2_p2_p2_1d (p2_p2_p2_1d, edge_length)

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

          do ivs = 1,nval_EM_s
            t_ival_s = v_ival_s(ivs)
            if (t_ival_s .eq. 0) then
               exit
            endif


            ! Nodes of the edge
            do j_1=1,3
               ! (x,y,z)-components of the electric field
               evec_sc(1) = conjg(soln_EM_S(1,ls_inod(j_1),t_ival_s,iel))
               evec_sc(2) = conjg(soln_EM_S(2,ls_inod(j_1),t_ival_s,iel))
               evec_sc(3) = conjg(soln_EM_S(3,ls_inod(j_1),t_ival_s,iel))
               ! ls_n_dot(1): Normal component of vec(:,1)
               ! ls_n_dot(1) = evec_sc(1) * edge_perp(1) + evec_sc(2) * edge_perp(2)
               ! ls_n_dot(1) = v2_dot_v3(edge_perp, evec_sc)
               ! ls_n_cross(1,1) = evec_sc(3) * edge_perp(2)
               ! ls_n_cross(2,1) = -1*evec_sc(3) * edge_perp(1)
               ! ls_n_cross(3,1) = evec_sc(2) * edge_perp(1) - evec_sc(1) * edge_perp(2)

               n_dot_ev_sc = v2_dot_v3(edge_perp, evec_sc)
               call v2_cross_v3(edge_perp, evec_sc, n_cross_ev_sc)


               do ivp = 1,nval_EM_p
                  t_ival_p = v_ival_p(ivp)
                  if (t_ival_p .eq. 0) then
                     exit
                  endif

                  do j_2=1,3
                     ! (x,y,z)-components of the electric field
                     evec_p(1)=soln_EM_p(1,ls_inod(j_2),t_ival_p,iel)
                     evec_p(2)=soln_EM_p(2,ls_inod(j_2),t_ival_p,iel)
                     evec_p(3)=soln_EM_p(3,ls_inod(j_2),t_ival_p,iel)
                     ! ls_n_dot(2): Normal component of vec(:,2)
                     ! ls_n_dot(2) = v2_dot_v3(edge_perp, evec_p)
                     ! ls_n_cross(1,2) = evec_p(3) * edge_perp(2)
                     ! ls_n_cross(2,2) = -1*evec_p(3) * edge_perp(1)
                     ! ls_n_cross(3,2) = evec_p(2) * edge_perp(1) - evec_p(1) * edge_perp(2)

                     n_dot_ev_p = v2_dot_v3(edge_perp, evec_p)
                     call v2_cross_v3(edge_perp, evec_p, n_cross_ev_p)


                     do ivac = 1,nval_AC
                        t_ival_ac = v_ival_ac(ivac)
                        if (t_ival_ac .eq. 0) then
                           exit
                        endif

                        do j_3=1,3
                           ! (x,y,z)-components of the acoustic field
                           uvec_ac(1) = soln_AC(1,ls_inod(j_3),t_ival_ac,iel)
                           uvec_ac(2) = soln_AC(2,ls_inod(j_3),t_ival_ac,iel)
                           uvec_ac(3) = soln_AC(3,ls_inod(j_3),t_ival_ac,iel)

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

                           tmp2 = (1.0d0/eps_b-1.0d0/eps_a)*(1.0d0/SI_EPS_0)
                           tmp2 = tmp2*(n_dot_d(1))*n_dot_d(2)
                           r_tmp = p2_p2_p2_1d(j_1, j_2, j_3)

!                            !Q_MB(ivals_p,ivals_s,ivals_ac) = Q_MB(ivals_p,ivals_s,ivals_ac)+ r_tmp*conjg(ls_n_dot(3))*(tmp1 + tmp2)
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
