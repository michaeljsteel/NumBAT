#include "numbat_decl.h"

! Calculate the unnormalised Q_PE integral of two EM modes and an AC mode using
! either analytic expressions (linear elements) or quadrature integration (curvilinear elements)

! Q_PE_unnorm = -eps_0  \int dx dy  eps_r^2  p_ijkl (E_s^*)_i (E_p)_j \partial k (u_l)^*



subroutine photoelastic_int_common (is_curvilinear, nval_em_p, nval_em_s, nval_ac_u, ival_p, &
   ival_s, ival_ac, &
   n_msh_el, n_msh_pts, elnd_to_mshpt, v_nd_xy, &
   n_elt_mats, el_material, p_tensor, beta_ac, soln_em_p, soln_em_s, soln_ac_u,&
   v_eps_rel, Q_PE, errco, emsg)

   use numbatmod
   use alloc
   use class_TriangleIntegrators

   integer(8) is_curvilinear
   integer(8) nval_em_p, nval_em_s, nval_ac_u, ival_p, ival_s, ival_ac
   integer(8) n_msh_el, n_msh_pts, n_elt_mats
   integer(8) el_material(n_msh_el), debug
   integer(8) elnd_to_mshpt(P2_NODES_PER_EL,n_msh_el)
   double precision v_nd_xy(2,n_msh_pts)
   complex(8) soln_em_p(3,P2_NODES_PER_EL,nval_em_p,n_msh_el)
   complex(8) soln_em_s(3,P2_NODES_PER_EL,nval_em_s,n_msh_el)
   complex(8) soln_ac_u(3,P2_NODES_PER_EL,nval_ac_u,n_msh_el)
   complex(8) p_tensor(3,3,3,3,n_elt_mats)

   complex(8) beta_ac

   complex(8), intent(out) :: Q_PE(nval_em_s, nval_em_p, nval_ac_u)
   integer(8), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) ::  emsg

   !---------------------------


   double precision nds_xy(2,P2_NODES_PER_EL)

   complex(8), dimension(:,:,:,:), allocatable :: basis_overlap


   complex(8) E_s_i_star, E_p_j, Ustar_l
   integer(8) j, typ_e
   integer(8) i_el, nd_i, nd_j, nd_l, xyz_i, xyz_j, xyz_k, xyz_l
   integer(8) ind_ip, ind_jp, ind_lp
   integer(8) ui
   complex(8) v_eps_rel(n_elt_mats), eps
   complex(8) zt1

   integer(8) v_ival_p(nval_em_p), v_ival_s(nval_EM_s), v_ival_u(nval_ac_u)
   integer(8) ivs, ivp, ivu, t_ival_s, t_ival_p, t_ival_u

   integer(8) n_curved
   logical is_curved
   type(QuadIntegrator) quadint

   debug = 0
   errco = 0
   emsg = ""

   !fo2py intent(in) nval_em_p, nval_em_s, nval_ac_u
   !fo2py intent(in) ival_p, ival_s, ival_ac, n_elt_mats
   !fo2py intent(in) n_msh_el, n_msh_pts, P2_NODES_PER_EL, elnd_to_mshpt, p_tensor, beta_ac, debug
   !fo2py intent(in) el_material, x, soln_em_p, soln_em_s, soln_ac_u, v_eps_rel
   !
   ! Need these dependencies to get f2py calling to work
   !f2py depend(elnd_to_mshpt) P2_NODES_PER_EL, n_msh_el
   !f2py depend(el_material) n_msh_pts
   !f2py depend(x) n_msh_pts
   !f2py depend(soln_em_p) P2_NODES_PER_EL, nval_em_p, n_msh_el
   !f2py depend(soln_em_s) P2_NODES_PER_EL, nval_em_s, n_msh_el
   !f2py depend(soln_ac_u) P2_NODES_PER_EL, nval_ac_u, n_msh_el
   !f2py depend(p_tensor) n_elt_mats
   !f2py depend(v_eps_rel) n_elt_mats
   !
   !fo2py intent(out) Q_PE
   !fo2py intent(out) errco
   !fo2py intent(out) emsg



   ui = stdout

   if (is_curvilinear .ne. 0) then
      call quadint%setup_reference_quadratures()
   endif

   ! build arrays holding which modes to be calculated
   call fill_ival_arrays(v_ival_p, nval_em_p, ival_p)
   call fill_ival_arrays(v_ival_s, nval_EM_s, ival_s)
   call fill_ival_arrays(v_ival_u, nval_ac_u, ival_ac)

   call complex_alloc_4d(basis_overlap, 3*P2_NODES_PER_EL, 3*P2_NODES_PER_EL, 3_8, 3*P2_NODES_PER_EL, &
      'basis_overlap', errco, emsg)

   Q_PE = D_ZERO


   do i_el=1,n_msh_el
      typ_e = el_material(i_el)
      eps = v_eps_rel(typ_e)

      ! find positions of all the P2 nodes for this elt
      do j=1,P2_NODES_PER_EL
         nds_xy(:, j) = v_nd_xy(:,  elnd_to_mshpt(j,i_el))
      enddo

      if (is_curvilinear .ne. 0) then
         is_curved = log_is_curved_elem_tri (P2_NODES_PER_EL, nds_xy)
         if (is_curved) then
            n_curved = n_curved + 1
         endif


         call make_P2_overlaps_i_j_dk_l_quadrature(i_el, beta_ac, typ_e, is_curved, &
            eps, p_tensor, n_elt_mats, nds_xy, quadint, basis_overlap, errco, emsg)
         RETONERROR(errco)
      else

         call make_P2_overlaps_i_j_dk_l_analytic(i_el, beta_ac, typ_e, &
            eps, p_tensor, n_elt_mats, nds_xy, basis_overlap)
      endif

      do nd_i=1,P2_NODES_PER_EL          ! nodes and components of the Stokes field
         do xyz_i=1,3
            ind_ip = xyz_i + 3*(nd_i-1)

            do ivs = 1,nval_em_s          ! for the requested Stokes modes
               t_ival_s = v_ival_s(ivs)
               if (t_ival_s .eq. 0) exit

               ! (E_S^*)_i
               E_s_i_star = conjg(soln_em_s(xyz_i, nd_i, t_ival_s, i_el))

               do nd_j=1,P2_NODES_PER_EL  ! nodes and components of the pump field
                  do xyz_j=1,3
                     ind_jp = xyz_j + 3*(nd_j-1)

                     do ivp = 1,nval_em_p   ! for the requested pump modes
                        t_ival_p = v_ival_p(ivp)
                        if (t_ival_p .eq. 0) exit

                        ! (E_P^)_j
                        E_p_j = soln_em_p(xyz_j, nd_j, t_ival_p, i_el)


                        do nd_l=1,P2_NODES_PER_EL  ! nodes and components of the AC field
                           do xyz_l=1,3
                              ind_lp = xyz_l + 3*(nd_l-1)

                              do ivu = 1,nval_ac_u
                                 t_ival_u = v_ival_u(ivu)
                                 if (t_ival_u .eq. 0) exit

                                 Ustar_l = conjg(soln_ac_u(xyz_l,nd_l,t_ival_u,i_el))
                                 zt1 = E_s_i_star * E_p_j * Ustar_l

                                 do xyz_k=1,3
                                    Q_PE(t_ival_s,t_ival_p,t_ival_u) = Q_PE(t_ival_s,t_ival_p,t_ival_u) &
                                       + zt1 * basis_overlap(ind_ip,ind_jp,xyz_k,ind_lp)
                                 enddo

                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   Q_PE = - Q_PE * SI_EPS_0


end subroutine photoelastic_int_common
