 ! Calculate the overlap integral of two EM modes and an AC mode
 ! using numerical quadrature.
 !
subroutine photoelastic_int_curvilinear_elts_old(nval_em_p, nval_em_s, nval_AC, ival_p, ival_s, ival_ac, &
   nel, npt, elnd_to_mesh, x,   &
   nb_typ_el, type_el, p_tensor, beta_AC, soln_em_p, soln_em_s, soln_AC,&
   eps_lst, overlap, errco, emsg)

   use numbatmod
   use alloc

   integer(8) nval_em_p, nval_em_s, nval_AC, ival_p, ival_s, ival_ac
   integer(8) nel, npt,  nb_typ_el
   integer(8) type_el(nel), debug
   integer(8) elnd_to_mesh(P2_NODES_PER_EL,nel)
   double precision x(2,npt)
   complex(8) soln_em_p(3,P2_NODES_PER_EL,nval_em_p,nel)
   complex(8) soln_em_s(3,P2_NODES_PER_EL,nval_em_s,nel)
   complex(8) soln_AC(3,P2_NODES_PER_EL,nval_AC,nel)
   complex(8) p_tensor(3,3,3,3,nb_typ_el)
   complex(8) beta_AC

   complex(8), intent(out) :: overlap(nval_em_s, nval_em_p, nval_AC)
   integer(8), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg

   !     Local variables

   double precision xel(2,P2_NODES_PER_EL)
   !complex(8) basis_overlap(3*P2_NODES_PER_EL0,3*P2_NODES_PER_EL0,3,3*P2_NODES_PER_EL0)
   complex(8), dimension(:,:,:,:), allocatable :: basis_overlap

   complex(8) E1star, E2, Ustar, eps
   integer(8) i, j, k, l, j1, typ_e
   integer(8) iel, ind_ip, i_eq
   integer(8) jtest, ind_jp, j_eq, k_eq
   integer(8) ltest, ind_lp, l_eq
   integer(8) itrial, ui, ival_ps, ival_ss, ival_acs
   complex(8) eps_lst(nb_typ_el)
   complex(8) zt1
   double precision mat_B(2,2), mat_T(2,2)
   !
   !     NQUAD: The number of quadrature points used in each element.
   integer(8) nquad, nquad_max, iq
   !       !  Limit to P2 polynomials
   parameter (nquad_max = 16)
   double precision wq(nquad_max)
   double precision xq(nquad_max), yq(nquad_max)
   double precision xx(2), xx_g(2), ww, det
   integer(8) n_curved
   logical is_curved

   complex(8) coeff_1, coeff_2
   double precision phi2_list(6), grad2_mat0(2,6)
   double precision grad2_mat(2,6)


   integer(8) v_ival_p(nval_em_p), v_ival_s(nval_EM_s), v_ival_ac(nval_AC)
   integer(8) ivs, ivp, ivac, t_ival_s, t_ival_p, t_ival_ac



   !fo2py intent(in) nval_em_p, nval_em_s, nval_AC
   !fo2py intent(in) ival_p, ival_s, ival_ac, nb_typ_el
   !fo2py intent(in) nel, npt, P2_NODES_PER_EL, elnd_to_mesh, p_tensor, beta_AC , debug
   !fo2py intent(in) type_el, x, soln_em_p, soln_em_s, soln_AC, eps_lst

   !f2py depend(elnd_to_mesh) P2_NODES_PER_EL, nel
   !f2py depend(type_el) npt
   !f2py depend(x) npt
   !f2py depend(soln_em_p) P2_NODES_PER_EL, nval_em_p, nel
   !f2py depend(soln_em_s) P2_NODES_PER_EL, nval_em_s, nel
   !f2py depend(soln_AC) P2_NODES_PER_EL, nval_AC, nel
   !f2py depend(p_tensor) nb_typ_el
   !f2py depend(eps_lst) nb_typ_el

   !fo2py intent(out) overlap
   !fo2py intent(out) errco
   !fo2py intent(out) emsg


   !!!!!!!!!!!!!!!!!!!!!!!!  Start Program  !!!!!!!!!!!!!!!!!!!!!!!!

   ui = stdout
debug = 0

   ! build arrays holding which modes to be calculated
   call fill_ival_arrays(v_ival_p, nval_em_p, ival_p)
   call fill_ival_arrays(v_ival_s, nval_EM_s, ival_s)
   call fill_ival_arrays(v_ival_ac, nval_ac, ival_ac)



   overlap = 0.0d0

   call quad_triangle (nquad, nquad_max, wq, xq, yq)

   if (debug .eq. 1) then
      write(ui,*) "photoelastic_int: nquad, nquad_max = ",&
      &nquad, nquad_max
   endif

   call complex_alloc_4d(basis_overlap, 3*P2_NODES_PER_EL, 3*P2_NODES_PER_EL, 3_8, 3*P2_NODES_PER_EL, &
      'basis_overlap', errco, emsg)

   ! do i=1,nval_em_s
   !    do j=1,nval_em_p
   !       do k=1,nval_AC
   !          overlap(i,j,k) = 0.0d0
   !       enddo
   !    enddo
   ! enddo
   overlap = D_ZERO


   ! Loop over elements - start

   do iel=1,nel
      typ_e = type_el(iel)
      do j=1,P2_NODES_PER_EL
         j1 = elnd_to_mesh(j,iel)
         xel(1,j) = x(1,j1)
         xel(2,j) = x(2,j1)
      enddo
      is_curved = log_is_curved_elem_tri (P2_NODES_PER_EL, xel)
      if (is_curved) then
         n_curved = n_curved + 1
      endif
      !ccccccccc
      do i=1,3*P2_NODES_PER_EL
         do j=1,3*P2_NODES_PER_EL
            do k=1,3
               do l=1,3*P2_NODES_PER_EL
                  basis_overlap(i,j,k,l) = 0.0d0
               enddo
            enddo
         enddo
      enddo
      !ccccccccc
      ! For each quadrature point evaluate overlap of Lagrange polynomials
      ! or derivative of Lagrange polynomials
      do iq=1,nquad
         xx(1) = xq(iq)
         xx(2) = yq(iq)
         ww = wq(iq)
         !         xx   = coordinate on the reference triangle
         !         xx_g = coordinate on the actual triangle
         !         phi2_list = values of Lagrange polynomials (1-6) at each local node.
         !         grad2_mat0 = gradient on the reference triangle (P2 element)
         call phi2_2d_mat(xx, phi2_list, grad2_mat0)
         !
         if (.not. is_curved) then
            !           Rectilinear element
            call jacobian_p1_2d(xx, xel, P2_NODES_PER_EL, xx_g, det, mat_B, mat_T, errco, emsg)
         else
            !           Isoparametric element !  fixed 2024/6/12
            call jacobian_p2_2d(xel, P2_NODES_PER_EL, phi2_list, grad2_mat0, xx_g, det, mat_B, mat_T, errco, emsg)
         endif
         if(abs(det) .lt. 1.0d-20) then
            write(*,*)
            write(*,*) "   ???"
            write(*,*) "PE_int: det = 0 : iel, det = ", iel, det
            write(*,*) "PE_int: Aborting..."
            stop
         endif
         !          grad_i  = gradient on the actual triangle
         !          grad_i  = Transpose(mat_T)*grad_i0
         !          Calculation of the matrix-matrix product:
         call DGEMM('Transpose','N', 2, 6, 2, D_ONE, mat_T, 2,&
         &grad2_mat0, 2, D_ZERO, grad2_mat, 2)
         coeff_1 = ww * abs(det)
         ! Calculate overlap of basis functions at quadrature point,
         ! which is a superposition of P2 polynomials for each function (field).
         do itrial=1,P2_NODES_PER_EL
            do i_eq=1,3
               ind_ip = i_eq + 3*(itrial-1)
               do jtest=1,P2_NODES_PER_EL
                  do j_eq=1,3
                     ind_jp = j_eq + 3*(jtest-1)
                     !                 Gradient of transverse components of basis function
                     do k_eq=1,2
                        do ltest=1,P2_NODES_PER_EL
                           do l_eq=1,3
                              ind_lp = l_eq + 3*(ltest-1)
                              zt1 = phi2_list(itrial) * phi2_list(jtest)&
                              &* grad2_mat(k_eq,ltest)
                              coeff_2 = p_tensor(i_eq,j_eq,k_eq,l_eq,typ_e)
                              eps = eps_lst(typ_e)
                              zt1 = coeff_1 * coeff_2 * eps**2 * zt1
                              basis_overlap(ind_ip,ind_jp,k_eq,ind_lp) =&
                              &basis_overlap(ind_ip,ind_jp,k_eq,ind_lp)&
                              &+ zt1
                           enddo
                        enddo
                     enddo
                     !                 Gradient of longitudinal components of basis function,
                     !                 which is i*beta*phi because field is assumed to be of
                     !                 form e^{i*beta*z} phi.
                     k_eq=3
                     do ltest=1,P2_NODES_PER_EL
                        do l_eq=1,3
                           ind_lp = l_eq + 3*(ltest-1)
                           zt1 = phi2_list(itrial) * phi2_list(jtest)&
                           &* phi2_list(ltest) * (-C_IM_ONE* beta_AC)
                           coeff_2 = p_tensor(i_eq,j_eq,k_eq,l_eq,typ_e)
                           eps = eps_lst(typ_e)
                           zt1 = coeff_1 * coeff_2 * eps**2 * zt1
                           basis_overlap(ind_ip,ind_jp,k_eq,ind_lp) =&
                           &basis_overlap(ind_ip,ind_jp,k_eq,ind_lp)&
                           &+ zt1
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      !
      !ccccccccc
      ! Having calculated overlap of basis functions on element
      ! now multiply by specific field values for modes of interest.
      !
      ! If only want overlap of one given combination of EM modes and AC mode.
      if (ival_p .ge. 0 .and. ival_s .ge. 0 .and. ival_ac .ge. 0) then
         do itrial=1,P2_NODES_PER_EL
            do i_eq=1,3
               ind_ip = i_eq + 3*(itrial-1)
               E1star = conjg(soln_em_s(i_eq,itrial,ival_p,iel))
               do jtest=1,P2_NODES_PER_EL
                  do j_eq=1,3
                     ind_jp = j_eq + 3*(jtest-1)
                     E2 = soln_em_p(j_eq,jtest,ival_s,iel)
                     do ltest=1,P2_NODES_PER_EL
                        do l_eq=1,3
                           ind_lp = l_eq + 3*(ltest-1)
                           Ustar = conjg(soln_AC(l_eq,ltest,ival_ac,iel))
                           do k_eq=1,3
                              zt1 = basis_overlap(ind_ip,ind_jp,k_eq,ind_lp)
                              zt1 = E1star * E2 * Ustar * zt1
                              overlap(ival_p,ival_s,ival_ac) = zt1 +&
                              &overlap(ival_p,ival_s,ival_ac)
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
         !
         ! If want overlap of given EM mode 1 and 2 and all AC modes.
      else if (ival_p .ge. 0 .and. ival_s .ge. 0 .and.&
      &ival_ac .eq. -1) then
         do itrial=1,P2_NODES_PER_EL
            do i_eq=1,3
               ind_ip = i_eq + 3*(itrial-1)
               E1star = conjg(soln_em_s(i_eq,itrial,ival_p,iel))
               do jtest=1,P2_NODES_PER_EL
                  do j_eq=1,3
                     ind_jp = j_eq + 3*(jtest-1)
                     E2 = soln_em_p(j_eq,jtest,ival_s,iel)
                     do ltest=1,P2_NODES_PER_EL
                        do l_eq=1,3
                           ind_lp = l_eq + 3*(ltest-1)
                           do ival_acs = 1,nval_AC
                              Ustar = conjg(soln_AC(l_eq,ltest,ival_acs,iel))
                              do k_eq=1,3
                                 zt1 = basis_overlap(ind_ip,ind_jp,k_eq,ind_lp)
                                 zt1 = E1star * E2 * Ustar * zt1
                                 overlap(ival_p,ival_s,ival_acs) = zt1 +&
                                 &overlap(ival_p,ival_s,ival_acs)
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
         !
         ! If want overlap of given EM mode 1 and all EM modes 2 and all AC modes.
      else if (ival_p .ge. 0 .and. ival_s .eq. -1 .and.&
      &ival_ac .eq. -1) then
         do itrial=1,P2_NODES_PER_EL
            do i_eq=1,3
               ind_ip = i_eq + 3*(itrial-1)
               E1star = conjg(soln_em_s(i_eq,itrial,ival_p,iel))
               do jtest=1,P2_NODES_PER_EL
                  do j_eq=1,3
                     ind_jp = j_eq + 3*(jtest-1)
                     do ival_ss = 1,nval_em_p
                        E2 = soln_em_p(j_eq,jtest,ival_ss,iel)
                        do ltest=1,P2_NODES_PER_EL
                           do l_eq=1,3
                              ind_lp = l_eq + 3*(ltest-1)
                              do ival_acs = 1,nval_AC
                                 Ustar = conjg(soln_AC(l_eq,ltest,ival_acs,iel))
                                 do k_eq=1,3
                                    zt1=basis_overlap(ind_ip,ind_jp,k_eq,ind_lp)
                                    zt1 = E1star * E2 * Ustar * zt1
                                    overlap(ival_p,ival_ss,ival_acs) = zt1 +&
                                    &overlap(ival_p,ival_ss,ival_acs)
                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
         !
         ! If want overlap of given EM mode 2 and all EM modes 1 and all AC modes.
      else if (ival_p .eq. -1 .and. ival_s .ge. 0 .and.&
      &ival_ac .eq. -1) then
         do itrial=1,P2_NODES_PER_EL
            do i_eq=1,3
               ind_ip = i_eq + 3*(itrial-1)
               do ival_ps = 1,nval_em_s
                  E1star = conjg(soln_em_s(i_eq,itrial,ival_ps,iel))
                  do jtest=1,P2_NODES_PER_EL
                     do j_eq=1,3
                        ind_jp = j_eq + 3*(jtest-1)
                        E2 = soln_em_p(j_eq,jtest,ival_s,iel)
                        do ltest=1,P2_NODES_PER_EL
                           do l_eq=1,3
                              ind_lp = l_eq + 3*(ltest-1)
                              do ival_acs = 1,nval_AC
                                 Ustar = conjg(soln_AC(l_eq,ltest,ival_acs,iel))
                                 do k_eq=1,3
                                    zt1=basis_overlap(ind_ip,ind_jp,k_eq,ind_lp)
                                    zt1 = E1star * E2 * Ustar * zt1
                                    overlap(ival_ps,ival_s,ival_acs) = zt1 +&
                                    &overlap(ival_ps,ival_s,ival_acs)
                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
         !
         ! If want overlap of all EM mode 1, all EM modes 2 and all AC modes.
      else if (ival_p .eq. -1 .and. ival_s .eq. -1 .and.&
      &ival_ac .eq. -1) then
         do itrial=1,P2_NODES_PER_EL
            do i_eq=1,3
               ind_ip = i_eq + 3*(itrial-1)
               do ival_ps = 1,nval_em_s
                  E1star = conjg(soln_em_s(i_eq,itrial,ival_ps,iel))
                  do jtest=1,P2_NODES_PER_EL
                     do j_eq=1,3
                        ind_jp = j_eq + 3*(jtest-1)
                        do ival_ss = 1,nval_em_p
                           E2 = soln_em_p(j_eq,jtest,ival_ss,iel)
                           do ltest=1,P2_NODES_PER_EL
                              do l_eq=1,3
                                 ind_lp = l_eq + 3*(ltest-1)
                                 do ival_acs = 1,nval_AC
                                    Ustar = conjg(soln_AC(l_eq,ltest,ival_acs,iel))
                                    do k_eq=1,3
                                       zt1=basis_overlap(ind_ip,ind_jp,k_eq,ind_lp)
                                       zt1 = E1star * E2 * Ustar * zt1
                                       overlap(ival_ps,ival_ss,ival_acs) = zt1 +&
                                       &overlap(ival_ps,ival_ss,ival_acs)
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
         !
         ! If want overlap of all EM mode 1, all EM modes 2 and one AC mode.
      else if (ival_p .eq. -1 .and. ival_s .eq. -1 .and.&
      &ival_ac .ge. 0) then
         do itrial=1,P2_NODES_PER_EL
            do i_eq=1,3
               ind_ip = i_eq + 3*(itrial-1)
               do ival_ps = 1,nval_em_s
                  E1star = conjg(soln_em_s(i_eq,itrial,ival_ps,iel))
                  do jtest=1,P2_NODES_PER_EL
                     do j_eq=1,3
                        ind_jp = j_eq + 3*(jtest-1)
                        do ival_ss = 1,nval_em_p
                           E2 = soln_em_p(j_eq,jtest,ival_ss,iel)
                           do ltest=1,P2_NODES_PER_EL
                              do l_eq=1,3
                                 ind_lp = l_eq + 3*(ltest-1)
                                 Ustar = conjg(soln_AC(l_eq,ltest,ival_ac,iel))
                                 do k_eq=1,3
                                    zt1=basis_overlap(ind_ip,ind_jp,k_eq,ind_lp)
                                    zt1 = E1star * E2 * Ustar * zt1
                                    overlap(ival_ps,ival_ss,ival_ac) = zt1 +&
                                    &overlap(ival_ps,ival_ss,ival_ac)
                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
      !ccccccccccc
      ! Loop over elements - end
      !ccccccccccc
   enddo
   ! Apply scaling that sits outside of integration.
   do i=1,nval_em_s
      do j=1,nval_em_p
         do k=1,nval_AC
            overlap(i,j,k) = overlap(i,j,k) * (-1.0d0) * SI_EPS_0
         enddo
      enddo
   enddo
   if(debug .eq. 1) then
      write(*,*) "PE_int: overlap"
      write(*,*) overlap
   endif
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
end subroutine photoelastic_int_curvilinear_elts_old
