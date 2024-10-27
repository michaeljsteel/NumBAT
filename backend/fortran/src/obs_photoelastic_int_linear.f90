 ! Calculate the overlap integral of two EM modes and an AC mode using
 ! analytic expressions for basis function overlaps on linear elements.
 !
subroutine photoelastic_int_linear_elts (nval_em_p, nval_em_s, nval_ac, ival_p, ival_s, ival_ac, &
    nel, npt, elnd_to_mshpt, type_el, x,&
    &nb_typ_el, p_tensor, beta_ac, soln_em_p, soln_em_s, soln_ac,&
    &eps_lst, overlap, errco, emsg)

       use numbatmod
       use alloc
       integer(8) nval_em_p, nval_em_s, nval_ac, ival_p, ival_s, ival_ac
       integer(8) nel, npt, nb_typ_el
       integer(8) type_el(nel), debug
       integer(8) elnd_to_mshpt(P2_NODES_PER_EL,nel)
       double precision x(2,npt)
       complex(8) soln_em_p(3,P2_NODES_PER_EL,nval_em_p,nel)
       complex(8) soln_em_s(3,P2_NODES_PER_EL,nval_em_s,nel)
       complex(8) soln_ac(3,P2_NODES_PER_EL,nval_ac,nel)
       complex(8) p_tensor(3,3,3,3,nb_typ_el)

       complex(8) beta_ac

       complex(8), intent(out) :: overlap(nval_em_s, nval_em_p, nval_ac)
       integer(8), intent(out) :: errco
       character(len=EMSG_LENGTH), intent(out) ::  emsg

       !---------------------------


       double precision xel(2,P2_NODES_PER_EL)

       !complex(8) basis_overlap(3*P2_NODES_PER_EL,3*P2_NODES_PER_EL,3,3*P2_NODES_PER_EL)
       complex(8), dimension(:,:,:,:), allocatable :: basis_overlap


       complex(8) E1star, E2, Ustar
       integer(8) i, j, k, j1, typ_e
       integer(8) iel, ind_ip, i_eq
       integer(8) jtest, ind_jp, j_eq, k_eq
       integer(8) ltest, ind_lp, l_eq
       integer(8) itrial, ui, ival_ps, ival_ss, ival_acs
       complex(8) eps_lst(nb_typ_el)
       complex(8) zt1
       double precision mat_B(2,2), mat_T(2,2), mat_T_tr(2,2)
       double precision det_b

       !     NQUAD: The number of quadrature points used in each element.
       integer(8) nquad, nquad_max
       !  Limit to P2 polynomials
       parameter (nquad_max = 16)
       double precision wq(nquad_max)
       double precision xq(nquad_max), yq(nquad_max)
       !c      integer(8) is_curved, n_curved
       complex(8) coeff

       double precision p2_p2_p2(6,6,6)
       double precision p2_p2_p2x(6,6,6), p2_p2_p2y(6,6,6)


       debug = 0
       !fo2py intent(in) nval_em_p, nval_em_s, nval_ac
       !fo2py intent(in) ival_p, ival_s, ival_ac, nb_typ_el
       !fo2py intent(in) nel, npt, P2_NODES_PER_EL, elnd_to_mshpt, p_tensor, beta_ac, debug
       !fo2py intent(in) type_el, x, soln_em_p, soln_em_s, soln_ac, eps_lst
       !
       ! Need these dependencies to get f2py calling to work
       !f2py depend(elnd_to_mshpt) P2_NODES_PER_EL, nel
       !f2py depend(type_el) npt
       !f2py depend(x) npt
       !f2py depend(soln_em_p) P2_NODES_PER_EL, nval_em_p, nel
       !f2py depend(soln_em_s) P2_NODES_PER_EL, nval_em_s, nel
       !f2py depend(soln_ac) P2_NODES_PER_EL, nval_ac, nel
       !f2py depend(p_tensor) nb_typ_el
       !f2py depend(eps_lst) nb_typ_el
       !
       !fo2py intent(out) overlap
       !fo2py intent(out) errco
       !fo2py intent(out) emsg



       ui = stdout


       !
       if ( P2_NODES_PER_EL .ne. 6 ) then
          write(ui,*) "photoelastic_int_v2: problem P2_NODES_PER_EL = ", P2_NODES_PER_EL
          write(ui,*) "photoelastic_int_v2: P2_NODES_PER_EL should be equal to 6 !"
          write(ui,*) "photoelastic_int_v2: Aborting..."
          stop
       endif
       !
       call quad_triangle (nquad, nquad_max, wq, xq, yq)
       if (debug .eq. 1) then
          write(ui,*) "photoelastic_int_v2: nquad, nquad_max = ",&
          &nquad, nquad_max
       endif

       call complex_alloc_4d(basis_overlap, 3*P2_NODES_PER_EL, 3*P2_NODES_PER_EL, 3_8, 3*P2_NODES_PER_EL, &
       'basis_overlap', errco, emsg)


       ! do i=1,nval_em_s
       !    do j=1,nval_em_p
       !       do k=1,nval_ac
       !          overlap(i,j,k) = 0.0d0
       !       enddo
       !    enddo
       ! enddo
       overlap = D_ZERO


       ! Loop over elements - start

       do iel=1,nel
          typ_e = type_el(iel)
          do j=1,P2_NODES_PER_EL
             j1 = elnd_to_mshpt(j,iel)
             xel(1,j) = x(1,j1)
             xel(2,j) = x(2,j1)
          enddo

          !       The geometric transformation (x,y) -> (x_g,y_g) = mat_B*(x,y)^t + (x_0, y_0, z_0)^t
          !       maps the current triangle to the reference triangle.
          do i=1,2
             do j=1,2
                mat_B(j,i) = xel(j,i+1) - xel(j,1)
             enddo
          enddo
          det_b = mat_B(1,1) * mat_B(2,2) - mat_B(1,2) * mat_B(2,1)
          if (abs(det_b) .le. 1.0d-22) then
             write(*,*) '?? PE_int_v2: Determinant = 0 :', det_b
             write(*,*) "xel = ", xel
             write(*,*) 'Aborting...'
             stop
          endif

          !       We also need, is the matrix mat_T of the reverse transformation
          !                (from reference to current triangle):
          !       mat_T = inverse matrix of de mat_B
          mat_T(1,1) =  mat_B(2,2) / det_b
          mat_T(2,2) =  mat_B(1,1) / det_b
          mat_T(1,2) = -mat_B(1,2) / det_b
          mat_T(2,1) = -mat_B(2,1) / det_b

          !       Note that if grad_i_0 is the gradient on the reference triangle,
          !       then the gradient on the actual triangle is:
          !       grad_i  = Transpose(mat_T)*grad_i0
          !
          !       mat_T_tr = Transpose(mat_T)
          mat_T_tr(1,1) = mat_T(1,1)
          mat_T_tr(2,2) = mat_T(2,2)
          mat_T_tr(1,2) = mat_T(2,1)
          mat_T_tr(2,1) = mat_T(1,2)

          call find_overlaps_p2_p2_p2 (p2_p2_p2, det_b)
          call find_overlaps_p2_p2_p2x (p2_p2_p2x, mat_T_tr, det_b)
          call find_overlaps_p2_p2_p2y (p2_p2_p2y, mat_T_tr, det_b)

          ! do i=1,3*P2_NODES_PER_EL
          !    do j=1,3*P2_NODES_PER_EL
          !       do k=1,3
          !          do l=1,3*P2_NODES_PER_EL
          !             basis_overlap(i,j,k,l) = 0.0d0
          !          enddo
          !       enddo
          !    enddo
          ! enddo

          basis_overlap = D_ZERO

          do itrial=1,P2_NODES_PER_EL
             do i_eq=1,3
                ind_ip = i_eq + 3*(itrial-1)
                do jtest=1,P2_NODES_PER_EL
                   do j_eq=1,3
                      ind_jp = j_eq + 3*(jtest-1)
                      !               Gradient of transverse components of basis function
                      do k_eq=1,3
                         do ltest=1,P2_NODES_PER_EL
                            do l_eq=1,3
                               ind_lp = l_eq + 3*(ltest-1)
                               if ( k_eq .eq. 1) then
                                  zt1 = p2_p2_p2x(itrial,jtest,ltest)
                               elseif ( k_eq .eq. 2) then
                                  zt1 = p2_p2_p2y(itrial,jtest,ltest)
                               elseif ( k_eq .eq. 3) then
                                  zt1 = p2_p2_p2(itrial,jtest,ltest)
                                  zt1 = zt1 * (-C_IM_ONE* beta_ac)
                               else
                                  write(*,*) "--- photoelastic_int_v2: "
                                  write(*,*) "k_eq has illegal value:"
                                  write(*,*) "k_eq = ", k_eq
                                  write(*,*) "Aborting..."
                                  stop
                               endif
                               coeff = p_tensor(i_eq,j_eq,k_eq,l_eq,typ_e)
                               zt1 = coeff * eps_lst(typ_e)**2 * zt1
                               basis_overlap(ind_ip,ind_jp,k_eq,ind_lp) = zt1
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
                               Ustar = conjg(soln_ac(l_eq,ltest,ival_ac,iel))
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
                               do ival_acs = 1,nval_ac
                                  Ustar = conjg(soln_ac(l_eq,ltest,ival_acs,iel))
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
                                  do ival_acs = 1,nval_ac
                                     Ustar = conjg(soln_ac(l_eq,ltest,ival_acs,iel))
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
                                  do ival_acs = 1,nval_ac
                                     Ustar = conjg(soln_ac(l_eq,ltest,ival_acs,iel))
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
                                     do ival_acs = 1,nval_ac
                                        Ustar = conjg(soln_ac(l_eq,ltest,ival_acs,iel))
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
                                     Ustar = conjg(soln_ac(l_eq,ltest,ival_ac,iel))
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
             do k=1,nval_ac
                overlap(i,j,k) = overlap(i,j,k) * (-1.0d0) * SI_EPS_0
             enddo
          enddo
       enddo
       !ccccccccccc
       if(debug .eq. 1) then
          write(*,*) "PE_int_v2: overlap"
          write(*,*) overlap
       endif
       !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !
       !       open (unit=26,file="Output/overlap_v2.txt")
       !         write(26,*) "overlap, SI_EPS_0 = "
       !         write(26,*) overlap, SI_EPS_0
       !       close (unit=26)
       !         open(4,file="Output/basis_overlap_v2.txt",status='unknown')
       !           iel = nel
       !           do itrial=1,P2_NODES_PER_EL
       !             do i_eq=1,3
       !               ind_ip = i_eq + 3*(itrial-1)
       !               do j_eq=1,3
       !                 do k_eq=1,3
       !                   do ltest=1,P2_NODES_PER_EL
       !                     do l_eq=1,3
       !                       ind_lp = l_eq + 3*(ltest-1)
       !                       zt1  = basis_overlap(ind_ip,j_eq,k_eq,ind_lp)
       !                       if (zt1 .ne. 0) then
       !                         write(4,*) ind_ip,j_eq,k_eq,ind_lp,
       !      *                  abs(zt1), zt1
       !                       endif
       !                     enddo
       !                   enddo
       !                 enddo
       !               enddo
       !             enddo
       !           enddo
       !         close(4)
       !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !
    end subroutine photoelastic_int_linear_elts
