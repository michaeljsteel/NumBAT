 ! Calculate the overlap integral of two EM modes and an AC mode using
 ! analytic expressions for basis function overlaps on linear elements.
 !
subroutine photoelastic_int_linear_elts (nval_EM_p, nval_EM_S, nval_AC,&
&ival1, ival2, ival3, nel, npt, nnodes, m_elnd_to_mshpt, type_el, x,&
&nb_typ_el, p_tensor, beta_AC, soln_EM_p, soln_EM_S, soln_AC,&
&eps_lst, debug, overlap, errco, emsg)

   use numbatmod
   use alloc
   integer(8) nval_EM_p, nval_EM_S, nval_AC, ival1, ival2, ival3
   integer(8) nel, npt, nnodes, nb_typ_el
   integer(8) type_el(nel), debug
   integer(8) m_elnd_to_mshpt(nnodes,nel)
   double precision x(2,npt)
   complex(8) soln_EM_p(3,nnodes,nval_EM_p,nel)
   complex(8) soln_EM_S(3,nnodes,nval_EM_S,nel)
   complex(8) soln_AC(3,nnodes,nval_AC,nel)
   complex(8) p_tensor(3,3,3,3,nb_typ_el)

   complex(8) beta_AC

   complex(8), intent(out) :: overlap(nval_EM_S, nval_EM_p, nval_AC)
   integer(8), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) ::  emsg

   !---------------------------

   integer(8) nnodes0
   parameter (nnodes0 = 6)
   double precision xel(2,nnodes0)

   !complex(8) basis_overlap(3*nnodes0,3*nnodes0,3,3*nnodes0)
   complex(8), dimension(:,:,:,:), allocatable :: basis_overlap


   complex(8) E1star, E2, Ustar
   integer(8) i, j, k, j1, typ_e
   integer(8) iel, ind_ip, i_eq
   integer(8) jtest, ind_jp, j_eq, k_eq
   integer(8) ltest, ind_lp, l_eq
   integer(8) itrial, ui, ival1s, ival2s, ival3s
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


   !fo2py intent(in) nval_EM_p, nval_EM_S, nval_AC
   !fo2py intent(in) ival1, ival2, ival3, nb_typ_el
   !fo2py intent(in) nel, npt, nnodes, m_elnd_to_mshpt, p_tensor, beta_AC, debug
   !fo2py intent(in) type_el, x, soln_EM_p, soln_EM_S, soln_AC, eps_lst
   !
   ! Need these dependencies to get f2py calling to work
   !f2py depend(m_elnd_to_mshpt) nnodes, nel
   !f2py depend(type_el) npt
   !f2py depend(x) npt
   !f2py depend(soln_EM_p) nnodes, nval_EM_p, nel
   !f2py depend(soln_EM_S) nnodes, nval_EM_S, nel
   !f2py depend(soln_AC) nnodes, nval_AC, nel
   !f2py depend(p_tensor) nb_typ_el
   !f2py depend(eps_lst) nb_typ_el
   !
   !fo2py intent(out) overlap
   !fo2py intent(out) errco
   !fo2py intent(out) emsg



   ui = stdout


   !
   if ( nnodes .ne. 6 ) then
      write(ui,*) "photoelastic_int_v2: problem nnodes = ", nnodes
      write(ui,*) "photoelastic_int_v2: nnodes should be equal to 6 !"
      write(ui,*) "photoelastic_int_v2: Aborting..."
      stop
   endif
   !
   call quad_triangle (nquad, nquad_max, wq, xq, yq)
   if (debug .eq. 1) then
      write(ui,*) "photoelastic_int_v2: nquad, nquad_max = ",&
      &nquad, nquad_max
   endif

   call complex_alloc_4d(basis_overlap, 3*nnodes0, 3*nnodes0, 3_8, 3*nnodes0, &
   'basis_overlap', errco, emsg)


   ! do i=1,nval_EM_S
   !    do j=1,nval_EM_p
   !       do k=1,nval_AC
   !          overlap(i,j,k) = 0.0d0
   !       enddo
   !    enddo
   ! enddo
   overlap = D_ZERO


   ! Loop over elements - start

   do iel=1,nel
      typ_e = type_el(iel)
      do j=1,nnodes
         j1 = m_elnd_to_mshpt(j,iel)
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

      ! do i=1,3*nnodes
      !    do j=1,3*nnodes
      !       do k=1,3
      !          do l=1,3*nnodes
      !             basis_overlap(i,j,k,l) = 0.0d0
      !          enddo
      !       enddo
      !    enddo
      ! enddo

      basis_overlap = D_ZERO

      do itrial=1,nnodes0
         do i_eq=1,3
            ind_ip = i_eq + 3*(itrial-1)
            do jtest=1,nnodes0
               do j_eq=1,3
                  ind_jp = j_eq + 3*(jtest-1)
                  !               Gradient of transverse components of basis function
                  do k_eq=1,3
                     do ltest=1,nnodes0
                        do l_eq=1,3
                           ind_lp = l_eq + 3*(ltest-1)
                           if ( k_eq .eq. 1) then
                              zt1 = p2_p2_p2x(itrial,jtest,ltest)
                           elseif ( k_eq .eq. 2) then
                              zt1 = p2_p2_p2y(itrial,jtest,ltest)
                           elseif ( k_eq .eq. 3) then
                              zt1 = p2_p2_p2(itrial,jtest,ltest)
                              zt1 = zt1 * (-C_IM_ONE* beta_AC)
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
      if (ival1 .ge. 0 .and. ival2 .ge. 0 .and. ival3 .ge. 0) then
         do itrial=1,nnodes0
            do i_eq=1,3
               ind_ip = i_eq + 3*(itrial-1)
               E1star = conjg(soln_EM_S(i_eq,itrial,ival1,iel))
               do jtest=1,nnodes0
                  do j_eq=1,3
                     ind_jp = j_eq + 3*(jtest-1)
                     E2 = soln_EM_p(j_eq,jtest,ival2,iel)
                     do ltest=1,nnodes0
                        do l_eq=1,3
                           ind_lp = l_eq + 3*(ltest-1)
                           Ustar = conjg(soln_AC(l_eq,ltest,ival3,iel))
                           do k_eq=1,3
                              zt1 = basis_overlap(ind_ip,ind_jp,k_eq,ind_lp)
                              zt1 = E1star * E2 * Ustar * zt1
                              overlap(ival1,ival2,ival3) = zt1 +&
                              &overlap(ival1,ival2,ival3)
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
         !
         ! If want overlap of given EM mode 1 and 2 and all AC modes.
      else if (ival1 .ge. 0 .and. ival2 .ge. 0 .and.&
      &ival3 .eq. -1) then
         do itrial=1,nnodes0
            do i_eq=1,3
               ind_ip = i_eq + 3*(itrial-1)
               E1star = conjg(soln_EM_S(i_eq,itrial,ival1,iel))
               do jtest=1,nnodes0
                  do j_eq=1,3
                     ind_jp = j_eq + 3*(jtest-1)
                     E2 = soln_EM_p(j_eq,jtest,ival2,iel)
                     do ltest=1,nnodes0
                        do l_eq=1,3
                           ind_lp = l_eq + 3*(ltest-1)
                           do ival3s = 1,nval_AC
                              Ustar = conjg(soln_AC(l_eq,ltest,ival3s,iel))
                              do k_eq=1,3
                                 zt1 = basis_overlap(ind_ip,ind_jp,k_eq,ind_lp)
                                 zt1 = E1star * E2 * Ustar * zt1
                                 overlap(ival1,ival2,ival3s) = zt1 +&
                                 &overlap(ival1,ival2,ival3s)
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
      else if (ival1 .ge. 0 .and. ival2 .eq. -1 .and.&
      &ival3 .eq. -1) then
         do itrial=1,nnodes0
            do i_eq=1,3
               ind_ip = i_eq + 3*(itrial-1)
               E1star = conjg(soln_EM_S(i_eq,itrial,ival1,iel))
               do jtest=1,nnodes0
                  do j_eq=1,3
                     ind_jp = j_eq + 3*(jtest-1)
                     do ival2s = 1,nval_EM_p
                        E2 = soln_EM_p(j_eq,jtest,ival2s,iel)
                        do ltest=1,nnodes0
                           do l_eq=1,3
                              ind_lp = l_eq + 3*(ltest-1)
                              do ival3s = 1,nval_AC
                                 Ustar = conjg(soln_AC(l_eq,ltest,ival3s,iel))
                                 do k_eq=1,3
                                    zt1=basis_overlap(ind_ip,ind_jp,k_eq,ind_lp)
                                    zt1 = E1star * E2 * Ustar * zt1
                                    overlap(ival1,ival2s,ival3s) = zt1 +&
                                    &overlap(ival1,ival2s,ival3s)
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
      else if (ival1 .eq. -1 .and. ival2 .ge. 0 .and.&
      &ival3 .eq. -1) then
         do itrial=1,nnodes0
            do i_eq=1,3
               ind_ip = i_eq + 3*(itrial-1)
               do ival1s = 1,nval_EM_S
                  E1star = conjg(soln_EM_S(i_eq,itrial,ival1s,iel))
                  do jtest=1,nnodes0
                     do j_eq=1,3
                        ind_jp = j_eq + 3*(jtest-1)
                        E2 = soln_EM_p(j_eq,jtest,ival2,iel)
                        do ltest=1,nnodes0
                           do l_eq=1,3
                              ind_lp = l_eq + 3*(ltest-1)
                              do ival3s = 1,nval_AC
                                 Ustar = conjg(soln_AC(l_eq,ltest,ival3s,iel))
                                 do k_eq=1,3
                                    zt1=basis_overlap(ind_ip,ind_jp,k_eq,ind_lp)
                                    zt1 = E1star * E2 * Ustar * zt1
                                    overlap(ival1s,ival2,ival3s) = zt1 +&
                                    &overlap(ival1s,ival2,ival3s)
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
      else if (ival1 .eq. -1 .and. ival2 .eq. -1 .and.&
      &ival3 .eq. -1) then
         do itrial=1,nnodes0
            do i_eq=1,3
               ind_ip = i_eq + 3*(itrial-1)
               do ival1s = 1,nval_EM_S
                  E1star = conjg(soln_EM_S(i_eq,itrial,ival1s,iel))
                  do jtest=1,nnodes0
                     do j_eq=1,3
                        ind_jp = j_eq + 3*(jtest-1)
                        do ival2s = 1,nval_EM_p
                           E2 = soln_EM_p(j_eq,jtest,ival2s,iel)
                           do ltest=1,nnodes0
                              do l_eq=1,3
                                 ind_lp = l_eq + 3*(ltest-1)
                                 do ival3s = 1,nval_AC
                                    Ustar = conjg(soln_AC(l_eq,ltest,ival3s,iel))
                                    do k_eq=1,3
                                       zt1=basis_overlap(ind_ip,ind_jp,k_eq,ind_lp)
                                       zt1 = E1star * E2 * Ustar * zt1
                                       overlap(ival1s,ival2s,ival3s) = zt1 +&
                                       &overlap(ival1s,ival2s,ival3s)
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
      else if (ival1 .eq. -1 .and. ival2 .eq. -1 .and.&
      &ival3 .ge. 0) then
         do itrial=1,nnodes0
            do i_eq=1,3
               ind_ip = i_eq + 3*(itrial-1)
               do ival1s = 1,nval_EM_S
                  E1star = conjg(soln_EM_S(i_eq,itrial,ival1s,iel))
                  do jtest=1,nnodes0
                     do j_eq=1,3
                        ind_jp = j_eq + 3*(jtest-1)
                        do ival2s = 1,nval_EM_p
                           E2 = soln_EM_p(j_eq,jtest,ival2s,iel)
                           do ltest=1,nnodes0
                              do l_eq=1,3
                                 ind_lp = l_eq + 3*(ltest-1)
                                 Ustar = conjg(soln_AC(l_eq,ltest,ival3,iel))
                                 do k_eq=1,3
                                    zt1=basis_overlap(ind_ip,ind_jp,k_eq,ind_lp)
                                    zt1 = E1star * E2 * Ustar * zt1
                                    overlap(ival1s,ival2s,ival3) = zt1 +&
                                    &overlap(ival1s,ival2s,ival3)
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
   do i=1,nval_EM_S
      do j=1,nval_EM_p
         do k=1,nval_AC
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
   !           do itrial=1,nnodes0
   !             do i_eq=1,3
   !               ind_ip = i_eq + 3*(itrial-1)
   !               do j_eq=1,3
   !                 do k_eq=1,3
   !                   do ltest=1,nnodes0
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
end subroutine photoelastic_int_v2
