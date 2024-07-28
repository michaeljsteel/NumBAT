 ! Calculate the overlap integral of two EM modes and an AC mode
 ! using numerical quadrature.
 !
subroutine photoelastic_int (nval_EM_p, nval_EM_S, nval_AC, ival1,&
&ival2, ival3, nel, npt, nnodes, table_nod, type_el, x,&
&nb_typ_el, p_tensor, beta_AC, soln_EM_p, soln_EM_S, soln_AC,&
&eps_lst, debug, overlap, errco, emsg)

   use numbatmod
   use alloc

   integer(8) nval_EM_p, nval_EM_S, nval_AC, ival1, ival2, ival3
   integer(8) nel, npt, nnodes, nb_typ_el
   integer(8) type_el(nel), debug
   integer(8) table_nod(nnodes,nel)
   double precision x(2,npt)
   complex(8) soln_EM_p(3,nnodes,nval_EM_p,nel)
   complex(8) soln_EM_S(3,nnodes,nval_EM_S,nel)
   complex(8) soln_AC(3,nnodes,nval_AC,nel)
   complex(8) p_tensor(3,3,3,3,nb_typ_el)
   complex(8) beta_AC

   complex(8), intent(out) :: overlap(nval_EM_S, nval_EM_p, nval_AC)
   integer(8), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg

   !     Local variables

   double precision xel(2,nnodes_0)
   !complex(8) basis_overlap(3*nnodes0,3*nnodes0,3,3*nnodes0)
   complex(8), dimension(:,:,:,:), allocatable :: basis_overlap

   complex(8) E1star, E2, Ustar, eps
   integer(8) i, j, k, l, j1, typ_e
   integer(8) iel, ind_ip, i_eq
   integer(8) jtest, ind_jp, j_eq, k_eq
   integer(8) ltest, ind_lp, l_eq
   integer(8) itrial, ui, ival1s, ival2s, ival3s
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
    
    
   !fo2py intent(in) nval_EM_p, nval_EM_S, nval_AC
   !fo2py intent(in) ival1, ival2, ival3, nb_typ_el
   !fo2py intent(in) nel, npt, nnodes, table_nod, p_tensor, beta_AC , debug
   !fo2py intent(in) type_el, x, soln_EM_p, soln_EM_S, soln_AC, eps_lst
    
   !f2py depend(table_nod) nnodes, nel
   !f2py depend(type_el) npt
   !f2py depend(x) npt
   !f2py depend(soln_EM_p) nnodes, nval_EM_p, nel
   !f2py depend(soln_EM_S) nnodes, nval_EM_S, nel
   !f2py depend(soln_AC) nnodes, nval_AC, nel
   !f2py depend(p_tensor) nb_typ_el
   !f2py depend(eps_lst) nb_typ_el
    
   !fo2py intent(out) overlap
   !fo2py intent(out) errco
   !fo2py intent(out) emsg
    
    
   !!!!!!!!!!!!!!!!!!!!!!!!  Start Program  !!!!!!!!!!!!!!!!!!!!!!!!
    
   ui = stdout


   !
   if ( nnodes .ne. 6 ) then
      write(ui,*) "photoelastic_int: problem nnodes = ", nnodes
      write(ui,*) "photoelastic_int: nnodes should be equal to 6 !"
      write(ui,*) "photoelastic_int: Aborting..."
      stop
   endif

   overlap = 0.0d0
   call quad_triangle (nquad, nquad_max, wq, xq, yq)
   if (debug .eq. 1) then
      write(ui,*) "photoelastic_int: nquad, nquad_max = ",&
      &nquad, nquad_max
   endif

   call complex_alloc_4d(basis_overlap, 3*nnodes, 3*nnodes, 3_8, 3*nnodes, &
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
         j1 = table_nod(j,iel)
         xel(1,j) = x(1,j1)
         xel(2,j) = x(2,j1)
      enddo
      is_curved = log_is_curved_elem_tri (nnodes, xel)
      if (is_curved) then
         n_curved = n_curved + 1
      endif
      !ccccccccc
      do i=1,3*nnodes
         do j=1,3*nnodes
            do k=1,3
               do l=1,3*nnodes
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
            call jacobian_p1_2d(xx, xel, nnodes,&
            &xx_g, det, mat_B, mat_T)
         else
            !           Isoparametric element !  fixed 2024/6/12
            call jacobian_p2_2d(xel, nnodes, phi2_list,&
            &grad2_mat0, xx_g, det, mat_B, mat_T)
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
         do itrial=1,nnodes_0
            do i_eq=1,3
               ind_ip = i_eq + 3*(itrial-1)
               do jtest=1,nnodes_0
                  do j_eq=1,3
                     ind_jp = j_eq + 3*(jtest-1)
                     !                 Gradient of transverse components of basis function
                     do k_eq=1,2
                        do ltest=1,nnodes_0
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
                     do ltest=1,nnodes_0
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
      if (ival1 .ge. 0 .and. ival2 .ge. 0 .and. ival3 .ge. 0) then
         do itrial=1,nnodes_0
            do i_eq=1,3
               ind_ip = i_eq + 3*(itrial-1)
               E1star = conjg(soln_EM_S(i_eq,itrial,ival1,iel))
               do jtest=1,nnodes_0
                  do j_eq=1,3
                     ind_jp = j_eq + 3*(jtest-1)
                     E2 = soln_EM_p(j_eq,jtest,ival2,iel)
                     do ltest=1,nnodes_0
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
         do itrial=1,nnodes_0
            do i_eq=1,3
               ind_ip = i_eq + 3*(itrial-1)
               E1star = conjg(soln_EM_S(i_eq,itrial,ival1,iel))
               do jtest=1,nnodes_0
                  do j_eq=1,3
                     ind_jp = j_eq + 3*(jtest-1)
                     E2 = soln_EM_p(j_eq,jtest,ival2,iel)
                     do ltest=1,nnodes_0
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
         do itrial=1,nnodes_0
            do i_eq=1,3
               ind_ip = i_eq + 3*(itrial-1)
               E1star = conjg(soln_EM_S(i_eq,itrial,ival1,iel))
               do jtest=1,nnodes_0
                  do j_eq=1,3
                     ind_jp = j_eq + 3*(jtest-1)
                     do ival2s = 1,nval_EM_p
                        E2 = soln_EM_p(j_eq,jtest,ival2s,iel)
                        do ltest=1,nnodes_0
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
         do itrial=1,nnodes_0
            do i_eq=1,3
               ind_ip = i_eq + 3*(itrial-1)
               do ival1s = 1,nval_EM_S
                  E1star = conjg(soln_EM_S(i_eq,itrial,ival1s,iel))
                  do jtest=1,nnodes_0
                     do j_eq=1,3
                        ind_jp = j_eq + 3*(jtest-1)
                        E2 = soln_EM_p(j_eq,jtest,ival2,iel)
                        do ltest=1,nnodes_0
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
         do itrial=1,nnodes_0
            do i_eq=1,3
               ind_ip = i_eq + 3*(itrial-1)
               do ival1s = 1,nval_EM_S
                  E1star = conjg(soln_EM_S(i_eq,itrial,ival1s,iel))
                  do jtest=1,nnodes_0
                     do j_eq=1,3
                        ind_jp = j_eq + 3*(jtest-1)
                        do ival2s = 1,nval_EM_p
                           E2 = soln_EM_p(j_eq,jtest,ival2s,iel)
                           do ltest=1,nnodes_0
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
         do itrial=1,nnodes_0
            do i_eq=1,3
               ind_ip = i_eq + 3*(itrial-1)
               do ival1s = 1,nval_EM_S
                  E1star = conjg(soln_EM_S(i_eq,itrial,ival1s,iel))
                  do jtest=1,nnodes_0
                     do j_eq=1,3
                        ind_jp = j_eq + 3*(jtest-1)
                        do ival2s = 1,nval_EM_p
                           E2 = soln_EM_p(j_eq,jtest,ival2s,iel)
                           do ltest=1,nnodes_0
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
   if(debug .eq. 1) then
      write(*,*) "PE_int: overlap"
      write(*,*) overlap
   endif
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
end subroutine photoelastic_int
