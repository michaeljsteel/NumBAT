! Calculate the overlap integral of an AC mode with itself using
! Direct integration
!
subroutine AC_alpha_int_v2 (nval,&
&nel, npt, nnodes, table_nod, type_el, x,&
&nb_typ_el, eta_tensor, beta_AC, Omega_AC, soln_AC,&
&AC_mode_energy_elastic, overlap)
!
   use numbatmod
   integer(8) nval, ival
   integer(8) nel, npt, nnodes, nb_typ_el
   integer(8) type_el(nel)
   integer(8) table_nod(nnodes,nel)
   double precision x(2,npt)
!       complex(8) x(2,npt)
   complex(8) soln_AC(3,nnodes,nval,nel)
   complex(8) Omega_AC(nval)
   complex(8) beta_AC, AC_mode_energy_elastic(nval)
   complex(8), dimension(nval) :: overlap
   complex(8) eta_tensor(3,3,3,3,nb_typ_el)
   ! integer(8) errco
   ! character(len=EMSG_LENGTH) emsg

!     Local variables
   integer(8) nnodes0
   parameter (nnodes0 = 6)
   integer(8) nod_el_p(nnodes0)
   double precision xel(2,nnodes0)
   complex(8) basis_overlap(3*nnodes0,3,3,3*nnodes0)
   complex(8) U, Ustar
   integer(8) i, j, j1, typ_e
   integer(8) iel, ind_ip, i_eq, k_eq
   integer(8) ltest, ind_lp, l_eq, j_eq
   integer(8) itrial, ui
   complex(8) z_tmp1

   double precision p2_p2(6,6), p2_p2x(6,6), p2_p2y(6,6)
   double precision p2x_p2x(6,6), p2y_p2y(6,6), p2x_p2y(6,6)
   double precision mat_B(2,2), mat_T(2,2), mat_T_tr(2,2)
   double precision det_b
!
   complex(8) coeff
!
!
!f2py intent(in) nval, nel, npt, nnodes, table_nod
!f2py intent(in) type_el, x, nb_typ_el, eta_tensor, beta_AC
!f2py intent(in) soln_AC, debug, Omega_AC, AC_mode_energy_elastic
!
!f2py depend(table_nod) nnodes, nel
!f2py depend(type_el) npt
!f2py depend(x) npt
!f2py depend(soln_AC) nnodes, nval, nel
!f2py depend(eta_tensor) nb_typ_el
!f2py depend(Omega_AC) nval
!f2py depend(AC_mode_energy_elastic) nval
!
!f2py intent(out) overlap

!
!
   !!!!!!!!!!!!!!!!!!!!!!!!  Start Program  !!!!!!!!!!!!!!!!!!!!!!!!
!
   z_tmp1=0.d0
   ui = stdout

!
   if ( nnodes .ne. 6 ) then
      write(ui,*) "AC_alpha_int_v2: problem nnodes = ", nnodes
      write(ui,*) "AC_alpha_int_v2: nnodes should be equal to 6 !"
      write(ui,*) "AC_alpha_int_v2: Aborting..."
      stop
   endif

   overlap = D_ZERO

   do iel=1,nel
      typ_e = type_el(iel)
      do j=1,nnodes
         j1 = table_nod(j,iel)
         nod_el_p(j) = j1
         xel(1,j) = x(1,j1)
         xel(2,j) = x(2,j1)
      enddo

      do i=1,2
         do j=1,2
            mat_B(j,i) = xel(j,i+1) - xel(j,1)
         enddo
      enddo
      det_b = mat_B(1,1) * mat_B(2,2) - mat_B(1,2) * mat_B(2,1)
      if (abs(det_b) .le. 1.0d-22) then
         write(*,*) '?? AC_alpha_int_v2: Determinant = 0 :', det_b
         write(*,*) "xel = ", xel
         write(*,*) 'Aborting...'
         stop
      endif
!
!     mat_T = Inverse of mat_B
      mat_T(1,1) = mat_B(2,2) / det_b
      mat_T(2,2) = mat_B(1,1) / det_b
      mat_T(1,2) = -mat_B(1,2) / det_b
      mat_T(2,1) = -mat_B(2,1) / det_b
!
!	mat_T_tr = Tanspose(mat_T)
      mat_T_tr(1,1) = mat_T(1,1)
      mat_T_tr(1,2) = mat_T(2,1)
      mat_T_tr(2,1) = mat_T(1,2)
      mat_T_tr(2,2) = mat_T(2,2)
!
      call mat_p2_p2(p2_p2, det_b)
      call mat_p2_p2x (p2_p2x, mat_T_tr, det_b)
      call mat_p2_p2y (p2_p2y, mat_T_tr, det_b)
      call mat_p2x_p2x (p2x_p2x, mat_T_tr, det_b)
      call mat_p2x_p2y (p2x_p2y, mat_T_tr, det_b)
      call mat_p2y_p2y (p2y_p2y, mat_T_tr, det_b)


! Calculate overlap of basis functions
! which is a superposition of P2 polynomials for each function (field).
      do itrial=1,nnodes0
         do i_eq=1,3
            ind_ip = i_eq + 3*(itrial-1)
            do j_eq=1,3
               do k_eq=1,3
                  do ltest=1,nnodes0
                     do l_eq=1,3
                        ind_lp = l_eq + 3*(ltest-1)
!                     See Eq. (45) of C. Wolff et al. PRB (2015)
                        if(j_eq == 1 .and. k_eq == 1) then
                           z_tmp1 = p2x_p2x(itrial,ltest)
                        elseif(j_eq == 1 .and. k_eq == 2) then
                           z_tmp1 = p2x_p2y(itrial,ltest)
                        elseif(j_eq == 1 .and. k_eq == 3) then
                           z_tmp1 = p2_p2x(ltest,itrial)
                           z_tmp1 = z_tmp1 * (C_IM_ONE* beta_AC)
                           !!!!!!!!!!!!!!!!!!!!!!!!
                        elseif(j_eq == 2 .and. k_eq == 1) then
                           z_tmp1 = p2x_p2y(ltest,itrial)
                        elseif(j_eq == 2 .and. k_eq == 2) then
                           z_tmp1 = p2y_p2y(itrial,ltest)
                        elseif(j_eq == 2 .and. k_eq == 3) then
                           z_tmp1 = p2_p2y(ltest,itrial)
                           z_tmp1 = z_tmp1 * (C_IM_ONE* beta_AC)
                           !!!!!!!!!!!!!!!!!!!!!!!!
                        elseif(j_eq == 3 .and. k_eq == 1) then
                           z_tmp1 = p2_p2x(itrial,ltest)
                           z_tmp1 = z_tmp1 * (-C_IM_ONE* beta_AC)
                        elseif(j_eq == 3 .and. k_eq == 2) then
                           z_tmp1 = p2_p2y(itrial,ltest)
                           z_tmp1 = z_tmp1 * (-C_IM_ONE* beta_AC)
                        elseif(j_eq == 3 .and. k_eq == 3) then
                           z_tmp1 = p2_p2(itrial,ltest)
                           z_tmp1 = z_tmp1 *  beta_AC**2
                        endif
                        coeff = eta_tensor(i_eq,j_eq,k_eq,l_eq,typ_e)
                        basis_overlap(ind_ip,j_eq,k_eq,ind_lp) =&
                        &coeff * z_tmp1
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

! Having calculated overlap of basis functions on element
! now multiply by specific field values for modes of interest.
      do ival=1,nval
         do itrial=1,nnodes0
            do i_eq=1,3
               ind_ip = i_eq + 3*(itrial-1)
               Ustar = conjg(soln_AC(i_eq,itrial,ival,iel))
               do ltest=1,nnodes0
                  do l_eq=1,3
                     ind_lp = l_eq + 3*(ltest-1)
                     U = soln_AC(l_eq,ltest,ival,iel)
                     do j_eq=1,3
                        do k_eq=1,3
                           z_tmp1 = basis_overlap(ind_ip,j_eq,k_eq,ind_lp)
                           z_tmp1 = Ustar * U * z_tmp1
                           overlap(ival) = overlap(ival) + z_tmp1
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo

! Multiply through prefactor
   do i=1,nval
!         z_tmp1 = -1.0 * Omega_AC(i)**2 / AC_mode_energy_elastic(i)
!       Flipped sign as assuming did not do integration by parts - going off CW advice.
      z_tmp1 = Omega_AC(i)**2 / AC_mode_energy_elastic(i)
      overlap(i) = z_tmp1 * overlap(i)
   enddo

end subroutine AC_alpha_int_v2
