! Calculate the overlap integral of an AC mode with itself using
! analytic expressions for basis function overlaps on linear elements.
!

! P_z = Re \int_A \dxdy (-2 i \Omega) c_zjkl u_j^* d_k u_l
subroutine AC_mode_power_int_v2 (n_modes, n_msh_el, n_msh_pts, &
   elnd_to_mesh, v_el_material, v_nd_xy, n_elt_mats, &
   stiffC_zjkl, beta_AC, Omega_AC, soln_ac_u, overlap)

   use numbatmod
   integer(8) n_modes, ival
   integer(8) n_msh_el, n_msh_pts, n_elt_mats
   integer(8) v_el_material(n_msh_el)
   integer(8) elnd_to_mesh(P2_NODES_PER_EL,n_msh_el)
   double precision v_nd_xy(2,n_msh_pts)
   complex(8) soln_ac_u(3,P2_NODES_PER_EL,n_modes,n_msh_el)
   complex(8) Omega_AC(n_modes)
   complex(8) beta_AC
   complex(8), dimension(n_modes) :: overlap
   complex(8) stiffC_zjkl(3,3,3,n_elt_mats)

   ! Locals

   integer(8) nod_el_p(P2_NODES_PER_EL)
   double precision xel(2,P2_NODES_PER_EL)
   complex(8) basis_overlap(3*P2_NODES_PER_EL,3,3*P2_NODES_PER_EL)
   complex(8) U, Ustar
   integer(8) i, j, j1, typ_e
   integer(8) iel, ind_ip, i_eq, k_eq
   integer(8) ltest, ind_lp, l_eq
   integer(8) itrial, ui
   double precision p2_p2(6,6), p2_p2x(6,6), p2_p2y(6,6)
   double precision mat_B(2,2), mat_T(2,2), mat_T_tr(2,2)
   double precision det_b

   complex(8) z_tmp1
   complex(8) coeff

!
!
!f2py intent(in) n_modes, n_msh_el, n_msh_pts, P2_NODES_PER_EL, elnd_to_mesh
!f2py intent(in) v_el_material, x, n_elt_mats, stiffC_zjkl, beta_AC
!f2py intent(in) soln_ac_u, debug, Omega_AC
!
!f2py depend(elnd_to_mesh) P2_NODES_PER_EL, n_msh_el
!f2py depend(v_el_material) n_msh_pts
!f2py depend(v_nd_xy) n_msh_pts
!f2py depend(soln_ac_u) P2_NODES_PER_EL, n_modes, n_msh_el
!f2py depend(stiffC_zjkl) n_elt_mats
!f2py depend(Omega_AC) n_modes
!
!f2py intent(out) overlap
!
!
   !!!!!!!!!!!!!!!!!!!!!!!!  Start Program  !!!!!!!!!!!!!!!!!!!!!!!!
!
   ui = stdout
!
   if ( P2_NODES_PER_EL .ne. 6 ) then
      write(ui,*) "AC_mode_power_int_v2: problem P2_NODES_PER_EL = ",&
      &P2_NODES_PER_EL
      write(ui,*) " --------- P2_NODES_PER_EL should be equal to 6 !"
      write(ui,*) "AC_mode_power_int_v2: Aborting..."
      stop
   endif
!
   do i=1,n_modes
      overlap(i) = 0.0d0
   enddo
!
!ccccccccccc
! Loop over elements - start
!ccccccccccc
   do iel=1,n_msh_el
      typ_e = v_el_material(iel)
      do j=1,P2_NODES_PER_EL
         j1 = elnd_to_mesh(j,iel)
         nod_el_p(j) = j1
         xel(1,j) = v_nd_xy(1,j1)
         xel(2,j) = v_nd_xy(2,j1)
      enddo
      do i=1,2
         do j=1,2
            mat_B(j,i) = xel(j,i+1) - xel(j,1)
         enddo
      enddo
      det_b = mat_B(1,1) * mat_B(2,2) - mat_B(1,2) * mat_B(2,1)
!       !  TEMPORARY CHANGE
      if (abs(det_b) .le. 1.0d-22) then
!c      if (abs(det_b) .le. 1.0d-8) then
         write(*,*) '?? AC_alpha_int_v2: Determinant = 0 :', det_b
         write(*,*) "xel = ", xel
         write(*,*) 'Aborting...'
         stop
      endif
!     mat_T = Inverse of mat_B
      mat_T(1,1) = mat_B(2,2) / det_b
      mat_T(2,2) = mat_B(1,1) / det_b
      mat_T(1,2) = -mat_B(1,2) / det_b
      mat_T(2,1) = -mat_B(2,1) / det_b
!
!
!	mat_T_tr = Tanspose(mat_T)
      mat_T_tr(1,1) = mat_T(1,1)
      mat_T_tr(1,2) = mat_T(2,1)
      mat_T_tr(2,1) = mat_T(1,2)
      mat_T_tr(2,2) = mat_T(2,2)

      call find_overlaps_p2_p2(p2_p2, det_b)
      call find_overlaps_p2_p2x (p2_p2x, mat_T_tr, det_b)
      call find_overlaps_p2_p2y (p2_p2y, mat_T_tr, det_b)
      do itrial=1,P2_NODES_PER_EL
         do i_eq=1,3
            ind_ip = i_eq + 3*(itrial-1)
!         Gradient of transverse components of basis function
            do k_eq=1,3
               do ltest=1,P2_NODES_PER_EL
                  do l_eq=1,3
                     ind_lp = l_eq + 3*(ltest-1)
                     if(k_eq == 1) then
                        z_tmp1 = p2_p2x(itrial,ltest)
                     elseif(k_eq == 2) then
                        z_tmp1 = p2_p2y(itrial,ltest)
                     elseif(k_eq == 3) then
                        z_tmp1 = p2_p2(itrial,ltest) * C_IM_ONE* beta_AC
                     else
                        write(ui,*) "AC_mode_power_int_v2: in_modesid value "
                        write(ui,*) "AC_mode_power_int_v2: k_eq = ", k_eq
                        write(ui,*) "AC_mode_power_int_v2: Aborting..."
                        stop
                     endif
                     coeff = stiffC_zjkl(i_eq,k_eq,l_eq,typ_e)
                     z_tmp1 = coeff * z_tmp1
                     basis_overlap(ind_ip,k_eq,ind_lp) = z_tmp1
                  enddo
               enddo
            enddo
         enddo
      enddo
!ccccccccc
! Having calculated overlap of basis functions on element
! now multiply by specific field values for modes of interest.
      do ival=1,n_modes
         do itrial=1,P2_NODES_PER_EL
            do i_eq=1,3
               ind_ip = i_eq + 3*(itrial-1)
               Ustar = conjg(soln_ac_u(i_eq,itrial,ival,iel))
               do ltest=1,P2_NODES_PER_EL
                  do l_eq=1,3
                     ind_lp = l_eq + 3*(ltest-1)
                     U = soln_ac_u(l_eq,ltest,ival,iel)
                     do k_eq=1,3
                        z_tmp1 = basis_overlap(ind_ip,k_eq,ind_lp)
                        overlap(ival) = overlap(ival) + Ustar * U * z_tmp1
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!ccccccccccc
! Loop over elements - end
!ccccccccccc
   enddo
! Multiply through prefactor
   do i=1,n_modes
      overlap(i) = -2.0 * C_IM_ONE* Omega_AC(i) * overlap(i)
   enddo

!       open (unit=26,file="Output/overlap.txt")
!       do i=1,n_modes
!         write(26,*) i, Omega_AC(i), abs(overlap(i)),
!      *              overlap(i)
!       enddo
!       close (unit=26)
!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end subroutine AC_mode_power_int_v2
