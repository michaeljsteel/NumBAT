#include "numbat_decl.h"
! Calculate the v_alpha integral of an AC mode with itself using
! Direct integration

! \alpha = \Omega^2/Energy_aC \int  eta_ijkl d_i u_j^* d_k u_l

subroutine ac_alpha_analytic (n_modes, n_msh_el, n_msh_pts, &
   elnd_to_mshpt, v_nd_xy, n_elt_mats, v_el_material,  &
    eta_ijkl, q_AC, Omega_AC, soln_ac_u, &
   AC_mode_energy_elastic, v_alpha)

   use numbatmod
   use class_TriangleIntegrators
   use class_BasisFunctions

   integer(8) n_modes, ival
   integer(8) n_msh_el, n_msh_pts, n_elt_mats
   integer(8) v_el_material(n_msh_el)
   integer(8) elnd_to_mshpt(P2_NODES_PER_EL,n_msh_el)
   double precision v_nd_xy(2,n_msh_pts)

   complex(8) soln_ac_u(3,P2_NODES_PER_EL,n_modes,n_msh_el)
   complex(8) Omega_AC(n_modes)
   complex(8) q_AC, AC_mode_energy_elastic(n_modes)
   complex(8), dimension(n_modes) :: v_alpha
   complex(8) eta_ijkl(3,3,3,3,n_elt_mats)

   ! Locals
   integer(8) errco
   character(len=EMSG_LENGTH) emsg

   integer(8) nod_el_p(P2_NODES_PER_EL)
   double precision xel(2,P2_NODES_PER_EL)
   double precision nds_xy(2,P2_NODES_PER_EL)

   complex(8) basis_v_alpha(3*P2_NODES_PER_EL,3,3,3*P2_NODES_PER_EL)
   complex(8) U, Ustar
   integer(8) i, j, j1, typ_e
   integer(8) i_el, ind_i, xyz_i, xyz_k
   integer(8) bf_l, ind_l, xyz_l, xyz_j
   integer(8) bf_i, ui
   complex(8) z_tmp1

   double precision p2_p2(6,6), p2_p2x(6,6), p2_p2y(6,6)
   double precision p2x_p2x(6,6), p2y_p2y(6,6), p2x_p2y(6,6)
   double precision mat_B(2,2), mat_T(2,2), mat_T_tr(2,2)
   double precision det_b

   complex(8) t_eta


   type(NBError) nberr
   type(PyFrontEnd) frontend
   type(BasisFunctions) basfuncs


!f2py intent(in) n_modes, n_msh_el, n_msh_pts, P2_NODES_PER_EL, elnd_to_mshpt
!f2py intent(in) v_el_material, x, n_elt_mats, eta_ijkl, q_AC
!f2py intent(in) soln_ac_u, debug, Omega_AC, AC_mode_energy_elastic

!f2py depend(elnd_to_mshpt) P2_NODES_PER_EL, n_msh_el
!f2py depend(v_el_material) n_msh_pts
!f2py depend(x) n_msh_pts
!f2py depend(soln_ac_u) P2_NODES_PER_EL, n_modes, n_msh_el
!f2py depend(eta_ijkl) n_elt_mats
!f2py depend(Omega_AC) n_modes
!f2py depend(AC_mode_energy_elastic) n_modes
!
!f2py intent(out) v_alpha


   z_tmp1=0.d0
   ui = stdout

   errco = 0
   call nberr%reset()

   call frontend%init_from_py(n_msh_el, n_msh_pts, elnd_to_mshpt, v_nd_xy, nberr)
   RET_ON_NBERR_UNFOLD(nberr)


   v_alpha = D_ZERO

   do i_el=1,n_msh_el
      typ_e = v_el_material(i_el)

      call frontend%nodes_at_el(i_el, nds_xy)

      call basfuncs%set_affine_for_elt(nds_xy, nberr)
      RET_ON_NBERR(nberr)

      do j=1,P2_NODES_PER_EL
         j1 = elnd_to_mshpt(j,i_el)
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
      call find_overlaps_p2_p2(p2_p2, det_b)
      call find_overlaps_p2_p2x (p2_p2x, mat_T_tr, det_b)
      call find_overlaps_p2_p2y (p2_p2y, mat_T_tr, det_b)
      call find_overlaps_p2x_p2x (p2x_p2x, mat_T_tr, det_b)
      call find_overlaps_p2x_p2y (p2x_p2y, mat_T_tr, det_b)
      call find_overlaps_p2y_p2y (p2y_p2y, mat_T_tr, det_b)


! Calculate v_alpha of basis functions
! which is a superposition of P2 polynomials for each function (fi_eld).
      do bf_i=1,P2_NODES_PER_EL
         do xyz_i=1,3
            ind_i = xyz_i + 3*(bf_i-1)
            do xyz_j=1,3
               do xyz_k=1,3
                  do bf_l=1,P2_NODES_PER_EL
                     do xyz_l=1,3
                        ind_l = xyz_l + 3*(bf_l-1)
!                     See Eq. (45) of C. Wolff et al. PRB (2015)
                        if(xyz_j == 1 .and. xyz_k == 1) then
                           z_tmp1 = p2x_p2x(bf_i,bf_l)
                        elseif(xyz_j == 1 .and. xyz_k == 2) then
                           z_tmp1 = p2x_p2y(bf_i,bf_l)
                        elseif(xyz_j == 1 .and. xyz_k == 3) then
                           z_tmp1 = p2_p2x(bf_l,bf_i)
                           z_tmp1 = z_tmp1 * (C_IM_ONE* q_AC)

                        elseif(xyz_j == 2 .and. xyz_k == 1) then
                           z_tmp1 = p2x_p2y(bf_l,bf_i)
                        elseif(xyz_j == 2 .and. xyz_k == 2) then
                           z_tmp1 = p2y_p2y(bf_i,bf_l)
                        elseif(xyz_j == 2 .and. xyz_k == 3) then
                           z_tmp1 = p2_p2y(bf_l,bf_i)
                           z_tmp1 = z_tmp1 * (C_IM_ONE* q_AC)
                           !!!!!!!!!!!!!!!!!!!!!!!!
                        elseif(xyz_j == 3 .and. xyz_k == 1) then
                           z_tmp1 = p2_p2x(bf_i,bf_l)
                           z_tmp1 = z_tmp1 * (-C_IM_ONE* q_AC)
                        elseif(xyz_j == 3 .and. xyz_k == 2) then
                           z_tmp1 = p2_p2y(bf_i,bf_l)
                           z_tmp1 = z_tmp1 * (-C_IM_ONE* q_AC)
                        elseif(xyz_j == 3 .and. xyz_k == 3) then
                           z_tmp1 = p2_p2(bf_i,bf_l)
                           z_tmp1 = z_tmp1 *  q_AC**2
                        endif
                        t_eta = eta_ijkl(xyz_i,xyz_j,xyz_k,xyz_l,typ_e)
                        basis_v_alpha(ind_i,xyz_j,xyz_k,ind_l) = t_eta * z_tmp1
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

! Having calculated v_alpha of basis functions on element
! now multiply by specific fi_eld values for modes of interest.
      do ival=1,n_modes
         do bf_i=1,P2_NODES_PER_EL
            do xyz_i=1,3
               ind_i = xyz_i + 3*(bf_i-1)
               Ustar = conjg(soln_ac_u(xyz_i,bf_i,ival,i_el))
               do bf_l=1,P2_NODES_PER_EL
                  do xyz_l=1,3
                     ind_l = xyz_l + 3*(bf_l-1)
                     U = soln_ac_u(xyz_l,bf_l,ival,i_el)
                     do xyz_j=1,3
                        do xyz_k=1,3
                           z_tmp1 = Ustar * U * basis_v_alpha(ind_i,xyz_j,xyz_k,ind_l)
                           v_alpha(ival) = v_alpha(ival) + z_tmp1
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo

! Multiply through prefactor
   do i=1,n_modes
!         z_tmp1 = -1.0 * Omega_AC(i)**2 / AC_mode_energy_elastic(i)
!       Flipped sign as assuming did not do integration by parts - going off CW advice.
      z_tmp1 = Omega_AC(i)**2 / AC_mode_energy_elastic(i)
      v_alpha(i) = z_tmp1 * v_alpha(i)
   enddo

end subroutine ac_alpha_analytic
