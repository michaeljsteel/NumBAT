#include "numbat_decl.h"
! Calculate the v_alpha integral of an AC mode with itself using
! Direct integration

! \alpha = \Omega^2/Energy_aC \int  eta_ijkl d_i u_j^* d_k u_l

subroutine ac_alpha_analytic (n_modes, n_msh_elts, n_msh_pts, &
   elnd_to_mshpt, v_mshpt_xy, n_elt_mats, v_elt_material,  &
   eta_ijkl, q_AC, Omega_AC, soln_ac_u, &
   v_ac_mode_energy, v_alpha_r, errco, emsg)

   use numbatmod
   use class_TriangleIntegrators
   use class_BasisFunctions

   integer(8) n_modes, md_i
   integer(8) n_msh_elts, n_msh_pts, n_elt_mats
   integer(8) v_elt_material(n_msh_elts)
   integer(8) elnd_to_mshpt(P2_NODES_PER_EL,n_msh_elts)
   double precision v_mshpt_xy(2,n_msh_pts)

   complex(8) soln_ac_u(3,P2_NODES_PER_EL,n_modes,n_msh_elts)
   complex(8) Omega_AC(n_modes)
   complex(8) q_AC, v_ac_mode_energy(n_modes)
   complex(8) eta_ijkl(3,3,3,3,n_elt_mats)
   double precision, dimension(n_modes), intent(out) :: v_alpha_r
   integer(8), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg

   ! Locals
   complex(8), dimension(n_modes) :: v_alpha
   double precision nds_xy(2,P2_NODES_PER_EL)

   complex(8) bas_ovrlp(3,3*P2_NODES_PER_EL,3,3*P2_NODES_PER_EL)
   complex(8) U, Ustar, Umod
   integer(8) typ_e
   integer(8) i_el, ind_j, xyz_i, xyz_k
   integer(8) bf_l, ind_l, xyz_l, xyz_j
   integer(8) bf_j
   complex(8) z_tmp1

   double precision p2_p2(6,6), p2_p2x(6,6), p2_p2y(6,6)
   double precision p2x_p2x(6,6), p2y_p2y(6,6), p2x_p2y(6,6)
   double precision det_b

   complex(8) t_eta


   type(NBError) nberr
   type(PyFrontEnd) frontend
   type(BasisFunctions) basfuncs


!f2py intent(in) n_modes, n_msh_elts, n_msh_pts, P2_NODES_PER_EL, elnd_to_mshpt
!f2py intent(in) v_elt_material, x, n_elt_mats, eta_ijkl, q_AC
!f2py intent(in) soln_ac_u, debug, Omega_AC, v_ac_mode_energy

!f2py depend(elnd_to_mshpt) P2_NODES_PER_EL, n_msh_elts
!f2py depend(v_elt_material) n_msh_pts
!f2py depend(v_mshpt_xy) n_msh_pts
!f2py depend(soln_ac_u) P2_NODES_PER_EL, n_modes, n_msh_elts
!f2py depend(eta_ijkl) n_elt_mats
!f2py depend(Omega_AC) n_modes
!f2py depend(v_ac_mode_energy) n_modes


   errco = 0
   emsg = ""
   call nberr%reset()

   call frontend%init_from_py(n_msh_elts, n_msh_pts, elnd_to_mshpt, v_mshpt_xy, nberr)
   RET_ON_NBERR_UNFOLD(nberr)

   v_alpha = D_ZERO
   z_tmp1 = C_ZERO

   do i_el=1,n_msh_elts
      typ_e = v_elt_material(i_el)

      call frontend%nodes_at_el(i_el, nds_xy)

      call basfuncs%set_affine_for_elt(nds_xy, nberr)
      RET_ON_NBERR(nberr)

      det_b = basfuncs%det
      call find_overlaps_p2_p2(p2_p2, det_b)
      call find_overlaps_p2_p2x (p2_p2x, basfuncs%mat_T_tr, det_b)
      call find_overlaps_p2_p2y (p2_p2y, basfuncs%mat_T_tr, det_b)
      call find_overlaps_p2x_p2x (p2x_p2x, basfuncs%mat_T_tr, det_b)
      call find_overlaps_p2x_p2y (p2x_p2y, basfuncs%mat_T_tr, det_b)
      call find_overlaps_p2y_p2y (p2y_p2y, basfuncs%mat_T_tr, det_b)


      ! Calculate v_alpha of basis functions
      ! which is a superposition of P2 polynomials for each function ().
      do xyz_i=1,3

         do bf_j=1,P2_NODES_PER_EL
            do xyz_j=1,3
               ind_j = xyz_j + 3*(bf_j-1)

               do xyz_k=1,3

                  do bf_l=1,P2_NODES_PER_EL
                     do xyz_l=1,3
                        ind_l = xyz_l + 3*(bf_l-1)

                        ! See Eq. (45) of C. Wolff et al. PRB (2015)
                        if (xyz_i == 1 .and. xyz_k == 1) then
                           z_tmp1 = p2x_p2x(bf_j,bf_l)
                        elseif (xyz_i == 1 .and. xyz_k == 2) then
                           z_tmp1 = p2x_p2y(bf_j,bf_l)
                        elseif (xyz_i == 1 .and. xyz_k == 3) then
                           z_tmp1 = C_IM_ONE* q_AC * p2_p2x(bf_l,bf_j)
                        elseif (xyz_i == 2 .and. xyz_k == 1) then
                           z_tmp1 = p2x_p2y(bf_l,bf_j)
                        elseif (xyz_i == 2 .and. xyz_k == 2) then
                           z_tmp1 = p2y_p2y(bf_j,bf_l)
                        elseif (xyz_i == 2 .and. xyz_k == 3) then
                           z_tmp1 = C_IM_ONE* q_AC * p2_p2y(bf_l,bf_j)
                        elseif (xyz_i == 3 .and. xyz_k == 1) then
                           z_tmp1 = -C_IM_ONE* q_AC * p2_p2x(bf_j,bf_l)  ! - sign because of u_j^*
                        elseif (xyz_i == 3 .and. xyz_k == 2) then
                           z_tmp1 = -C_IM_ONE* q_AC * p2_p2y(bf_j,bf_l)
                        elseif (xyz_i == 3 .and. xyz_k == 3) then
                           z_tmp1 = q_AC**2 * p2_p2(bf_j,bf_l)
                        endif
                        t_eta = eta_ijkl(xyz_i,xyz_j,xyz_k,xyz_l,typ_e)
                        bas_ovrlp(xyz_i,ind_j,xyz_k,ind_l) = t_eta * z_tmp1
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      ! Having calculated v_alpha of basis functions on element
      ! now multiply by specific field values for modes of interest.
      do md_i=1,n_modes

         do bf_j=1,P2_NODES_PER_EL
            do xyz_j=1,3
               ind_j = xyz_j + 3*(bf_j-1)
               Ustar = conjg(soln_ac_u(xyz_j,bf_j,md_i,i_el))

               do bf_l=1,P2_NODES_PER_EL
                  do xyz_l=1,3
                     ind_l = xyz_l + 3*(bf_l-1)
                     U = soln_ac_u(xyz_l,bf_l,md_i,i_el)

                     Umod = U*Ustar
                     do xyz_i=1,3
                        do xyz_k=1,3
                           v_alpha(md_i) = v_alpha(md_i) + Umod * bas_ovrlp(xyz_i,ind_j,xyz_k,ind_l)
                        enddo
                     enddo

                  enddo
               enddo

            enddo
         enddo
      enddo
   enddo

   v_alpha_r = real(v_alpha * Omega_AC**2 / v_ac_mode_energy)

end subroutine ac_alpha_analytic

