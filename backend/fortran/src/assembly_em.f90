#include "numbat_decl.h"
 !  Construct the left hand and right hand matrices  mOp_stiff_re/im and mat_2
 !  for the main linear equations

subroutine assembly_em (bdy_cdn, i_base, shift_ksqr, bloch_vec, &
   perm_pp, perm_qq, mesh_raw, entities, cscmat, pbcs, nberr)


   use numbatmod
   use alloc
   use class_MeshRaw
   use class_PeriodicBCs
   use class_SparseCSC_EM
   use class_BasisFunctions
   use class_TriangleIntegrators

   type(MeshRawEM) :: mesh_raw
   type(MeshEntities) :: entities
   type(PeriodicBCs) :: pbcs
   type(SparseCSC_EM) :: cscmat
   type(NBError) :: nberr


   integer(8) bdy_cdn, i_base

   !if(E_H_field .eq. FEM_FORMULATION_E) then
   !   perm_qq = eps_eff*vacwavenum_k0**2
   !   pp = 1.0d0
   !elseif(E_H_field .eq. FEM_FORMULATION_H) then
   !   perm_qq = vacwavenum_k0**2
   !   pp = 1.0d0/eps_eff
   !endif


   complex(8) perm_pp(mesh_raw%n_elt_mats), perm_qq(mesh_raw%n_elt_mats), shift_ksqr
   double precision bloch_vec(2)


   ! -----------------
   integer(8) errco
   character(len=EMSG_LENGTH) emsg

   type(QuadIntegrator) :: quadint
   type(BasisFunctions) :: basfuncs

   integer(8), dimension(:), allocatable :: i_work



   integer(8) i_base2

   double precision xy_ref(2), wt

   integer(8) el_nds_i(P2_NODES_PER_EL)
   double precision el_nds_xy(2, P2_NODES_PER_EL)
   integer(8) n_curved
   logical is_curved

   double precision vec_phi2_i(2), curlt_phi2_i
   double precision vec_phi2_j(2), curlt_phi2_j
   double precision phi3_z_i, gradt_phi3_i(2)
   double precision phi3_z_j, gradt_phi3_j(2)


   integer(8) i,  k,  i_el, iq, typ_e
   integer(8) ety_i, msh_pt_i, eqn_i, dof_i
   integer(8) ety_j, msh_pt_j, eqn_j, dof_j
   integer(8) col_start, col_end

   double precision ddot, r_tmp1, r_tmp2
   complex(8) M_tt, M_zz, M_tz, M_zt
   complex(8) K_tt, K_zz, K_tz, K_zt
   complex(8) tperm_pp, tperm_qq
   complex(8) K_elt, M_elt

   !complex(8) val_exp(N_ENTITY_PER_EL),
   complex(8) z_phase_fact
!   double precision delta_xx(2)


   errco=0

   call integer_alloc_1d(i_work, 3*entities%n_entities, 'i_work', nberr); RET_ON_NBERR(nberr)

   !  The CSC indexing, i.e., cscmat%v_col_ptr, is 1-based
   !  But valpr.f may require the CSC indexing to be 0-based
   if (i_base .eq. 0) then ! This is the active case
      i_base2 = 1
   else
      i_base2 = 0
   endif


   ! These are the K and M matrices in Kokou's paper expressed in CSC format
   cscmat%mOp_stiff = C_ZERO
   cscmat%mOp_mass = C_ZERO


   call quadint%setup_reference_quadratures()
   n_curved = 0
   z_phase_fact = 1.0

   do i_el=1,mesh_raw%n_msh_elts              ! For each element
      typ_e = mesh_raw%v_elt_material(i_el)    ! Find the material and local material properties

      tperm_pp = perm_pp(typ_e)             !  1 (E-mode), 1/eps_r (H-mode)
      tperm_qq = perm_qq(typ_e)             !  eps_r * k0^2 (E-mode), k0^2 (H-mode)

      call mesh_raw%find_nodes_for_elt(i_el, el_nds_i, el_nds_xy, is_curved)

      if (is_curved) then
         n_curved = n_curved + 1
      endif

      ! if (bdy_cdn .eq. BCS_PERIODIC) then
      !    do j=1,P2_NODES_PER_EL
      !       mesh_pt = pbcs%iperiod_N(el_nds_i(j))
      !       if (mesh_pt .ne. 0) el_nds_i(j) = mesh_pt
      !    enddo
      ! endif

      call basfuncs%build_vector_elt_map(el_nds_i)


      ! ! TODO:  move to function
      ! if (bdy_cdn .eq. BCS_PERIODIC) then
      ! val_exp = C_ONE
      !    !  val_exp: Bloch mod ephase factor between the origin point and destination point
      !    !  For a pair of periodic points, one is chosen as origin and the other is the destination
      !    do j=1,N_ENTITY_PER_EL
      !       ip = entities%v_tags(j,i_el)
      !       mesh_pt = pbcs%iperiod_N_E_F(ip)
      !       if (mesh_pt .ne. 0) then
      !          delta_xx(:) = entities%v_xy(:,ip) - entities%v_xy(:,mesh_pt)
      !          r_tmp1 = ddot(2, bloch_vec, 1, delta_xx, 1)
      !          val_exp(j) = exp(C_IM_ONE*r_tmp1)
      !       endif
      !    enddo
      ! endif

      do iq=1,quadint%n_quad ! for each quadrature point in reference triangle

         ! find quad point in reference triangle
         call quadint%get_quad_point(iq, xy_ref, wt)
         ! Evaluate the basis functions and gradients at the quadrature point
         ! Gradients are evaluated in the actual triangle
         call basfuncs%evaluate_at_position(i_el, xy_ref, is_curved, el_nds_xy, nberr)
         RET_ON_NBERR(nberr)


         ! N_ENTITY_PER_EL is number of field sites per elt, each associated with one P2 or P3 node
         ! N_ETY_TRANSVERSE is number of transverse ones (face and 3 P2 edges)

         ! iterating columns
         do ety_j=1,N_ENTITY_PER_EL
            msh_pt_j = entities%v_tags(ety_j, i_el) ! global name for the corresopnding mesh point

            do dof_j=1,3  ! max of 3 basis functions associated with any node
               eqn_j = cscmat%m_eqs(dof_j,msh_pt_j)

               if (eqn_j .gt. 0) then
                  col_start = cscmat%v_col_ptr(eqn_j) + i_base2
                  col_end = cscmat%v_col_ptr(eqn_j+1) - 1 + i_base2
                  !  unpack row into i_work
                  do i=col_start,col_end
                     i_work(cscmat%v_row_ind(i) + i_base2) = i
                  enddo

                  ! evaluate this entity's P2 and P3 functions and derivatives
                  ! at the current quad point
                  call basfuncs%find_derivatives(dof_j, ety_j, vec_phi2_j, curlt_phi2_j, phi3_z_j, gradt_phi3_j)


                  ! iterating rows
                  do ety_i=1,N_ENTITY_PER_EL
                     !z_phase_fact = val_exp(ety_j) * conjg(val_exp(ety_i))
                     msh_pt_i = entities%v_tags(ety_i,i_el)

                     do dof_i=1,3
                        eqn_i = cscmat%m_eqs(dof_i,msh_pt_i)

                        if (eqn_i .gt. 0) then

                           call basfuncs%find_derivatives(dof_i, ety_i, vec_phi2_i, curlt_phi2_i, phi3_z_i, gradt_phi3_i)


                           !!!!!!!!!!!!!!!!!!!!!!!!!!
                           !  Reference; see Eq. (40) of the FEM paper:
                           !  K. Dossou and M. Fontaine
                           !  "A high order isoparametric finite element method for the computation of waveguide modes"
                           !  Computer Methods in Applied Mechanics and Engineering, vol. 194, no. 6-8, pp. 837-858, 2005.
                           !!!!!!!!!!!!!!!!!!!!!!!!!!

                           ! Here we are building Eqs. 13
                           ! We use capital E, F for transverse vector parts
                           !        lower    e,f for hatted longitudinal part
                           ! F_i, f_i are the test functions indexed by rows
                           ! E_i, e_i are the soln functions indexed by cols
                           if (ety_i .le. N_ETY_TRANSVERSE .and. ety_j .le. N_ETY_TRANSVERSE) then    ! [tt] part
                              ! K_tt =   (curlt E_j).(curlt F_i)-k^2 eps (E_j, F_i)
                              ! M_tt = - (E_j, F_i)
                              r_tmp1 = curlt_phi2_j * curlt_phi2_i
                              r_tmp2 = ddot(2, vec_phi2_j, 1, vec_phi2_i, 1)
                              K_tt = r_tmp1 * tperm_pp - r_tmp2 * tperm_qq
                              M_tt = - r_tmp2 * tperm_pp
                              K_elt = K_tt
                              M_elt = M_tt

                           elseif (ety_i .le. N_ETY_TRANSVERSE .and. ety_j .gt. N_ETY_TRANSVERSE) then  ! [tz] part
                              ! K_tz =  0
                              ! M_tz = (gradt e_j, F_i)
                              r_tmp1 = ddot(2, gradt_phi3_j, 1, vec_phi2_i, 1)
                              K_tz = 0.0d0
                              M_tz = r_tmp1 * tperm_pp
                              K_elt = K_tz
                              M_elt = M_tz

                           elseif (ety_i .gt. N_ETY_TRANSVERSE .and. ety_j .le. N_ETY_TRANSVERSE) then ! [zt] part
                              ! K_tz =  (E_j, gradt f_i)
                              ! M_tz =  0
                              r_tmp1 = ddot(2, vec_phi2_j, 1, gradt_phi3_i, 1)
                              K_zt = r_tmp1 * tperm_pp
                              M_zt = 0.0d0
                              K_elt = K_zt
                              M_elt = M_zt

                           elseif (ety_i .gt. N_ETY_TRANSVERSE .and. ety_j .gt. N_ETY_TRANSVERSE) then ! [zz] part
                              ! K_zz =  - (gradt e_j, gradt f_i) + eps_r k^2 (e_j, f_i )
                              ! M_zz =  0
                              r_tmp1 = ddot(2, gradt_phi3_j, 1, gradt_phi3_i, 1)
                              r_tmp2 = phi3_z_j * phi3_z_i
                              K_zz = - r_tmp1 * tperm_pp + r_tmp2 * tperm_qq
                              M_zz = 0.0d0
                              K_elt = K_zz
                              M_elt = M_zz

                           else
                              write(emsg,*) "ety_i or ety_j has an ", "invalid value",  &
                                 "ety_i ety_j, = ", ety_i, ety_j
                                 call nberr%set(NBERR_BAD_ASSEMBLY, emsg);
                              return
                           endif

                           ! without periodic bcs, z_phase_fact = 0 always
                           ! add term to sum weighted by quadrature factor and area det
                           K_elt = K_elt * wt * abs(basfuncs%det) * z_phase_fact
                           M_elt = M_elt * wt * abs(basfuncs%det) * z_phase_fact

                           ! TODO:
                           !  This is an arpack solver thing.
                           !   No reason to do this deep in loop.
                           !  Should do once the matrices have been constructed
                           !  And closer to the solver
                           !
                           K_elt = K_elt - shift_ksqr*M_elt

                           k = i_work(eqn_i)
                           if (k .gt. 0 .and. k .le. cscmat%n_nonz) then   !is this test necessary?
                              cscmat%mOp_stiff(k) = cscmat%mOp_stiff(k) + K_elt
                              cscmat%mOp_mass(k) = cscmat%mOp_mass(k) + M_elt
                           else
                              write(emsg,*) "asmbly: problem with row_ind: k, nonz = ", k, cscmat%n_nonz
                              call nberr%set(NBERR_BAD_ASSEMBLY, emsg);
                              return

                           endif

                        endif
                     enddo
                  enddo
               endif
            enddo
         enddo
      enddo
   enddo

end
