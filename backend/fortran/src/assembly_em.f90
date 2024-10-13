#include "numbat_decl.h"
 !  Construct the left hand and right hand matrices  mOp_stiff_re/im and mat_2
 !  for the main linear equations

!call find_basis_derivatives(j_eq, jcol, phi_vec_map, phi_P2_ref, phi_P3_ref, &
!                     gradt_P1_act, gradt_P2_act, gradt_P3_act, vec_phi_j, curlt_phi_j, gradt_j, phi_z_j)



subroutine find_basis_derivatives(idof, ifunc, phi_vec_map, phi_P2_ref, phi_P3_ref, &
   gradt_P1_act, gradt_P2_act, gradt_P3_act, vec_phi_i, curlt_phi_i, gradt_i, phi_z_i)

   use numbatmod

   integer(8) idof, ifunc
   integer(8) phi_vec_map(4,3,N_DDL_T)
   double precision phi_P2_ref(P2_NODES_PER_EL), phi_P3_ref(P3_NODES_PER_EL)
   double precision gradt_P1_act(2,P1_NODES_PER_EL), gradt_P2_act(2,P2_NODES_PER_EL), gradt_P3_act(2,P3_NODES_PER_EL)

   double precision vec_phi_i(2), curlt_phi_i
   double precision gradt_i(2)
   double precision phi_z_i

   if (ifunc .le. N_DDL_T) then ! A transverse dof (edge or face)
      ! Uses P2 vector elements so determine the basis vector
      call make_phi_vector_basis(idof, ifunc, phi_vec_map, phi_P2_ref, &
         gradt_P1_act, gradt_P2_act, vec_phi_i, curlt_phi_i)
      gradt_i = D_ZERO
      phi_z_i = D_ZERO
   else   ! a longitudinal dof, use P3 scalar element
      vec_phi_i = D_ZERO
      curlt_phi_i = D_ZERO
      gradt_i(:) = gradt_P3_act(:,ifunc-N_DDL_T)
      phi_z_i = phi_P3_ref(ifunc-N_DDL_T)
   endif


end subroutine


subroutine assembly  (bdy_cdn, i_base, n_msh_el, n_msh_pts, n_ddl, neq, nnodes, &
   shift_ksqr, bloch_vec, nb_typ_el, perm_pp, perm_qq, &
   mesh_raw, entities, &
   m_eqs, ip_period_N, ip_period_E_F, &
   nonz, row_ind, col_ptr, &
   mOp_stiff, mOp_mass, errco, emsg)

   !  NQUAD: The number of quadrature points used in each element.

   use numbatmod
   use alloc
   use class_MeshRaw

   type(MeshRaw) :: mesh_raw
   type(MeshEntities) :: entities


   integer(8) bdy_cdn, i_base,  nb_typ_el, nonz
   integer(8) n_msh_el, n_msh_pts, n_ddl, neq, nnodes

   !if(E_H_field .eq. FEM_FORMULATION_E) then
   !   perm_qq = eps_eff*vacwavenum_k0**2
   !   pp = 1.0d0
   !elseif(E_H_field .eq. FEM_FORMULATION_H) then
   !   perm_qq = vacwavenum_k0**2
   !   pp = 1.0d0/eps_eff
   !endif


   complex(8) perm_pp(nb_typ_el), perm_qq(nb_typ_el), shift_ksqr
   double precision bloch_vec(2)

   integer(8) ip_period_N(n_msh_pts), ip_period_E_F(n_ddl)

   integer(8) m_eqs(3,n_ddl)

   integer(8) row_ind(nonz), col_ptr(neq+1)

   complex(8), intent(out) :: mOp_stiff(nonz), mOp_mass(nonz)

   integer(8) errco
   character(len=EMSG_LENGTH) emsg

   ! -----------------
   integer(8), dimension(:), allocatable :: i_work



   integer(8) i_base2

   integer(8), parameter :: nquad_max = 25
   integer(8) nquad

   double precision wq(nquad_max)
   double precision xq(nquad_max), yq(nquad_max)
   double precision x_act(2), x_ref(2), ww, det
   double precision mat_B(2,2), mat_T(2,2)
   double precision gradt_i(2), gradt_j(2)
   double precision phi_z_i, phi_z_j

   integer(8)   ui_stdout

   integer(8), parameter :: nddl_0 = 14


   integer(8) nod_el_p(P2_NODES_PER_EL)
   integer(8) phi_vec_map(4,3,N_DDL_T)
   double precision el_xy(2, P2_NODES_PER_EL)

   ! values of basis functions and gradients at given point in reference and actual triangles
   double precision phi1_ref(P1_NODES_PER_EL)
   double precision gradt_P1_ref(2,P1_NODES_PER_EL), gradt_P1_act(2,P1_NODES_PER_EL)

   double precision phi_P2_ref(P2_NODES_PER_EL)
   double precision gradt_P2_ref(2,P2_NODES_PER_EL), gradt_P2_act(2,P2_NODES_PER_EL)

   double precision phi_P3_ref(P3_NODES_PER_EL)
   double precision gradt_P3_ref(2,P3_NODES_PER_EL), gradt_P3_act(2,P3_NODES_PER_EL)

   double precision vec_phi_i(2), curlt_phi_i
   double precision vec_phi_j(2), curlt_phi_j

   complex(8) val_exp(nddl_0), z_phase_fact

   integer(8) i, j, k, j_mshpt, iel, iq, typ_e
   integer(8) jcol, jp, ind_jp, j_eq
   integer(8) irow, ip, ind_ip, i_eq
   integer(8) n_curved, debug, col_start, col_end
   complex(8) K_elt, M_elt
   logical is_curved

   double precision delta_xx(2)
   double precision ddot, r_tmp1, r_tmp2
   complex(8) M_tt, M_zz, M_tz, M_zt
   complex(8) K_tt, K_zz, K_tz, K_zt
   complex(8) tperm_pp, tperm_qq


   ui_stdout = stdout
   debug = 0

   call integer_alloc_1d(i_work, 3*n_ddl, 'i_work', errco, emsg); RETONERROR(errco)

   !  The CSC indexing, i.e., col_ptr, is 1-based
   !  But valpr.f may have changed the CSC indexing to 0-based indexing)
   if (i_base .eq. 0) then
      i_base2 = 1
   else
      i_base2 = 0
   endif

   ! Determine quadrature weights for triangle integrations at 16 points (seets nquad)
   call quad_triangle (nquad, nquad_max, wq, xq, yq)

   ! These are the K and M matrices in Kokou's paper expressed in CSR format
   mOp_stiff = C_ZERO
   mOp_mass = C_ZERO


   n_curved = 0

   do iel=1,n_msh_el                     ! For each element
      typ_e = mesh_raw%el_material(iel)       ! Find the material

      tperm_pp = perm_pp(typ_e)             !  1 (E-mode), 1/eps_r (H-mode)
      tperm_qq = perm_qq(typ_e)             !  eps_r * k0^2 (E-mode), k0^2 (H-mode)

      do j=1,nnodes                              ! For each of the 6 P2 nodes
         j_mshpt = mesh_raw%elnd_to_mesh(j,iel)        !    find the index of the mesh point
         nod_el_p(j) = j_mshpt                        !    store the mesh point indices for this element
         el_xy(:,j) = mesh_raw%v_nd_xy(:,j_mshpt)  !    find their physical positions
      enddo

      is_curved = log_is_curved_elem_tri (nnodes, el_xy)

      if (is_curved) then
         n_curved = n_curved + 1
      endif

      if (bdy_cdn .eq. BCS_PERIODIC) then
         do j=1,nnodes
            j_mshpt = ip_period_N(nod_el_p(j))
            if (j_mshpt .ne. 0) nod_el_p(j) = j_mshpt
         enddo
      endif

      call make_phi_vector_map(nod_el_p, phi_vec_map)
   !    if (iel .eq. 1) then
   !       write(*,*) 'phimap 1 1', (phi_vec_map(k,1, 1), k=1,4)
   !       write(*,*) 'phimap 1 2', (phi_vec_map(k,1, 2), k=1,4)
   !       write(*,*) 'phimap 1 3', (phi_vec_map(k,1, 3), k=1,4)
   !       write(*,*) 'phimap 1 4', (phi_vec_map(k,1, 4), k=1,4)

   !       write(*,*) 'phimap 2 1', (phi_vec_map(k,2, 1), k=1,4)
   !       write(*,*) 'phimap 2 2', (phi_vec_map(k,2, 2), k=1,4)
   !       write(*,*) 'phimap 2 3', (phi_vec_map(k,2, 3), k=1,4)
   !       write(*,*) 'phimap 2 4', (phi_vec_map(k,2, 4), k=1,4)

   !       write(*,*) 'phimap 2 1', (phi_vec_map(k,3, 1), k=1,4)
   !       write(*,*) 'phimap 2 2', (phi_vec_map(k,3, 2), k=1,4)
   !       write(*,*) 'phimap 2 3', (phi_vec_map(k,3, 3), k=1,4)
   !       write(*,*) 'phimap 2 4', (phi_vec_map(k,3, 4), k=1,4)


   ! endif

      val_exp = C_ONE

      if (bdy_cdn .eq. BCS_PERIODIC) then
         !  val_exp: Bloch mod ephase factor between the origin point and destination point
         !  For a pair of periodic points, one is chosen as origin and the other is the destination
         do j=1,nddl_0
            ip = entities%v_tags(j,iel)
            j_mshpt = ip_period_E_F(ip)
            if (j_mshpt .ne. 0) then
               delta_xx(:) = entities%v_xy(:,ip) - entities%v_xy(:,j_mshpt)
               r_tmp1 = ddot(2, bloch_vec, 1, delta_xx, 1)
               val_exp(j) = exp(C_IM_ONE*r_tmp1)
            endif
         enddo
      endif

      do iq=1,nquad         ! for each quadrature point in reference triangle
         x_ref(1) = xq(iq)     !    find local points and weighting
         x_ref(2) = yq(iq)
         ww = wq(iq)
         !  xx   = coordinate on the reference triangle
         !  xx_g = coordinate on the actual triangle

         ! Evaluate the basis functions and gradients at the quadrature point

         call phi1_2d_mat(x_ref, phi1_ref, gradt_P1_ref)  ! P1 elements
         call phi2_2d_mat(x_ref, phi_P2_ref, gradt_P2_ref)  ! P2 elements
         call phi3_2d_mat(x_ref, phi_P3_ref, gradt_P3_ref)  ! P3 elements

         if (.not. is_curved ) then
            !  Rectilinear element
            call jacobian_p1_2d(x_ref, el_xy, nnodes, x_act, det, mat_B, mat_T, errco, emsg)

            if (det .le. 0 .and. debug .eq. 1 .and. iq .eq. 1) then
               write(ui_stdout,*) "   !!!"
               write(ui_stdout,*) "asmbly: det <= 0: iel, det ", iel, det
               write(ui_stdout,*) "x : ", (nod_el_p(j),j=1,nnodes)
               write(ui_stdout,*) "x : ", (el_xy(1,j),j=1,3)
               write(ui_stdout,*) "y : ", (el_xy(2,j),j=1,3)
               write(ui_stdout,*)
            endif

         else !  Isoparametric element, 2024-06 fix
            call jacobian_p2_2d(el_xy, nnodes, phi_P2_ref, gradt_P2_ref, x_act, det, mat_B, mat_T)
         endif

         !  gradt_i_mat  = gradtient on the actual triangle
         !  gradt_i_act  = Transpose(mat_T)*gradt_i_ref

         call DGEMM('Transpose','N', 2, 3, 2,  D_ONE, mat_T, 2, gradt_P1_ref, 2, D_ZERO, gradt_P1_act, 2)
         call DGEMM('Transpose','N', 2, 6, 2,  D_ONE, mat_T, 2, gradt_P2_ref, 2, D_ZERO, gradt_P2_act, 2)
         call DGEMM('Transpose','N', 2, 10, 2, D_ONE, mat_T, 2, gradt_P3_ref, 2, D_ZERO, gradt_P3_act, 2)

         ! nddl_0 is number of field dof per elt, each associated with one mesh point
         ! N_DDL_T is number of transverse ones
         do jcol=1,nddl_0
            jp = entities%v_tags(jcol,iel)

            do j_eq=1,3
               !  jp = entities%v_tags(jcol,iel)
               ind_jp = m_eqs(j_eq,jp)
               if (ind_jp .gt. 0) then
                  col_start = col_ptr(ind_jp) + i_base2
                  col_end = col_ptr(ind_jp+1) - 1 + i_base2
                  !  unpack row into i_work
                  do i=col_start,col_end
                     i_work(row_ind(i) + i_base2) = i
                  enddo


                  call find_basis_derivatives(j_eq, jcol, phi_vec_map, phi_P2_ref, phi_P3_ref, &
                     gradt_P1_act, gradt_P2_act, gradt_P3_act, vec_phi_j, curlt_phi_j, gradt_j, phi_z_j)

                  ! if (jcol .le. N_DDL_T) then ! A transverse dof (edge or face)
                  !    ! Uses P2 vector elements so determine the basis vector
                  !    call make_phi_vector_basis(j_eq, jcol, phi_vec_map, phi_P2_ref, &
                  !       gradt_P1_act, gradt_P2_act, vec_phi_j, curlt_phi_j)
                  !    gradt_j = D_ZERO
                  !    phi_z_j = D_ZERO
                  ! else   ! a longitudinal dof, use P3 scalar element
                  !    vec_phi_j = D_ZERO
                  !    curlt_phi_j = D_ZERO
                  !    gradt_j(:) = gradt_P3_act(:,jcol-N_DDL_T)
                  !    phi_z_j = phi_P3_ref(jcol-N_DDL_T)
                  ! endif

                  do irow=1,nddl_0
                     z_phase_fact = val_exp(jcol) * conjg(val_exp(irow))
                     ip = entities%v_tags(irow,iel)
                     do i_eq=1,3
                        ind_ip = m_eqs(i_eq,ip)
                        if (ind_ip .gt. 0) then

                           !if (ind_jp .eq. ind_ip .and. &
                           !   abs(imag(z_phase_fact)) .gt. 1.0d-15) then
                           !   write(ui_stdout,*) "phase_fact: ", ind_jp, ind_ip, &
                           !      z_phase_fact, val_exp(jcol), val_exp(irow)
                           !endif

                           call find_basis_derivatives(i_eq, irow, phi_vec_map, phi_P2_ref, phi_P3_ref, &
                              gradt_P1_act, gradt_P2_act, gradt_P3_act, vec_phi_i, curlt_phi_i, gradt_i, phi_z_i)

                           ! if (irow .le. N_DDL_T) then  !  edge or face element
                           !    call make_phi_vector_basis(i_eq, irow, phi_vec_map, &
                           !       phi_P2_ref, gradt_P1_act, gradt_P2_act, vec_phi_i, &
                           !       curlt_phi_i)
                           !    gradt_i = D_ZERO
                           !    phi_z_i = D_ZERO
                           ! else
                           !    vec_phi_i = D_ZERO
                           !    curlt_phi_i = D_ZERO
                           !    gradt_i(:) = gradt_P3_act(:,irow-N_DDL_T)
                           !    phi_z_i = phi_P3_ref(irow-N_DDL_T)
                           ! endif

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
                           if (irow .le. N_DDL_T .and. jcol .le. N_DDL_T) then    ! [tt] part
                              ! K_tt =   (curl E_j).(curl F_i)-k^2 eps (E_j, F_i)
                              ! M_tt = - (E_j, F_i)
                              r_tmp1 = curlt_phi_j * curlt_phi_i
                              r_tmp2 = ddot(2, vec_phi_j, 1, vec_phi_i, 1)
                              K_tt = r_tmp1 * tperm_pp - r_tmp2 * tperm_qq
                              M_tt = - r_tmp2 * tperm_pp
                              K_elt = K_tt
                              M_elt = M_tt

                           elseif (irow .le. N_DDL_T .and. jcol .gt. N_DDL_T) then  ! [tz] part
                              ! K_tz =  0
                              ! M_tz = (gradt e_j, F_i)
                              r_tmp1 = ddot(2, gradt_j, 1, vec_phi_i, 1)
                              K_tz = 0.0d0
                              M_tz = r_tmp1 * tperm_pp
                              K_elt = K_tz
                              M_elt = M_tz

                           elseif (irow .gt. N_DDL_T .and. jcol .le. N_DDL_T) then ! [zt] part
                              ! K_tz =  (E_j, gradt f_i)
                              ! M_tz =  0
                              r_tmp1 = ddot(2, vec_phi_j, 1, gradt_i, 1)
                              K_zt = r_tmp1 * tperm_pp
                              M_zt = 0.0d0
                              K_elt = K_zt
                              M_elt = M_zt

                           elseif (irow .gt. N_DDL_T .and. jcol .gt. N_DDL_T) then ! [zz] part
                              ! K_zz =  - (gradt e_j, gradt f_i) + eps_r k^2 (e_j, f_i )
                              ! M_zz =  0
                              r_tmp1 = ddot(2, gradt_j, 1, gradt_i, 1)
                              r_tmp2 = phi_z_j * phi_z_i
                              K_zz = - r_tmp1 * tperm_pp + r_tmp2 * tperm_qq
                              M_zz = 0.0d0
                              K_elt = K_zz
                              M_elt = M_zz

                           else
                              errco = NBERR_BAD_ASSEMBLY
                              write(emsg,*) "irow or jcol has an ", "invalid value",  &
                                 "irow jcol, = ", irow, jcol
                              return
                           endif

                           ! without periodic bcs, z_phase_fact = 0 always
                           ! add term to sum weighted by quadrature factor and area det
                           K_elt = K_elt * ww * abs(det) * z_phase_fact
                           M_elt = M_elt * ww * abs(det) * z_phase_fact

                           ! TODO:
                           !  This is an arpack solver thing.
                           !   No reason to do this deep in loop.
                           !  Should do once the matrices have been constructed
                           !  And closer to the solver
                           !
                           K_elt = K_elt - shift_ksqr*M_elt

                           k = i_work(ind_ip)
                           if (k .gt. 0 .and. k .le. nonz) then   !is this test necessary?
                              mOp_stiff(k) = mOp_stiff(k) + K_elt
                              mOp_mass(k) = mOp_mass(k) + M_elt
                           else
                              errco = NBERR_BAD_ASSEMBLY
                              write(emsg,*) "asmbly: problem with row_ind: k, nonz = ", k, nonz
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

   if (debug .eq. 1) then
      write(ui_stdout,*) "asmbly: shift_ksqr = ", shift_ksqr
      write(ui_stdout,*) "asmbly: number of curved elements = ", n_curved
      write(ui_stdout,*) "asmbly: n_msh_el, (n_msh_el-n_curved) = ", n_msh_el, &
         (n_msh_el-n_curved)
   endif

   if (debug .eq. 1) then
      write(ui_stdout,*)
      write(ui_stdout,*) "  Re perm_pp = ", dble(perm_pp)
      write(ui_stdout,*) "imag perm_pp = ", imag(perm_pp)
      write(ui_stdout,*)
      write(ui_stdout,*) "  Re perm_qq = ", dble(perm_qq)
      write(ui_stdout,*) "imag perm_qq = ", imag(perm_qq)
   endif


   return
end
