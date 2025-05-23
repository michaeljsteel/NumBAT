 !  Calculate the Overlap integral of the prime and adjoint Plane Waves

 !  TODO: this does much more than calculate the overlap. It is building fields
 !  Needs understanding and renaming

subroutine orthogonal (n_modes, n_msh_elts, n_msh_pts, &
   nnodes, nb_typ_el, pp, m_elnd_to_mshpt, &
   type_el, x, beta1, soln_k1, &
   mat_overlap, overlap_file, PrintAll, &
   pair_warning, k_0)

   use numbatmod

   integer(8) :: n_modes
   integer(8) n_msh_elts, n_msh_pts, nnodes, nb_typ_el
   integer(8) type_el(n_msh_elts)
   integer(8) m_elnd_to_mshpt(nnodes,n_msh_elts)
   double precision x(2,n_msh_pts)
   complex(8) soln_k1(3,nnodes+7,n_modes,n_msh_elts)
   complex(8) soln_k2(3,nnodes+7,n_modes,n_msh_elts)
   complex(8) pp(nb_typ_el)
   complex(8) beta1(n_modes)
   complex(8) beta2(n_modes)
   integer(8) errco
   character(len=EMSG_LENGTH) emsg

   !  complex(8) mat_overlap(n_modes,n_modes)
   complex(8), dimension(n_modes,n_modes) :: mat_overlap
   character overlap_file*100
   double precision k_0
   !  Local variables

   integer(8) nod_el_p(P2_NODES_PER_EL)
   complex(8) sol_el_1(2*P2_NODES_PER_EL+10) , sol_el_2(2*P2_NODES_PER_EL)
   complex(8) vec_1(2*P2_NODES_PER_EL)
   complex(8) mat_scal(2*P2_NODES_PER_EL,2*P2_NODES_PER_EL+10)
   integer(8) i, j, j1, typ_e
   integer(8) iel, ival, jval
   integer(8) jtest, ind_jp, j_eq
   integer(8) itrial, ind_ip, i_eq
   integer(8) n_curved, debug, ui
   logical is_curved
   double precision xel(2,P2_NODES_PER_EL)
   double precision phi2_list(6), grad2_mat0(2,6)
   double precision grad2_mat(2,6)
   double precision phi3_list(10), grad3_mat0(2,10)
   double precision grad3_mat(2,10)
   double precision vec_phi_j(2), vec_phi_i(2)
   double precision  r_tmp1
   complex(8) z_tmp1, z_tmp2, z_beta_1, coeff_1

   !  NQUAD: The number of quadrature points used in each element.
   integer(8) nquad, nquad_max, iq
   parameter (nquad_max = 25)
   double precision wq(nquad_max)
   double precision xq(nquad_max), yq(nquad_max)
   double precision xx(2), xx_g(2), ww, det
   double precision mat_B(2,2), mat_T(2,2)
   !  Mode ordering
   integer(8) skip, PrintAll
   logical pair_warning
   complex(8) betatmp1(1), betatmp2(1)
   complex(8) soltmp1(3,nnodes+7,n_msh_elts,1)
   complex(8) soltmp2(3,nnodes+7,n_msh_elts,1)
   integer(8) compcount, elcount, nodecount, redo, j2
   double precision val_max_diag, val_max_off


   !  Start Program  !!!!!!!!!!!!!!!!!!!!!!!!

   ui = stdout
   debug = 0
   pair_warning = .false.

   call quad_triangle (nquad, nquad_max, wq, xq, yq)
   if (debug .eq. 1) then
      write(ui,*) "orthogonal: nquad, nquad_max = ", nquad, nquad_max
   endif

   redo = 0
   !second rearranged overlap
122 continue

   !do jval=1,n_modes
   !  do ival=1,n_modes
   !  mat_overlap(ival,jval) = 0.0d0
   !  enddo
   !enddo
   mat_overlap  = C_ZERO

   n_curved = 0

   do iel=1,n_msh_elts
      typ_e = type_el(iel)
      do j=1,nnodes
         j1 = m_elnd_to_mshpt(j,iel)
         nod_el_p(j) = j1
         xel(1,j) = x(1,j1)
         xel(2,j) = x(2,j1)
      enddo
      is_curved =  log_is_curved_elem_tri (nnodes, xel)
      if (is_curved) then
         n_curved = n_curved + 1
      endif

      !  do i=1,2*nnodes
      !  do j=1,2*nnodes+10
      !  mat_scal(i,j) = 0.0d0
      !  enddo
      !  enddo
      mat_scal = C_ZERO

      do iq=1,nquad
         xx(1) = xq(iq)
         xx(2) = yq(iq)
         ww = wq(iq)
         !  xx   = coordinate on the reference triangle
         !  xx_g = coordinate on the actual triangle
         !  We will also need the gradients of the P1 element
         !  grad2_mat0 = gradient on the reference triangle (P2 element)
         call phi2_2d_mat(xx, phi2_list, grad2_mat0)
         !  grad3_mat0 = gradient on the reference triangle (P3 element)
         call phi3_2d_mat(xx, phi3_list, grad3_mat0)

         if (.not. is_curved ) then
            !  Rectilinear element
            call jacobian_p1_2d(xx, xel, nnodes, xx_g, det, mat_B, mat_T, errco, emsg)
            !  if (det .le. 0) then
            if (det .le. 0 .and. debug .eq. 1 .and. iq .eq. 1) then
               write(*,*) "   !!!"
               write(*,*) "orthogonal: det <= 0: iel, det ", iel, det
            endif
         else
            !  Isoparametric element
            !  24/6/12 Deleting first broken argument xx:
            !  p2_2d is diff to p1_2d.
            call jacobian_p2_2d(xel, nnodes, phi2_list, grad2_mat0, xx_g, det, mat_B, mat_T, errco, emsg)
         endif
         !  if(abs(det) .lt. 1.0d-10) then
         if(abs(det) .lt. 1.0d-20) then
            write(*,*)
            write(*,*) "   ???"
            write(*,*) "orthogonal: det = 0 : iel, det = ", iel, det
            write(*,*) "orthogonal: Aborting..."
            stop
         endif
         !  grad_i  = gradient on the actual triangle
         !  grad_i  = Transpose(mat_T)*grad_i0
         !  Calculation of the matrix-matrix product:
         call DGEMM('Transpose','N', 2, 6, 2, D_ONE, mat_T, 2, &
            grad2_mat0, 2, D_ZERO, grad2_mat, 2)
         call DGEMM('Transpose','N', 2, 10, 2, D_ONE, mat_T, 2, &
            grad3_mat0, 2, D_ZERO, grad3_mat, 2)
         coeff_1 = ww * abs(det) * pp(typ_e)
         do itrial=1,P2_NODES_PER_EL
            do i_eq=1,2
               ind_ip = i_eq + 2*(itrial-1)
               !  Determine the basis vector
               do i=1,2
                  vec_phi_i(i) = 0.0d0
               enddo
               vec_phi_i(i_eq) = phi2_list(itrial)
               do jtest=1,P2_NODES_PER_EL
                  do j_eq=1,2
                     ind_jp = j_eq + 2*(jtest-1)
                     !  Determine the basis vector
                     do i=1,2
                        vec_phi_j(i) = 0.0d0
                     enddo
                     vec_phi_j(j_eq) = phi2_list(jtest)
                     !  z_tmp1 = ddot(2, vec_phi_i, 1, vec_phi_j, 1)
                     z_tmp1 = vec_phi_i(1)*vec_phi_j(1) + vec_phi_i(2)*vec_phi_j(2)
                     z_tmp1 = coeff_1 * z_tmp1
                     z_tmp1 = z_tmp1/k_0
                     mat_scal(ind_ip,ind_jp) = mat_scal(ind_ip,ind_jp) + z_tmp1
                  enddo
               enddo
               do jtest=1,10
                  j_eq = 3
                  ind_jp = jtest + 2*P2_NODES_PER_EL
                  !  Determine the basis vector
                  do i=1,2
                     vec_phi_j(i) = -grad3_mat(i,jtest)
                  enddo
                  !  z_tmp1 = ddot(2, vec_phi_i, 1, vec_phi_j, 1)
                  z_tmp1 = vec_phi_i(1)*vec_phi_j(1) + vec_phi_i(2)*vec_phi_j(2)
                  z_tmp1 = coeff_1 * z_tmp1
                  z_tmp1 = z_tmp1/k_0
                  mat_scal(ind_ip,ind_jp) = mat_scal(ind_ip,ind_jp) + z_tmp1
               enddo
            enddo
         enddo
      enddo

      do ival=1,n_modes
         !  do i=1,nnodes
         !  do j=1,2
         !  The 2 transverse components of the mode ival
         !  ind_ip = j + 2*(i-1)
         !  z_tmp1 = soln_k2(j,i,ival,iel)
         !  sol_el_2(ind_ip) = z_tmp1
         !  enddo
         !  enddo

         do jval=1,n_modes
            z_beta_1 = beta1(jval)
            do i=1,nnodes
               do j=1,2
                  !  The 2 transverse components of the mode jval
                  ind_jp = j + 2*(i-1)
                  z_tmp1 = soln_k1(j,i,jval,iel)
                  sol_el_1(ind_jp) = z_tmp1 * z_beta_1
               enddo
            enddo

            do i=1,3
               !  The longitudinal component at the vertices (P3 elements)
               ind_jp = i + 2*nnodes
               z_tmp1 = soln_k1(3,i,jval,iel)
               sol_el_1(ind_jp) = z_tmp1 * z_beta_1
            enddo
            do i=nnodes+1,13
               !  The longitudinal component at the edge nodes and interior node (P3 elements)
               ind_jp = i + 2*nnodes - nnodes + 3
               z_tmp1 = soln_k1(3,i,jval,iel)
               sol_el_1(ind_jp) = z_tmp1 * z_beta_1
            enddo

            !  Matrix-Vector product
            do i=1,2*nnodes
               vec_1(i) = 0.0d0
               do j=1,2*nnodes+10
                  z_tmp1 = sol_el_1(j)
                  z_tmp2 = mat_scal(i,j)
                  vec_1(i) = vec_1(i) + z_tmp1 * z_tmp2
               enddo
            enddo

            !  Scalar product
            z_tmp1 = 0.0d0
            do i=1,2*nnodes
               z_tmp1 = vec_1(i) * sol_el_2(i)
               mat_overlap(ival,jval) = mat_overlap(ival,jval) + z_tmp1
            enddo
         enddo
      enddo
   enddo



   !  reorder complex conjugate douplets
   j = 1
   skip = 0
   redo = redo + 1
   if (redo .eq. 2) goto 123

121 continue
   if (j .gt. n_modes) then
      goto 122
      !  if all is well - correct orthogonality with its self
   elseif (abs(mat_overlap(j,j)) .gt. 1.0d-7) then
      j = j+1
      goto 121
      !  first of a wrongly ordered complex conjugate pair (save values)
   elseif (skip .eq. 0) then
      if (j .eq. n_modes) then
         pair_warning = .true.
      endif
      !  find jvals (j and j2) of swaped pair and switch
      do j2 = j+1,n_modes
         if (abs(mat_overlap(j2,j)) .gt. 1.0d-7) then
            betatmp1(1) = beta2(j)
            betatmp2(1) = beta2(j2)
            beta2(j)  = betatmp2(1)
            beta2(j2) = betatmp1(1)
            do compcount = 1,3
               do elcount = 1,n_msh_elts
                  do nodecount = 1,nnodes+7
                     soltmp1(compcount,nodecount,elcount,1) = soln_k2(compcount,nodecount,j,elcount)
                     soltmp2(compcount,nodecount,elcount,1) = soln_k2(compcount,nodecount,j2,elcount)
                     soln_k2(compcount,nodecount,j,elcount) = soltmp2(compcount,nodecount,elcount,1)
                     soln_k2(compcount,nodecount,j2,elcount) = soltmp1(compcount,nodecount,elcount,1)
                  enddo
               enddo
            enddo
         endif
      enddo
      skip = 1
      j = j+1
      goto 121
      !  dont touch second half of pair (as already switched)
   elseif (skip .gt. 0) then
      skip = 0
      j = j+1
      goto 121
   endif



123 continue
   if (PrintAll .eq. 1) then
      open(3,file=overlap_file)
      val_max_diag = 0.0d0
      val_max_off = 0.0d0
      do jval=1,n_modes
         do ival=1,n_modes
            r_tmp1 = abs(mat_overlap(ival,jval))
            if (ival .eq. jval) then
               if (val_max_diag .lt. r_tmp1) then
                  val_max_diag = r_tmp1
               endif
            else
               if (val_max_off .lt. r_tmp1) then
                  val_max_off = r_tmp1
               endif
            endif
            write(3,12) ival, jval, mat_overlap(ival,jval),&
               abs(mat_overlap(ival,jval)),&
               abs(beta1(ival)-beta2(jval)), beta1(ival), beta2(jval),&
               abs(mat_overlap(ival,jval) -&
               conjg(mat_overlap(jval,ival)))
         enddo
      enddo
      write(3,*) "val_max_diag, val_max_off = ", val_max_diag, val_max_off
      close(3)
12    format(2(I4),2(g25.17), 2(g16.8), 8(g18.10))
   endif


   return
end
