
subroutine asmbly_AC (i_base, nel, npt, n_dof, nnodes, shift, beta, nb_typ_el, rho, c_tensor, &
   elnd_to_mshpt, type_el, in_dof, x, nonz, row_ind, col_ptr, &
   mat1_re, mat1_im, mat2, i_work, symmetry_flag, debug)


   use numbatmod

   integer(8) nnodes
   integer(8) nel, npt, n_dof
   integer(8) i_base, nb_typ_el, nonz
   integer(8) row_ind(nonz), col_ptr(n_dof+1)
   integer(8) type_el(nel)
   integer(8) elnd_to_mshpt(nnodes,nel), in_dof(3,npt)
   integer(8) i_work(3*npt), symmetry_flag
   complex(8) shift, beta
   double precision x(2,npt)

   complex(8) rho(nb_typ_el), c_tensor(6,6,nb_typ_el)
   complex(8) mat2(nonz)
   double precision mat1_re(nonz), mat1_im(nonz)


   integer(8) nod_el(P2_NODES_PER_EL)
   double precision xel(2,P2_NODES_PER_EL)

   integer(8) i_base2, debug, typ_e
   integer(8) i, j, k, iel, j1
   integer(8) jtest, jp, ind_jp, j_eq
   integer(8) itrial, ip, ind_ip, i_eq
   integer(8) col_start, col_end

   complex(8) mat_K(18,18), mat_M(18,18)
   complex(8) c_tensor_el(6,6), rho_el
   complex(8) z_tmp1, z_tmp2

!  The CSC indexing, i.e., col_ptr, is 1-based
!  But valpr.f may have changed the CSC indexing to 0-based indexing)
   debug = 0

   if (i_base .eq. 0) then
      i_base2 = 1
   else
      i_base2 = 0
   endif

   if (nnodes .ne. 6 ) then
      write(*,*) "asmbly_AC: problem nnodes = ", nnodes
      write(*,*) "asmbly_AC: nnodes should be equal to 6 !"
      write(*,*) "asmbly_AC: Aborting..."
      stop
   endif


   mat1_re  = D_ZERO
   mat1_im  = D_ZERO
   mat2  = D_ZERO

   do iel=1,nel


      typ_e = type_el(iel)
      do j=1,nnodes
         j1 = elnd_to_mshpt(j,iel)
         nod_el(j) = j1
         xel(:,j) = x(:,j1)
      enddo

      rho_el = rho(typ_e)

      !  do j=1,6
      !  do i=1,6
      !  c_tensor_el(i,j) = c_tensor(i,j,typ_e)
      !  enddo
      !  enddo

      c_tensor_el = c_tensor(:,:,typ_e)

!  If c_tensor has regular symmetries use more efficient formulation

      if (symmetry_flag .eq. 1) then
         call mat_el_v2 (xel,beta,c_tensor_el,rho_el,mat_K,mat_M)
      elseif (symmetry_flag .eq. 0) then
         call mat_el_v3 (xel,beta,c_tensor_el,rho_el,mat_K,mat_M)
      else
         write(*,*) "asmbly_AC: problem with symmetry_flag "
         write(*,*) "asmbly_AC: c_tensor = ", symmetry_flag
         write(*,*) "asmbly_AC: Aborting..."
         stop
      endif

      do jtest=1,nnodes
         jp = elnd_to_mshpt(jtest,iel)
         do j_eq=1,3
            ind_jp = in_dof(j_eq,jp)
            if (ind_jp .gt. 0) then
               col_start = col_ptr(ind_jp) + i_base2
               col_end = col_ptr(ind_jp+1) - 1 + i_base2

!  unpack row into i_work
               do i=col_start,col_end
                  i_work(row_ind(i) + i_base2) = i
               enddo

               do itrial=1,nnodes
                  do i_eq=1,3
                     ip = elnd_to_mshpt(itrial,iel)
                     ind_ip = in_dof(i_eq,ip)
                     if (ind_ip .gt. 0) then
                        z_tmp1 = mat_K(3*(itrial-1) + i_eq, 3*(jtest-1) + j_eq)
                        z_tmp2 = mat_M(3*(itrial-1) + i_eq, 3*(jtest-1) + j_eq)
                        z_tmp1 = z_tmp1 - shift*z_tmp2
                        k = i_work(ind_ip)
                        if (k .gt. 0 .and. k .le. nonz) then
                           mat1_re(k) = mat1_re(k) + real(z_tmp1)
                           mat1_im(k) = mat1_im(k) + imag(z_tmp1)
                           mat2(k) = mat2(k) + z_tmp2
                        else
                           write(*,*) "asmbly_AC: problem with row_ind !!"
                           write(*,*) "asmbly_AC: k, nonz = ", k, nonz
                           write(*,*) "asmbly_AC: Aborting..."
                           stop
                        endif
                     endif
                  enddo
               enddo
            endif
         enddo
      enddo
   enddo

   return
end
