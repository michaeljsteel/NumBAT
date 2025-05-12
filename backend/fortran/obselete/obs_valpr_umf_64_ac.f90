#include "numbat_decl.h"


subroutine valpr_64_ac (i_base, dim_krylov, n_modes, itermax, &
   arp_tol, cscmat, v_evals_nu, v_evecs_arp, nberr)

!  ------------------------------------------------------------------

   use numbatmod
   use alloc
   use class_SparseCSC_AC
   use class_ValprVecs



   integer(8) i_base, dim_krylov, n_modes

   type(SparseCSC_AC) :: cscmat

   integer(8) itermax
   double precision arp_tol

   complex(8) v_evals_nu(n_modes),  v_evecs_arp(cscmat%n_dof,n_modes)

   type(NBError) nberr

   !---------------------------------------------

   type(ValprVecs) vecs

   double precision umf_control (20), umf_info (90)
   integer(8) umf_numeric, symbolic
   !, sys

   integer(8) i, j, n_dof, n_nonz


   complex(8) arp_shift

   character(len=EMSG_LENGTH) emsg


   double precision, allocatable, dimension(:) :: mOp_stiff_re, mOp_stiff_im

   integer(8) n_conv

   !  32-bit integers for ARPACK
   integer(4) n_modes_32, dim_krylov_32, n_dof_32
   integer(4) arp_ido, arp_info, arp_iparam(11)
   integer(4) ipntr_32(14), lworkl_32

   logical arp_rvec
   character arp_bmat
   character(2) arp_which
   logical arp_active


   integer(8) ui

   !  ------------------------------------------------------------------

   ui = stdout
   emsg = ""

   if (i_base .ne. 0) then
      write(emsg,*) "valpr_64: i_base != 0 : ", i_base,&
      &"valpr_64: UMFPACK requires 0-based indexing"
      call nberr%set(-102_8, emsg);
   endif

   n_dof = cscmat%n_dof
   n_nonz = cscmat%n_nonz

   call vecs%init(n_modes, dim_krylov, n_dof, nberr); RET_ON_NBERR(nberr)

   call double_alloc_1d(mOp_stiff_re, n_nonz, 'mOp_stiff_re', nberr); RET_ON_NBERR(nberr)
   call double_alloc_1d(mOp_stiff_im, n_nonz, 'mOp_stiff_im', nberr); RET_ON_NBERR(nberr)

   mOp_stiff_re = dble(cscmat%mOp_stiff)
   mOp_stiff_im = dimag(cscmat%mOp_stiff)


!  ----------------------------------------------------------------
!  factor the matrix as A = LU and save to a file
!  ----------------------------------------------------------------



!  set default parameters

   !  umfpack * report status (print level = umf_control(1)) :
   !  print level = 0 or less : No output, even when an error occurs.
   !  print level = 1 (default value) : then error messages are printed,
   !  and nothing is printed if the status is UMFPACK OK.
   !  print level = 2 or more : then the status is always printed.

   call umf4zdef (umf_control)
   umf_control (1) = 1
   call umf4zpcon (umf_control)

!  pre-order and symbolic analysis
   call umf4zsym (n_dof, n_dof, cscmat%v_col_ptr, cscmat%v_row_ind, &
      mOp_stiff_re, mOp_stiff_im, symbolic, umf_control, umf_info)

   if (umf_info (1) .lt. 0) then
      write(emsg,'(A,i4)') 'Error occurred in sparse matrix symbolic factorization umf4zsym:', int(umf_info (1))
      call nberr%set(NBERR_BAD_UMF4ZSYM, emsg)
      return
   endif


   !  Complete the umf_numeric factorization
   !  TODO: This call does not appear to be thread safe!  Breaks tutorial 3 B in thread mode
   call umf4znum (cscmat%v_col_ptr, cscmat%v_row_ind, mOp_stiff_re, mOp_stiff_im, symbolic, umf_numeric, umf_control, umf_info)

   if (umf_info (1) .lt. 0) then
      write(emsg,*) 'Error occurred in sparse matrix umf_numeric factorization umf4znum: ', umf_info (1)
      call nberr%set(NBERR_BAD_UMF4ZNUM, emsg)
      return
   endif

   !  free the symbolic analysis
   call umf4zfsym (symbolic)

   v_evals_nu = D_ZERO



!  ##################################################################

!  On commence le workail avec znaupd
!  ----------------------------------


   !  The arp_ido parameter is used for reverse communication.
   !  Initially, it should be set to 0.

   !  Setting arp_info to 0 instructs ARPACK to construct an initial vector with random components.
   !  Setting arp_iparam(1) to 1 indicates that ARPACK should calculate translations
   !  based on the projected matrix and according to the "arp_which criterion.


   n_dof_32 = int(n_dof, 4)
   n_modes_32 = int(n_modes, 4)
   dim_krylov_32 = int(dim_krylov, 4)
   lworkl_32 = int(vecs%lworkl, 4)

   arp_ido = 0
   arp_iparam(1) = 1
   arp_iparam(3) = int(itermax, 4)
!  arp_iparam(7) = 3
   arp_iparam(7) = 1
   arp_info = 0


   !----------------------------------------------------
   !  Main loop in inverse communication mode
   !----------------------------------------------------
   arp_bmat = 'I'    !  plain (not generalised) eigenvalue problem
   arp_which = 'LM'  !  seek largest magnitude eigs

   arp_active = .true.
   arp_shift = C_ONE   !  Is ignored, as long as iparam(7)=1


   !  Test for dimesnion conditions in znaupd (err code = -3)
   !  Test for N=n_dof_32, NEV=n_modes_32, NCV=dim_krylov_32
   !  Need 0<n_modes_32<n_dof_32-1, 1<= dim_krylov_32-n_modes_32, dim_krylov_32<=n_dof_32
   if ((n_dof_32-1 .le. n_modes_32) .or. (dim_krylov_32-n_modes_32 .lt. 1) &
      .or.  dim_krylov_32 > n_dof_32) then
      write(emsg,'(A,A)') 'ARPACK eigensolver dimensional'//&
         ' conditions failed (would generate ARPACK znaupd error code of -3).' // NEW_LINE('A'),&
         'You should probably increase the grid resolution.'
      call nberr%set(NBERR_BAD_ZNAUPD, emsg);
      return
   endif

   do while (arp_active)

      call znaupd (arp_ido, arp_bmat, n_dof_32, arp_which, n_modes_32, arp_tol,&
         vecs%resid, dim_krylov_32, vecs%v_schur, n_dof_32, arp_iparam,&
         ipntr_32, vecs%workd, vecs%workl, lworkl_32, vecs%rwork, arp_info)

      if (arp_ido .eq. -1 .or. arp_ido .eq. 1) then

         !------------------------------------------------------
         !  Apply  y <--- OP*x = inv[A-SIGMA*M]*M*x
         !  with x at x = vecs%workd(ipntr_32(1))
         !  and place the result at  y = vecs%workd(ipntr_32(2))
         !------------------------------------------------------

         call zcopy(n_dof_32, vecs%workd(ipntr_32(1)), 1,vecs%vect1_z, 1)
         call z_mxv_csc (n_dof, vecs%vect1_z, vecs%vect2_z, cscmat%n_nonz, cscmat%v_row_ind,&
            cscmat%v_col_ptr, cscmat%mOp_mass)

         call vecs%vect2_complex_to_real()

         ! solve Ax=b, without iterative refinement (UMFPACK_A)
         ! x = inv[K-sigma M] * vect2
         !sys = 0
         call umf4zsol (UMFPACK_A,  vecs%vect1_re, vecs%vect1_im, vecs%vect2_re, vecs%vect2_im, &
            umf_numeric, umf_control, umf_info)

         if (umf_info (1) .lt. 0) then
            write(emsg,*) 'Error occurred in umf4zsol: ', umf_info (1)
            call nberr%set(NBERR_BAD_UMF4ZSOL, emsg)
            return
         endif

         call vecs%vect1_real_to_complex()

         call zcopy(n_dof_32, vecs%vect1_z, 1, vecs%workd(ipntr_32(2)), 1)

      else if (arp_ido .eq. 2) then  ! never happens for numbat

         write(ui,*) 'VALPR_64: ATTENTION arp_ido = ', arp_ido
         write(ui,*) 'check the results...'

         !  ----------------------------------------------
         !  | Apply  y <--- M*x                       |
         !  | x = vecs%workd(ipntr_32(1))  et  y = vecs%workd(ipntr_32(2)) |
         !  ----------------------------------------------

         call zcopy(n_dof_32, vecs%workd(ipntr_32(1)), 1, vecs%vect1_z, 1)
         call z_mxv_csc (n_dof, vecs%vect1_z, vecs%vect2_z, cscmat%n_nonz, cscmat%v_row_ind,&
            cscmat%v_col_ptr, cscmat%mOp_mass)


         ! rhs_re = dble(vecs%vect2_z)
         ! rhs_im = imag(vecs%vect2_z)

         ! !  solve Ax=b, without iterative refinement
         ! sys = 0
         ! call umf4zsol (sys, lhs_re, lhs_im, rhs_re, rhs_im,&
         !    umf_numeric, umf_control, umf_info)

         ! if (umf_info (1) .lt. 0) then
         !    write(emsg, *) 'Error occurred in umf4zsol: ', umf_info (1)
         !    call nberr%set(NBERR_BAD_UMF4ZSOL, emsg);
         !    return
         ! endif

         ! vecs%vect2_z  = dcmplx (lhs_re, lhs_im )

         call zcopy(n_dof_32, vecs%vect2_z,1, vecs%workd(ipntr_32(2)), 1)

      else
         arp_active = .false.
      end if

   enddo

   call umf4zfnum (umf_numeric) !  free the umf_numeric factorization

   !  ---------------------------------------------------
   !  | Either we have convergence, or there is an error. |
   !  ---------------------------------------------------

   n_conv = arp_iparam(5)

   call check_arpack_results(n_modes,  n_conv, arp_info, nberr)
   RET_ON_NBERR(nberr)



   !  Get the final eigenvectors
   !'A' means get the actual eigenvectors, not just schur/arnolid vectors

   arp_rvec = .true.
   arp_shift = C_ZERO  !  Is ignored, as long as iparam(7)=1

   ! Output goes into eval_ritz and evecs_arp

   call zneupd (arp_rvec, 'A', vecs%arp_select, &
      vecs%eval_ritz, v_evecs_arp, n_dof_32, arp_shift,&
      vecs%workev, arp_bmat, n_dof_32, arp_which, n_modes_32, arp_tol,&
      vecs%resid, dim_krylov_32, vecs%v_schur, n_dof_32, arp_iparam, ipntr_32,&
      vecs%workd, vecs%workl, lworkl_32, vecs%rwork, arp_info)

   if (arp_info .ne. 0) then
      write(emsg,*) 'VALPR_64: Error with _neupd, arp_info = ', arp_info, &
         ' This error can occur if the mesh is too coarse.'
      call nberr%set(NBERROR_109, emsg)
      return
   endif


   v_evals_nu = vecs%eval_ritz(1:n_modes) ! eval_ritz is one longer due to zneupd requirements

   ! do i = 1, n_modes
   !    do j = 1, n_dof
   !       v_evecs_arp(j,i) = vecs%v_schur(j,i)
   !    enddo
   ! enddo


end

