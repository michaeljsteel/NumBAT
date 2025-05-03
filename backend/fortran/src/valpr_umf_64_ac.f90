#include "numbat_decl.h"


subroutine valpr_64_ac (i_base, dim_krylov, n_modes,  itermax,&
   arp_tol, cscmat, v_evals_nu, v_evecs_arp, nberr)

!  ------------------------------------------------------------------

   use numbatmod
   use alloc
   use class_SparseCSC_AC

   use class_ValprVecs

   type(NBError) nberr

   type(SparseCSC_AC) :: cscmat

   integer(8) :: n_modes
   integer(8) n_conv, i_base, dim_krylov, lwork


   !complex(8) mat2(cscmat%n_nonz)


   complex(8) v_evals_nu(n_modes),  v_evecs_arp(cscmat%n_dof,n_modes)

   double precision time1_fact, time2_fact

   double precision umf_control (20), umf_info (90)
   integer(8) umf_numeric, symbolic, sys

   integer(8) itermax, i, j, n_dof, n_nonz


   type(ValprVecs) :: vecs

   double precision arp_tol
   complex(8) arp_shift

   integer(8) errco
   character(len=EMSG_LENGTH) emsg

   integer(8) alloc_stat
   !complex(8), dimension(:), allocatable :: workev
   !double precision, dimension(:), allocatable :: rwork
   !logical, dimension(:), allocatable :: selecto


   !complex(8), dimension(:), allocatable :: vect1, vect2, workd, resid, work
   !complex(8), dimension(:,:), allocatable :: vschur


   double precision, allocatable, dimension(:) :: lhs_re, lhs_im
   double precision, allocatable, dimension(:) :: rhs_re, rhs_im
   double precision, allocatable, dimension(:) :: mOp_stiff_re, mOp_stiff_im


   !  Local variables
   !  32-bit integers for ARPACK
   integer(4) n_modes_32, dim_krylov_32, n_dof_32
   integer(4) arp_ido, arp_info, ierr_32, arp_iparam(11)
   integer(4) ipntr_32(14), lworkl_32

   logical arp_rvec
   character arp_bmat
   character(2) arp_which
   logical arp_active

   integer(8) ui, debug, show_mem_est

   !  ------------------------------------------------------------------

   debug = 0
   show_mem_est = 0

   ui = stdout
   errco = 0
   emsg = ""

   if (i_base .ne. 0) then
      write(emsg,*) "valpr_64: i_base != 0 : ", i_base,&
      &"valpr_64: UMFPACK requires 0-based indexing"
      errco = -102
      call nberr%set(errco, emsg);
   endif

   lwork = 3*dim_krylov*(dim_krylov+2)

   n_dof = cscmat%n_dof
   n_nonz = cscmat%n_nonz

   call vecs%init(n_modes, dim_krylov, n_dof, nberr); RET_ON_NBERR(nberr)

   !call complex_nalloc_1d(vect1, n_dof, 'vect1_ac', nberr); RET_ON_NBERR(nberr)
   !call complex_nalloc_1d(vect2, n_dof, 'vect1_ac', nberr); RET_ON_NBERR(nberr)
   !call complex_nalloc_1d(workd, 3*n_dof, 'vect1_ac', nberr); RET_ON_NBERR(nberr)
   !call complex_nalloc_1d(work, lwork, 'vect1_ac', nberr); RET_ON_NBERR(nberr)

   !call complex_nalloc_1d(resid, n_dof, 'vect1_ac', nberr); RET_ON_NBERR(nberr)

   !call complex_nalloc_2d(vschur, n_dof, dim_krylov, 'vect1_ac', nberr); RET_ON_NBERR(nberr)


   call double_nalloc_1d(mOp_stiff_re, n_nonz, 'mOp_stiff_re', nberr); RET_ON_NBERR(nberr)
   call double_nalloc_1d(mOp_stiff_im, n_nonz, 'mOp_stiff_im', nberr); RET_ON_NBERR(nberr)

   call double_nalloc_1d(lhs_re, n_dof, 'lhs_re', nberr); RET_ON_NBERR(nberr)
   call double_nalloc_1d(lhs_im, n_dof, 'lhs_im', nberr); RET_ON_NBERR(nberr)
   call double_nalloc_1d(rhs_re, n_dof, 'rhs_re', nberr); RET_ON_NBERR(nberr)
   call double_nalloc_1d(rhs_im, n_dof, 'rhs_im', nberr); RET_ON_NBERR(nberr)

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



   !call complex_nalloc_1d(workev, 3*dim_krylov, 'workev', nberr);
   !call double_nalloc_1d(rwork, dim_krylov, 'rwork', nberr);
   !call logical_nalloc_1d(selecto, dim_krylov, 'selecto', nberr);

   ! alloc

   ! do i=1,n_modes
   !    v_evals_nu(i) = 0.0d0
   ! enddo

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
   lworkl_32 = int(lwork, 4)

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

   do while (arp_active)

      !  Test for dimesnion conditions in znaupd (err code = -3)
      !  Test for N=n_dof_32, NEV=n_modes_32, NCV=dim_krylov_32
      !  Need 0<n_modes_32<n_dof_32-1, 1<= dim_krylov_32-n_modes_32, dim_krylov_32<=n_dof_32
      if ((n_dof_32-1 .le. n_modes_32) .or. (dim_krylov_32-n_modes_32 .lt. 1) &
         .or.  dim_krylov_32 > n_dof_32) then
         write(emsg,'(A,A)') 'ARPACK eigensolver dimensional'//&
            ' conditions failed (would generate ARPACK znaupd error ' //&
            'code of -3).' // NEW_LINE('A'),&
            'You should probably increase the grid resolution.'
         call nberr%set(NBERR_BAD_ZNAUPD, emsg);
         return
      endif

      call znaupd (arp_ido, arp_bmat, n_dof_32, arp_which, n_modes_32, arp_tol,&
         vecs%resid, dim_krylov_32, vecs%v_schur, n_dof_32, arp_iparam,&
         ipntr_32, vecs%workd, vecs%workl, lworkl_32, vecs%rwork, arp_info)

      if (arp_ido .eq. -1 .or. arp_ido .eq. 1) then

         !------------------------------------------------------
         !  Apply  y <--- OP*x = inv[A-SIGMA*M]*M*x
         !  with x at x = vecs%workd(ipntr_32(1))
         !  and place the result at  y = vecs%workd(ipntr_32(2))           |
         !------------------------------------------------------

         call zcopy(n_dof_32, vecs%workd(ipntr_32(1)), 1,vecs%vect1_z, 1)
         call z_mxv_csc (n_dof, vecs%vect1_z, vecs%vect2_z, cscmat%n_nonz, cscmat%v_row_ind,&
            cscmat%v_col_ptr, cscmat%mOp_mass)

         rhs_re = dble(vecs%vect2_z)
         rhs_im = imag(vecs%vect2_z)


         !  solve Ax=b, without iterative refinement
         sys = 0
         call umf4zsol (UMFPACK_A, lhs_re, lhs_im, rhs_re, rhs_im,&
            umf_numeric, umf_control, umf_info)

         if (umf_info (1) .lt. 0) then
            write(ui,*) 'Error occurred in umf4zsol: ', umf_info (1)
            emsg = 'Error occurred in umf4zsol: '
            call nberr%set(NBERR_BAD_UMF4ZSOL, emsg);
            return
         endif

         vecs%vect2_z  = dcmplx (lhs_re, lhs_im )

         call zcopy(n_dof_32, vecs%vect2_z, 1, vecs%workd(ipntr_32(2)), 1)

      else if (arp_ido .eq. 2) then

         write(ui,*) 'VALPR_64: ATTENTION arp_ido = ', arp_ido
         write(ui,*) 'check the results...'

         !  ----------------------------------------------
         !  | Apply  y <--- M*x                       |
         !  | x = vecs%workd(ipntr_32(1))  et  y = vecs%workd(ipntr_32(2)) |
         !  ----------------------------------------------

         call zcopy(n_dof_32, vecs%workd(ipntr_32(1)), 1, vecs%vect1_z, 1)
         call z_mxv_csc (n_dof, vecs%vect1_z, vecs%vect2_z, cscmat%n_nonz, cscmat%v_row_ind,&
            cscmat%v_col_ptr, cscmat%mOp_mass)


         rhs_re = dble(vecs%vect2_z)
         rhs_im = imag(vecs%vect2_z)

         !  solve Ax=b, without iterative refinement
         sys = 0
         call umf4zsol (sys, lhs_re, lhs_im, rhs_re, rhs_im,&
            umf_numeric, umf_control, umf_info)

         if (umf_info (1) .lt. 0) then
            write(emsg, *) 'Error occurred in umf4zsol: ', umf_info (1)
            call nberr%set(NBERR_BAD_UMF4ZSOL, emsg);
            return
         endif

         vecs%vect2_z  = dcmplx (lhs_re, lhs_im )

         call zcopy(n_dof_32, vecs%vect2_z,1, vecs%workd(ipntr_32(2)), 1)

      else
         arp_active = .false.
      end if

   enddo

   !  ---------------------------------------------------
   !  | Either we have convergence, or there is an error. |
   !  ---------------------------------------------------

   n_conv = arp_iparam(5)

   ! if (arp_info .gt. 0) then
   !    write(ui,*)
   !    write(ui,*) "VALPR_64: arp_info != 0 : ", arp_info
   !    write(ui,*) "VALPR_64: arp_iparam(5)=", arp_iparam(5), n_modes_32
   !    write(ui,*) "VALPR_64: number of converged values = ",&
   !    &arp_iparam(5)
   !    write(ui,*)
   ! endif



   if (n_conv .ne. n_modes) then
      write(emsg, '(A,I5,I5)') &
         "py_calc_modes_AC: convergence problem " // &
         "in valpr_64: n_conv != n_modes  ", n_conv, n_modes
      errco = -7
      call nberr%set(errco, emsg); RET_ON_NBERR(nberr)

   endif






   if (arp_info.lt.0) then

!  ---------------------------------------------------
!  | Error message, check the documentation in DNAUPD. |

!  INFO    Integer.  (INPUT/OUTPUT)
!  If INFO .EQ. 0, a randomly initial residual vector is used.
!  If INFO .NE. 0, RESID contains the initial residual vector,
!  possibly from a previous run.
!  Error flag on output.
!  =  0: Normal exit.
!  =  1: Maximum number of iterations taken.
!  All possible eigen_modesues of OP has been found. IPARAM(5)
!  returns the number of wanted converged Ritz values.
!  =  2: No longer an informational error. Deprecated starting
!  with release 2 of ARPACK.
!  =  3: No shifts could be applied during a cycle of the
!  Implicitly restarted Arnoldi iteration. One possibility
!  is to increase the size of NCV relative to NEV.
!  See remark 4 below.
!  = -1: N must be positive.
!  = -2: NEV must be positive.
!  = -3: NCV-NEV >= 2 and less than or equal to N.
!  = -4: The maximum number of Arnoldi update iteration
!  must be greater than zero.
!  = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
!  = -6: BMAT must be one of 'I' or 'G'.
!  = -7: Length of private work array is not sufficient.
!  = -8: Error return from LAPACK eigen_modesue calculation;
!  = -9: Starting vector is zero.
!  = -10: IPARAM(7) must be 1,2,3,4.
!  = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
!  = -12: IPARAM(1) must be equal to 0 or 1.
!  = -9999: Could not build an Arnoldi factorization.
!  IPARAM(5) returns the size of the current Arnoldi
!  factorization.

!  ---------------------------------------------------

      write(emsg, '(A,I5,/,A)') 'Error occurred in _naupd:'//&
         ' ARPACK error code = ', arp_info,&
         ' You should probably increase the grid resolution.'

      write(*,*)
      return

   else


      !  Get the final eigenvectors
      !'A' means get the actual eigenvectors, not just schur/arnolid vectors

      arp_rvec = .true.
      arp_shift = (0.0d0,0.0d0)

      call zneupd (arp_rvec, 'A', vecs%arp_select, v_evals_nu, vecs%v_schur, n_dof_32, arp_shift,&
         vecs%workev, arp_bmat, n_dof_32, arp_which, n_modes_32, arp_tol,&
         vecs%resid, dim_krylov_32, vecs%v_schur, n_dof_32, arp_iparam, ipntr_32,&
         vecs%workd, vecs%workl, lworkl_32, vecs%rwork, ierr_32)

      !  Eigenvalues and eigenvectors:
      !  The real part of an eigenvalue is listed in the first column of the D table.
      !  The imaginary part of an eigenvalue is listed in the second column of the D table.
      !  The eigenvectors are stored in the first n_modes_32 columns of the V table
      !  when the arp_rvec option is set to true.
      !  Otherwise, the V table contains an orthogonal basis of the eigenspace.

      if (ierr_32.ne.0) then
         write(emsg,*) 'VALPR_64:' ,&
            ' Error with _neupd, arp_info = ', ierr_32,&
            ' Check the documentation of zneupd. ',&
            ' This error can occur if the mesh is too coarse.'
         errco = -109
         call nberr%set(errco, emsg);
         return

      else
         do i = 1, n_modes
            do j = 1, n_dof
               v_evecs_arp(j,i) = vecs%v_schur(j,i)
            enddo
         enddo
      endif
   endif

!  free the umf_numeric factorization
   call umf4zfnum (umf_numeric)


   !deallocate(rwork, STAT=alloc_stat)
   !deallocate(selecto)

   return
end

!------------------------------------------------------------------------
!------------------------------------------------------------------------

subroutine report_stats_umf4zsym(debug, show_mem_est, umf_info)

   use numbatmod

   double precision umf_info (90)
   integer(8) ui, debug, show_mem_est

   ui = stdout

   if (debug .eq. 1 .or. show_mem_est .eq. 1) then
      write(ui,80) umf_info (1), umf_info (16),&
      &(umf_info (21) * umf_info (4)) / 2**20,&
      &(umf_info (22) * umf_info (4)) / 2**20,&
      &umf_info (23), umf_info (24), umf_info (25)
80    format ('  symbolic analysis:',/,&
      &'   status:  ', f5.0, /,&
      &'   time:    ', e10.2, ' (sec)'/,&
      &'   estimates (upper bound) for umf_numeric LU:', /,&
      &'   size of LU:    ', f12.2, ' (MB)', /,&
      &'   memory needed: ', f12.2, ' (MB)', /,&
      &'   flop count:    ', e12.2, /&
      &'   nnz (L):       ', f12.0, /&
      &'   nnz (U):       ', f12.0)
   endif

end

!------------------------------------------------------------------------
!------------------------------------------------------------------------

subroutine report_stats_umf4znum(debug, show_mem_est, umf_info)

   use numbatmod


   double precision umf_info (90)
   integer(8) ui, debug, show_mem_est

   ui = stdout

   if (debug .eq. 1 .or. show_mem_est .eq. 1) then

      write(ui,90) umf_info (1), umf_info (66),&
      &(umf_info (41) * umf_info (4)) / 2**20,&
      &(umf_info (42) * umf_info (4)) / 2**20,&
      &umf_info (43), umf_info (44), umf_info (45)
90    format ('  umf_numeric factorization:',/,&
      &'   status:  ', f5.0, /,&
      &'   time:    ', e10.2, /,&
      &'   actual umf_numeric LU statistics:', /,&
      &'   size of LU:    ', f12.2, ' (MB)', /,&
      &'   memory needed: ', f12.2, ' (MB)', /,&
      &'   flop count:    ', e12.2, /&
      &'   nnz (L):       ', f12.0, /&
      &'   nnz (U):       ', f12.0)

   endif

end
