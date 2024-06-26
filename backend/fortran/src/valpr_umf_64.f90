#include "numbat_decl.h"


! UMFPACK docs:
! https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/dev/UMFPACK/Doc/UMFPACK_UserGuide.pdf
! https://fossies.org/linux/SuiteSparse/UMFPACK/Doc/UMFPACK_UserGuide.pdf
!
! Calls in here are for the double complex precision, long int format, because compiled with ZLONG
! And definitions in sswrap/umf4_f77_wrapper
! This corresponds to integer(8) arrays
! Yet calling is done with int4 scalar control parameters. TODO: necessary?

! Arpack docs:
! http://li.mit.edu/Archive/Activities/Archive/CourseWork/Ju_Li/MITCourses/18.335/Doc/ARPACK/Lehoucq97.pdf
!


!      ---------------------------------------------------
!     | Arpack Error message, check the documentation in ZNAUPD. |

!  INFO    Integer.  (INPUT/OUTPUT)
!          If INFO .EQ. 0, a randomly initial residual vector is used.
!          If INFO .NE. 0, RESID contains the initial residual vector,
!                          possibly from a previous run.
!          Error flag on output.
!          =  0: Normal exit.
!          =  1: Maximum number of iterations taken.
!                All possible eigen_modesues of OP has been found. IPARAM(5)
!                returns the number of wanted converged Ritz values.
!          =  2: No longer an informational error. Deprecated starting
!                with release 2 of ARPACK.
!          =  3: No shifts could be applied during a cycle of the
!                Implicitly restarted Arnoldi iteration. One possibility
!                is to increase the size of NCV relative to NEV.
!                See remark 4 below.
!          = -1: N must be positive.
!          = -2: NEV must be positive.
!          = -3: NCV-NEV >= 2 and less than or equal to N.
!          = -4: The maximum number of Arnoldi update iteration
!                must be greater than zero.
!          = -5: arp_which must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
!          = -6: arp_bmat must be one of 'I' or 'G'.
!          = -7: Length of private work array is not sufficient.
!          = -8: Error return from LAPACK eigen_modesue calculation;
!          = -9: Starting vector is zero.
!          = -10: IPARAM(7) must be 1,2,3,4.
!          = -11: IPARAM(7) = 1 and arp_bmat = 'G' are incompatable.
!          = -12: IPARAM(1) must be equal to 0 or 1.
!          = -9999: Could not build an Arnoldi factorization.
!                   IPARAM(5) returns the size of the current Arnoldi
!                   factorization.

!      ---------------------------------------------------



module valpr_help
   use numbatmod

   ! interface
   !    subroutine valpr_allocs(dim_krylov, workev, rwork, select, errco, emsg)
   !       use numbatmod
   !       integer(8) dim_krylov
   !       integer alloc_stat

   !       complex(8), dimension(:), allocatable, intent(inout) :: workev
   !       double precision, dimension(:), allocatable, intent(inout) :: rwork
   !       logical, dimension(:), allocatable, intent(inout) :: select

   !       integer, intent(out) :: errco
   !       character(len=EMSG_LENGTH), intent(out) :: emsg
   !    end subroutine
   ! end interface

contains

   subroutine valpr_allocs(dim_krylov, workev, rwork, select, errco, emsg)

      use numbatmod
      integer(8) dim_krylov
      integer alloc_stat

      complex(8), dimension(:), allocatable, intent(inout) :: workev
      double precision, dimension(:), allocatable, intent(inout) :: rwork
      logical, dimension(:), allocatable, intent(inout) :: select

      integer, intent(out) :: errco
      character(len=EMSG_LENGTH), intent(out) :: emsg

      alloc_stat = 0
      allocate(workev(3*dim_krylov), rwork(dim_krylov), STAT=alloc_stat)
      if (alloc_stat /= 0) then
         write(emsg,*) "VALPR_64: Mem. allocation is unsuccessfull", &
            " for the arrays workev, rwork. alloc_stat, dim_krylov = ", alloc_stat, dim_krylov
         errco = NBERROR_101
         return
      endif

      allocate(select(dim_krylov), STAT=alloc_stat)
      if (alloc_stat /= 0) then
         write(emsg,*) "VALPR_64: Mem. allocation is unsuccessfull", &
            " for the array select. alloc_stat, dim_krylov = ", alloc_stat, dim_krylov
         errco = NBERROR_102
         return
      endif

   end subroutine

end module valpr_help


subroutine apply_arpack_OPx(neq, x, y, nonz, row_ind, col_ptr, mat2, vect1, vect2, &
   sys, lhs_re, lhs_im,  numeric, umf_control, umf_info, errco, emsg)

   use numbatmod
   integer(8) neq
   complex(8) :: x(neq), y(neq)
   integer(8) nonz
   integer(8) row_ind(neq), col_ptr(neq)
   complex(8) mat2(nonz), vect1(neq), vect2(neq)
   double precision lhs_re(neq), lhs_im(neq)
   integer(8) numeric, symbolic, sys
   double precision umf_control(UMFPACK_CONTROL)
   double precision umf_info(UMFPACK_INFO)
   integer errco
   character(len=EMSG_LENGTH) emsg


   double precision rhs_re(neq), rhs_im(neq)
   integer(4) neq_32

   neq_32 = int(neq, 4)

   call zcopy(neq_32, x, 1, vect1, 1)                ! LAPACK routine
   call z_mxv_csc (neq, vect1, vect2, nonz, row_ind, col_ptr, mat2)   ! Local routine: vec2 = mat2.vect1

   rhs_re = realpart(vect2)
   rhs_im = imagpart(vect2)

!       solve Ax=b, without iterative refinement
   sys = 0
   call umf4zsol (sys, lhs_re, lhs_im, rhs_re, rhs_im, numeric, umf_control, umf_info)

   if (umf_info (1) .lt. 0) then
      write(emsg,*) 'Error occurred in umf4zsol: ', umf_info (1)
      errco = NBERROR_106
      return
   endif

   vect2 = lhs_re + C_IM_ONE * lhs_im

   call zcopy(neq_32, vect2, 1, y, 1)

end subroutine





!     ------------------------------------------------------------------
!
!                    sous-routine VALPR_64.
!                                 ------
!     gere l'utilisation de ARPACK
!
!     ------------------------------------------------------------------
!

!     Input Parameters:
!     -----------------
!
!     dim_krylov (I)              : dimension of the Krylov space
!     n_modes (I)            : number of desired eigenvalues
!     neq (I)                : number of equations
!     workd (DP)             : working matrix for dnaupd,
!                              size 3 * neq
!     resid (DP)             : working vector for dnaupd,
!                              size neq
!     v (DP)                 : matrix of Schur vectors,
!                              size neq * dim_krylov
!     asup, ainf, adiag (DP) : storage of the tangent matrix
!     msup, minf, mdiag (DP) : storage of the mass matrix
!     kld (I)                : vector of column start locations
!     vect1, vect2, vect3 (DP) : working vectors,
!                              size neq
!     long (I)               : length of super-vectors for
!                              the matrices
!     ddot (DP)              : function called to compute
!                              the dot product
!     workl (DP)              : working space for dnaupd,
!                              size lworkl >= 3 * dim_krylov^2 + 6 * dim_krylov
!     action                 : 0: restart with a random vector
!                              1: read an initial vector
!
!
!     Output Parameters:
!     ------------------
!
!     reel (DP)              : real parts, size dim_krylov
!     imag (DP)              : imaginary parts, size dim_krylov
!
!     ------------------------------------------------------------------
!  v = sovled eigvecs
!  d = solved eigvals
!
subroutine valpr_64 (i_base, dim_krylov, n_modes, neq, itermax, lworkl,&
   tol, nonz, row_ind, col_ptr, mat1_re, mat1_im, mat2,&
   del_vect1, del_vect2, workd, resid, v, d, del_workl, vp,&
   rhs_re, rhs_im, lhs_re, lhs_im, &
   n_conv, time_fact, time_arpack, debug, errco, emsg)

   ! mat1 is the inverse shift operator  Op = inv[A-SIGMA*M]*M, where M=Idenity
   ! mat2 is the identity and hopefully never used ?

   use numbatmod
   use class_stopwatch
   use valpr_help


   integer(8) neq, nonz, n_conv, i_base, n_modes
   integer(8) row_ind(nonz), col_ptr(neq+1)
   complex(8) mat2(nonz)
   double precision mat1_re(nonz), mat1_im(nonz)
   double precision rhs_re(neq), rhs_im(neq)
   double precision lhs_re(neq), lhs_im(neq)

   double precision time_fact, time_arpack
   integer errco
   character(len=EMSG_LENGTH) emsg


   integer(8) itermax, dim_krylov,lworkl
   complex(8) resid(neq), v(neq,dim_krylov), workd(3*neq)
   complex(8) d(n_modes+1),  vp(neq,n_modes)



   complex(8) shift2
   integer(8) i, j

   integer(8) counter
   complex(8) del_vect1(neq), del_vect2(neq), del_workl(lworkl)

   double precision umf_control(UMFPACK_CONTROL)
   double precision umf_info(UMFPACK_INFO)
   complex(8) :: vect1(neq), vect2(neq), workl(lworkl)

   integer(8) numeric, symbolic, sys

   double precision tol
   !
   !      integer(8) max_dim_krylov
   !      parameter(max_dim_krylov=3000) ! previously 1500



   integer alloc_stat
   double precision, dimension(:), allocatable :: rwork

   complex(8), dimension(:), allocatable :: workev
   logical, dimension(:), allocatable :: select



   type(Stopwatch) :: clock_main

!     Local variables
   !     32-bit integers for ARPACK
   integer(4) neq_32, n_modes_32, dim_krylov_32
   integer(4) arp_ido, arp_info, arp_iparam(11)
   integer(4) ipntr_32(14), lworkl_32
   !

   logical rvec

   integer(8) ui, debug

   character arp_bmat
   character(2) arp_which
   logical arp_active


   !      common/imp/ui, debug
   !
   !     ------------------------------------------------------------------
   !
   ui = stdout
   errco = 0
   emsg = ""

   if (i_base .ne. 0) then
      write(emsg,*) "valpr_64: i_base != 0 : ", i_base, &
         "valpr_64: UMFPACK requires 0-based indexing"
      errco = NBERROR_103
      return

   endif

   ! TODO: check this works on intel
   call valpr_allocs(dim_krylov, workev, rwork, select, errco, emsg)

   ! alloc_stat = 0
   ! allocate(workev(3*dim_krylov), rwork(dim_krylov), STAT=alloc_stat)
   ! if (alloc_stat /= 0) then
   !    write(emsg,*) "VALPR_64: Mem. allocation is unsuccessfull", &
   !       " for the arrays workev, rwork. alloc_stat, dim_krylov = ", alloc_stat, dim_krylov
   !    errco = NBERROR_101
   !    return
   ! endif

   ! allocate(select(dim_krylov), STAT=alloc_stat)
   ! if (alloc_stat /= 0) then
   !    write(emsg,*) "VALPR_64: Mem. allocation is unsuccessfull", &
   !       " for the array select. alloc_stat, dim_krylov = ", alloc_stat, dim_krylov
   !    errco = NBERROR_102
   !    return
   ! endif


   shift2 = C_ONE

   ! TODO: why is this not zeroed to n_modes+1?
   do i=1,n_modes
      d(i) = 0.0d0
   enddo

!       ----------------------------------------------------------------
!       factor the matrix and save to a file
!       ----------------------------------------------------------------


   call clock_main%reset()


   !     umfpack * report status (print level = umf_control(1)) :
   !     print level = 0 or less : No output, even when an error occurs.
   !     print level = 1 (default value) : then error messages are printed,
   !                      and nothing is printed if the status is UMFPACK OK.
   !     print level = 2 or more : then the status is always printed.
   !
   !       print umf_control parameters.  set umf_control (1) to 1 to print
   !       error messages only
   call umf4zdef (umf_control)    ! load default UMFPACK umf_control parameters
   umf_control (1) = 1            ! print only on errors (2 for everything)
   call umf4zpcon (umf_control)   ! print umf_control parameters (null op if umf_control(1)=1)


   ! Pre-order and symbolic analysis
   ! factors neq x neq matrix  in CSR format with col and row arrays col_ptr, row_ind
   ! complex entries are in mat1_re and mat1_im
   call umf4zsym (neq, neq, col_ptr, row_ind, mat1_re, mat1_im, &
      symbolic, umf_control, umf_info)


   if (umf_info (1) .lt. 0) then
      write(emsg,*) 'Error occurred in sparse matrix symbolic factorization umf4zsym: ', umf_info (1)
      errco = NBERROR_104
      return
   endif

   ! print statistics computed so far
   ! call umf4zpinf (umf_control, umf_info) could also be done.
   if (debug .eq. 1) then
      write(ui,80) umf_info (1), umf_info (16),&
      &(umf_info (21) * umf_info (4)) / 2**20,&
      &(umf_info (22) * umf_info (4)) / 2**20,&
      &umf_info (23), umf_info (24), umf_info (25)
80    format ('  symbolic analysis:',/,&
         '   status:  ', f5.0, /,&
         '   time:    ', e10.2, ' (sec)'/,&
         '   estimates (upper bound) for numeric LU:', /,&
         '   size of LU:    ', f12.2, ' (MB)', /,&
         '   memory needed: ', f12.2, ' (MB)', /,&
         '   flop count:    ', e12.2, /&
         '   nnz (L):       ', f12.0, /&
         '   nnz (U):       ', f12.0)
   endif


   !  Complete the numeric factorization
   call umf4znum (col_ptr, row_ind, mat1_re, mat1_im, symbolic, &
      numeric, umf_control, umf_info)

   if (umf_info (1) .lt. 0) then
      write(emsg,*) 'Error occurred in sparse matrix numeric factorization umf4znum: ', umf_info (1)
      errco = NBERROR_105
      return
   endif

!  print statistics for the numeric factorization
!  call umf4zpinf (umf_control, umf_info) could also be done.
   if (debug .eq. 1) then
      write(ui,90) umf_info (1), umf_info (66),&
      &(umf_info (41) * umf_info (4)) / 2**20,&
      &(umf_info (42) * umf_info (4)) / 2**20,&
      &umf_info (43), umf_info (44), umf_info (45)
90    format ('  numeric factorization:',/,&
      &'   status:  ', f5.0, /,&
      &'   time:    ', e10.2, /,&
      &'   actual numeric LU statistics:', /,&
      &'   size of LU:    ', f12.2, ' (MB)', /,&
      &'   memory needed: ', f12.2, ' (MB)', /,&
      &'   flop count:    ', e12.2, /&
      &'   nnz (L):       ', f12.0, /&
      &'   nnz (U):       ', f12.0)
   endif



!       save the symbolic analysis to the file s42.umf
!       note that this is not needed until another matrix is
!       factorized, below.
!	filenum = 42
!        call umf4zssym (symbolic, filenum, status)
!        if (status .lt. 0) then
!            write(ui,*) 'Error occurred in umf4zssym: ', status
!            stop
!        endif
!
!       save the LU factors to the file n0.umf
!        call umf4zsnum (numeric, filenum, status)
!        if (status .lt. 0) then
!            write(ui,*) 'Error occurred in umf4zsnum: ', status
!            stop
!        endif

!       free the symbolic analysis
   call umf4zfsym (symbolic)

!       free the numeric factorization
!        call umf4zfnum (numeric)

   call clock_main%stop()
   time_fact = clock_main%cpu_time()

!
!       No LU factors (symbolic or numeric) are in memory at this point.
!
!c       ----------------------------------------------------------------
!c       load the LU factors back in, and solve the system
!c       ----------------------------------------------------------------
!
!c       At this point the program could terminate and load the LU
!C       factors (numeric) from the n0.umf file, and solve the
!c       system (see below).  Note that the symbolic object is not
!c       required.
!
!c       load the numeric factorization back in (filename: n0.umf)
!        call umf4zlnum (numeric, filenum, status)
!        if (status .lt. 0) then
!            write(ui,*) 'Error occurred in umf4zlnum: ', status
!            stop
!        endif
!
!       ##################################################################
!





! The arp_ido parameter is used for reverse communication.
! Initially, it should be set to 0.

! Setting arp_info to 0 instructs ARPACK to construct an initial vector with random components.
! Setting arp_iparam(1) to 1 indicates that ARPACK should calculate translations based
!  on the projected matrix and according to the "arp_which" criterion.                                           |


   !Not sure that casting of these scalar vairables to int4 is really required
   neq_32 = int(neq, 4)
   n_modes_32 = int(n_modes, 4)
   dim_krylov_32 = int(dim_krylov, 4)
   lworkl_32 = int(lworkl, 4)

   arp_ido = 0
   arp_iparam(1) = 1                 ! exact shift strategy
   arp_iparam(3) = int(itermax, 4)   ! max iterations
!      arp_iparam(7) = 3
   arp_iparam(7) = 1                 ! comp mode: matrix-vec products only, no shift?
   arp_info = 0                      ! Random initial candidate vector



!----------------------------------------------------
!    Main loop in inverse communication mode
!----------------------------------------------------

   call clock_main%reset()
   counter = 0

   arp_bmat = 'I'    ! plain (not generalised) eigenvalue problem
   arp_which = 'LM'  ! seek largest magnitude eigs

   arp_active = .true.

   do while (arp_active)


   ! zn: double complex non symmetric

   call znaupd (arp_ido, arp_bmat, neq_32, arp_which, n_modes_32, tol, &
    resid, dim_krylov_32, v, neq_32, arp_iparam,&
    ipntr_32, workd, workl, lworkl_32, rwork, arp_info)
!
   counter = counter + 1








   if (arp_ido .eq. -1 .or. arp_ido .eq. 1) then   ! Request for y = OP*x = inv[A-SIGMA*M]*M*x

      call apply_arpack_OPx(neq, workd(ipntr_32(1)), workd(ipntr_32(2)), &
      nonz, row_ind, col_ptr, mat2, vect1, vect2, &
      sys, lhs_re, lhs_im,  numeric, umf_control, umf_info, errco, emsg)
!      ------------------------------------------------------
!     | Apply  y <--- OP*x = inv[A-SIGMA*M]*M*x         |
!     | with x at x = workd(ipntr_32(1))
      ! and place the result at  y = workd(ipntr_32(2))           |
!      ------------------------------------------------------


      !TODO: delete all this as functdion is working
!       call zcopy(neq_32, workd(ipntr_32(1)), 1, vect1, 1)                ! LAPACK routine
!       call z_mxv_csc (neq, vect1, vect2, nonz, row_ind, col_ptr, mat2)   ! Local routine: vec2 = mat2.vect1

!       !do i=1,neq
!       !   rhs_re(i) = dble(vect2(i))
!       !   rhs_im(i) = dimag(vect2(i))
!       !enddo
!       rhs_re = realpart(vect2)
!       rhs_im = imagpart(vect2)

! !       solve Ax=b, without iterative refinement
!       sys = 0
!       call umf4zsol (sys, lhs_re, lhs_im, rhs_re, rhs_im, numeric, umf_control, umf_info)

!       if (umf_info (1) .lt. 0) then
!          write(emsg,*) 'Error occurred in umf4zsol: ', umf_info (1)
!          errco = NBERROR_106
!          return
!       endif

!       ! TODO: how to make this a one liner?
!       ! do i=1,neq
!       !    vect2 (i) = dcmplx (lhs_re (i), lhs_im (i))
!       ! enddo
!       vect2 = lhs_re + C_IM_ONE * lhs_im

!       call zcopy(neq_32, vect2, 1, workd(ipntr_32(2)), 1)


   else if (arp_ido .eq. 2) then  ! Request for y = M*x    !TODO:  IO don't think this ever happens for bmat=I, ie M=I

      write(ui,*) 'VALPR_64: ATTENTION arp_ido = ', arp_ido
      write(ui,*) 'check the results...'

!           ----------------------------------------------
!          | On execute  y <--- M*x                       |
!          | x = workd(ipntr_32(1))  et  y = workd(ipntr_32(2)) |
!           ----------------------------------------------

      call zcopy(neq_32, workd(ipntr_32(1)), 1, vect1, 1)
      call z_mxv_csc (neq, vect1, vect2, nonz, row_ind, col_ptr, mat2)
!
      ! do i=1,neq
      !    rhs_re(i) = dble(vect2(i))
      !    rhs_im(i) = imag(vect2(i))
      ! enddo
      rhs_re = realpart(vect2)
      rhs_im = imagpart(vect2)
!
!       solve Ax=b, without iterative refinement
      sys = 0
      call umf4zsol (sys, lhs_re, lhs_im, rhs_re, rhs_im, numeric, umf_control, umf_info)
      if (umf_info (1) .lt. 0) then
         write(emsg,*) 'Error occurred in umf4zsol: ', umf_info (1)
         errco = NBERROR_107
         return
      endif

      !do i=1,neq
      !   vect2 (i) = dcmplx (lhs_re (i), lhs_im (i))
      !enddo

      vect2 = lhs_re + C_IM_ONE * lhs_im

      call zcopy(neq_32, vect2, 1, workd(ipntr_32(2)), 1)

   else ! we are done, for better or rose
      arp_active = .false.
   end if
end do
!      ---------------------------------------------------
!     | Either we have convergence, or there is an error. |
!      ---------------------------------------------------

   n_conv = arp_iparam(5)

   if (arp_info .gt. 0) then
      write(ui,*)
      write(ui,*) "VALPR_64: The Arnoldi iteration scheme has failed"
      write(ui,*) "VALPR_64: The znaupd error flag has the value ",&
      &"arp_info=", arp_info
      if (arp_info .eq. 1) then
         write(ui,*) "VALPR_64: Max iterations exceeded."
         write(ui,*) " Requested eigen_modes = ", n_modes_32
         write(ui,*) " Converged eigen_modes = ", n_conv
         write(ui,*) " You might try:"
         write(ui,*) "  1) Increasing the requested number",&
         &" of eigen_modesues"
         write(ui,*) "  2) Increasing the grid resolution"
      endif
      if (arp_info .eq. 3) then
         write(ui,*) "VALPR_64: Shift could not be applied."
      endif
      write(ui,*) "VALPR_64: For details on znaupd errors see",&
      &" https://www.caam.rice.edu/software/ARPACK/UG/node138.html"
!        write(ui,*) "VALPR_64: arp_iparam(5) = ", arp_iparam(5), n_modes_32
!        write(ui,*) "VALPR_64: number of converged values = ",
!     *                arp_iparam(5)
      write(ui,*)
   endif

   if (arp_info.lt.0) then
      write(emsg,*) "VALPR_64: The Arnoldi iteration scheme has failed.", &
      "VALPR_64: The znaupd error flag has the value ", &
      "arp_info=", arp_info

      errco = NBERROR_108
      return
   else

!      -------------------------------------
!     | Ici on recupere les valeurs propres |
!      -------------------------------------

      rvec = .true.

      call zneupd (rvec, 'A', select, d, v, neq_32, shift2,&
      &workev, arp_bmat, neq_32, arp_which, n_modes_32, tol,&
      &resid, dim_krylov_32, v, neq_32, arp_iparam, ipntr_32,&
      &workd, workl, lworkl_32, rwork, arp_info)
!      ------------------------------------------------------------
!     | La partie reelle d'une valeur propre se trouve dans la     |
!     | premiere colonne du tableau D, la partie imaginaire est    |
!     | dans la seconde.                                           |
!     | Les vecteurs propres sont dans les premieres n_modes_32 colonnes |
!     | du tableau V, lorsque demande (rvec=.true.). Sinon, on y   |
!     | trouve une base orthogonale de l'espace propre.            |
!      ------------------------------------------------------------

      if (arp_info.ne.0) then

!      -----------------------------------------------------
!     | Error condition: Check the documentation of DNEUPD. |
!      -----------------------------------------------------

         write(ui,*) 'VALPR_64:'
         write(ui,*) ' Error with _neupd, arp_info = ', arp_info
         write(ui,*) ' Check the documentation of _neupd. '
         write(ui,*) 'Aborting...'
         stop

      else
         do i = 1, n_modes
            do j = 1, neq
               vp(j,i) = v(j,i)
            enddo
         enddo
      endif
   endif

!       free the numeric factorization
   call umf4zfnum (numeric)

   call clock_main%stop()
   time_arpack = clock_main%cpu_time()

   deallocate(workev, rwork, STAT=alloc_stat)
   deallocate(select)

   return
end
