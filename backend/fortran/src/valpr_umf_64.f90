#include "numbat_decl.h"

module class_ValprVecs
   use numbatmod
   use alloc

   implicit none
   private

   type, public  :: ValprVecs

      integer(8) :: lworkl

      complex(8), dimension(:,:), allocatable:: v_schur

      complex(8), dimension(:), allocatable :: vect1
      complex(8), dimension(:), allocatable :: vect2
      complex(8), dimension(:), allocatable :: workd
      complex(8), dimension(:), allocatable :: workl
      complex(8), dimension(:), allocatable :: resid
      complex(8), dimension(:), allocatable :: eval_ritz

      complex(8), dimension(:), allocatable :: workev
      double precision, dimension(:), allocatable :: rwork
      logical, dimension(:), allocatable :: arp_select

   contains

      procedure :: init => ValprVecs_init
      final :: destructor

   end type ValprVecs

contains

   subroutine ValprVecs_init(this, n_modes, dim_krylov, neq, errco, emsg)
      class(ValprVecs) :: this
      integer(8) :: n_modes, dim_krylov, neq
      integer(8) errco
      character(len=EMSG_LENGTH) emsg

      errco = 0

      this%lworkl = 3 * dim_krylov**2 + 5 * dim_krylov

      call complex_alloc_2d(this%v_schur, neq, dim_krylov, 'v_schur', errco, emsg); RETONERROR(errco)
      call complex_alloc_1d(this%vect1, neq, 'vect1', errco, emsg); RETONERROR(errco)
      call complex_alloc_1d(this%vect2, neq, 'vect2', errco, emsg); RETONERROR(errco)
      call complex_alloc_1d(this%workd, 3*neq, 'workd', errco, emsg); RETONERROR(errco)
      call complex_alloc_1d(this%workl, this%lworkl, 'workl', errco, emsg); RETONERROR(errco)
      call complex_alloc_1d(this%resid, neq, 'resid', errco, emsg); RETONERROR(errco)
      call complex_alloc_1d(this%eval_ritz, n_modes+1, 'eval_ritz', errco, emsg); RETONERROR(errco)
      call complex_alloc_1d(this%workev, 3*dim_krylov, 'workl', errco, emsg); RETONERROR(errco)
      call double_alloc_1d(this%rwork, dim_krylov, 'rwork', errco, emsg); RETONERROR(errco)
      call logical_alloc_1d(this%arp_select, dim_krylov, 'arp_select', errco, emsg); RETONERROR(errco)



   end subroutine

   subroutine destructor(this)
      type(ValprVecs) :: this

      if (allocated(this%v_schur)) then

         deallocate(this%v_schur)
         deallocate(this%vect1,  this%vect2,    this%workd, this%workl)
         deallocate(this%resid,  this%eval_ritz)
         deallocate(this%workev, this%rwork,    this%arp_select)

      endif

   end subroutine


end module


!  UMFPACK docs:
!  https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/dev/UMFPACK/Doc/UMFPACK_UserGuide.pdf
!  https://fossies.org/linux/SuiteSparse/UMFPACK/Doc/UMFPACK_UserGuide.pdf

!  Calls in here are for the double complex precision, long int format, because compiled with ZLONG
!  And definitions in sswrap/umf4_f77_wrapper
!  This corresponds to integer(8) arrays
!  Yet calling is done with int4 scalar control parameters. TODO: necessary?

!  Arpack docs:
!  http://li.mit.edu/Archive/Activities/Archive/CourseWork/Ju_Li/MITCourses/18.335/Doc/ARPACK/Lehoucq97.pdf



!  ---------------------------------------------------
!  | Arpack Error message, check the documentation in ZNAUPD. |

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
!  = -5: arp_which must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
!  = -6: arp_bmat must be one of 'I' or 'G'.
!  = -7: Length of private work array is not sufficient.
!  = -8: Error return from LAPACK eigen_modesue calculation;
!  = -9: Starting vector is zero.
!  = -10: IPARAM(7) must be 1,2,3,4.
!  = -11: IPARAM(7) = 1 and arp_bmat = 'G' are incompatable.
!  = -12: IPARAM(1) must be equal to 0 or 1.
!  = -9999: Could not build an Arnoldi factorization.
!  IPARAM(5) returns the size of the current Arnoldi
!  factorization.

!  ---------------------------------------------------



subroutine apply_arpack_OPx(neq, x, y, nonz, row_ind, col_ptr, mat2, vect1, vect2, &
   lhs_re, lhs_im,  umf_numeric, umf_control, umf_info, errco, emsg)

   use numbatmod

   integer(8) neq
   complex(8) :: x(neq), y(neq)
   integer(8) nonz
   integer(8) row_ind(neq), col_ptr(neq)
   complex(8) mat2(nonz), vect1(neq), vect2(neq)
   double precision lhs_re(neq), lhs_im(neq)
   integer(8) umf_numeric
   double precision umf_control(UMFPACK_CONTROL)
   double precision umf_info(UMFPACK_INFO)
   integer(8) errco
   character(len=EMSG_LENGTH) emsg


   double precision rhs_re(neq), rhs_im(neq)
   integer(4) neq_32
   !integer(8) sys

   neq_32 = int(neq, 4)

   call zcopy(neq_32, x, 1, vect1, 1)                !  LAPACK routine
   call z_mxv_csc (neq, vect1, vect2, nonz, row_ind, col_ptr, mat2)   !  Local routine: vec2 = mat2.vect1

   !rhs_re = realpart(vect2)
   !rhs_im = imagpart(vect2)
   rhs_re = dble(vect2)
   rhs_im = dimag(vect2)


   !  solve Ax=b, without iterative refinement (UMFPACK_A)
   !sys = 0
   call umf4zsol (UMFPACK_A, lhs_re, lhs_im, rhs_re, rhs_im, &
      umf_numeric, umf_control, umf_info)

   if (umf_info (1) .lt. 0) then
      write(emsg,*) 'Error occurred in umf4zsol: ', umf_info (1)
      errco = NBERROR_106
      return
   endif

   vect2 = lhs_re + C_IM_ONE * lhs_im

   call zcopy(neq_32, vect2, 1, y, 1)

end subroutine





!  ------------------------------------------------------------------

!  sous-routine VALPR_64.
!  ------
!  gere l'utilisation de ARPACK

!  ------------------------------------------------------------------


!  Input Parameters:
!  -----------------

!  dim_krylov (I)              : dimension of the Krylov space
!  n_modes (I)            : number of desired eigenvalues
!  neq (I)                : number of equations
!  workd (DP)             : working matrix for dnaupd,
!  size 3 * neq
!  resid (DP)             : working vector for dnaupd,
!  size neq
!  v (DP)                 : matrix of Schur vectors,
!  size neq * dim_krylov
!  asup, ainf, adiag (DP) : storage of the tangent matrix
!  msup, minf, mdiag (DP) : storage of the mass matrix
!  kld (I)                : vector of column start locations
!  vect1, vect2, vect3 (DP) : working vectors,
!  size neq
!  long (I)               : length of super-vectors for
!  the matrices
!  ddot (DP)              : function called to compute
!  the dot product
!  workl (DP)              : working space for dnaupd,
!  size lworkl >= 3 * dim_krylov^2 + 6 * dim_krylov
!  action                 : 0: restart with a random vector
!  1: read an initial vector


!  Output Parameters:
!  ------------------

!  reel (DP)              : real parts, size dim_krylov
!  imag (DP)              : imaginary parts, size dim_krylov

!  ------------------------------------------------------------------
!  v = sovled eigvecs
!  d = solved eigvals

! del_vect1, del_vect2, del_workl, &
! ext_workd, ext_resid, ext_lworkl, &

subroutine valpr_64 (&
   i_base, dim_krylov, n_modes, neq, itermax, &
   arp_tol, nonz, &
   errco, emsg, &
   row_ind, col_ptr, &
   mat1, mat2, &
   v_evals, v_evecs )

   !  mat1 is the inverse shift operator  Op = inv[A-SIGMA*M]*M, where M=Idenity
   !  mat2 is the identity and hopefully never used ?

   use numbatmod
   use class_ValprVecs
   use alloc

   integer(8), intent(in) :: itermax, dim_krylov
   integer(8) neq, nonz, i_base, n_modes
   integer(8) row_ind(nonz), col_ptr(neq+1)

   complex(8), intent(in) :: mat1(nonz)
   complex(8), intent(in) :: mat2(nonz)

   complex(8), intent(out) :: v_evals(n_modes)

   complex(8), intent(out) :: v_evecs(neq, n_modes)

   integer(8) errco
   character(len=EMSG_LENGTH) emsg

   ! ----------------------------------------------------------

   ! UMFPACK requires complex arrays as pairs of doubles
   double precision, allocatable, dimension(:) :: mat1_re, mat1_im
   double precision, allocatable, dimension(:) :: lhx_re, lhx_im
   double precision, allocatable, dimension(:) :: rhs_re, rhs_im


   double precision umf_control(UMFPACK_CONTROL)
   double precision umf_info(UMFPACK_INFO)
   integer(8) umf_numeric, umf_symbolic


   type(ValprVecs) :: vecs

   integer(8) :: lworkl, n_conv

   !  32-bit integers for ARPACK
   integer(4) neq_32, n_modes_32, dim_krylov_32
   integer(4) arp_ido, arp_info, arp_iparam(11)
   integer(8),  parameter :: ARP_IPNTR_DIM = 14
   integer(4) ipntr_32(ARP_IPNTR_DIM), lworkl_32
   double precision arp_tol
   complex(8) arp_shift

   logical arp_rvec
   character arp_bmat
   character(2) arp_which
   logical arp_active

   integer(8) ui


   ui = stdout
   errco = 0
   emsg = ""

   if (i_base .ne. 0) then
      write(emsg,*) "valpr_64: i_base != 0 : ", i_base, &
         "valpr_64: UMFPACK requires 0-based indexing"
      errco = NBERROR_103
      return

   endif

   call vecs%init(n_modes, dim_krylov, neq, errco, emsg)
   RETONERROR(errco)

   lworkl = 3 * dim_krylov**2 + 5 * dim_krylov   ! length specified in znaupd.f source file

   v_evals = C_ZERO

   !  ----------------------------------------------------------------
   !  factor the matrix and save to a file
   !  ----------------------------------------------------------------


   !  umfpack * report status (print level = umf_control(1)) :
   !  print level = 0 or less : No output, even when an error occurs.
   !  print level = 1 (default value) : then error messages are printed,
   !  and nothing is printed if the status is UMFPACK OK.
   !  print level = 2 or more : then the status is always printed.
   !
   !  print umf_control parameters.  set umf_control (1) to 1 to print
   !  error messages only
   call umf4zdef (umf_control)    !  load default UMFPACK umf_control parameters
   umf_control (1) = 1            !  print only on errors (2 for everything)
   call umf4zpcon (umf_control)   !  print umf_control parameters (null op if umf_control(1)=1)


   !  Pre-order and symbolic analysis
   !  factors neq x neq matrix  in CSR format with col and row arrays col_ptr, row_ind
   !  complex entries are in mat1_re and mat1_im

   call double_alloc_1d(mat1_re, nonz, 'mat1_re', errco, emsg); RETONERROR(errco)
   call double_alloc_1d(mat1_im, nonz, 'mat1_im', errco, emsg); RETONERROR(errco)

   call double_alloc_1d(lhx_re, neq, 'lhx_re', errco, emsg); RETONERROR(errco)
   call double_alloc_1d(lhx_im, neq, 'lhx_im', errco, emsg); RETONERROR(errco)
   call double_alloc_1d(rhs_re, neq, 'rhs_re', errco, emsg); RETONERROR(errco)
   call double_alloc_1d(rhs_im, neq, 'rhs_im', errco, emsg); RETONERROR(errco)




   mat1_re = dble(mat1)
   mat1_im = dimag(mat1)



   ! do jj=1,neq+1
   !    write(*,*) jj, col_ptr(jj)
   ! end do

   ! do jj=1,nonz
   !    write(*,*) jj, row_ind(jj), mat1_re(jj), mat1_im(jj)
   ! end do

   umf_info(1) = 0
    call umf4zsym (neq, neq, col_ptr, row_ind, mat1_re, mat1_im, &
      umf_symbolic, umf_control, umf_info)


   if (umf_info (1) .lt. 0) then
      write(emsg,'(A,i4)') 'Error occurred in sparse matrix symbolic factorization umf4zsym:', int(umf_info (1))
      errco = NBERROR_104
      return
   endif


   !  Complete the umf_numeric factorization
   call umf4znum (col_ptr, row_ind, mat1_re, mat1_im, umf_symbolic, &
      umf_numeric, umf_control, umf_info)

   if (umf_info (1) .lt. 0) then
      write(emsg,*) 'Error occurred in sparse matrix umf_numeric factorization umf4znum: ', umf_info (1)
      errco = NBERROR_105
      return
   endif


   call umf4zfsym (umf_symbolic)   !  free the symbolic analysis

   !  The arp_ido parameter is used for reverse communication.
   !  Initially, it should be set to 0.

   !  Setting arp_info to 0 instructs ARPACK to construct an initial vector with random components.
   !  Setting arp_iparam(1) to 1 indicates that ARPACK should calculate translations
   !  based on the projected matrix and according to the "arp_which criterion.


   !Not sure that casting of these scalar vairables to int4 is really required
   neq_32 = int(neq, 4)
   n_modes_32 = int(n_modes, 4)
   dim_krylov_32 = int(dim_krylov, 4)
   lworkl_32 = int(lworkl, 4)

   arp_ido = 0
   arp_iparam(1) = 1                 !  exact shift strategy
   arp_iparam(3) = int(itermax, 4)   !  max iterations
   !  arp_iparam(7) = 3
   arp_iparam(7) = 1                 !  comp mode: matrix-vec products only, no shift?
   arp_info = 0                      !  Random initial candidate vector



   !----------------------------------------------------
   !  Main loop in inverse communication mode
   !----------------------------------------------------
   arp_bmat = 'I'    !  plain (not generalised) eigenvalue problem
   arp_which = 'LM'  !  seek largest magnitude eigs

   arp_active = .true.
   arp_shift = C_ONE   !  Is ignored, as long as iparam(7)=1


   do while (arp_active)
      !call znaupd (arp_ido, arp_bmat, neq_32, arp_which, n_modes_32, arp_tol, &
      !   resid, dim_krylov_32, v_schur, neq_32, arp_iparam,&
      !   ipntr_32, workd, workl, lworkl_32, rwork, arp_info)

      call znaupd (arp_ido, arp_bmat, neq_32, arp_which, n_modes_32, arp_tol, &
         vecs%resid, dim_krylov_32, vecs%v_schur, neq_32, arp_iparam,&
         ipntr_32, vecs%workd, vecs%workl, lworkl_32, vecs%rwork, arp_info)


      if (arp_ido .eq. -1 .or. arp_ido .eq. 1) then   !  Request for y = OP*x = inv[A-SIGMA*M]*M*x

         !------------------------------------------------------
         !  Apply  y <--- OP*x = inv[A-SIGMA*M]*M*x
         !  with x at x = workd(ipntr_32(1))
         !  and place the result at  y = workd(ipntr_32(2))           |
         !------------------------------------------------------

         call apply_arpack_OPx(neq, vecs%workd(ipntr_32(1)), vecs%workd(ipntr_32(2)), &
            nonz, row_ind, col_ptr, mat2, vecs%vect1, vecs%vect2, &
            lhx_re, lhx_im,  umf_numeric, umf_control, umf_info, errco, emsg)


      else if (arp_ido .eq. 2) then  !  Request for y = M*x    !TODO:  IO don't think this ever happens for bmat=I, ie M=I

         write(ui,*) 'VALPR_64: ATTENTION arp_ido = ', arp_ido
         write(ui,*) 'check the results...'

         !----------------------------------------------
         !  On execute  y <--- M*x
         !  x = workd(ipntr_32(1))  et  y = workd(ipntr_32(2))
         !----------------------------------------------

         call zcopy(neq_32, vecs%workd(ipntr_32(1)), 1, vecs%vect1, 1)
         call z_mxv_csc (neq, vecs%vect1, vecs%vect2, nonz, row_ind, col_ptr, mat2)

         rhs_re = dble(vecs%vect2)
         rhs_im = dimag(vecs%vect2)

         !  solve Ax=b, without iterative refinement
         call umf4zsol (UMFPACK_A, lhx_re, lhx_im, rhs_re, rhs_im, umf_numeric, umf_control, umf_info)
         if (umf_info (1) .lt. 0) then
            write(emsg,*) 'Error occurred in umf4zsol: ', umf_info (1)
            errco = NBERROR_107
            return
         endif

         vecs%vect2 = lhx_re + C_IM_ONE * lhx_im

         call zcopy(neq_32, vecs%vect2, 1, vecs%workd(ipntr_32(2)), 1)

      else !  we are done, for better or worse
         arp_active = .false.
      end if
   end do

   !--------------------------------------------------
   !  Either we have convergence, or there is an error. |
   !---------------------------------------------------

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
      !  write(ui,*) "VALPR_64: arp_iparam(5) = ", arp_iparam(5), n_modes_32
      !  write(ui,*) "VALPR_64: number of converged values = ",
      !  *                arp_iparam(5)
      write(ui,*)
   endif

   if (arp_info.lt.0) then
      write(emsg,*) "VALPR_64: The Arnoldi iteration scheme has failed.", &
         "VALPR_64: The znaupd error flag has the value ", &
         "arp_info=", arp_info

      errco = NBERROR_108
      return
   endif

   if (n_conv .ne. n_modes) then
      write(emsg,*) "Convergence problem in valpr_64: n_conv != n_modes : ", &
         n_conv, n_modes ,"You should probably increase resolution of mesh!"
      errco = -19
      return
   endif


   !  Get the final eigenvectors
   !'A' means get the actual eigenvectors, not just schur/arnolid vectors
   arp_rvec = .true. !  get the full set of vectors

   call zneupd (arp_rvec, 'A', vecs%arp_select, vecs%eval_ritz, &
      v_evecs, &
      neq_32, arp_shift, &
      vecs%workev, arp_bmat, neq_32, arp_which, n_modes_32, arp_tol, &
      vecs%resid, dim_krylov_32, vecs%v_schur, neq_32, arp_iparam, ipntr_32, &
      vecs%workd, vecs%workl, lworkl_32, vecs%rwork, arp_info)


   !  Eigenvalues and eigenvectors:
   !  The real part of an eigenvalue is listed in the first column of the D table.
   !  The imaginary part of an eigenvalue is listed in the second column of the D table.
   !  The eigenvectors are stored in the first n_modes_32 columns of the V table
   !  when the arp_rvec option is set to true.
   !  Otherwise, the V table contains an orthogonal basis of the eigenspace.

   if (arp_info .ne. 0) then
      write(emsg,*) 'VALPR_64: Error with _neupd, arp_info = ', arp_info
      errco = NBERROR_109
   else
      v_evals = vecs%eval_ritz(1:n_modes) ! eval_ritz is one longer due to zneupd requirements
   endif



   call umf4zfnum (umf_numeric)   !  free the umf_numeric factorization
   deallocate(mat1_re, mat1_im)
   deallocate(lhx_re, lhx_im)

   return
end
