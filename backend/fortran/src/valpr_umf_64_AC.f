c     ------------------------------------------------------------------
c
c     VALPR_64.
c                                 ------
c     Manages the use of ARPACK
c
c     ------------------------------------------------------------------

      subroutine valpr_64_AC (i_base, nvect, n_modes, neq, itermax,
     *  ltrav, tol, nonz, row_ind, col_ptr, mat1_re, mat1_im, mat2,
     *  vect1, vect2, workd, resid, vschur, nu_out, trav, vp,
     *  rhs_re, rhs_im, lhs_re, lhs_im, n_conv,
     *  debug, show_mem_est, errno, emsg)

c     ------------------------------------------------------------------

      use numbatmod
c
      integer(8) :: n_modes
      integer(8) neq, nonz, n_conv, i_base, nvect, ltrav
      integer(8) row_ind(nonz), col_ptr(neq+1)

      integer errno
      character(len=EMSG_LENGTH) emsg


      double precision mat1_re(nonz), mat1_im(nonz)
      double precision rhs_re(neq), rhs_im(neq)
      double precision lhs_re(neq), lhs_im(neq)

      complex(8) mat2(nonz)
      complex(8) resid(neq), vschur(neq,nvect), workd(3*neq)
      complex(8) vect1(neq), vect2(neq), trav(ltrav)
      complex(8) nu_out(n_modes+1), shift2, vp(neq,n_modes)
c
      double precision time1_fact, time2_fact
c
      double precision control (20), info_umf (90)
      integer(8) numeric, symbolic, sys
c
      integer(8) itermax, i, j
      integer(8) compteur



c
      double precision tol
c
      integer alloc_stat
      complex(8), dimension(:), allocatable :: workev
      double precision, dimension(:), allocatable :: rwork
      logical, dimension(:), allocatable :: selecto


c     Local variables
c     32-bit integers for ARPACK
      integer*4 neq_32, n_modes_32, nvect_32
      integer*4 ido_32, info_32, ierr_32, iparam_32(11)
      integer*4 ipntr_32(14), ltrav_32
c
      logical rvec
      character bmat*1, which*2

c      data bmat/'G'/
C       data which/'SR'/
c      data which/'SM'/
      data bmat/'I'/
      data which/'LM'/
c
      integer(8) ui, debug, show_mem_est
c      common/imp/ui, debug
c
c     ------------------------------------------------------------------
c
      ui = 6
      errno = 0
      emsg = ""

      call cpu_time(time1_fact)

c       ----------------------------------------------------------------
c       factor the matrix as A = LU and save to a file
c       ----------------------------------------------------------------
c
      if (debug .eq. 1) then
        write(ui,*) "valpr_64: factorisation (UMFPACK)"
      endif

c     set default parameters
      call umf4zdef (control)

c     umfpack * report status (print level = control(1)) :
c     print level = 0 or less : No output, even when an error occurs.
c     print level = 1 (default value) : then error messages are printed,
c                      and nothing is printed if the status is UMFPACK OK.
c     print level = 2 or more : then the status is always printed.
      control (1) = 1
      call umf4zpcon (control)

c     pre-order and symbolic analysis
      call umf4zsym (neq, neq, col_ptr, row_ind,
     *       mat1_re, mat1_im, symbolic, control, info_umf)

c     print statistics computed so far
c     call umf4zpinf (control, info_umf) could also be done.
      call report_stats_umf4zsym(debug, show_mem_est, info_umf)
      if (info_umf (1) .lt. 0) then
          write(emsg,*) 'Error occurred in umf4zsym: ', info_umf (1)
          errno = -104
          return
      endif


c     numeric factorization
c     TODO: This call does not appear to be thread safe! Breaks tutorial 3 B in thread mode
      call umf4znum (col_ptr, row_ind, mat1_re,
     *         mat1_im, symbolic, numeric, control, info_umf)

c     call umf4zpinf (control, info_umf) could also be done.
      call report_stats_umf4znum(debug, show_mem_est, info_umf)

      if (info_umf (1) .lt. 0) then
        write(emsg,*) 'Error occurred in umf4znum: ', info_umf(1)
        errno = -105
        return
      endif

      alloc_stat = 0
      allocate(workev(3*nvect), rwork(nvect), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(emsg,*) "VALPR_64: Mem. allocation is unsuccessfull ",
     *  "for the arrays workev, rwork",
     *  "alloc_stat, nvect = ", alloc_stat, nvect
        errno = -100
        return
      endif

      allocate(selecto(nvect), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(emsg,*) "VALPR_64: Mem. allocation is unsuccessfull ",
     *  "for the array selecto",
     *  "alloc_stat, nvect = ", alloc_stat, nvect
        errno = -101
        return
      endif

c
c
      do i=1,n_modes
        nu_out(i) = 0.0d0
      enddo
c
c       ##################################################################
c
      if (i_base .ne. 0) then
        write(emsg,*) "valpr_64: i_base != 0 : ", i_base,
     *   "valpr_64: UMFPACK requires 0-based indexing"
        errno = -102
      endif
c

c       save the symbolic analysis to the file s42.umf
c       note that this is not needed until another matrix is
c       factorized, below.
c	filenum = 42
c        call umf4zssym (symbolic, filenum, status)
c        if (status .lt. 0) then
c            write(ui,*) 'Error occurred in umf4zssym: ', status
c            stop
c        endif
c
c       save the LU factors to the file n0.umf
c        call umf4zsnum (numeric, filenum, status)
c        if (status .lt. 0) then
c            write(ui,*) 'Error occurred in umf4zsnum: ', status
c            stop
c        endif

c       free the symbolic analysis
        call umf4zfsym (symbolic)

cc       free the numeric factorization
c        call umf4zfnum (numeric)
c
      call cpu_time(time2_fact)
      if (debug .eq. 1) then
        write(ui,*) "valpr_64: factorisation completed"
        write(ui,*) "LU factorisation : CPU time = ",
     *         (time2_fact-time1_fact)
c ,
c     *         100*(time2_fact-time1_fact)/(time2-time1),"%"
      endif
c
c       No LU factors (symbolic or numeric) are in memory at this point.
c
cc       ----------------------------------------------------------------
cc       load the LU factors back in, and solve the system
cc       ----------------------------------------------------------------
c
cc       At this point the program could terminate and load the LU
cC       factors (numeric) from the n0.umf file, and solve the
cc       system (see below).  Note that the symbolic object is not
cc       required.
c
cc       load the numeric factorization back in (filename: n0.umf)
c        call umf4zlnum (numeric, filenum, status)
c        if (status .lt. 0) then
c            write(ui,*) 'Error occurred in umf4zlnum: ', status
c            stop
c        endif
c
c       ##################################################################
c
c     On commence le travail avec znaupd
c     ----------------------------------
c
c      ------------------------------------------------------------
c     | Le parametre IDO_32 est utilise pour la communication.        |
c     | A l'etape initiale il doit valoir 0.                       |
c     | Le choix INFO_32=0 correspond a la construction par           |
c     | Arpack d'un vecteur initial a composantes aleatoires       |
c     | iparam_32(1)=1 signifie que Arpack calcule les translations   |
c     | a partir de la matrice projetee et conformement au critere |
c     | "which".                                                   |
c      ------------------------------------------------------------

      neq_32 = int(neq, 4)
      n_modes_32 = int(n_modes, 4)
      nvect_32 = int(nvect, 4)
      ltrav_32 = int(ltrav, 4)

      ido_32 = 0
      iparam_32(1) = 1
      iparam_32(3) = int(itermax, 4)
c      iparam_32(7) = 3
      iparam_32(7) = 1
      info_32 = 0
ccccccccccccccccccc
c
      compteur = 0

c      ----------------------------------------------------
c     | Boucle principale en mode de communication inverse |
c      ----------------------------------------------------
c
c
20    continue
c
c     Test for dimesnion conditions in znaupd (err code = -3)
c     Test for N=neq_32, NEV=n_modes_32, NCV=nvect_32
c     Need 0<n_modes_32<neq_32-1, 1<= nvect_32-n_modes_32, nvect_32<=neq_32
      if ((neq_32-1 .le. n_modes_32) .or. (nvect_32-n_modes_32 .lt. 1)
     *    .or.  nvect_32 > neq_32) then
        write(emsg,'(A,A)') 'ARPACK eigensolver dimensional'//
     *   ' conditions failed (would generate ARPACK znaupd error ' //
     *   'code of -3).' // NEW_LINE('A'),
     *   'You should probably increase the grid resolution.'
         errno = -106
          return
      endif

      call znaupd (ido_32, bmat, neq_32, which, n_modes_32, tol,
     *             resid, nvect_32, vschur, neq_32, iparam_32,
     *             ipntr_32, workd, trav, ltrav_32, rwork, info_32)
c
      compteur = compteur + 1

c      if (ido_32.eq.-1) then
         if (ido_32 .eq. -1 .or. ido_32 .eq. 1) then

c      ------------------------------------------------------
c     | On execute  y <--- OP*x = inv[A-SIGMA*M]*M*x         |
c     | pour obtenir un vecteur de depart dans l'image de OP |
c     | x = workd(ipntr_32(1)) et y = workd(ipntr_32(2))           |
c      ------------------------------------------------------

         call zcopy(neq_32, workd(ipntr_32(1)), 1,vect1, 1)
         call z_mxv_csc (neq, vect1, vect2, nonz, row_ind,
     *     col_ptr, mat2)
c
         do i=1,neq
           rhs_re(i) = dble(vect2(i))
           rhs_im(i) = imag(vect2(i))
         enddo
c
c       solve Ax=b, without iterative refinement
        sys = 0
        call umf4zsol (sys, lhs_re, lhs_im, rhs_re, rhs_im,
     *     numeric, control, info_umf)
        if (info_umf (1) .lt. 0) then
            write(ui,*) 'Error occurred in umf4zsol: ', info_umf (1)
            emsg = 'Error occurred in umf4zsol: '
            errno = -107
            return
        endif
        do i=1,neq
          vect2 (i) = dcmplx (lhs_re (i), lhs_im (i))
        enddo
c
         call zcopy(neq_32, vect2, 1, workd(ipntr_32(2)), 1)
         go to 20
c
         else if (ido_32.eq.2) then
c
         write(ui,*) 'VALPR_64: ATTENTION ido_32 = ', ido_32
         write(ui,*) 'check the results...'

c           ----------------------------------------------
c          | On execute  y <--- M*x                       |
c          | x = workd(ipntr_32(1))  et  y = workd(ipntr_32(2)) |
c           ----------------------------------------------

            call zcopy(neq_32, workd(ipntr_32(1)), 1, vect1, 1)
            call z_mxv_csc (neq, vect1, vect2, nonz, row_ind,
     *        col_ptr, mat2)
c
         do i=1,neq
           rhs_re(i) = dble(vect2(i))
           rhs_im(i) = imag(vect2(i))
         enddo
c
c       solve Ax=b, without iterative refinement
        sys = 0
        call umf4zsol (sys, lhs_re, lhs_im, rhs_re, rhs_im,
     *     numeric, control, info_umf)
        if (info_umf (1) .lt. 0) then
            write(ui,*) 'Error occurred in umf4zsol: ', info_umf (1)
            emsg = 'Error occurred in umf4zsol: '
            errno = -108
            return
        endif
        do i=1,neq
          vect2 (i) = dcmplx (lhs_re (i), lhs_im (i))
        enddo
            call zcopy(neq_32, vect2,1, workd(ipntr_32(2)), 1)
            go to 20

         end if

c      ---------------------------------------------------
c     | Either we have convergence, or there is an error. |
c      ---------------------------------------------------

      n_conv = iparam_32(5)

      if (info_32 .gt. 0) then
        write(ui,*)
        write(ui,*) "VALPR_64: info_32 != 0 : ", info_32
        write(ui,*) "VALPR_64: iparam_32(5)=", iparam_32(5), n_modes_32
        write(ui,*) "VALPR_64: number of converged values = ",
     *                iparam_32(5)
        write(ui,*)
      endif

      if (info_32.lt.0) then

c      ---------------------------------------------------
c     | Error message, check the documentation in DNAUPD. |

c  INFO    Integer.  (INPUT/OUTPUT)
c          If INFO .EQ. 0, a randomly initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          Error flag on output.
c          =  0: Normal exit.
c          =  1: Maximum number of iterations taken.
c                All possible eigen_modesues of OP has been found. IPARAM(5)
c                returns the number of wanted converged Ritz values.
c          =  2: No longer an informational error. Deprecated starting
c                with release 2 of ARPACK.
c          =  3: No shifts could be applied during a cycle of the
c                Implicitly restarted Arnoldi iteration. One possibility
c                is to increase the size of NCV relative to NEV.
c                See remark 4 below.
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV-NEV >= 2 and less than or equal to N.
c          = -4: The maximum number of Arnoldi update iteration
c                must be greater than zero.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work array is not sufficient.
c          = -8: Error return from LAPACK eigen_modesue calculation;
c          = -9: Starting vector is zero.
c          = -10: IPARAM(7) must be 1,2,3,4.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
c          = -12: IPARAM(1) must be equal to 0 or 1.
c          = -9999: Could not build an Arnoldi factorization.
c                   IPARAM(5) returns the size of the current Arnoldi
c                   factorization.

c      ---------------------------------------------------

         write(emsg, '(A,I5,/,A)') 'Error occurred in _naupd:'//
     *       ' ARPACK error code = ', info_32,
     *       ' You should probably increase the grid resolution.'

         write(*,*)
         return

      else

c      -------------------------------------
c     | Ici on recupere les valeurs propres |
c      -------------------------------------

      rvec = .true.
      shift2 = (0.0d0,0.0d0)

      call zneupd (rvec, 'A', selecto, nu_out, vschur, neq_32, shift2,
     *  workev, bmat, neq_32, which, n_modes_32, tol,
     *  resid, nvect_32, vschur, neq_32, iparam_32, ipntr_32,
     *  workd, trav, ltrav_32, rwork, ierr_32)
c      ------------------------------------------------------------
c     | La partie reelle d'une valeur propre se trouve dans la     |
c     | premiere colonne du tableau D, la partie imaginaire est    |
c     | dans la seconde.                                           |
c     | Les vecteurs propres sont dans les premieres n_modes_32 colonnes |
c     | du tableau V, lorsque demande (rvec=.true.). Sinon, on y   |
c     | trouve une base orthogonale de l'espace propre.            |
c      ------------------------------------------------------------

      if (ierr_32.ne.0) then
         write(emsg,*) 'VALPR_64:' ,
     *     ' Error with _neupd, info_32 = ', ierr_32,
     *     ' Check the documentation of zneupd. ',
     *     ' This error can occur if the mesh is too coarse.'
         errno = -109
         return

      else
           do i = 1, n_modes
             do j = 1, neq
               vp(j,i) = vschur(j,i)
             enddo
           enddo
         endif
      endif
c
c     free the numeric factorization
      call umf4zfnum (numeric)
c
c
c      if (debug .eq. 1) then
c        do i=1,n_modes
c          write (*,*) "i, nu_out(i) = ", i, nu_out(i)
c        enddo
c      endif
c
      deallocate(workev, rwork, STAT=alloc_stat)
      deallocate(selecto)
c
      return
      end

c------------------------------------------------------------------------
c------------------------------------------------------------------------

      subroutine report_stats_umf4zsym(debug, show_mem_est, info_umf)

      implicit none

      double precision info_umf (90)
      integer(8) ui, debug, show_mem_est

      ui = 6

      if (debug .eq. 1 .or. show_mem_est .eq. 1) then
        write(ui,80) info_umf (1), info_umf (16),
     $      (info_umf (21) * info_umf (4)) / 2**20,
     $      (info_umf (22) * info_umf (4)) / 2**20,
     $      info_umf (23), info_umf (24), info_umf (25)
80      format ('  symbolic analysis:',/,
     $      '   status:  ', f5.0, /,
     $      '   time:    ', e10.2, ' (sec)'/,
     $      '   estimates (upper bound) for numeric LU:', /,
     $      '   size of LU:    ', f12.2, ' (MB)', /,
     $      '   memory needed: ', f12.2, ' (MB)', /,
     $      '   flop count:    ', e12.2, /
     $      '   nnz (L):       ', f12.0, /
     $      '   nnz (U):       ', f12.0)
      endif

      end

c------------------------------------------------------------------------
c------------------------------------------------------------------------

      subroutine report_stats_umf4znum(debug, show_mem_est, info_umf)

      implicit none

      double precision info_umf (90)
      integer(8) ui, debug, show_mem_est

      ui = 6

      if (debug .eq. 1 .or. show_mem_est .eq. 1) then

        write(ui,90) info_umf (1), info_umf (66),
     $      (info_umf (41) * info_umf (4)) / 2**20,
     $      (info_umf (42) * info_umf (4)) / 2**20,
     $      info_umf (43), info_umf (44), info_umf (45)
90      format ('  numeric factorization:',/,
     $      '   status:  ', f5.0, /,
     $      '   time:    ', e10.2, /,
     $      '   actual numeric LU statistics:', /,
     $      '   size of LU:    ', f12.2, ' (MB)', /,
     $      '   memory needed: ', f12.2, ' (MB)', /,
     $      '   flop count:    ', e12.2, /
     $      '   nnz (L):       ', f12.0, /
     $      '   nnz (U):       ', f12.0)

      endif

      end
