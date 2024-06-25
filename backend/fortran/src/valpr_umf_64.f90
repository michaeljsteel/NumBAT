!     ------------------------------------------------------------------
!
!                    sous-routine VALPR_64.
!                                 ------
!     gere l'utilisation de ARPACK
!
!     ------------------------------------------------------------------
!
!     Parametres d'entree:
!     --------------------
!
!     nvect (I)              : dimension de l'espace de Krylov
!     n_modes (I)               : nombre de valeurs propres desirees.
!     neq (I)                : nombre d'equations
!     workd (DP)             : matrice de travail pour dnaupd,
!                              taille  3 neq
!     resid   (DP)           : vecteur de trvail pour dnaupd,
!                              taille neq
!     v (DP)                 : matrice des vecteurs de Schur,
!                              taille neq*nvect
!     asup,ainf,adiag (DP)   : stockage de la matrice tangente.
!     msup,minf,mdiag (DP)   : stockage de la matrice masse.
!     kld (I)                : vecteur de localisation des debuts de
!                              colonnes
!     vect1,vect2,vect3 (DP) : vecteurs de travail,
!                              taille neq
!     long (I)               : longueur des super-vecteurs pour
!                              les matrices.
!     ddot (DP)              : fonction appelee pour calculer
!                              le produit scalaire
!     trav(DP)               : espace de travail pour dnaupd,
!                              taille ltrav>= 3*nvect^2+6*nvect
!     action                 : 0: redemarrage avec un vecteur aleatoire
!                              1: lire une vecteur initial

!    ls_data:  just timing info

!     Parametres de sortie:
!     ---------------------
!
!     reel (DP)              : parties reelles, taille nvect
!     imag (DP)              : parties imaginaires, taille nvect
!
!     ------------------------------------------------------------------
!
subroutine valpr_64 (i_base, nvect, n_modes, neq, itermax, ltrav,&
&tol, nonz, row_ind, col_ptr, mat1_re, mat1_im, mat2,&
&vect1, vect2, workd, resid, v, d, trav, vp,&
&rhs_re, rhs_im, lhs_re, lhs_im, n_conv, ls_data,&
&numeric, control, info_umf, debug, errco, emsg)
!
!     ------------------------------------------------------------------
!
   use numbatmod
!
   integer(8) neq, nonz, n_conv, i_base
   integer(8) row_ind(nonz), col_ptr(neq+1)
   complex(8) mat2(nonz)
   double precision mat1_re(nonz), mat1_im(nonz)
   double precision rhs_re(neq), rhs_im(neq)
   double precision lhs_re(neq), lhs_im(neq)
   double precision ls_data(10)
!
   double precision time1_fact, time2_fact
   double precision time1, time2
!
   double precision control (20), info_umf (90)
   integer(8) numeric, symbolic, sys
!
   integer(8) itermax, nvect, i, j, ltrav
   integer(8) :: n_modes
   integer(8) compteur
   complex(8) resid(neq), v(neq,nvect), workd(3*neq)
   complex(8) vect1(neq), vect2(neq), trav(ltrav)
   complex(8) d(n_modes+1), shift2, vp(neq,n_modes)
!
   double precision tol
!
!      integer(8) max_nvect
!      parameter(max_nvect=3000) ! previously 1500
!
   integer alloc_stat
!       !  (3*max_nvect),
   complex(8), dimension(:), allocatable :: workev
!       !  (max_nvect)
   double precision, dimension(:), allocatable :: rwork
!       !  (max_nvect)
   logical, dimension(:), allocatable :: select

   integer errco
   character(len=EMSG_LENGTH) emsg


!     Local variables
!     32-bit integers for ARPACK
   integer(4) neq_32, n_modes_32, nvect_32
   integer(4) ido_32, info_32, ierr_32, iparam_32(11)
   integer(4) ipntr_32(14), ltrav_32
!
   logical rvec
   character bmat*1, which*2

!      data bmat/'G'/
!      data which/'SR'/
!      data which/'SM'/
   data bmat/'I'/
   data which/'LM'/
!
   integer(8) ui, debug
!      common/imp/ui, debug
!
!     ------------------------------------------------------------------
!
   ui = stdout

   errco = 0
   emsg = ""

   alloc_stat = 0
   allocate(workev(3*nvect), rwork(nvect), STAT=alloc_stat)
   if (alloc_stat /= 0) then
      write(*,*) "VALPR_64: Mem. allocation is unseccesfull"
      write(*,*) "for the arrays workev, rwork"
      write(*,*) "alloc_stat, nvect = ", alloc_stat, nvect
      write(*,*) "Aborting..."
      stop
   endif

   allocate(select(nvect), STAT=alloc_stat)
   if (alloc_stat /= 0) then
      write(*,*) "VALPR_64: Mem. allocation is unseccesfull"
      write(*,*) "for the array select"
      write(*,*) "alloc_stat, nvect = ", alloc_stat, nvect
      write(*,*) "Aborting..."
      stop
   endif

!
   shift2 = (0.0d0,0.0d0)
!
   do i=1,n_modes
      d(i) = 0.0d0
   enddo
!
!       ##################################################################
!
   if (i_base .ne. 0) then
      write(ui,*) "valpr_64: i_base != 0 : ", i_base
      write(ui,*) "valpr_64: UMFPACK requires 0-based indexing"
      write(ui,*) "valpr_64l: Aborting..."
      stop
   endif
!
!       ----------------------------------------------------------------
!       factor the matrix and save to a file
!       ----------------------------------------------------------------
!
   if (debug .eq. 1) then
      write(ui,*) "valpr_64: factorisation (UMFPACK)"
   endif

   !TODO: replace with stopwatch
   call cpu_time(time1_fact)
   ls_data(1) = time1_fact


   !     umfpack * report status (print level = control(1)) :
   !     print level = 0 or less : No output, even when an error occurs.
   !     print level = 1 (default value) : then error messages are printed,
   !                      and nothing is printed if the status is UMFPACK OK.
   !     print level = 2 or more : then the status is always printed.
   !
   !       print control parameters.  set control (1) to 1 to print
   !       error messages only
   call umf4zdef (control)
   control (1) = 1
   call umf4zpcon (control)

!       pre-order and symbolic analysis
   call umf4zsym (neq, neq, col_ptr, row_ind,&
   &mat1_re, mat1_im, symbolic, control, info_umf)
!
!       print statistics computed so far
!       call umf4zpinf (control, info_umf) could also be done.
   if (debug .eq. 1) then
      write(ui,80) info_umf (1), info_umf (16),&
      &(info_umf (21) * info_umf (4)) / 2**20,&
      &(info_umf (22) * info_umf (4)) / 2**20,&
      &info_umf (23), info_umf (24), info_umf (25)
80    format ('  symbolic analysis:',/,&
      &'   status:  ', f5.0, /,&
      &'   time:    ', e10.2, ' (sec)'/,&
      &'   estimates (upper bound) for numeric LU:', /,&
      &'   size of LU:    ', f12.2, ' (MB)', /,&
      &'   memory needed: ', f12.2, ' (MB)', /,&
      &'   flop count:    ', e12.2, /&
      &'   nnz (L):       ', f12.0, /&
      &'   nnz (U):       ', f12.0)
   endif

!       check umf4zsym error condition
   if (info_umf (1) .lt. 0) then
      write(ui,*) 'Error occurred in umf4zsym: ', info_umf (1)
      stop
   endif

!       numeric factorization
   call umf4znum (col_ptr, row_ind, mat1_re,&
   &mat1_im, symbolic, numeric, control, info_umf)

!       print statistics for the numeric factorization
!       call umf4zpinf (control, info_umf) could also be done.
   if (debug .eq. 1) then
      write(ui,90) info_umf (1), info_umf (66),&
      &(info_umf (41) * info_umf (4)) / 2**20,&
      &(info_umf (42) * info_umf (4)) / 2**20,&
      &info_umf (43), info_umf (44), info_umf (45)
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

!       check umf4znum error condition
   if (info_umf (1) .lt. 0) then
      write(ui,*) 'Error occurred in umf4znum: ', info_umf (1)
      stop
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

!c       free the numeric factorization
!        call umf4zfnum (numeric)
!
   call cpu_time(time2_fact)
   ls_data(2) = time2_fact
   if (debug .eq. 1) then
      write(ui,*) "valpr_64: factorisation completed"
      write(ui,*) "LU factorisation : CPU time = ",&
      &(time2_fact-time1_fact)
! ,
!     *         100*(time2_fact-time1_fact)/(time2-time1),"%"
   endif
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
!     On commence le travail avec znaupd
!     ----------------------------------
!
!      ------------------------------------------------------------
!     | Le parametre IDO_32 est utilise pour la communication.        |
!     | A l'etape initiale il doit valoir 0.                       |
!     | Le choix INFO_32=0 correspond a la construction par           |
!     | Arpack d'un vecteur initial a composantes aleatoires       |
!     | iparam_32(1)=1 signifie que Arpack calcule les translations   |
!     | a partir de la matrice projetee et conformement au critere |
!     | "which".                                                   |
!      ------------------------------------------------------------

   neq_32 = int(neq, 4)
   n_modes_32 = int(n_modes, 4)
   nvect_32 = int(nvect, 4)
   ltrav_32 = int(ltrav, 4)

   ido_32 = 0
   iparam_32(1) = 1
   iparam_32(3) = int(itermax, 4)
!      iparam_32(7) = 3
   iparam_32(7) = 1
   info_32 = 0
!cccccccccccccccccc
!
   compteur = 0

!      ----------------------------------------------------
!     | Boucle principale en mode de communication inverse |
!      ----------------------------------------------------
!
   call cpu_time(time1)
   ls_data(3) = time1
!
20 continue
!
   call znaupd (ido_32, bmat, neq_32, which, n_modes_32, tol,&
   &resid, nvect_32, v, neq_32, iparam_32,&
   &ipntr_32, workd, trav, ltrav_32, rwork, info_32)
!
   compteur = compteur + 1

!      if (ido_32.eq.-1) then
   if (ido_32 .eq. -1 .or. ido_32 .eq. 1) then

!      ------------------------------------------------------
!     | On execute  y <--- OP*x = inv[A-SIGMA*M]*M*x         |
!     | pour obtenir un vecteur de depart dans l'image de OP |
!     | x = workd(ipntr_32(1)) et y = workd(ipntr_32(2))           |
!      ------------------------------------------------------

      call zcopy(neq_32, workd(ipntr_32(1)), 1,vect1, 1)
      call z_mxv_csc (neq, vect1, vect2, nonz, row_ind,&
      &col_ptr, mat2)

      do i=1,neq
         rhs_re(i) = dble(vect2(i))
         rhs_im(i) = imag(vect2(i))
      enddo

!       solve Ax=b, without iterative refinement
      sys = 0
      call umf4zsol (sys, lhs_re, lhs_im, rhs_re, rhs_im,&
      &numeric, control, info_umf)

      if (info_umf (1) .lt. 0) then
         write(ui,*) 'Error occurred in umf4zsol: ', info_umf (1)
         stop
      endif

      do i=1,neq
         vect2 (i) = dcmplx (lhs_re (i), lhs_im (i))
      enddo

      call zcopy(neq_32, vect2, 1, workd(ipntr_32(2)), 1)
      go to 20

   else if (ido_32.eq.2) then

      write(ui,*) 'VALPR_64: ATTENTION ido_32 = ', ido_32
      write(ui,*) 'check the results...'

!           ----------------------------------------------
!          | On execute  y <--- M*x                       |
!          | x = workd(ipntr_32(1))  et  y = workd(ipntr_32(2)) |
!           ----------------------------------------------

      call zcopy(neq_32, workd(ipntr_32(1)), 1, vect1, 1)
      call z_mxv_csc (neq, vect1, vect2, nonz, row_ind,&
      &col_ptr, mat2)
!
      do i=1,neq
         rhs_re(i) = dble(vect2(i))
         rhs_im(i) = imag(vect2(i))
      enddo
!
!       solve Ax=b, without iterative refinement
      sys = 0
      call umf4zsol (sys, lhs_re, lhs_im, rhs_re, rhs_im,&
      &numeric, control, info_umf)
      if (info_umf (1) .lt. 0) then
         write(ui,*) 'Error occurred in umf4zsol: ', info_umf (1)
         stop
      endif
      do i=1,neq
         vect2 (i) = dcmplx (lhs_re (i), lhs_im (i))
      enddo
      call zcopy(neq_32, vect2,1, workd(ipntr_32(2)), 1)
      go to 20

   end if

!      ---------------------------------------------------
!     | Either we have convergence, or there is an error. |
!      ---------------------------------------------------

   n_conv = iparam_32(5)

   if (info_32 .gt. 0) then
      write(ui,*)
      write(ui,*) "VALPR_64: The Arnoldi iteration scheme has failed"
      write(ui,*) "VALPR_64: The znaupd error flag has the value ",&
      &"info_32=", info_32
      if (info_32 .eq. 1) then
         write(ui,*) "VALPR_64: Max iterations exceeded."
         write(ui,*) " Requested eigen_modesues = ", iparam_32(5)
         write(ui,*) " Converged eigen_modesues = ", n_modes_32
         write(ui,*) " You might try:"
         write(ui,*) "  1) Increasing the requested number",&
         &" of eigen_modesues"
         write(ui,*) "  2) Increasing the grid resolution"
      endif
      if (info_32 .eq. 3) then
         write(ui,*) "VALPR_64: Shift could not be applied."
      endif
      write(ui,*) "VALPR_64: For details on znaupd errors see",&
      &" https://www.caam.rice.edu/software/ARPACK/UG/node138.html"
!        write(ui,*) "VALPR_64: iparam_32(5) = ", iparam_32(5), n_modes_32
!        write(ui,*) "VALPR_64: number of converged values = ",
!     *                iparam_32(5)
      write(ui,*)
   endif

   if (info_32.lt.0) then

!      ---------------------------------------------------
!     | Error message, check the documentation in DNAUPD. |

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
!          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
!          = -6: BMAT must be one of 'I' or 'G'.
!          = -7: Length of private work array is not sufficient.
!          = -8: Error return from LAPACK eigen_modesue calculation;
!          = -9: Starting vector is zero.
!          = -10: IPARAM(7) must be 1,2,3,4.
!          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
!          = -12: IPARAM(1) must be equal to 0 or 1.
!          = -9999: Could not build an Arnoldi factorization.
!                   IPARAM(5) returns the size of the current Arnoldi
!                   factorization.

!      ---------------------------------------------------

      write(ui,*) 'VALPR_64:'
      write(ui,*) "VALPR_64: The Arnoldi iteration scheme has failed."
      write(ui,*) "VALPR_64: The znaupd error flag has the value ",&
      &"info_32=", info_32
      write(ui,*) "VALPR_64: For details on znaupd errors see",&
      &" https://www.caam.rice.edu/software/ARPACK/UG/node138.html"
      write(ui,*) 'Aborting...'
      stop

   else

!      -------------------------------------
!     | Ici on recupere les valeurs propres |
!      -------------------------------------

      rvec = .true.

      call zneupd (rvec, 'A', select, d, v, neq_32, shift2,&
      &workev, bmat, neq_32, which, n_modes_32, tol,&
      &resid, nvect_32, v, neq_32, iparam_32, ipntr_32,&
      &workd, trav, ltrav_32, rwork, ierr_32)
!      ------------------------------------------------------------
!     | La partie reelle d'une valeur propre se trouve dans la     |
!     | premiere colonne du tableau D, la partie imaginaire est    |
!     | dans la seconde.                                           |
!     | Les vecteurs propres sont dans les premieres n_modes_32 colonnes |
!     | du tableau V, lorsque demande (rvec=.true.). Sinon, on y   |
!     | trouve une base orthogonale de l'espace propre.            |
!      ------------------------------------------------------------

      if (ierr_32.ne.0) then

!      -----------------------------------------------------
!     | Error condition: Check the documentation of DNEUPD. |
!      -----------------------------------------------------

         write(ui,*) 'VALPR_64:'
         write(ui,*) ' Error with _neupd, info_32 = ', ierr_32
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

   call cpu_time(time2)
   ls_data(4) = time2


   deallocate(workev, rwork, STAT=alloc_stat)
   deallocate(select)

   return
end
