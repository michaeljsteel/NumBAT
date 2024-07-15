#include "numbat_decl.h"

module conv_gmsh_interfaces

!     interface
!         subroutine parse_msh_file(fname_msh, gmsh_version, n_pts, n_elts, v_ipts, &
!                 v_ielts, vx, vy, v_gmsh_elt_type, errco, emsg)

!             use numbatmod
!             implicit none

!             character(len=*), intent(in) :: fname_msh
!             integer(8) gmsh_version

!             integer(8) n_pts, n_elts
!             integer(8) errco
!             character(len=EMSG_LENGTH) emsg
!             integer(8) v_ipts(MAX_N_PTS)
!             integer(8) v_ielts(MAX_N_ELTS), v_gmsh_elt_type(MAX_N_ELTS)

!             !double precision vx(MAX_N_PTS), vy(MAX_N_PTS)
!             double precision, dimension(:), allocatable, intent(inout) :: vx, vy

!         end subroutine

! end interface

contains


   subroutine parse_msh_file(fname_msh, gmsh_version, n_pts, n_elts, &
      vx, vy, v_ipts, v_ielts,  v_gmsh_elt_type, errco, emsg)

      use numbatmod
      use alloc


      character(len=*), intent(in) :: fname_msh
      integer(8) gmsh_version

      integer(8) n_pts, n_elts
      integer(8) errco
      character(len=EMSG_LENGTH) emsg

      !double precision vx(MAX_N_PTS), vy(MAX_N_PTS)
      double precision, dimension(:), allocatable, intent(inout) :: vx, vy
      integer(8), dimension(:), allocatable, intent(inout) :: v_ipts, v_ielts, v_gmsh_elt_type
      !integer(8) v_ipts(MAX_N_PTS)
      !integer(8) v_ielts(MAX_N_ELTS), v_gmsh_elt_type(MAX_N_ELTS)

      character str_in*(FNAME_LENGTH)
      integer(8) ui_in, i, j, tmp1

      !  Check size of .msh file
      ui_in = 24
      open (unit=ui_in,file=fname_msh)

      if (gmsh_version .eq. 2) then  !  what is alternative to v2 ?
         read(ui_in,'(a1)') str_in   !  $MeshFormat
         read(ui_in,'(a1)') str_in   !  $Version string
         read(ui_in,'(a1)') str_in   !  $EndMeshFormat
      endif

      read(ui_in,'(a5)') str_in      !  $Nodes
      read(ui_in,*) n_pts            !  $Number of nodes

      if (MAX_N_PTS .lt. n_pts) then
         write(emsg,*) 'Grid too large: The requested number of ', &
            'nodes ', n_pts, &
            'exceeds the maximum allowed of ', MAX_N_PTS, ' nodes.', &
            'Try reducing the lc_* grid resolution parameters.'
         errco = NBERROR_111
         close(ui_in)
         return
      endif


      call double_alloc_1d(vx, n_pts, 'vx', errco, emsg); RETONERROR(errco)
      call double_alloc_1d(vy, n_pts, 'vy', errco, emsg); RETONERROR(errco)
      call integer_alloc_1d(v_ipts, n_pts, 'v_ipts', errco, emsg); RETONERROR(errco)

      !  Seems like v_ipts just ends up as trivial mapping 1..n_pts,
      !  but perhaps this is not guaranteed by gmsh?
      do i=1,n_pts
         read(ui_in,*) j, vx(i), vy(i), tmp1 !  Node number, xval, yval
         v_ipts(j) = i

         if (v_ipts(j) .ne. j)  then !  REMOVE ME
            write(emsg,*) 'v_ipts misalignment in conv_gmsh'
            errco = NBERROR_112
            close(ui_in)
            return
         endif

      enddo

      read(ui_in,'(a5)') str_in      !  $EndNodes


      !  Read elements
      read(ui_in,'(a5)') str_in      !  $Elements
      read(ui_in,*) n_elts           !  Number of elements

      if(MAX_N_ELTS .lt. n_elts) then
         write(emsg,*) 'Too many FEM elts: MAX_N_ELTSlts < n_elts', &
            MAX_N_ELTS, n_elts
         errco = NBERROR_114
         close(ui_in)
         return
      endif


      call integer_alloc_1d(v_gmsh_elt_type, n_elts, 'v_gmsh_elt_type', errco, emsg); RETONERROR(errco)
      call integer_alloc_1d(v_ielts, n_elts, 'v_ielts', errco, emsg); RETONERROR(errco)


      !  Read array of elements:
      !  Index EltType(8 or 9) Number_of_tags <tag>  node-number-list (3 or 6 nodes)
      do i=1,n_elts
         read(ui_in,*) j, v_gmsh_elt_type(i)
         v_ielts(j) = i   !  map gmsh index to our index (if they are ever different?!) TODO: this is never used, REMOVE
      enddo

      close(ui_in)

      !  confirm that v_ipts and v_ielts are trivial
      do i=1,n_pts
         if (v_ipts(i) .ne. i) then
            write(emsg,*) 'v_ipts(i) mismatch at i=', i
            errco = NBERROR_115
            return
         endif
      enddo

      do i=1,n_elts
         if (v_ielts(i) .ne. i) then
            write(emsg,*) 'v_ielts(i) mismatch at i=', i
            errco = NBERROR_116
            return
         endif
      enddo

   end subroutine





end module conv_gmsh_interfaces
