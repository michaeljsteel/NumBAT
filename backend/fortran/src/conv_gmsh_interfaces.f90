#include "numbat_decl.h"

module conv_gmsh_interfaces


contains


   subroutine parse_msh_file(fname_msh, gmsh_version, n_msh_pts, n_msh_elts, &
      vx, vy, v_i_mshpts,  v_gmsh_elt_type, nberr)

      use numbatmod
      use alloc

      type(NBError) nberr

      character(len=*), intent(in) :: fname_msh
      integer(8) gmsh_version

      integer(8) n_msh_pts, n_msh_elts
      character(len=EMSG_LENGTH) emsg

      double precision, dimension(:), allocatable, intent(inout) :: vx, vy
      integer(8), dimension(:), allocatable, intent(inout) :: v_i_mshpts, v_gmsh_elt_type

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
      read(ui_in,*) n_msh_pts            !  $Number of nodes

      if (MAX_N_PTS .lt. n_msh_pts) then
         write(emsg,*) 'Grid too large: The requested number of ', &
            'nodes ', n_msh_pts, &
            'exceeds the maximum allowed of ', MAX_N_PTS, ' nodes.', &
            'Try reducing the lc_* grid resolution parameters.'
            close(ui_in)
            call nberr%set(NBERROR_111, emsg)
         return
      endif


      call double_alloc_1d(vx, n_msh_pts, 'vx', nberr); RET_ON_NBERR(nberr)
      call double_alloc_1d(vy, n_msh_pts, 'vy', nberr); RET_ON_NBERR(nberr)
      call integer_alloc_1d(v_i_mshpts, n_msh_pts, 'v_i_mshpts', nberr); RET_ON_NBERR(nberr)

      !  Seems like v_i_mshpts just ends up as trivial mapping 1..n_msh_pts,
      !  but perhaps this is not guaranteed by gmsh?
      do i=1,n_msh_pts
         read(ui_in,*) j, vx(i), vy(i), tmp1 !  Node number, xval, yval
         v_i_mshpts(j) = i

         if (v_i_mshpts(j) .ne. j)  then !  REMOVE ME
            write(emsg,*) 'v_i_mshpts misalignment in conv_gmsh'
            close(ui_in)
            call nberr%set(NBERROR_112, emsg)
            return
         endif

      enddo

      read(ui_in,'(a5)') str_in      !  $EndNodes

      !  Read elements
      read(ui_in,'(a5)') str_in      !  $Elements
      read(ui_in,*) n_msh_elts           !  Number of elements

      if(MAX_N_ELTS .lt. n_msh_elts) then
         write(emsg,*) 'Too many FEM elts: MAX_N_ELTSlts < n_msh_elts', &
            MAX_N_ELTS, n_msh_elts
         close(ui_in)
         call nberr%set(NBERROR_114, emsg)
         return
      endif


      call integer_alloc_1d(v_gmsh_elt_type, n_msh_elts, 'v_gmsh_elt_type', nberr); RET_ON_NBERR(nberr)
      !call integer_alloc_1d(v_i_elts, n_msh_elts, 'v_i_elts', nberr); RET_ON_NBERR(nberr)


      !  Read array of elements:
      !  Index EltType(8 or 9) Number_of_tags <tag>  node-number-list (3 or 6 nodes)
      do i=1,n_msh_elts
         read(ui_in,*) j, v_gmsh_elt_type(i)
         !v_i_elts(j) = i   !  map gmsh index to our index (if they are ever different?!) TODO: this is never used, REMOVE
      enddo

      close(ui_in)

      ! !  confirm that v_i_mshpts and v_i_elts are trivial
      ! do i=1,n_msh_pts
      !    if (v_i_mshpts(i) .ne. i) then
      !       write(emsg,*) 'v_i_mshpts(i) mismatch at i=', i
      !       errco = NBERROR_115
      !       return
      !    endif
      ! enddo

      ! do i=1,n_msh_elts
      !    if (v_i_elts(i) .ne. i) then
      !       write(emsg,*) 'v_i_elts(i) mismatch at i=', i
      !       errco = NBERROR_116
      !       return
      !    endif
      ! enddo

   end subroutine





end module conv_gmsh_interfaces
