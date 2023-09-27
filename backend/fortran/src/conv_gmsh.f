#include "numbat_decl.h"
c*******************************************************
c
c     conv_gmsh: covert the Gmsh .geo to mesh format and 
c     then to the NumBAT .mail mesh format
c
c*******************************************************
 
      subroutine conv_gmsh(geoname, errco, emsg)
 
      implicit none
 
      character geoname*1024
      integer errco
      character emsg*1024

      integer i_sym
      integer gmsh_version

      character file0_mesh*1024
      character file1_mesh*1024, file2_mesh*1024

      ! purely for debugging
      integer i_mesh(3)   

      integer n_elts, n_pts
      integer max_n_elts, max_n_pts

      parameter (max_n_pts=250000)
      parameter (max_n_elts=120000)

      ! d1 vars contain 3 node lines (node, node, edge)
      ! d2 vars contain 6 node triangs (node, node, node, edge, edge, edge)
! number of elts found
      integer n_gmsh_lines, n_gmsh_tri
! node/edge codes for each elt
      integer v_gmsh_line_nodes(3,max_n_elts)
      integer v_elt_nodes(6,max_n_elts)     
! material index for each elt
      integer typ_el_d1(max_n_elts), typ_el_d2(max_n_elts) 

      integer v_nd_imat(max_n_pts)


! Individual nodes, number and position
      integer v_ipts(max_n_pts)       
      double precision vx(max_n_pts),  vy(max_n_pts)

!Elements, number and material index
      integer v_elt_type(max_n_elts)   
      integer v_ielts(max_n_elts)
c
      integer n_tags, material_tag
      integer i, j, k
      integer dummy(10)
      integer GMSH_TYPE_LINE2NODE, GMSH_TYPE_TRIANG6NODE
      parameter(GMSH_TYPE_LINE2NODE=8)
      parameter(GMSH_TYPE_TRIANG6NODE=9)
c
      double precision tmp1, tmp2, tmp3
c
      integer debug, ui, namelength2
      double precision stime1, stime2
      double precision ctime1, ctime2
      character str_in*1024
      character file_ui*10100
      integer ui_in, ui_out 
      integer sysret
      integer  imat, nd
C
Cf2py intent(in) geoname
Cf2py intent(out) errco, emsg
c
ccccccccccccccccccccccccc
c
      debug = 0
      ui = 66
      gmsh_version = 2

      errco = 0


      namelength2 = len_trim(geoname)
      if (namelength2 .ge. 1024) then
        write(emsg, *) "Name of .geo file is too long extend in ",
     *  "conv_gmsh_py.f"
        errco = 1
      endif
C
      file0_mesh = geoname(1:namelength2)//".geo"
      file1_mesh = geoname(1:namelength2)//".msh"
      file2_mesh = geoname(1:namelength2)//".mail"
      file_ui    = geoname(1:namelength2)//".log"
C

      ! Stage 1 conversion:  .geo to .msh
      ! Run gmsh by system call
      call make_msh_file(file0_mesh, file1_mesh, sysret)
      if (sysret .gt. 0) then
        errco = sysret
        return
      endif
 

      !Second stage conversion:  .msh to .mail
      i_sym = 0
      call get_clocks(stime1, ctime1)


c   Initialisation
      i_mesh(1) = -1
      i_mesh(2) = -1
      i_mesh(3) = -1
 
c     Check size of .msh file
      ui_in = 24
      open (unit=ui_in,file=file1_mesh)
        if(gmsh_version .eq. 2) then  ! what is alternative to v2 ?
          read(ui_in,'(a1)') str_in       ! $MeshFormat
          read(ui_in,'(a1)') str_in       ! $Version string
          read(ui_in,'(a1)') str_in       ! $EndMeshFormat
        endif
        read(ui_in,'(a5)') str_in         ! $Nodes
        read(ui_in,*) n_pts                ! $Number of nodes
c     close(ui_in)
c
      if(max_n_pts .lt. n_pts) then ! 
        write(emsg,*) 'Grid too large: The requested number of ',
     *   'nodes ', n_pts, 
     *   'exceeds the maximum allowed of ', max_n_pts, ' nodes.',
     *   'Try reducing the lc_* grid resolution parameters.'
        errco = -1
        close(ui_in)
        return 
      endif

c
c      et nodes and positions
c     open (unit=24,file=file1_mesh)
c       if(gmsh_version .eq. 2) then  ! what is alternative to v2 ?
c         read(24,'(a1)') str_in       ! $MeshFormat
c         read(24,'(a1)') str_in       ! $Version string
c         read(24,'(a1)') str_in       ! $EndMeshFormat
c       endif
c       read(24,'(a5)') str_in         ! $Nodes
c       read(24,*) j                   ! $Number of nodes (already in n_pts)

c     Maps gmsh node number to our sequence of node numbers 1..n_pts
c     Seems like v_ipts just ends up as trivial mapping 1..n_pts, but perhaps this is not guaranteed by gmsh
      do i=1,n_pts
        read(ui_in,*) j, vx(i), vy(i), tmp1 ! Node number, xval, yval
        v_ipts(j) = i
        if (v_ipts(j) .ne. j)  then ! REMOVE ME
        write(emsg,*) 'v_ipts misalignment in conv_gmsh'
        errco = -22
        close(ui_in)
        endif 
      enddo
 
      read(ui_in,'(a5)') str_in           ! $EndNodes


      ! Read elements 
      read(ui_in,'(a5)') str_in           ! $Elements
      read(ui_in,*) n_elts                ! Number of elements
 
      if(max_n_elts .lt. n_elts) then
        write(emsg,*) 'Too many FEM elts: max_n_eltslts < n_elts',
     *    max_n_elts, n_elts
        errco = -2
        close(ui_in)
        return
      endif
 
c     Read array of elements:  Index EltType(8 or 9) Number_of_tags <tag>  node-number-list (3 or 6 nodes)
      do i=1,n_elts
        read(ui_in,*) j, v_elt_type(i)
        v_ielts(j) = i   ! map gmsh index to our index (if they are ever different?!) TODO: this is never used, REMOVE
      enddo
c
      close(ui_in)


 
      ! these elements are edges at the outer boundary and so are associated with one material
      ! line format is  index "8" "2" PhysCurve_j Line_j Node_index x 3
      ! GMSH_TYPE_LINE2NODE = 8    

      ! these are interior triangles and enclose one material
      ! line format is  index "9" "2" PhysSurfac_j PlaneSurface_j Node_index x 6
      ! GMSH_TYPE_TRIANG6NODE = 9  

      open (unit=ui_in,file=file1_mesh)
 
 
c     Skip the lines already treated
      if(gmsh_version .eq. 2) then
          read(ui_in,'(a1)') str_in    !5 line header
          read(ui_in,'(a1)') str_in
          read(ui_in,'(a1)') str_in
      endif
      read(ui_in,'(a5)') str_in
      read(ui_in,*) tmp1                    !Num of points
      do i=1,n_pts                       !Nodes and positions
        read(ui_in,*) j, tmp1, tmp2, tmp3
      enddo
      read(ui_in,'(a5)') str_in             !3 line Elements header
      read(ui_in,'(a5)') str_in
      read(ui_in,*) tmp1
 
c     Decode element tags     
c     Details at: https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format

      ! Set number of expected tags in Element lines
      if(gmsh_version .eq. 2) then
c       formerly 6 on windows gmsh 2.5.0
        n_tags = 5   ! Number of tags before fields defining node numbers in an element
        material_tag = 4  ! Which field identifies physical curve or physical surface of an element (associated with a particular material)
      else
        n_tags = 5
        material_tag = 3
      endif
c

      n_gmsh_lines = 0    !number of line2nodes
      n_gmsh_tri = 0    !number of triang6nodes

      do i=1,n_elts
        if(v_elt_type(i) .eq. GMSH_TYPE_LINE2NODE) then   ! 2nd field is 8 = 3-node second order line
          
          n_gmsh_lines = n_gmsh_lines + 1

          read(ui_in,*) (dummy(k), k=1,n_tags), 
     *      (v_gmsh_line_nodes(k,n_gmsh_lines), k=1,3) ! Get gmsh node numbers for this element

          do k=1,3
            j = v_gmsh_line_nodes(k,n_gmsh_lines)
            v_gmsh_line_nodes(k,n_gmsh_lines) = v_ipts(j)   ! Map to our node numbers (actually the same)
          enddo
          ! TODO: Would expect this to all be the same outer material number but seems not
          typ_el_d1(n_gmsh_lines) = dummy(material_tag)   ! Get phys_curve index for this element. 

        elseif(v_elt_type(i) .eq. GMSH_TYPE_TRIANG6NODE) then  !2nd field is 9 = 6-node second order triangle

          n_gmsh_tri = n_gmsh_tri + 1
          read(ui_in,*) (dummy(k), k=1,n_tags), ! Get gmsh node numbers for this element
     *      (v_elt_nodes(k,n_gmsh_tri), k=1,6)

          do k=1,6  ! this loop seems to be a no-op in that v_ipts(j)=j afaict
            j = v_elt_nodes(k,n_gmsh_tri)
            v_elt_nodes(k,n_gmsh_tri) = v_ipts(j)  ! Map to our node numbers (actually the same)
          enddo
          typ_el_d2(n_gmsh_tri) = dummy(material_tag)  ! Get phys_surface index for this element. 

        else

          write(emsg,*) 'Unknown gmsh elt type:',
     *       'v_elt_type(i), i = ', v_elt_type(i), i
          errco = -5
          return 

        endif
      enddo
 
      close(ui_in)
 
      ! Now we have:
      !   v_gmsh_line_nodes[:, j] = (node#, node#, edge# for 3-node line number j)
      !   v_elt_nodes[:, j] = (node# x3,  node# x 3 for 6-node triangle number ! j) , or is it 3 nodes, 3 edges?
      !  typ_el_d1[] = material next to boundary line 
      !  typ_el_d2[] = material inside triangle

      ! Next, associate the boundary edges with their physical_curve stored in v_nd_imat
      ! If v_nd_imat(j) = pc_k != 0,  node j lies on an outer boundary on physical_curve pc_k

      do i=1,n_pts
        v_nd_imat(i) = 0
      enddo

      do i=1,n_gmsh_lines  ! for each line element
        imat = typ_el_d1(i)   ! adjacent material 

C         v_ipts(v_gmsh_line_nodes(1,i))
        nd = v_gmsh_line_nodes(1,i)      
        v_nd_imat(nd) = imat                ! assign this material to a node on the line

C         v_ipts(v_gmsh_line_nodes(2,i))
        nd = v_gmsh_line_nodes(2,i) 
        v_nd_imat(nd) = imat

        nd = v_gmsh_line_nodes(3,i)
        v_nd_imat(nd) = imat
      enddo


c     i_sym = 0 always, so this is switched off
      if(i_sym .ne. 0) then
        call symmetry(n_pts, n_gmsh_tri, 
     *    max_n_elts, max_n_pts, v_nd_imat, v_elt_nodes, 
     *    typ_el_d2, vx, vy, i_sym)
      endif
c
      ! TODO: Figure out what this is doing. 
      ! Some kind of useful reordering?
      call renumerote(n_pts, n_gmsh_tri, v_nd_imat, 
     * v_elt_nodes, vx, vy, errco, emsg)
      RETONERROR(errco)
c
c   

      ! Write the NumBAT format .mail file
      ! Format:
      ! Number_of_nodes  Number_of_6node_triangles
      ! Node_number  x_j   y_j  v_nd_imat_j                    ! Node
      ! locations and material for those which are edge points (v_nd_imat(i) != 0)
      ! Triangle_number  6x node indices  element_type    ! only for 6-node triangles, not 3-node lines

      ui_out = 26
      open (unit=ui_out,file=file2_mesh)
 
      write(ui_out,*) n_pts, n_gmsh_tri  
 
      do i=1,n_pts
        write(ui_out, '(2x, i7, 2(g25.15), i6)') i, vx(i), vy(i), 
     *        v_nd_imat(i)
      enddo
 
      do i=1,n_gmsh_tri
        write(ui_out,*) i, (v_elt_nodes(k,i), k=1,6), 
     *         typ_el_d2(i)
      enddo
      close(ui_out)
 

      if(debug .eq. 1) then
C         temps initial (sert pour la durree des calcul)
        call get_clocks(stime2, ctime2)
        open (unit=ui,file=file_ui)
        write(ui,*) "conv_gmsh_m: debug = ", debug
        write(ui,*) "gmsh_version = ", gmsh_version
        write(ui,*) "i_mesh = ", i_mesh
        write(ui,*) " file1_mesh = ", file1_mesh
        write(ui,*) " file2_mesh = ", file2_mesh
        write(ui,*) " i_sym = ", i_sym
        write(ui,*) 'Number of points = ', n_pts
        write(ui,*) 'Number of elements = ', n_elts
        write(ui,*) 'The program terminates normally'
        write(ui,*) 'Symmetry code = ', i_sym
        write(ui,*) 'CPU time (sec.)           = ', ctime2-ctime1
        write(ui,*) 'Wall time (sec.)          = ', stime2-stime1
c     *                                   dble(time2-ctime1)/100.0
        close(ui)
      endif

      continue  !! TODO: where is this from 

      i_mesh(1) = n_pts
      i_mesh(2) = n_elts
      i_mesh(3) = gmsh_version

      end


c##################################################################################
      subroutine make_msh_file(file0_mesh, file1_mesh, sysret)

      character file0_mesh*1024, file1_mesh*1024
      character com_line*1024
      integer sysret

      com_line = "gmsh -0 -2  -order 2 -v 0 -o " // 
     *    trim(file1_mesh) // " " // trim(file0_mesh)

      call system(com_line, sysret)

      return 
      end
