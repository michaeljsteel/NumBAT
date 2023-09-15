c*******************************************************
c
c     conv_gmsh: covert the Gmsh .geo to mesh format and 
c     then to the NumBAT .mail mesh format
c
c*******************************************************
c
      subroutine conv_gmsh(geoname, err_no, err_msg)
c
      implicit none
c
      character geoname*1024
      integer err_no
      character err_msg*1024

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
      integer ne_d1, ne_d2                                 
! node/edge codes for each elt
      integer nu_d1(3,max_n_elts), nu_d2(6,max_n_elts)     
! material index for each elt
      integer typ_el_d1(max_n_elts), typ_el_d2(max_n_elts) 

      integer idfn(max_n_pts)


! Individual nodes, number and position
      integer v_ipts(max_n_pts)       
      double precision x(max_n_pts),  y(max_n_pts)

!Elements, number and material index
      integer v_elt_type(max_n_elts)   
      integer v_ielts(max_n_elts)
c
      integer number_tags, material_tag
      integer i, j, k
      integer dummy(10)
      integer gmsh_type_line2node, gmsh_type_triang6node
c
      double precision tmp1, tmp2, tmp3
c
      integer debug, ui, namelength2
      double precision stime1, stime2
      double precision ctime1, ctime2
      character str_in*1024
      character file_ui*10100
      integer sysret
C
Cf2py intent(in) geoname
Cf2py intent(out) err_no, err_msg
c
ccccccccccccccccccccccccc
c
      debug = 0
      ui = 66
      gmsh_version = 2

      err_no = 0


      namelength2 = len_trim(geoname)
      if (namelength2 .ge. 1024) then
        write(err_msg, *) "Name of .geo file is too long extend in ",
     *  "conv_gmsh_py.f"
        err_no = 1
      endif
C
      file0_mesh = geoname(1:namelength2)//".geo"
      file1_mesh = geoname(1:namelength2)//".msh"
      file2_mesh = geoname(1:namelength2)//".mail"
      file_ui    = geoname(1:namelength2)//".log"
C


C     Run Gmsh to get first stage conversion:  .geo to .msh
      call make_msh_file(file0_mesh, file1_mesh, sysret)
      if (sysret .gt. 0) then
        err_no = sysret
        return
      endif
 


c     Second stage conversion:  .msh to .mail
      i_sym = 0
      call get_clocks(stime1, ctime1)


c   Initialisation
      i_mesh(1) = -1
      i_mesh(2) = -1
      i_mesh(3) = -1
c
c     Check size of .msh file
      open (unit=24,file=file1_mesh)
        if(gmsh_version .eq. 2) then  ! what is alternative to v2 ?
          read(24,'(a1)') str_in       ! $MeshFormat
          read(24,'(a1)') str_in       ! $Version string
          read(24,'(a1)') str_in       ! $EndMeshFormat
        endif
        read(24,'(a5)') str_in         ! $Nodes
        read(24,*) n_pts                ! $Number of nodes
      close(24)
c
      if(max_n_pts .lt. n_pts) then ! TODO: pass this back as message
c       open (unit=ui,file=file_ui)
        write(err_msg,*) 'Grid too large: The requested number of ',
     *   'nodes ', n_pts, 
     *   'exceeds the maximum allowed of ', max_n_pts, ' nodes.',
     *   'Try reducing the lc_* grid resolution parameters.'
        err_no = -1
c       close(ui)
        return 
      endif
c
c      write(*,*) 'n_pts = ', n_pts
c
c     Get nodes and positions
      open (unit=24,file=file1_mesh)
        if(gmsh_version .eq. 2) then  ! what is alternative to v2 ?
          read(24,'(a1)') str_in       ! $MeshFormat
          read(24,'(a1)') str_in       ! $Version string
          read(24,'(a1)') str_in       ! $EndMeshFormat
        endif
        read(24,'(a5)') str_in         ! $Nodes
        read(24,*) j                   ! $Number of nodes (already in n_pts)

c     Maps gmsh node number to our sequence of node numbers 1..n_pts
c     Seems like v_ipts just ends up as trivial mapping 1..n_pts, but perhaps this is not guaranteed by gmsh
      do i=1,n_pts
        read(24,*) j, x(i), y(i), tmp1 ! Node number, xval, yval
        v_ipts(j) = i
        if (v_ipts(j) .ne. j)  then ! REMOVE ME
        write(*,*) 'v_ipts misalignment'
        endif 
      enddo
c
      read(24,'(a5)') str_in           ! $EndNodes
      read(24,'(a5)') str_in           ! $Elements
      read(24,*) n_elts                ! Number of elements
c
      if(max_n_elts .lt. n_elts) then

        write(err_msg,*) 'Too many FEM elts: max_n_eltslts < n_elts',
     *    max_n_elts, n_elts
        err_no = -2
        return
      endif
c
c     Read array of elements:  Index EltType(Material code) Number_of_tags <tag>  node-number-list (3 or 6 nodes)
      do i=1,n_elts
        read(24,*) j, v_elt_type(i)
        v_ielts(j) = i
      enddo
c
      close(24)
c
      open (unit=25,file=file1_mesh)
c
c
c     Skip the lines already treated
        if(gmsh_version .eq. 2) then
          read(25,'(a1)') str_in    !5 line header
          read(25,'(a1)') str_in
          read(25,'(a1)') str_in
        endif
        read(25,'(a5)') str_in
        read(25,*) j
      do i=1,n_pts                       !Nodes and positions
        read(25,*) j, tmp1, tmp2, tmp3
      enddo
      read(25,'(a5)') str_in             !3 line Elements header
      read(25,'(a5)') str_in
      read(25,*) j
c
c     Decode element tags     
c     Details at: https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format

      ! Set number of expected tags in Element lines
      if(gmsh_version .eq. 2) then
c       formerly 6 on windows gmsh 2.5.0
        number_tags = 5   ! Number of tags before fields defining node numbers in an element
        material_tag = 4  ! Which field identifies physical curve or physical surface of an element (associated with a particular material)
      else
        number_tags = 5
        material_tag = 3
      endif
c

      ! these elements are edges at the outer boundary and so are associated with one material
      ! line format is  index "8" "2" PhysCurve_j Line_j Node_index x 3
      gmsh_type_line2node = 8    

      ! these are interior triangles and enclose one material
      ! line format is  index "9" "2" PhysSurfac_j PlaneSurface_j Node_index x 6
      gmsh_type_triang6node = 9  

      ne_d1 = 0    !number of line2nodes
      ne_d2 = 0    !number of triang6nodes

      do i=1,n_elts
        if(v_elt_type(i) .eq. gmsh_type_line2node) then   ! 2nd field is 8 = 3-node second order line
          
          ne_d1 = ne_d1 + 1
          read(25,*) (dummy(k), k=1,number_tags), 
     *      (nu_d1(k,ne_d1), k=1,3) ! Get gmsh node numbers for this element
          do k=1,3
            j = nu_d1(k,ne_d1)
            nu_d1(k,ne_d1) = v_ipts(j)   ! Map to our node numbers (actually the same)
          enddo
          ! TODO: Would expect this to all be the same outer material number but seems not
          typ_el_d1(ne_d1) = dummy(material_tag)   ! Get phys_curve index for this element. 

        elseif(v_elt_type(i) .eq. gmsh_type_triang6node) then  !2nd field is 9 = 6-node second order triangle

          ne_d2 = ne_d2 + 1
          read(25,*) (dummy(k), k=1,number_tags), ! Get gmsh node numbers for this element
     *      (nu_d2(k,ne_d2), k=1,6)
          do k=1,6
            j = nu_d2(k,ne_d2)
            nu_d2(k,ne_d2) = v_ipts(j)  ! Map to our node numbers (actually the same)
          enddo
          typ_el_d2(ne_d2) = dummy(material_tag)  ! Get phys_curve index for this element. 

        else

          write(err_msg,*) 'Unknown gmsh elt type:',
     *       'v_elt_type(i), i = ', v_elt_type(i), i
          err_no = -5
          return 

        endif
      enddo
c
      close(25)
c
      do i=1,n_pts
        idfn(i) = 0
      enddo
c
      ! nu_d1[:, j] = (node#, node#, edge# for 3-node line number j)
      ! nu_d2[:, j] = (node# x3,  node# x 3 for 6-node triangle number j)
      ! Associate the boundary edges with their physical_curve stored in idfn
      ! If idfn(j) = pc_k,  node j lies on an outer boundary on physical_curve pc_k
      do i=1,ne_d1
        j = typ_el_d1(i)

C         v_ipts(nu_d1(1,i))
        k = nu_d1(1,i) 
        idfn(k) = j

C         v_ipts(nu_d1(2,i))
        k = nu_d1(2,i) 
        idfn(k) = j

        k = nu_d1(3,i)
        idfn(k) = j
      enddo

c
c
c
c     i_sym = 0 always, so this is switched off
      if(i_sym .ne. 0) then
        call symmetry(n_pts, ne_d2, 
     *      max_n_elts, max_n_pts, idfn, nu_d2, typ_el_d2, 
     *      x, y, i_sym)
      endif
c
      ! TODO: Figure out what this is doing. 
      ! Some kind of useful reordering?
      call renumerote(n_pts, ne_d2, idfn, nu_d2, x, y, ui)
c
c   

      ! Write the NumBAT format .mail file
      ! Format:
      ! Number_of_nodes  Number_of_6node_triangles
      ! Node_number  x_j   y_j  idfn_j                    ! Node
      ! locations and material for those which are edge points (idfn(i) != 0)
      ! Triangle_number  6x node indices  element_type    ! only for 6-node triangles, not 3-node lines

      open (unit=26,file=file2_mesh)
c
      write(26,*) n_pts, ne_d2  
c
      do i=1,n_pts
        write(26,88) i, x(i), y(i), idfn(i)
      enddo
c
      do i=1,ne_d2
        write(26,*) i, (nu_d2(k,i), k=1,6), typ_el_d2(i)
      enddo
      close(26)
88    format(2x, i7, 2(g25.15), i6)
c

      If(debug .eq. 1) then
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
        write(ui,*) 'Number of elements = ', ne_d2
        write(ui,*) 'The program terminates normally'
        write(ui,*) 'Symmetry code = ', i_sym
        write(ui,*) 'CPU time (sec.)           = ', ctime2-ctime1
        write(ui,*) 'Wall time (sec.)          = ', stime2-stime1
c     *                                   dble(time2-ctime1)/100.0
        close(ui)
      EndIf

      continue

      i_mesh(1) = n_pts
      i_mesh(2) = ne_d2
      i_mesh(3) = gmsh_version

C      stop
      return
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
