#include "numbat_decl.h"

!*******************************************************

!  conv_gmsh: convert the Gmsh .geo to mesh format and
!  then to the NumBAT .mail mesh format

!*******************************************************

subroutine conv_gmsh_impl(geoname, assertions_on, errco, emsg)

   use numbatmod

   character(len=FNAME_LENGTH) geoname
   integer assertions_on, errco
   character(len=EMSG_LENGTH) emsg

   integer i_sym, gmsh_version

   character(len=FNAME_LENGTH) :: fname_geo, fname_msh, fname_mail
   integer n_elts, n_pts

   !  d1 vars contain 3 node lines (node, node, edge)
   !  d2 vars contain 6 node triangs (node, node, node, edge, edge, edge)
   !  number of elts found
   integer n_gelts_lines, n_gelts_triangs
   !  node/edge codes for each elt
   integer v_lines_nodes(3,MAX_N_ELTS)
   integer v_triang_nodes(6,MAX_N_ELTS)

   !  material index for each elt
   integer v_eltbdy_physcurve(MAX_N_ELTS)
   integer v_eltint_physsurf(MAX_N_ELTS)
   integer v_nd_iphyscurve(MAX_N_PTS)

!  Individual nodes, number and position
   integer v_ipts(MAX_N_PTS)
   double precision vx(MAX_N_PTS),  vy(MAX_N_PTS)

!Elements, number and material index
   integer v_gmsh_elt_type(MAX_N_ELTS)
   integer v_ielts(MAX_N_ELTS)
   integer i, j
   integer namelength2
   double precision stime1, ctime1
   integer sysret
   integer  iphyscurve, nd

   call get_clocks(stime1, ctime1)

   !debug = 0
   gmsh_version = 2
   errco = 0
   i_sym = 0   !  symmetry is off

   namelength2 = len_trim(geoname)
   if (namelength2 .ge. FNAME_LENGTH) then
      write(emsg, *) "Name of .geo file is too long extend in ",&
         "conv_gmsh_py.f"
      errco = 1
   endif

   fname_geo = geoname(1:namelength2)//".geo"
   fname_msh = geoname(1:namelength2)//".msh"
   fname_mail = geoname(1:namelength2)//".mail"
   !fname_log  = geoname(1:namelength2)//".log"

   !  Stage 1 conversion:  .geo to .msh
   call make_msh_file(fname_geo, fname_msh, sysret)
   if (sysret .gt. 0) then
      errco = sysret
      return
   endif

   !Second stage conversion:  .msh to .mail
   call parse_msh_file(fname_msh, gmsh_version, &
      n_pts, n_elts, v_ipts, v_ielts, &
      vx, vy, v_gmsh_elt_type, errco, emsg)
   RETONERROR(errco)

   if (assertions_on .ne. 0) then
      write(*,*) 'Node check fresh from gmsh:'
      call check_point_separations(n_pts, vx, vy, errco, emsg)
   endif

   !  Now we know the number of points and mappings (even if trivial)
   !  Next we load elt data according to the gmsh element types

   call decode_element_tags(fname_msh, gmsh_version, &
      n_pts, n_elts, v_gmsh_elt_type, v_ipts, &
      n_gelts_lines, n_gelts_triangs, v_lines_nodes,  &
      v_triang_nodes, v_eltbdy_physcurve, v_eltint_physsurf,  &
      errco, emsg)
   RETONERROR(errco)

   !  Now we have:
   !  v_lines_nodes[:, j] = (node#, node#, node# for 3-node line number j)
   !  v_triang_nodes[:, j] = (cornernode# x3,  edgenode# x 3 for 6-node triangle number j)
   !  v_eltbdy_physcurve[] = number of Physical_Curve on which elt lies (may or may not be outer boundary ?)
   !  v_eltint_physsurf[] = number of Physical_Surface inside the  domain

   !  Next, associate the three nodes on each boundary elt with their physical_curve stored in v_nd_iphyscurve
   !  If v_nd_iphyscurve(j) = pc_k != 0,  node j lies on PhysicalCurve pc_k

   do i=1,n_pts
      v_nd_iphyscurve(i) = 0
   enddo

   do i=1,n_gelts_lines           !  for each boundary line element
      iphyscurve = v_eltbdy_physcurve(i)   !  associate its 3 nodes with the PhysicalCurve it is on

      do j=1,3
         nd = v_lines_nodes(j,i)
         v_nd_iphyscurve(nd) = iphyscurve                !  assign this material to a node on the line
      enddo
   enddo


!  i_sym = 0 always, so this is switched off  TODO: the call to symmetry is very broken with parameteres
   !if(i_sym .ne. 0) then

      !call symmetry(n_pts, n_gelts_triangs, &
       !  MAX_N_ELTS, MAX_N_PTS, v_nd_iphyscurve, v_triang_nodes, &
        !  v_eltint_physsurf, vx, vy, i_sym)
   !endif

   call balance_fem_node_graph(n_pts, n_gelts_triangs, &
      v_triang_nodes, v_nd_iphyscurve, vx, vy,  &
      assertions_on, errco, emsg)
   RETONERROR(errco)

   call write_mail_file(fname_mail, n_pts, n_gelts_triangs, &
      vx, vy, v_nd_iphyscurve, v_triang_nodes, &
      v_eltint_physsurf)

      if (assertions_on .ne. 0) then
         write(*,*) 'Node check after FEM rebalance'
         call check_point_separations(n_pts, vx, vy, errco, emsg)
      endif

end


!##################################################################################
!TODO: why is this not done direct from python where system calling is easier?

subroutine make_msh_file(fname_geo, fname_msh, sysret)

   use numbatmod

   character(len=*), intent(in) :: fname_geo, fname_msh

   character(len=EMSG_LENGTH) :: com_line
   character(len=256) :: gmsh_app
   integer sysret


#ifdef __APPLE__
   gmsh_app="/Applications/Gmsh.app/Contents/MacOS/gmsh"
#else
   gmsh_app = "gmsh"
#endif !!  __APPLE__

   com_line = trim(gmsh_app) // " " // "-0 -2  -order 2 -v 0 -o " // &
       trim(fname_msh) // " " // trim(fname_geo)


   sysret = nb_system(com_line)


   end

   !##################################################################################

   subroutine parse_msh_file(fname_msh, gmsh_version, n_pts, n_elts, v_ipts, &
           v_ielts, vx, vy, v_gmsh_elt_type, errco, emsg)

       use numbatmod
       implicit none
       integer gmsh_version

       integer n_pts, n_elts
       integer errco
       character(len=EMSG_LENGTH) emsg
       character fname_msh*(FNAME_LENGTH)

       integer v_ipts(MAX_N_PTS)
       integer v_ielts(MAX_N_ELTS), v_gmsh_elt_type(MAX_N_ELTS)
       double precision vx(MAX_N_PTS), vy(MAX_N_PTS)

       character str_in*(FNAME_LENGTH)
       integer ui_in, i, j, tmp1

       !  Check size of .msh file
       ui_in = 24
       open (unit=ui_in,file=fname_msh)
       if(gmsh_version .eq. 2) then  !  what is alternative to v2 ?
           read(ui_in,'(a1)') str_in       !  $MeshFormat
           read(ui_in,'(a1)') str_in       !  $Version string
           read(ui_in,'(a1)') str_in       !  $EndMeshFormat
       endif
       read(ui_in,'(a5)') str_in         !  $Nodes
       read(ui_in,*) n_pts                !  $Number of nodes
       !  close(ui_in)
       !
       if(MAX_N_PTS .lt. n_pts) then !
           write(emsg,*) 'Grid too large: The requested number of ', &
               'nodes ', n_pts, &
               'exceeds the maximum allowed of ', MAX_N_PTS, ' nodes.', &
               'Try reducing the lc_* grid resolution parameters.'
           errco = -1
           close(ui_in)
           return
       endif

       !  Maps gmsh node number to our sequence of node numbers 1..n_pts
       !  Seems like v_ipts just ends up as trivial mapping 1..n_pts,
       !  but perhaps this is not guaranteed by gmsh?
       do i=1,n_pts
       read(ui_in,*) j, vx(i), vy(i), tmp1 !  Node number, xval, yval
       v_ipts(j) = i
       if (v_ipts(j) .ne. j)  then !  REMOVE ME
           write(emsg,*) 'v_ipts misalignment in conv_gmsh'
           errco = -22
           close(ui_in)
           return
       endif
       enddo

       read(ui_in,'(a5)') str_in           !  $EndNodes


       !  Read elements
       read(ui_in,'(a5)') str_in           !  $Elements
       read(ui_in,*) n_elts                !  Number of elements

       if(MAX_N_ELTS .lt. n_elts) then
           write(emsg,*) 'Too many FEM elts: MAX_N_ELTSlts < n_elts', &
               MAX_N_ELTS, n_elts
           errco = -2
           close(ui_in)
           return
       endif

       !  Read array of elements:  Index EltType(8 or 9) Number_of_tags <tag>  node-number-list (3 or 6 nodes)
       do i=1,n_elts
       read(ui_in,*) j, v_gmsh_elt_type(i)
       v_ielts(j) = i   !  map gmsh index to our index (if they are ever different?!) TODO: this is never used, REMOVE
       enddo
       !
       close(ui_in)

       !  confirm that v_ipts and v_ielts are trivial
       do i=1,n_pts
       if (v_ipts(i) .ne. i) then
           write(emsg,*) 'v_ipts(i) mismatch at i=', i
           errco=-3
           close(ui_in)
           return
       endif
       enddo
       do i=1,n_elts
       if (v_ielts(i) .ne. i) then
           write(emsg,*) 'v_ielts(i) mismatch at i=', i
           errco=-4
           close(ui_in)
           return
       endif
       enddo

       end

       !##################################################################################

       subroutine decode_element_tags(fname_msh, gmsh_version, &
               n_pts, n_elts, v_gmsh_elt_type, v_ipts, &
               n_gelts_lines, n_gelts_triangs, v_lines_nodes, &
               v_triang_nodes, v_eltbdy_physcurve, v_eltint_physsurf,  &
               errco, emsg)

           use numbatmod

           !  Format of Gmsh Element lines
           !  GMSH_TYPE_LINE2NODE = 8
           !  These elements are edges at the outer boundary and so are associated with one material
           !  line format:  index "8" "2" Physical_Line Line Node_index x 3
           !  Here Physical_Line and Line are exactly the values in the .geo file on which this elt lies

           !  GMSH_TYPE_TRIANG6NODE = 9
           !  These are interior triangles and enclose one material
           !  line format:  index "9" "2" Physical_Surface Plane_Surface VertexNode_index_x_3 EdgeNode_index_x_3
           !  Here Physical_Surface and Plane_Surface are exactly the values in the .geo file in which this elt lies

           integer, parameter :: GMSH_TYPE_LINE2NODE=8
           integer, parameter :: GMSH_TYPE_TRIANG6NODE=9

           integer errco
           character(len=EMSG_LENGTH) emsg
           character fname_msh*(FNAME_LENGTH)
           character str_in*(FNAME_LENGTH)
           integer ui_in
           integer dummy(10), n_pretags, physmat_tag, i,j,k
           double precision tmp1, tmp2, tmp3


           integer v_ipts(MAX_N_PTS)
           integer v_gmsh_elt_type(MAX_N_ELTS)

           integer n_elts, n_pts, gmsh_version
           integer n_gelts_triangs, n_gelts_lines

           integer v_lines_nodes(3, MAX_N_ELTS)
           integer v_triang_nodes(6,MAX_N_ELTS)

           integer v_eltbdy_physcurve(MAX_N_ELTS)
           integer v_eltint_physsurf(MAX_N_ELTS)


           !  Details at: https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format

           !  Set number of expected tags in Element lines
           if(gmsh_version .eq. 2) then
               !  formerly 6 on windows gmsh 2.5.0
               n_pretags = 5        !  Number of tags before fields defining node numbers in an element
               physmat_tag = 4  !  Which field identifies physical curve or physical surface of an element (associated with a particular material)
           else   !  TODO: This should never happen. remove
               n_pretags = 5
               physmat_tag = 3
           endif
           !



           !  Skip over lines before the elements
           ui_in = 24
           open (unit=ui_in,file=fname_msh)
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

           n_gelts_lines = 0    !number of line2nodes
           n_gelts_triangs = 0    !number of triang6nodes

           do i=1,n_elts
           !select(v_gmsh_elt_type(i))
           !  case(GMSH_TYPE_LINE2NODE)
           if(v_gmsh_elt_type(i) .eq. GMSH_TYPE_LINE2NODE) then   !  2nd field is 8 = 3-node second order line
               n_gelts_lines = n_gelts_lines + 1

               read(ui_in,*) (dummy(k), k=1,n_pretags), &
                   (v_lines_nodes(k,n_gelts_lines), k=1,3) !  Get gmsh node numbers for this element

               do k=1,3
               j = v_lines_nodes(k,n_gelts_lines)
               v_lines_nodes(k,n_gelts_lines) = v_ipts(j)   !  Map to our node numbers (actually the same)
               enddo

               !  TODO: Would expect this to all be the same outer material number but seems not
               v_eltbdy_physcurve(n_gelts_lines) = dummy(physmat_tag)   !  Get phys_curve index for this element.

               !case(GMSH_TYPE_TRIANG6NODE)
           elseif(v_gmsh_elt_type(i) .eq. GMSH_TYPE_TRIANG6NODE) then  !2nd field is 9 = 6-node second order triangle

               n_gelts_triangs = n_gelts_triangs + 1
               read(ui_in,*) (dummy(k), k=1,n_pretags), &
                   (v_triang_nodes(k,n_gelts_triangs), k=1,6)        !  Get gmsh node numbers for this element

               do k=1,6  !  this loop seems to be a no-op in that v_ipts(j)=j afaict
               j = v_triang_nodes(k,n_gelts_triangs)
               v_triang_nodes(k,n_gelts_triangs) = v_ipts(j)  !  Map to our node numbers (actually the same)
               enddo

               v_eltint_physsurf(n_gelts_triangs) = dummy(physmat_tag)  !  Get phys_surface index for this element.

           else
               !case default

               write(emsg,*) 'Unknown gmsh elt type:', &
                   'v_gmsh_elt_type(i), i = ', v_gmsh_elt_type(i), i
               errco = -5
               return

           endif
           !end select
           enddo

           close(ui_in)

           end

           !##################################################################################
           subroutine write_mail_file(fname_mail, n_pts, n_gelts_triangs, &
                   vx, vy, v_nd_iphyscurve, v_triang_nodes, &
                   v_eltint_physsurf)

               !  Write the NumBAT format .mail file
               !  Format:
               !  Number_of_nodes  Number_of_6node_triangles
               !  Node_number  x_j   y_j  v_nd_iphyscurve_j :   Node positions(x,y), and physcurve number for those nodes which are on bdy elts (v_nd_iphyscurve(i) != 0)
               !  Triangle_number  6x node indices  physsurf_index   !  only for 6-node triangles, not 3-node lines

               use numbatmod
               implicit none
               integer ui_out, i, k

               character fname_mail*(FNAME_LENGTH)
               integer n_pts,  n_gelts_triangs

               integer v_triang_nodes(6,MAX_N_ELTS)


               integer v_eltint_physsurf(MAX_N_ELTS)
               double precision vx(MAX_N_PTS), vy(MAX_N_PTS)
               integer v_nd_iphyscurve(MAX_N_PTS)


               ui_out = 26
               open (unit=ui_out,file=fname_mail)

               write(ui_out,*) n_pts, n_gelts_triangs

               do i=1,n_pts
               write(ui_out, '(2x, i7, 2(g25.15), i6)') i, vx(i), vy(i), &
                   v_nd_iphyscurve(i)
               enddo

               do i=1,n_gelts_triangs
               !write(ui_out,*) i, (v_triang_nodes(k,i), k=1,6), &
               !  v_eltint_physsurf(i)
               write(ui_out,'(2x, i9, i9, i9, i9, i9, i9, i9, i9)') &
                   i, (v_triang_nodes(k,i), k=1,6), v_eltint_physsurf(i)
               enddo
               close(ui_out)

               end
