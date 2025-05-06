#include "numbat_decl.h"




!*******************************************************

!  conv_gmsh: convert the Gmsh .geo to mesh format and
!  then to the NumBAT .mail mesh format

!*******************************************************

subroutine conv_gmsh_impl(geo_fname, assertions_on, nberr)

   use numbatmod
   use alloc
   use conv_gmsh_interfaces


   character(len=*), intent(in) :: geo_fname
   integer(8), intent(in) :: assertions_on
   type(NBError), intent(inout) :: nberr

   integer(8) errco
   character(len=EMSG_LENGTH) emsg

   integer(8) i_sym, gmsh_version

   character(len=FNAME_LENGTH) :: fname_geo, fname_msh, fname_mail
   integer(8) n_msh_elts, n_msh_pts

   !  d1 vars contain 3 node lines (node, node, edge)
   !  d2 vars contain 6 node triangs (node, node, node, edge, edge, edge)
   !  number of elts found
   integer(8) n_gelt_lines, n_gelt_triangs


   !  node/edge codes for each elt
   integer(8),  dimension(:,:), allocatable :: v_lineelt_nodes
   integer(8),  dimension(:,:), allocatable :: v_trielt_nodes

   !  material index for each elt
   integer(8),  dimension(:), allocatable :: v_lineelt_physcurve
   integer(8),  dimension(:), allocatable :: v_trielt_physsurf
   integer(8),  dimension(:), allocatable :: v_mshpt_iphyscurve

   !  Individual nodes, number and position
   integer(8),  dimension(:), allocatable :: v_i_mshpts
   double precision, dimension(:), allocatable :: vx,  vy

   !Elements, number and material index
   integer(8),  dimension(:), allocatable :: v_gmsh_elt_type

   integer(8) i, j
   integer(8) fnamelen
   integer(8) iphyscurve, nd

   gmsh_version = 2
   errco = 0
   i_sym = 0   !  symmetry is off

   fnamelen = len_trim(geo_fname)
   if (fnamelen .ge. FNAME_LENGTH) then
      write(emsg, *) "Name of .geo file is too long extend in ", &
         "conv_gmsh_py.f"
      call nberr%set(NBERROR_110, emsg)
      return
   endif

   fname_geo = geo_fname(1:fnamelen)//".geo"
   fname_msh = geo_fname(1:fnamelen)//".msh"
   fname_mail = geo_fname(1:fnamelen)//".mail"

   ! First stage conversion
   ! Find number of mesh points and their locations
   !      number of triangle elements or boundary lines and their Gmsh physical indices
   call parse_msh_file(fname_msh, gmsh_version, n_msh_pts, n_msh_elts, &
      vx, vy, v_i_mshpts, v_gmsh_elt_type, nberr)
   RET_ON_NBERR(nberr)

   if (assertions_on .ne. 0) then
      write(*,*) 'Node check fresh from gmsh:'
      call check_point_separations(n_msh_pts, vx, vy, errco, emsg)
   endif

   if (n_msh_pts .gt. MAX_N_ELTS) then
      write(emsg, '(A,I7,A,I7,A)') 'The generated mesh has ', n_msh_pts, &
         ' nodes, which exceeds the maximum of ', MAX_N_ELTS, '.'
      call nberr%set(NBERR_MESH_TOO_LARGE, emsg)
      return
   endif

   ! We could work out the exact numbers from parse_msh_file
   ! But these sizes are sufficient

   call integer_alloc_2d(v_lineelt_nodes, 3_8, n_msh_elts, 'v_lineelt_nodes', nberr); RET_ON_NBERR(nberr)
   call integer_alloc_2d(v_trielt_nodes, 6_8, n_msh_elts, 'v_trielt_nodes', nberr); RET_ON_NBERR(nberr)
   call integer_alloc_1d(v_lineelt_physcurve, n_msh_elts, 'v_lineelt_nodes', nberr); RET_ON_NBERR(nberr)
   call integer_alloc_1d(v_trielt_physsurf, n_msh_elts, 'v_lineelt_nodes', nberr); RET_ON_NBERR(nberr)
   call integer_alloc_1d(v_mshpt_iphyscurve, n_msh_pts, 'v_lineelt_nodes', nberr); RET_ON_NBERR(nberr)

   !  Now we know the number of points and mappings (even if trivial)
   !  Next we load elt data according to the gmsh element types

   call decode_element_tags(fname_msh, gmsh_version, &
      n_msh_pts, n_msh_elts, v_i_mshpts, &
      v_gmsh_elt_type, n_gelt_lines, n_gelt_triangs, &
      v_lineelt_nodes, v_trielt_nodes, v_lineelt_physcurve, v_trielt_physsurf, nberr)
   RET_ON_NBERR(nberr)

   !  Now we have:
   !  v_lineelt_nodes[:, j] = (node#, node#, node# for 3-node line number j)
   !  v_trielt_nodes[:, j] = (cornernode# x3,  edgenode# x 3 for 6-node triangle number j)
   !  v_lineelt_physcurve[] = label of Physical_Curve on which meshpt lies (may or may not be outer boundary ?)
   !  v_trielt_physsurf[] = label of Physical_Surface inside the  domain

   !  Next, associate the three nodes on each boundary elt with their physical_curve stored in
   !  v_lineelt_nodes gives the set of node triples on a boundary
   !  translate this to v_mshpt_iphyscurve, which specifies a physical boundary or zero for _every_ mesh elt.
   !  If v_mshpt_iphyscurve(j) = pc_k != 0,  node j lies on PhysicalCurve pc_k

   v_mshpt_iphyscurve = 0

   do i=1,n_gelt_lines           !  for each boundary line element
      iphyscurve = v_lineelt_physcurve(i)   !  associate its 3 nodes with the PhysicalCurve it is on

      do j=1,3
         nd = v_lineelt_nodes(j,i)
         v_mshpt_iphyscurve(nd) = iphyscurve                !  assign this material to a node on the line
      enddo
   enddo

   ! v_lineelt_nodes is now finished with

   !  i_sym = 0 always, so this is switched off  TODO: the call to symmetry is very broken with parameteres
   !if(i_sym .ne. 0) then

   !call symmetry(n_msh_pts, n_gelt_triangs, &
   !  MAX_N_ELTS, MAX_N_PTS, v_mshpt_iphyscurve, v_trielt_nodes, &
   !  v_trielt_physsurf, vx, vy, i_sym)
   !endif

   ! reorder the mesh points to make the mesh better somehow
   call balance_fem_node_graph(n_msh_pts, n_gelt_triangs, &
      v_trielt_nodes, v_mshpt_iphyscurve, vx, vy,  &
      assertions_on, nberr)
   RET_ON_NBERR(nberr)

   ! generate final numbat mesh file format
   call write_mail_file(fname_mail, n_msh_pts, n_msh_elts, n_gelt_triangs, &
      vx, vy, v_mshpt_iphyscurve, v_trielt_nodes, v_trielt_physsurf)

   if (assertions_on .ne. 0) then
      write(*,*) 'Node check after FEM rebalance'
      call check_point_separations(n_msh_pts, vx, vy, errco, emsg)
   endif

end




 !  Format of Gmsh Element lines
!  GMSH_TYPE_LINE2NODE = 8
!  These elements are edges at the outer boundary and so are associated with one material
!  line format:  index "8" "2" Physical_Line Line Node_index x 3
!  Here Physical_Line and Line are exactly the values in the .geo file on which this elt lies

!  GMSH_TYPE_TRIANG6NODE = 9
!  These are interior triangles and enclose one material
!  line format:  index "9" "2" Physical_Surface Plane_Surface VertexNode_index_x_3 EdgeNode_index_x_3
!  Here Physical_Surface and Plane_Surface are exactly the values in the .geo file in which this elt lies

subroutine decode_element_tags(fname_msh, gmsh_version, &
   n_msh_pts, n_msh_elts, v_i_mshpts, v_gmsh_elt_type,  &
   n_gelt_lines, n_gelt_triangs, v_lineelt_nodes, &
   v_trielt_nodes, v_lineelt_physcurve, v_trielt_physsurf, nberr)

   use numbatmod
   type(NBError) nberr

   integer(8),  parameter :: GMSH_TYPE_LINE2NODE=8
   integer(8),  parameter :: GMSH_TYPE_TRIANG6NODE=9

   character(len=*), intent(in) :: fname_msh

   integer(8) n_msh_elts, n_msh_pts, gmsh_version
   integer(8) n_gelt_triangs, n_gelt_lines



   integer(8) v_i_mshpts(n_msh_pts)
   integer(8) v_gmsh_elt_type(n_msh_elts)

   integer(8) v_lineelt_nodes(3, n_msh_elts)
   integer(8) v_trielt_nodes(6, n_msh_elts)

   integer(8) v_lineelt_physcurve(n_msh_elts)
   integer(8) v_trielt_physsurf(n_msh_elts)

   character(len=EMSG_LENGTH) emsg

   !-----------------------------------

   character str_in*(FNAME_LENGTH)
   integer(8) dummy(10), n_pretags, physmat_tag, i,j,k

   integer(8) ui_in
   double precision tmp1, tmp2, tmp3


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
   do i=1,n_msh_pts                       !Nodes and positions
      read(ui_in,*) j, tmp1, tmp2, tmp3
   enddo
   read(ui_in,'(a5)') str_in             !3 line Elements header
   read(ui_in,'(a5)') str_in
   read(ui_in,*) tmp1

   n_gelt_lines = 0    !number of line2nodes
   n_gelt_triangs = 0    !number of triang6nodes

   do i=1,n_msh_elts
      if(v_gmsh_elt_type(i) .eq. GMSH_TYPE_LINE2NODE) then   !  2nd field is 8 = 3-node second order line
         n_gelt_lines = n_gelt_lines + 1

         read(ui_in,*) (dummy(k), k=1,n_pretags), &
            (v_lineelt_nodes(k,n_gelt_lines), k=1,3) !  Get gmsh node numbers for this element

         do k=1,3
            j = v_lineelt_nodes(k,n_gelt_lines)
            v_lineelt_nodes(k,n_gelt_lines) = v_i_mshpts(j)   !  Map to our node numbers (actually the same)
         enddo

         !  TODO: Would expect this to all be the same outer material number but seems not
         v_lineelt_physcurve(n_gelt_lines) = dummy(physmat_tag)   !  Get phys_curve index for this element.

      elseif(v_gmsh_elt_type(i) .eq. GMSH_TYPE_TRIANG6NODE) then  !2nd field is 9 = 6-node second order triangle

         n_gelt_triangs = n_gelt_triangs + 1
         read(ui_in,*) (dummy(k), k=1,n_pretags), &
            (v_trielt_nodes(k,n_gelt_triangs), k=1,6)        !  Get gmsh node numbers for this element

         do k=1,6  !  this loop seems to be a no-op in that v_i_mshpts(j)=j afaict
            j = v_trielt_nodes(k,n_gelt_triangs)
            v_trielt_nodes(k,n_gelt_triangs) = v_i_mshpts(j)  !  Map to our node numbers (actually the same)
         enddo

         v_trielt_physsurf(n_gelt_triangs) = dummy(physmat_tag)  !  Get phys_surface index for this element.

      else

         write(emsg,*) 'Unknown gmsh elt type:', &
            'v_gmsh_elt_type(i), i = ', v_gmsh_elt_type(i), i
         call nberr%set(NBERROR_117, emsg)
         return

      endif
      !end select
   enddo

   close(ui_in)

end


subroutine write_mail_file(fname_mail, n_msh_pts, n_msh_elts, n_gelt_triangs, &
   vx, vy, v_mshpt_iphyscurve, v_trielt_nodes, &
   v_trielt_physsurf)

   !  Write the NumBAT format .mail file
   !  Format:
   !  Number_of_nodes  Number_of_6node_triangles
   !  Node_number  x_j   y_j  v_mshpt_iphyscurve_j :   Node positions(x,y), and physcurve number for those nodes which are on bdy elts (v_mshpt_iphyscurve(i) != 0)
   !  Triangle_number  6x node indices  physsurf_index   !  only for 6-node triangles, not 3-node lines

   use numbatmod
   implicit none
   integer(8) ui_out, i, k

   character(len=*) fname_mail
   integer(8) n_msh_pts,  n_msh_elts, n_gelt_triangs

   double precision vx(n_msh_pts), vy(n_msh_pts)
   integer(8) v_mshpt_iphyscurve(n_msh_pts)
   integer(8) v_trielt_nodes(6,n_msh_elts)
   integer(8) v_trielt_physsurf(n_msh_elts)

   ui_out = 26
   open (unit=ui_out,file=fname_mail)

   write(ui_out,*) n_msh_pts, n_gelt_triangs

   do i=1,n_msh_pts
      write(ui_out, '(2x, i7, 2(g25.15), i6)') i, vx(i), vy(i), v_mshpt_iphyscurve(i)
   enddo

   do i=1,n_gelt_triangs
      write(ui_out,'(2x, i9, i9, i9, i9, i9, i9, i9, i9)') &
         i, (v_trielt_nodes(k,i), k=1,6), v_trielt_physsurf(i)
   enddo

   close(ui_out)

end
