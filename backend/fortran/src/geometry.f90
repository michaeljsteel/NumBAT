#include "numbat_decl.h"

!  Construct the FEM mesh

!  type_nod = 0  => interior point
!  type_nod != 0 => boundary point

!  Reads .mail file to find
!  - x,y coords of mesh points  (xy_nodes)
!  - mesh points associated with each element (table_nod)
!  - whether number of material types read matches expected value n_elt_mats

!  -  Fills:  xy_nodes, type_nod, type_el, table_nod

subroutine construct_fem_node_tables_em(mesh_file, dim_x, dim_y, &
   n_msh_el, n_msh_pts, d_nodes_per_el, n_elt_mats, &
   mesh_props, errco, emsg)

   use numbatmod
   use class_MeshProps

   ! ins
   character(len=*) mesh_file
   double precision dim_x, dim_y

   ! outs

   type(MeshProps) :: mesh_props

   integer(8) n_msh_el, n_msh_pts, d_nodes_per_el, n_elt_mats

   integer(8) type_nod(n_msh_pts), type_el(n_msh_el)
   integer(8) table_nod(d_nodes_per_el, n_msh_el)
   double precision xy_nodes(2,n_msh_pts)

   integer errco
   character(len=EMSG_LENGTH) :: emsg


   ! locals
   double precision xx(2)

   integer(8) n_elt_mats2
   !integer, parameter :: max_typ_el=10
   integer(8) n_msh_pts2, n_msh_el2
   integer(8) i, j, k


   !  check the mesh file is consistent with what we expect
   open (unit=24,file=mesh_file, status='old')
   read(24,*) n_msh_pts2, n_msh_el2

   if(n_msh_pts .ne. n_msh_pts2) then
      write(emsg,*) "construct_fem_nodal_tables: n_msh_pts != n_msh_pts2 : ", n_msh_pts, n_msh_pts2
      errco=-101
      return
   endif

   if(n_msh_el .ne. n_msh_el2) then
      write(emsg,*) "construct_fem_nodal_tables: n_msh_el != n_msh_el2 : ", n_msh_el, n_msh_el2
      errco=-102
      return
   endif


!  Read coordinates of the FEM mesh points
   do i=1,n_msh_pts
      read(24,*) k, (xx(j),j=1,2), type_nod(i)
      xy_nodes(1,i) = xx(1)*dim_x
      xy_nodes(2,i) = xx(2)*dim_y
   enddo

!  Connectivity table
   n_elt_mats2 = 1   !  largest index of materials in the file
   do i=1,n_msh_el
      read(24,*) k, (table_nod(j,i),j=1,d_nodes_per_el), type_el(i)
      j = type_el(i)
      if(n_elt_mats2 .lt. j) n_elt_mats2 = j
      if(j .lt. 0) then
         write(emsg,*) "geometry: type_el(i) < 0 : ", i, type_el(i)
         errco=-9
         return
      endif
   enddo
   close(24)

   if(n_elt_mats2 .gt. n_elt_mats) then
      write(emsg,*) "geometry: n_elt_mats2 > n_elt_mats : ", n_elt_mats2, n_elt_mats
      errco=-10
      return
   endif

   mesh_props%type_el = type_el
   mesh_props%type_nod = type_nod
   mesh_props%table_nod = table_nod
   mesh_props%xy_nodes = xy_nodes

   return
end


subroutine construct_fem_node_tables(mesh_file, dim_x, dim_y, &
   n_msh_el, n_msh_pts, d_nodes_per_el, n_elt_mats, &
    xy_nodes, type_nod, type_el, table_nod, &
   errco, emsg)

   use numbatmod

   ! ins
   character(len=*) mesh_file
   double precision dim_x, dim_y

   ! outs
   integer(8) n_msh_el, n_msh_pts, d_nodes_per_el, n_elt_mats
   integer(8) type_nod(n_msh_pts), type_el(n_msh_el)
   integer(8) table_nod(d_nodes_per_el, n_msh_el)
   double precision xy_nodes(2,n_msh_pts)

   integer errco
   character(len=EMSG_LENGTH) :: emsg


   ! locals
   double precision xx(2)

   integer(8) n_elt_mats2
   !integer, parameter :: max_typ_el=10
   integer(8) n_msh_pts2, n_msh_el2
   integer(8) i, j, k


   !  check the mesh file is consistent with what we expect
   open (unit=24,file=mesh_file, status='old')
   read(24,*) n_msh_pts2, n_msh_el2

   if(n_msh_pts .ne. n_msh_pts2) then
      write(emsg,*) "construct_fem_nodal_tables: n_msh_pts != n_msh_pts2 : ", n_msh_pts, n_msh_pts2
      errco=-101
      return
   endif

   if(n_msh_el .ne. n_msh_el2) then
      write(emsg,*) "construct_fem_nodal_tables: n_msh_el != n_msh_el2 : ", n_msh_el, n_msh_el2
      errco=-102
      return
   endif


!  Read coordinates of the FEM mesh points
   do i=1,n_msh_pts
      read(24,*) k, (xx(j),j=1,2), type_nod(i)
      xy_nodes(1,i) = xx(1)*dim_x
      xy_nodes(2,i) = xx(2)*dim_y
   enddo

!  Connectivity table
   n_elt_mats2 = 1   !  largest index of materials in the file
   do i=1,n_msh_el
      read(24,*) k, (table_nod(j,i),j=1,d_nodes_per_el), type_el(i)
      j = type_el(i)
      if(n_elt_mats2 .lt. j) n_elt_mats2 = j
      if(j .lt. 0) then
         write(emsg,*) "geometry: type_el(i) < 0 : ", i, type_el(i)
         errco=-9
         return
      endif
   enddo
   close(24)

   if(n_elt_mats2 .gt. n_elt_mats) then
      write(emsg,*) "geometry: n_elt_mats2 > n_elt_mats : ", n_elt_mats2, n_elt_mats
      errco=-10
      return
   endif

   return
end

