#include "numbat_decl.h"



subroutine construct_fem_node_tables_ac(mesh_file, dimscale_in_m, dim_y, &
   n_msh_el, n_msh_pts, nodes_per_el, n_elt_mats, &
    v_nd_xy, type_nod, type_el, elnd_to_mesh, &
   errco, emsg)

   use numbatmod

   ! ins
   character(len=*) mesh_file
   double precision dimscale_in_m, dim_y

   ! outs
   integer(8) n_msh_el, n_msh_pts, nodes_per_el, n_elt_mats
   integer(8) type_nod(n_msh_pts), type_el(n_msh_el)
   integer(8) elnd_to_mesh(nodes_per_el, n_msh_el)
   double precision v_nd_xy(2,n_msh_pts)

   integer(8) errco
   character(len=EMSG_LENGTH) :: emsg


   ! locals
   double precision xx(2)

   integer(8) n_elt_mats2
   !integer(8),  parameter :: max_typ_el=10
   integer(8) n_msh_pts2, n_msh_el2
   integer(8) i, j, k
   integer(8) ui


   ui = 24

   !  check the mesh file is consistent with what we expect
   open (unit=ui,file=mesh_file, status='old')
   read(ui,*) n_msh_pts2, n_msh_el2

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
      read(ui,*) k, (xx(j),j=1,2), type_nod(i)
      v_nd_xy(1,i) = xx(1)*dimscale_in_m
      v_nd_xy(2,i) = xx(2)*dim_y
   enddo

!  Connectivity table
   n_elt_mats2 = 1   !  largest index of materials in the file
   do i=1,n_msh_el
      read(ui,*) k, (elnd_to_mesh(j,i),j=1,nodes_per_el), type_el(i)
      j = type_el(i)
      if(n_elt_mats2 .lt. j) n_elt_mats2 = j
      if(j .lt. 0) then
         write(emsg,*) "geometry: type_el(i) < 0 : ", i, type_el(i)
         errco=-9
         return
      endif
   enddo
   close(ui)

   if(n_elt_mats2 .gt. n_elt_mats) then
      write(emsg,*) "geometry: n_elt_mats2 > n_elt_mats : ", n_elt_mats2, n_elt_mats
      errco=-10
      return
   endif

   return
end

