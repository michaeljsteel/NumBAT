#include "numbat_decl.h"



subroutine construct_fem_node_tables_ac(mesh_file, dimscale_in_m,  &
   n_msh_elts, n_msh_pts, n_elt_mats, &
    v_mshpt_xy, type_nod, type_el, elnd_to_mshpt, nberr)

   use numbatmod

   ! ins
   character(len=*) mesh_file
   double precision dimscale_in_m
   ! outs
   integer(8) n_msh_elts, n_msh_pts, n_elt_mats
   integer(8) type_nod(n_msh_pts), type_el(n_msh_elts)
   integer(8) elnd_to_mshpt(P2_NODES_PER_EL, n_msh_elts)
   double precision v_mshpt_xy(2,n_msh_pts)

   type(NBError) nberr

   ! locals
   integer(8) errco
   character(len=EMSG_LENGTH) :: emsg
   double precision xx(2)

   integer(8) n_elt_mats2
   !integer(8),  parameter :: max_typ_el=10
   integer(8) n_msh_pts2, n_msh_elts2
   integer(8) i, j, k
   integer(8) ui


   ui = 24

   !  check the mesh file is consistent with what we expect
   open (unit=ui,file=mesh_file, status='old')
   read(ui,*) n_msh_pts2, n_msh_elts2

   if(n_msh_pts .ne. n_msh_pts2) then
      write(emsg,*) "construct_fem_nodal_tables: n_msh_pts != n_msh_pts2 : ", n_msh_pts, n_msh_pts2
      errco=-101
      call nberr%set(errco, emsg)
      return
   endif

   if(n_msh_elts .ne. n_msh_elts2) then
      write(emsg,*) "construct_fem_nodal_tables: n_msh_elts != n_msh_elts2 : ", n_msh_elts, n_msh_elts2
      errco=-102
      call nberr%set(errco, emsg)
      return
   endif


!  Read coordinates of the FEM mesh points
   do i=1,n_msh_pts
      read(ui,*) k, (xx(j),j=1,2), type_nod(i)
      v_mshpt_xy(1,i) = xx(1)*dimscale_in_m
      v_mshpt_xy(2,i) = xx(2)*dimscale_in_m
   enddo

!  Connectivity table
   n_elt_mats2 = 1   !  largest index of materials in the file
   do i=1,n_msh_elts
      read(ui,*) k, (elnd_to_mshpt(j,i),j=1,P2_NODES_PER_EL), type_el(i)
      j = type_el(i)
      if(n_elt_mats2 .lt. j) n_elt_mats2 = j
      if(j .lt. 0) then
         write(emsg,*) "geometry: type_el(i) < 0 : ", i, type_el(i)
         errco=-9
         call nberr%set(errco, emsg)
         return
      endif
   enddo
   close(ui)

   if(n_elt_mats2 .gt. n_elt_mats) then
      write(emsg,*) "geometry: n_elt_mats2 > n_elt_mats : ", n_elt_mats2, n_elt_mats
      errco=-10
      call nberr%set(errco, emsg)
      return
   endif

   return
end

