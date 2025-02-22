

subroutine MeshRawAC_allocate(this, n_msh_pts, n_msh_el, n_elt_mats, nberr)

   class(MeshRawAC) :: this
   integer(8) :: n_msh_el, n_msh_pts, n_elt_mats
   type(NBError) nberr

   this%n_msh_pts = n_msh_pts
   this%n_msh_el = n_msh_el
   this%n_elt_mats = n_elt_mats

   call double_nalloc_2d(this%v_nd_xy, 2_8, n_msh_pts, 'v_nd_xy', nberr); RET_ON_NBERR(nberr)

   call integer_nalloc_1d(this%v_nd_physindex, n_msh_pts, 'v_nd_physindex', nberr); RET_ON_NBERR(nberr)

   call integer_nalloc_2d(this%elnd_to_mshpt, P2_NODES_PER_EL, n_msh_el, 'elnd_to_mshpt', nberr);
   RET_ON_NBERR(nberr)

end subroutine


 ! Seems identical to the EM version.
subroutine MeshRawAC_load_node_tables_from_py(this, &
   v_nd_xy, v_nd_physindex, el_material, elnd_to_mshpt, nberr)

   class(MeshRawAC) :: this

   type(NBError) nberr


   integer(8) :: v_nd_physindex(this%n_msh_pts)
   double precision :: v_nd_xy(2, this%n_msh_pts)
   integer(8) :: el_material(this%n_msh_el)
   integer(8) :: elnd_to_mshpt(P2_NODES_PER_EL, this%n_msh_el)

   this%v_nd_physindex = v_nd_physindex
   this%v_nd_xy = v_nd_xy
   this%el_material = el_material
   this%elnd_to_mshpt = elnd_to_mshpt

end subroutine

 ! Seems identical to the EM version.
subroutine MeshRawAC_construct_node_tables(this, mesh_file, dimscale_in_m,nberr)

   class(MeshRawAC) :: this

   ! ins
   character(len=*) mesh_file
   double precision dimscale_in_m

   ! outs
   !integer(8) elnd_to_mshpt(P2_NODES_PER_EL, this%n_msh_el)

   type(NBError) nberr

   ! locals
   character(len=EMSG_LENGTH) :: emsg


   double precision xx(2)

   integer(8) n_elt_mats2, n_msh_pts2, n_msh_el2
   integer(8) i, j, k
   integer(8) ui


   ui = 24

   write(*,*) 'opening meshfile', mesh_file, '<<'

   !  check the mesh file is consistent with what we expect
   open (unit=ui,file=mesh_file, status='old')
   read(ui,*) n_msh_pts2, n_msh_el2

   if (this%n_msh_pts .ne. n_msh_pts2) then
      write(emsg,*) "construct_fem_nodal_tables: n_msh_pts != n_msh_pts2 : ", this%n_msh_pts, n_msh_pts2
      call nberr%set(-101_8, emsg)
      return
   endif

   if (this%n_msh_el .ne. n_msh_el2) then
      write(emsg,*) "construct_fem_nodal_tables: n_msh_el != n_msh_el2 : ", this%n_msh_el, n_msh_el2
      call nberr%set(-102_8, emsg)
      return
   endif


   !  Read coordinates of the FEM mesh points
   do i=1,this%n_msh_pts
      read(ui,*) k, (xx(j),j=1,2), this%v_nd_physindex(i)
      this%v_nd_xy(:,i) = xx*dimscale_in_m
   enddo

   !  Connectivity table
   n_elt_mats2 = 1   !  largest index of materials in the file
   do i=1,this%n_msh_el

      read(ui,*) k, (this%elnd_to_mshpt(j,i),j=1,P2_NODES_PER_EL), this%el_material(i)

      j = this%el_material(i)
      if(n_elt_mats2 .lt. j) n_elt_mats2 = j

      ! shouldn't happen, delete
      if(j .lt. 0) then
         write(emsg,*) "geometry: type_el(i) < 0 : ", i, this%el_material(i)
         call nberr%set(-9_8, emsg)
         return
      endif
   enddo

   close(ui)

   if(n_elt_mats2 .gt. this%n_elt_mats) then
      write(emsg,*) "geometry: n_elt_mats2 > n_elt_mats : ", n_elt_mats2, this%n_elt_mats
      call nberr%set(-10_8, emsg)
      return
   endif

end

