
c     Construct the FEM mesh
c
c   type_nod = 0  => interiour point
c   type_nod != 0 => boundary point
c
c
      subroutine geometry (n_msh_el, n_msh_pts, nnodes, n_typ_el,
     *    lx, ly, type_nod, type_el, table_nod, x, 
     *    mesh_file, errno_arpack, emsg_arpack)
c
      implicit none
      integer*8 n_msh_el, n_msh_pts, nnodes, n_typ_el, max_typ_el
      integer*8 type_nod(n_msh_pts), type_el(n_msh_el)
      integer*8 table_nod(nnodes,n_msh_el), n_typ_el2
      double precision lx, ly
      parameter (max_typ_el=10)
      double precision x(2,n_msh_pts), xx(2)

      character mesh_file*1000
      integer*8 errno_arpack
      character*2048 emsg_arpack


      integer*8 n_msh_pts2, n_msh_el2, ui
      integer*8 i, j, k
c
      ui = 6
c

      ! is the mesh file consistent with what we expect
      open (unit=24,file=mesh_file, status='old')
      read(24,*) n_msh_pts2, n_msh_el2
c
      if(n_msh_pts .ne. n_msh_pts2) then
        write(ui,*) "geometry: n_msh_pts != n_msh_pts2 : ", 
     *    n_msh_pts, n_msh_pts2
      endif
      if(n_msh_el .ne. n_msh_el2) then
        write(ui,*) "geometry: n_msh_el != n_msh_el2 : ", 
     *    n_msh_el, n_msh_el2
      endif

c    Coordinate of the FEM points
      do i=1,n_msh_pts
        read(24,*) k, (xx(j),j=1,2), type_nod(i)
        x(1,i) = xx(1)*lx
        x(2,i) = xx(2)*ly
      enddo

c     Connectivity table
      n_typ_el2 = 1   ! largest index of materials in the file
      do i=1,n_msh_el
        read(24,*) k, (table_nod(j,i),j=1,nnodes), type_el(i)
        j = type_el(i)
        if(n_typ_el2 .lt. j) n_typ_el2 = j
        if(j .lt. 0) then
          write(emsg_arpack,*) "geometry: type_el(i) < 0 : ", 
     *      i, type_el(i)
         errno_arpack=-9
         return 
        endif
      enddo
      close(24)

      if(n_typ_el2 .gt. n_typ_el) then
         write(emsg_arpack,*) "geometry: n_typ_el2 > n_typ_el : ", 
     *    n_typ_el2, n_typ_el
         errno_arpack=-10
         return 
      endif

      return
      end

