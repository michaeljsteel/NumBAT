
      subroutine list_face (nel, table_edge_face)


      implicit none
      integer(8) nel
      integer(8) table_edge_face(14,nel)
      integer(8) i

!     Table of connectivity for the face (for 2D FEM, face = triangle element)

      do i=1,nel
         ! each element is a face
        table_edge_face(1,i) = i
      enddo

      return
      end
