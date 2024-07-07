
subroutine label_faces (nel, table_node)


   implicit none
   integer(8) nel
   integer(8) table_node(14,nel)
   integer(8) i

   !  Table of connectivity for the face (for 2D FEM, face = triangle element)

   do i=1,nel !  each element is a face
      table_node(1,i) = i
   enddo

   return
end
