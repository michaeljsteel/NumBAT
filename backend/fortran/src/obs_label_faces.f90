
subroutine label_faces (nel, elnd_to_meshe)


   implicit none
   integer(8) nel
   integer(8) elnd_to_meshe(14,nel)
   integer(8) i

   !  Table of connectivity for the face (for 2D FEM, face = triangle element)

   do i=1,nel !  each element is a face
      elnd_to_meshe(1,i) = i
   enddo

   return
end
