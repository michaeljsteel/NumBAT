
subroutine label_faces (nel, m_elnd_to_mshpte)


   implicit none
   integer(8) nel
   integer(8) m_elnd_to_mshpte(14,nel)
   integer(8) i

   !  Table of connectivity for the face (for 2D FEM, face = triangle element)

   do i=1,nel !  each element is a face
      m_elnd_to_mshpte(1,i) = i
   enddo

   return
end
