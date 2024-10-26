

subroutine array_material_EM (nel, nb_typ_el, n_index, type_el, ls_material)

   implicit none
   integer(8) P2_NODES_PER_EL
   parameter (P2_NODES_PER_EL = 6)
   integer(8) nel, nb_typ_el
   integer(8) type_el(nel)
   complex(8) n_index(nb_typ_el)
   complex(8) ls_material(1,P2_NODES_PER_EL+7,nel)

   !     Local variables
   integer(8) k_typ
   integer(8) iel, inod
   complex(8) n_index_el

   do iel=1,nel
      k_typ = type_el(iel)
      n_index_el = n_index(k_typ)
      do inod=1,P2_NODES_PER_EL+7
         ls_material(1,inod,iel) = n_index_el
      enddo
   enddo


end
