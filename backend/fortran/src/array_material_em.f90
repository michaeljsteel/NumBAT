

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


subroutine array_material_AC (nel, nb_typ_el, type_el, rho, c_tensor, p_tensor, eta_tensor, ls_material)

      use numbatmod
   implicit none
   integer(8) nel, nb_typ_el
   integer(8) type_el(nel)
   complex(8) rho(nb_typ_el), c_tensor(6,6,nb_typ_el)
   complex(8) p_tensor(3,3,3,3,nb_typ_el)
   complex(8) eta_tensor(3,3,3,3,nb_typ_el)
   complex(8) ls_material(10,P2_NODES_PER_EL,nel)

   !     Locals
   integer(8) k_typ
   integer(8) iel, inod
   complex(8) rho_el, c_11, c_12, c_44
   complex(8) p_11, p_12, p_44, eta_11, eta_12, eta_44

   !f2py intent(in) nel, nb_typ_el, type_el
   !f2py intent(in) rho, c_tensor, p_tensor
   !f2py intent(in) eta_tensor

   !f2py depend(type_el) nel
   !f2py depend(rho) nb_typ_el
   !f2py depend(c_tensor) nb_typ_el
   !f2py depend(p_tensor) nb_typ_el
   !f2py depend(eta_tensor) nb_typ_el

   !f2py intent(out) ls_material



   do iel=1,nel
      k_typ = type_el(iel)
      rho_el = rho(k_typ)
      c_11 = c_tensor(1,1,k_typ)
      c_12 = c_tensor(1,2,k_typ)
      c_44 = c_tensor(4,4,k_typ)
      p_11 = p_tensor(1,1,1,1,k_typ)
      p_12 = p_tensor(1,1,2,2,k_typ)
      p_44 = p_tensor(2,3,2,3,k_typ)
      eta_11 = eta_tensor(1,1,1,1,k_typ)
      eta_12 = eta_tensor(1,1,2,2,k_typ)
      eta_44 = eta_tensor(2,3,2,3,k_typ)
      do inod=1,P2_NODES_PER_EL
         ls_material(1,inod,iel) = rho_el
         ls_material(2,inod,iel) = c_11
         ls_material(3,inod,iel) = c_12
         ls_material(4,inod,iel) = c_44
         ls_material(5,inod,iel) = p_11
         ls_material(6,inod,iel) = p_12
         ls_material(7,inod,iel) = p_44
         ls_material(8,inod,iel) = eta_11
         ls_material(9,inod,iel) = eta_12
         ls_material(10,inod,iel) = eta_44
      enddo
   enddo



end subroutine array_material_AC
