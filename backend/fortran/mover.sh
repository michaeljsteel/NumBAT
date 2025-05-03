for fb in csr_length_AC.f
do
    fn=`basename $fb .f`

    echo Handling:: $fb
    cp $fb spare
    echo "findent -ofree < ${fn}.f > ${fn}.f90"
    findent -ofree < ${fn}.f > ${fn}.f90
    mv ${fn}.f90 $fb
    git mv ${fn}.f ${fn}.f90

    echo ""
done

#ac_alpha_int.f ac_alpha_int_v2.f ac_mode_elastic_energy_int.f ac_mode_elastic_energy_int_v4.f
#ac_mode_power_int.f ac_mode_power_int_v2.f ac_mode_power_int_v4.f
#array_material_ac.f array_material_em.f
#csr_length_AC.f
#em_mode_e_energy_int.f em_mode_energy_int_ez.f em_mode_energy_int.f
#em_mode_energy_int_v2_ez.f em_mode_energy_int_v2.f em_mode_energy_int_v2_wg.f
#jacobian_p1_2d.f jacobian_p2_2d.f
#lattice_vec.f moving_boundary.f
#photoelastic_int.f photoelastic_int_v2.f quad_triangle.f
#sort_csr.f sort_int.f find_eigvals_order_AC.f find_eigvals_order.f
