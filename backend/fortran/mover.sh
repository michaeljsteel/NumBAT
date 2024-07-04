for fb in z_indexx_AC.f z_indexx.f
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
#get_H_field_p3.f
#gmsh_post_process_ac.f gmsh_post_process.f
#h_mode_field_ez.f interp_nod_2d.f
#jacobian_p1_2d.f jacobian_p2_2d.f
#lattice_vec.f moving_boundary.f periodic_cond.f periodic_N_E_F.f periodic_node.f
#photoelastic_int.f photoelastic_int_v2.f quad_triangle.f 
#sort_csr.f sort_int.f z_indexx_AC.f z_indexx.f
