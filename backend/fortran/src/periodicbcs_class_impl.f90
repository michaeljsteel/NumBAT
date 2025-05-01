subroutine PeriodicBCs_allocate(this, mesh_raw, entities, nberr)


   class(PeriodicBCs) :: this

   type(MeshRawEM) :: mesh_raw
   type(MeshEntities) :: entities

   type(NBError) nberr


   call integer_nalloc_1d(this%iperiod_N, mesh_raw%n_msh_pts, 'iperiod_N', nberr);
   RET_ON_NBERR(nberr)

   call integer_nalloc_1d(this%iperiod_N_E_F, entities%n_entities, 'iperiod_N_E_F', nberr);
   RET_ON_NBERR(nberr)

   call integer_nalloc_1d(this%inperiod_N, mesh_raw%n_msh_pts, 'inperiod_N', nberr);
   RET_ON_NBERR(nberr)

   call integer_nalloc_1d(this%inperiod_N_E_F, entities%n_entities, 'inperiod_N_E_F', nberr);
   RET_ON_NBERR(nberr)


end subroutine
