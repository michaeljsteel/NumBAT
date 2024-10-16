subroutine PeriodicBCs_allocate(this, mesh_raw, entities, errco, emsg)


   class(PeriodicBCs) :: this

   type(MeshRaw) :: mesh_raw
   type(MeshEntities) :: entities

   integer(8),  intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg


   call integer_alloc_1d(this%iperiod_N, mesh_raw%n_msh_pts, 'iperiod_N', errco, emsg);
   RETONERROR(errco)

   call integer_alloc_1d(this%iperiod_N_E_F, entities%n_entities, &
      'iperiod_N_E_F', errco, emsg); RETONERROR(errco)

   call integer_alloc_1d(this%inperiod_N, mesh_raw%n_msh_pts, 'inperiod_N', errco, emsg);
   RETONERROR(errco)

   call integer_alloc_1d(this%inperiod_N_E_F, entities%n_entities, &
      'inperiod_N_E_F', errco, emsg); RETONERROR(errco)



end subroutine
