        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:47:57 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GEOMETRY__genmod
          INTERFACE 
            SUBROUTINE GEOMETRY(N_MSH_EL,N_MSH_PTS,NODES_PER_EL,N_TYP_EL&
     &,DIM_X,DIM_Y,TYPE_NOD,TYPE_EL,TABLE_NOD,MESH_XY,MESH_FILE,ERRCO,  &
     &EMSG)
              INTEGER(KIND=8) :: NODES_PER_EL
              INTEGER(KIND=8) :: N_MSH_PTS
              INTEGER(KIND=8) :: N_MSH_EL
              INTEGER(KIND=8) :: N_TYP_EL
              REAL(KIND=8) :: DIM_X
              REAL(KIND=8) :: DIM_Y
              INTEGER(KIND=8) :: TYPE_NOD(N_MSH_PTS)
              INTEGER(KIND=8) :: TYPE_EL(N_MSH_EL)
              INTEGER(KIND=8) :: TABLE_NOD(NODES_PER_EL,N_MSH_EL)
              REAL(KIND=8) :: MESH_XY(2,N_MSH_PTS)
              CHARACTER(LEN=1000) :: MESH_FILE
              INTEGER(KIND=8) :: ERRCO
              CHARACTER(LEN=2048) :: EMSG
            END SUBROUTINE GEOMETRY
          END INTERFACE 
        END MODULE GEOMETRY__genmod
