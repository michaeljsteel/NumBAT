        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:47:58 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_COORD_P3__genmod
          INTERFACE 
            SUBROUTINE GET_COORD_P3(NEL,NPT,NNODES,N_DDL,TABLE_NOD,     &
     &TYPE_NOD,TABLE_N_E_F,TYPE_N_E_F,X,X_N_E_F,VISITE)
              INTEGER(KIND=8) :: N_DDL
              INTEGER(KIND=8) :: NNODES
              INTEGER(KIND=8) :: NPT
              INTEGER(KIND=8) :: NEL
              INTEGER(KIND=8) :: TABLE_NOD(NNODES,NEL)
              INTEGER(KIND=8) :: TYPE_NOD(NPT)
              INTEGER(KIND=8) :: TABLE_N_E_F(14,NEL)
              INTEGER(KIND=8) :: TYPE_N_E_F(2,N_DDL)
              REAL(KIND=8) :: X(2,NPT)
              REAL(KIND=8) :: X_N_E_F(2,N_DDL)
              INTEGER(KIND=8) :: VISITE(N_DDL)
            END SUBROUTINE GET_COORD_P3
          END INTERFACE 
        END MODULE GET_COORD_P3__genmod
