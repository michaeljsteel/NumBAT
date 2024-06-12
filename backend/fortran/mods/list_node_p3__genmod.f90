        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:48:01 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LIST_NODE_P3__genmod
          INTERFACE 
            SUBROUTINE LIST_NODE_P3(NEL,NPT,NNODES,N_EDGE,NPT_P3,       &
     &TABLE_NOD,TABLE_N_E_F,VISITE)
              INTEGER(KIND=8) :: NNODES
              INTEGER(KIND=8) :: NPT
              INTEGER(KIND=8) :: NEL
              INTEGER(KIND=8) :: N_EDGE
              INTEGER(KIND=8) :: NPT_P3
              INTEGER(KIND=8) :: TABLE_NOD(NNODES,NEL)
              INTEGER(KIND=8) :: TABLE_N_E_F(14,NEL)
              INTEGER(KIND=8) :: VISITE(NPT)
            END SUBROUTINE LIST_NODE_P3
          END INTERFACE 
        END MODULE LIST_NODE_P3__genmod
