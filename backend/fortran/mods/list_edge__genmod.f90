        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:48:01 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LIST_EDGE__genmod
          INTERFACE 
            SUBROUTINE LIST_EDGE(NEL,NPT,NNODES,N_EDGE,TYPE_NOD,        &
     &TABLE_NOD,TABLE_EDGE,TABLE_EDGE_FACE,VISITE)
              INTEGER(KIND=8) :: NNODES
              INTEGER(KIND=8) :: NPT
              INTEGER(KIND=8) :: NEL
              INTEGER(KIND=8) :: N_EDGE
              INTEGER(KIND=8) :: TYPE_NOD(NPT)
              INTEGER(KIND=8) :: TABLE_NOD(NNODES,NEL)
              INTEGER(KIND=8) :: TABLE_EDGE(4,NPT)
              INTEGER(KIND=8) :: TABLE_EDGE_FACE(14,NEL)
              INTEGER(KIND=8) :: VISITE(NPT)
            END SUBROUTINE LIST_EDGE
          END INTERFACE 
        END MODULE LIST_EDGE__genmod
