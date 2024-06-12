        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:53:06 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PERIODIC_NODE__genmod
          INTERFACE 
            SUBROUTINE PERIODIC_NODE(NEL,NPT,NNODES,TYPE_NOD,X,IP_PERIOD&
     &,N_PERIOD,TABLE_NOD,LAT_VECS)
              INTEGER(KIND=8) :: NNODES
              INTEGER(KIND=8) :: NPT
              INTEGER(KIND=8) :: NEL
              INTEGER(KIND=8) :: TYPE_NOD(NPT)
              REAL(KIND=8) :: X(2,NPT)
              INTEGER(KIND=8) :: IP_PERIOD(NPT)
              INTEGER(KIND=8) :: N_PERIOD(NPT)
              INTEGER(KIND=8) :: TABLE_NOD(NNODES,NEL)
              REAL(KIND=8) :: LAT_VECS(2,2)
            END SUBROUTINE PERIODIC_NODE
          END INTERFACE 
        END MODULE PERIODIC_NODE__genmod
