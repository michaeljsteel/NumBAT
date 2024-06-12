        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:47:54 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CURVED_ELEM_TRI__genmod
          INTERFACE 
            SUBROUTINE CURVED_ELEM_TRI(NNODES,XEL,INFO_CURVED,TMP)
              INTEGER(KIND=8) :: NNODES
              REAL(KIND=8) :: XEL(2,NNODES)
              INTEGER(KIND=8) :: INFO_CURVED
              REAL(KIND=8) :: TMP
            END SUBROUTINE CURVED_ELEM_TRI
          END INTERFACE 
        END MODULE CURVED_ELEM_TRI__genmod
