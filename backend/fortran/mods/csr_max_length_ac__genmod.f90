        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:47:54 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CSR_MAX_LENGTH_AC__genmod
          INTERFACE 
            SUBROUTINE CSR_MAX_LENGTH_AC(NEL,NPT,NEQ,NNODES,TABLE_NOD,  &
     &INEQ,LB,NONZ)
              INTEGER(KIND=8) :: NNODES
              INTEGER(KIND=8) :: NEQ
              INTEGER(KIND=8) :: NPT
              INTEGER(KIND=8) :: NEL
              INTEGER(KIND=8) :: TABLE_NOD(NNODES,NEL)
              INTEGER(KIND=8) :: INEQ(3,NPT)
              INTEGER(KIND=8) :: LB(NEQ+1)
              INTEGER(KIND=8) :: NONZ
            END SUBROUTINE CSR_MAX_LENGTH_AC
          END INTERFACE 
        END MODULE CSR_MAX_LENGTH_AC__genmod
