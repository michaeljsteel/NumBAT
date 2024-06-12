        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:47:54 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CSR_MAX_LENGTH__genmod
          INTERFACE 
            SUBROUTINE CSR_MAX_LENGTH(NEL,N_DDL,NEQ,NNODES,TABLE_N_E_F, &
     &INEQ,LB,NONZ)
              INTEGER(KIND=8) :: NEQ
              INTEGER(KIND=8) :: N_DDL
              INTEGER(KIND=8) :: NEL
              INTEGER(KIND=8) :: NNODES
              INTEGER(KIND=8) :: TABLE_N_E_F(14,NEL)
              INTEGER(KIND=8) :: INEQ(3,N_DDL)
              INTEGER(KIND=8) :: LB(NEQ+1)
              INTEGER(KIND=8) :: NONZ
            END SUBROUTINE CSR_MAX_LENGTH
          END INTERFACE 
        END MODULE CSR_MAX_LENGTH__genmod
