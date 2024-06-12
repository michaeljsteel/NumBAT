        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:47:53 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE BOUND_COND__genmod
          INTERFACE 
            SUBROUTINE BOUND_COND(I_COND,N_DDL,NEQ,TYPE_N_E_F,INEQ)
              INTEGER(KIND=8) :: N_DDL
              INTEGER(KIND=8) :: I_COND
              INTEGER(KIND=8) :: NEQ
              INTEGER(KIND=8) :: TYPE_N_E_F(2,N_DDL)
              INTEGER(KIND=8) :: INEQ(3,N_DDL)
            END SUBROUTINE BOUND_COND
          END INTERFACE 
        END MODULE BOUND_COND__genmod
