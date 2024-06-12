        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:53:05 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PERIODIC_COND__genmod
          INTERFACE 
            SUBROUTINE PERIODIC_COND(I_COND,N_DDL,NEQ,TYPE_N_E_F,       &
     &IP_PERIOD_E_F,INEQ,DEBUG)
              INTEGER(KIND=8) :: N_DDL
              INTEGER(KIND=8) :: I_COND
              INTEGER(KIND=8) :: NEQ
              INTEGER(KIND=8) :: TYPE_N_E_F(2,N_DDL)
              INTEGER(KIND=8) :: IP_PERIOD_E_F(N_DDL)
              INTEGER(KIND=8) :: INEQ(3,N_DDL)
              INTEGER(KIND=8) :: DEBUG
            END SUBROUTINE PERIODIC_COND
          END INTERFACE 
        END MODULE PERIODIC_COND__genmod
