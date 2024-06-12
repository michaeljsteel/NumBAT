        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:53:06 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PERIODIC_N_E_F__genmod
          INTERFACE 
            SUBROUTINE PERIODIC_N_E_F(N_DDL,TYPE_N_E_F,X_N_E_F,         &
     &IP_PERIOD_N_E_F,N_PERIOD_N_E_F,LAT_VECS)
              INTEGER(KIND=8) :: N_DDL
              INTEGER(KIND=8) :: TYPE_N_E_F(2,N_DDL)
              REAL(KIND=8) :: X_N_E_F(2,N_DDL)
              INTEGER(KIND=8) :: IP_PERIOD_N_E_F(N_DDL)
              INTEGER(KIND=8) :: N_PERIOD_N_E_F(N_DDL)
              REAL(KIND=8) :: LAT_VECS(2,2)
            END SUBROUTINE PERIODIC_N_E_F
          END INTERFACE 
        END MODULE PERIODIC_N_E_F__genmod
