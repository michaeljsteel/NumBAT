        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:47:52 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ASMBLY__genmod
          INTERFACE 
            SUBROUTINE ASMBLY(I_COND,I_BASE,NEL,NPT,N_DDL,NEQ,NNODES,   &
     &SHIFT,BLOCH_VEC,NB_TYP_EL,PP,QQ,TABLE_NOD,TABLE_N_E_F,TYPE_EL,INEQ&
     &,IP_PERIOD_N,IP_PERIOD_E_F,X,X_N_E_F,NONZ,ROW_IND,COL_PTR,MAT1_RE,&
     &MAT1_IM,MAT2,I_WORK)
              INTEGER(KIND=8) :: NONZ
              INTEGER(KIND=8) :: NB_TYP_EL
              INTEGER(KIND=8) :: NNODES
              INTEGER(KIND=8) :: NEQ
              INTEGER(KIND=8) :: N_DDL
              INTEGER(KIND=8) :: NPT
              INTEGER(KIND=8) :: NEL
              INTEGER(KIND=8) :: I_COND
              INTEGER(KIND=8) :: I_BASE
              COMPLEX(KIND=8) :: SHIFT
              REAL(KIND=8) :: BLOCH_VEC(2)
              COMPLEX(KIND=8) :: PP(NB_TYP_EL)
              COMPLEX(KIND=8) :: QQ(NB_TYP_EL)
              INTEGER(KIND=8) :: TABLE_NOD(NNODES,NEL)
              INTEGER(KIND=8) :: TABLE_N_E_F(14,NEL)
              INTEGER(KIND=8) :: TYPE_EL(NEL)
              INTEGER(KIND=8) :: INEQ(3,N_DDL)
              INTEGER(KIND=8) :: IP_PERIOD_N(NPT)
              INTEGER(KIND=8) :: IP_PERIOD_E_F(N_DDL)
              REAL(KIND=8) :: X(2,NPT)
              REAL(KIND=8) :: X_N_E_F(2,N_DDL)
              INTEGER(KIND=8) :: ROW_IND(NONZ)
              INTEGER(KIND=8) :: COL_PTR(NEQ+1)
              REAL(KIND=8) :: MAT1_RE(NONZ)
              REAL(KIND=8) :: MAT1_IM(NONZ)
              COMPLEX(KIND=8) :: MAT2(NONZ)
              INTEGER(KIND=8) :: I_WORK(3*N_DDL)
            END SUBROUTINE ASMBLY
          END INTERFACE 
        END MODULE ASMBLY__genmod
