        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 19:12:05 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PHOTOELASTIC_INT__genmod
          INTERFACE 
            SUBROUTINE PHOTOELASTIC_INT(NVAL_EM_P,NVAL_EM_S,NVAL_AC,    &
     &IVAL1,IVAL2,IVAL3,NEL,NPT,NNODES,TABLE_NOD,TYPE_EL,X,NB_TYP_EL,   &
     &P_TENSOR,BETA_AC,SOLN_EM_P,SOLN_EM_S,SOLN_AC,EPS_LST,DEBUG,OVERLAP&
     &)
              INTEGER(KIND=8) :: NB_TYP_EL
              INTEGER(KIND=8) :: NNODES
              INTEGER(KIND=8) :: NPT
              INTEGER(KIND=8) :: NEL
              INTEGER(KIND=8) :: NVAL_AC
              INTEGER(KIND=8) :: NVAL_EM_S
              INTEGER(KIND=8) :: NVAL_EM_P
              INTEGER(KIND=8) :: IVAL1
              INTEGER(KIND=8) :: IVAL2
              INTEGER(KIND=8) :: IVAL3
              INTEGER(KIND=8) :: TABLE_NOD(NNODES,NEL)
              INTEGER(KIND=8) :: TYPE_EL(NEL)
              REAL(KIND=8) :: X(2,NPT)
              COMPLEX(KIND=8) :: P_TENSOR(3,3,3,3,NB_TYP_EL)
              COMPLEX(KIND=8) :: BETA_AC
              COMPLEX(KIND=8) :: SOLN_EM_P(3,NNODES,NVAL_EM_P,NEL)
              COMPLEX(KIND=8) :: SOLN_EM_S(3,NNODES,NVAL_EM_S,NEL)
              COMPLEX(KIND=8) :: SOLN_AC(3,NNODES,NVAL_AC,NEL)
              COMPLEX(KIND=8) :: EPS_LST(NB_TYP_EL)
              INTEGER(KIND=8) :: DEBUG
              COMPLEX(KIND=8) :: OVERLAP(NVAL_EM_S,NVAL_EM_P,NVAL_AC)
            END SUBROUTINE PHOTOELASTIC_INT
          END INTERFACE 
        END MODULE PHOTOELASTIC_INT__genmod
