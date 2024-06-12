        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:48:07 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MOVING_BOUNDARY__genmod
          INTERFACE 
            SUBROUTINE MOVING_BOUNDARY(NVAL_EM_P,NVAL_EM_S,NVAL_AC,IVAL1&
     &,IVAL2,IVAL3,NEL,NPT,NNODES,TABLE_NOD,TYPE_EL,X,NB_TYP_EL,        &
     &TYP_SELECT_IN,TYP_SELECT_OUT,SOLN_EM_P,SOLN_EM_S,SOLN_AC,EPS_LST, &
     &DEBUG,OVERLAP)
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
              INTEGER(KIND=8) :: TABLE_NOD(6,NEL)
              INTEGER(KIND=8) :: TYPE_EL(NEL)
              REAL(KIND=8) :: X(2,NPT)
              INTEGER(KIND=8) :: TYP_SELECT_IN
              INTEGER(KIND=8) :: TYP_SELECT_OUT
              COMPLEX(KIND=8) :: SOLN_EM_P(3,NNODES,NVAL_EM_P,NEL)
              COMPLEX(KIND=8) :: SOLN_EM_S(3,NNODES,NVAL_EM_S,NEL)
              COMPLEX(KIND=8) :: SOLN_AC(3,NNODES,NVAL_AC,NEL)
              COMPLEX(KIND=8) :: EPS_LST(NB_TYP_EL)
              INTEGER(KIND=4) :: DEBUG
              COMPLEX(KIND=8) :: OVERLAP(NVAL_EM_S,NVAL_EM_P,NVAL_AC)
            END SUBROUTINE MOVING_BOUNDARY
          END INTERFACE 
        END MODULE MOVING_BOUNDARY__genmod
