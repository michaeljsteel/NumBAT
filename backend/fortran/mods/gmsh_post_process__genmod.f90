        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:47:58 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GMSH_POST_PROCESS__genmod
          INTERFACE 
            SUBROUTINE GMSH_POST_PROCESS(PLOT_VAL,E_H_FIELD,NVAL,NEL,NPT&
     &,NNODES,TABLE_NOD,TYPE_EL,NB_TYP_EL,N_EFF,X,VAL_CMPLX,SOL,VISITE, &
     &GMSH_FILE_POS,DIR_NAME,Q_AVERAGE,PLOT_REAL,PLOT_IMAG,PLOT_ABS)
              INTEGER(KIND=8) :: NB_TYP_EL
              INTEGER(KIND=8) :: NNODES
              INTEGER(KIND=8) :: NPT
              INTEGER(KIND=8) :: NEL
              INTEGER(KIND=8) :: NVAL
              INTEGER(KIND=8) :: PLOT_VAL
              INTEGER(KIND=8) :: E_H_FIELD
              INTEGER(KIND=8) :: TABLE_NOD(NNODES,NEL)
              INTEGER(KIND=8) :: TYPE_EL(NEL)
              COMPLEX(KIND=8) :: N_EFF(NB_TYP_EL)
              REAL(KIND=8) :: X(2,NPT)
              COMPLEX(KIND=8) :: VAL_CMPLX(NVAL)
              COMPLEX(KIND=8) :: SOL(3,NNODES+7,NVAL,NEL)
              INTEGER(KIND=8) :: VISITE(NPT)
              CHARACTER(*) :: GMSH_FILE_POS
              CHARACTER(*) :: DIR_NAME
              INTEGER(KIND=8) :: Q_AVERAGE
              INTEGER(KIND=8) :: PLOT_REAL
              INTEGER(KIND=8) :: PLOT_IMAG
              INTEGER(KIND=8) :: PLOT_ABS
            END SUBROUTINE GMSH_POST_PROCESS
          END INTERFACE 
        END MODULE GMSH_POST_PROCESS__genmod
