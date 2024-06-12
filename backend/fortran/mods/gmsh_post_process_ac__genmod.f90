        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:47:59 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GMSH_POST_PROCESS_AC__genmod
          INTERFACE 
            SUBROUTINE GMSH_POST_PROCESS_AC(PLOT_VAL,NVAL,NEL,NPT,NNODES&
     &,TABLE_NOD,TYPE_EL,X,VAL_CMPLX,SOL,SOL_AVG,VISITE,GMSH_FILE_POS,  &
     &DIR_NAME,D_IN_NM,DEBUG)
              INTEGER(KIND=8) :: NNODES
              INTEGER(KIND=8) :: NPT
              INTEGER(KIND=8) :: NEL
              INTEGER(KIND=8) :: NVAL
              INTEGER(KIND=8) :: PLOT_VAL
              INTEGER(KIND=8) :: TABLE_NOD(NNODES,NEL)
              INTEGER(KIND=8) :: TYPE_EL(NEL)
              REAL(KIND=8) :: X(2,NPT)
              COMPLEX(KIND=8) :: VAL_CMPLX(NVAL)
              COMPLEX(KIND=8) :: SOL(3,NNODES,NVAL,NEL)
              COMPLEX(KIND=8) :: SOL_AVG(3,NPT)
              INTEGER(KIND=8) :: VISITE(NPT)
              CHARACTER(*) :: GMSH_FILE_POS
              CHARACTER(*) :: DIR_NAME
              REAL(KIND=8) :: D_IN_NM
              INTEGER(KIND=8) :: DEBUG
            END SUBROUTINE GMSH_POST_PROCESS_AC
          END INTERFACE 
        END MODULE GMSH_POST_PROCESS_AC__genmod
