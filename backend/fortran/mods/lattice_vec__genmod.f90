        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:48:01 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LATTICE_VEC__genmod
          INTERFACE 
            SUBROUTINE LATTICE_VEC(NPT,X,LAT_VECS,DEBUG)
              INTEGER(KIND=8) :: NPT
              REAL(KIND=8) :: X(2,NPT)
              REAL(KIND=8) :: LAT_VECS(2,2)
              INTEGER(KIND=8) :: DEBUG
            END SUBROUTINE LATTICE_VEC
          END INTERFACE 
        END MODULE LATTICE_VEC__genmod
