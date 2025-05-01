#include "numbat_decl.h"

module class_PeriodicBCs

    use numbatmod
    use alloc
    use class_MeshRawEM

    private

    type, public :: PeriodicBCs

    integer(8), dimension(:), allocatable :: iperiod_N
    integer(8), dimension(:), allocatable :: iperiod_N_E_F
    integer(8), dimension(:), allocatable :: inperiod_N
    integer(8), dimension(:), allocatable :: inperiod_N_E_F

    contains

    procedure :: allocate => PeriodicBCs_allocate

    end type PeriodicBCs

    contains

#include "periodicbcs_class_impl.f90"


end module class_PeriodicBCs