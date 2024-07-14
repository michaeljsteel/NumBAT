
! This procedure will not be exported and is not
! directly callable by users of this library.

module modfoo

implicit none
private
public :: mes_func

contains

integer(8) function internal_function()
    internal_function = 0
end function internal_function

integer(8) function mes_func()
    mes_func = internal_function()
end function mes_func

end module modfoo
