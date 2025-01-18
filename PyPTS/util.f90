module util_m
    use iso_c_binding
    implicit none
    
contains

function fstring(c_strings) result(cha)
    !! convert c string to fortran character
    character(c_char),dimension(*),intent(in) :: c_strings

    integer nchar
    character(:),allocatable :: cha

    ! get length of c_strings
    nchar = 1
    do while(c_strings(nchar) /= c_null_char)
      nchar = nchar + 1
    end do

    ! convert array of characters to character variable
    cha = repeat(" ", nchar-1)
    cha = transfer(c_strings(1:nchar-1),cha)

end function

end module util_m