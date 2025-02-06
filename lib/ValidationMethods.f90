module validation_methods_m
    use kind_parameters_m
    implicit none
    private
    interface validate_array
        module procedure validate_int_array
        module procedure validate_real_array
    end interface

    interface validate_array_strict
        module procedure validate_int_array_strict
        module procedure validate_real_array_strict
    end interface

    public :: validate_array, validate_array_strict

    contains

logical function validate_real_array(arr, arr_shape) result(valid)
    real(DP),allocatable,intent(in) :: arr(..)
    integer(IP),intent(in) :: arr_shape(rank(arr))

    valid = .true.
    if ( all(shape(arr) /= arr_shape) ) then
        valid = .false.
    end if

end function

logical function validate_int_array(arr, arr_shape) result(valid)
    integer(IP),allocatable,intent(in) :: arr(..)
    integer(IP),intent(in) :: arr_shape(rank(arr))

    valid = .true.
    if ( all(shape(arr) /= arr_shape) ) then
        valid = .false.
    end if

end function

!==========================================

logical function validate_real_array_strict(arr, arr_shape) result(valid)
    real(DP),allocatable,intent(in) :: arr(..)
    integer(IP),intent(in) :: arr_shape(rank(arr))

    valid = .true.
    if ( .not. allocated(arr) ) then
        valid = .false.
    else
        valid = validate_real_array(arr, arr_shape)
    end if

end function

logical function validate_int_array_strict(arr, arr_shape) result(valid)
    integer(IP),allocatable,intent(in) :: arr(..)
    integer(IP),intent(in) :: arr_shape(rank(arr))

    valid = .true.
    if ( .not. allocated(arr) ) then
        valid = .false.
    else
        valid = validate_int_array(arr, arr_shape)
    end if

end function
    
end module validation_methods_m