module test_util_m
    use kind_parameters_m
    implicit none

    real(DP),parameter :: F64_EQUAL_TORRELANCE = 1.0e-12
    
    interface assert_equal
        module procedure assert_equal_i32
        module procedure assert_equal_f64
        module procedure assert_equal_f32_and_f64
        module procedure assert_equal_f32_and_i32
    end interface

contains

pure integer(IP) function assert_equal_i32(x, y) result(code)
    integer(IP),intent(in) :: x, y
    
    code = merge(0, 1, x == y)

end function


pure integer(IP) function assert_equal_f64(x, y) result(code)
    real(DP),intent(in) :: x, y
    
    code = merge(0, 1, abs(x - y) <= F64_EQUAL_TORRELANCE)

end function

pure integer(IP) function assert_equal_f32_and_f64(x, y) result(code)
    real(DP),intent(in) :: x
    real(SP),intent(in) :: y
    
    code = assert_equal_f64(x, real(y, kind=DP))

end function

pure integer(IP) function assert_equal_f32_and_i32(x, y) result(code)
    real(DP),intent(in) :: x
    integer(IP),intent(in) :: y
    
    code = assert_equal_f64(x, real(y, kind=DP))

end function


end module test_util_m