program test
    use kind_parameters_m
    use struct_array_m
    use test_util_m
    ! use struct_array_4_m
    implicit none


    type(scalar3_t) a, b, c
    type(scalar4_t) x
    real(DP) :: k = 2.0


    a = scalar3_t(1.0, 2.0, 3.0)
    b = scalar3_t(4.0, 5.0, 6.0)

    ! summation
    c = a + b
    if ( assert_equal(c%x, 5.0) /= 0) error stop  "1"
    if ( assert_equal(c%y, 7.0) /= 0) error stop  "2"
    if ( assert_equal(c%z, 9.0) /= 0) error stop  "3"

    ! subtraction
    c = a - b
    if ( assert_equal(c%x, -3.0) /= 0) error stop  "2-1"
    if ( assert_equal(c%y, -3.0) /= 0) error stop  "2-2"
    if ( assert_equal(c%z, -3.0) /= 0) error stop  "2-3"

    !multiply (scalar)
    c = a*k
    if ( assert_equal(c%x, 2) /= 0 ) error stop "3-1"
    if ( assert_equal(c%y, 4) /= 0 ) error stop "3-2"
    if ( assert_equal(c%z, 6) /= 0 ) error stop "3-3"

    ! dot product
    b = scalar3_t(1.0, 0.0, 0.0)

    if ( assert_equal(a.dot.b, 1) /= 0) error stop "4"

    !orthogonal vector
    c = a .orth. b
    if ( assert_equal(c%x, 0) /= 0 ) error stop "5-1"
    if ( assert_equal(c%y, 2) /= 0 ) error stop "5-2"
    if ( assert_equal(c%z, 3) /= 0 ) error stop "5-3"

    x = scalar4_t(4.0, 5.0, 6.0, 1.0)
    c = a + x
    if ( assert_equal(c%x, 5.0) /= 0) error stop  "1"
    if ( assert_equal(c%y, 7.0) /= 0) error stop  "2"
    if ( assert_equal(c%z, 9.0) /= 0) error stop  "3"
    
end program test