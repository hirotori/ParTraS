module struct_array_m
    !! module for array structure with 3 members
    use kind_parameters_m
    implicit none

    type scalar3_t
        !! struct with 3 components.
        !! Instances should be initialized via default constructer (e.g., "scalar3_t(1.0, 1.0, 1.0)")
        !! Basic operators (+, -, *) are overloaded to provide basic vector operations.
        !! 
        !! @note
        !!
        !! Note that operator "*" is only valid for multilpication by scalar.
        !! For cross product operation, use ".cross." operator.
        !!
        !! @endnote
        !!
        real(DP) x
            !! 1st component
        real(DP) y
            !! 2nd component
        real(DP) z
            !! 3rd component
    end type

    type,extends(scalar3_t) :: scalar4_t
        integer(IP) m
    end type

    !<<<<< constructer >>>>>
    interface scalar3_t
        module procedure scalar3_from_3_
        module procedure scalar3_from_array_
    end interface

    !<<<<< operator overload >>>>>>

    interface operator(+)
        module procedure add_
    end interface

    interface operator(-)
        module procedure subtract_
    end interface

    interface operator(*)
        module procedure multi_double_
    end interface

    interface operator(.cross.)
        module procedure multi_scalar3_
    end interface

    interface operator(.dot.)
        module procedure dot_
    end interface

    interface operator(.orth.)
        module procedure ortho_
    end interface
    
contains

    type(scalar3_t) function scalar3_from_3_(x, y, z) result(obj)
        real(DP),intent(in) :: x, y, z

        obj%x = x
        obj%y = y
        obj%z = z

    end function

    type(scalar3_t) function scalar3_from_array_(x3) result(obj)
        real(DP),intent(in) :: x3(3)

        obj%x = x3(1)
        obj%y = x3(2)
        obj%z = x3(3)

    end function


    pure type(scalar3_t) function add_(this, b) result(obj)
        class(scalar3_t),intent(in) :: this
        class(scalar3_t),intent(in) :: b
        
        obj%x = this%x + b%x
        obj%y = this%y + b%y
        obj%z = this%z + b%z

    end function

    pure type(scalar3_t) function subtract_(this, b) result(obj)
        class(scalar3_t),intent(in) :: this
        class(scalar3_t),intent(in) :: b
        
        obj%x = this%x - b%x
        obj%y = this%y - b%y
        obj%z = this%z - b%z

    end function

    pure type(scalar3_t) function multi_double_(this, b) result(obj)
        class(scalar3_t),intent(in) :: this
        real(DP),intent(in) :: b
        
        obj%x = b*this%x
        obj%y = b*this%y
        obj%z = b*this%z

    end function

    pure type(scalar3_t) function multi_scalar3_(this, b) result(obj)
        !! cross product operation
        class(scalar3_t),intent(in) :: this
        class(scalar3_t),intent(in) :: b
        
        obj%x = this%y*b%z - this%z*b%y
        obj%y = this%z*b%x - this%x*b%z
        obj%z = this%x*b%y - this%y*b%x

    end function


    pure real(DP) function dot_(this, b) result(res)
        class(scalar3_t),intent(in) :: this
        class(scalar3_t),intent(in) :: b

        res = this%x*b%x + this%y*b%y + this%z*b%z
    
    end function

    pure type(scalar3_t) function ortho_(this, b) result(obj)
        !! orthogonal vector with b
        class(scalar3_t),intent(in) :: this
        class(scalar3_t),intent(in) :: b
        
        obj = this - b*(this .dot. b)

    end function

    pure real(DP) function L2norm(a) result(obj)
        !! compute L2-norm.
        class(scalar3_t),intent(in) :: a

        obj = sqrt(a .dot. a)

    end function


    subroutine to_array(arr_of_struct, array_2d)
        !! convert array of structure to 2d array.
        class(scalar3_t),intent(inout) :: arr_of_struct(:)
            !! array of structure.
        real(DP),allocatable,intent(inout) :: array_2d(:,:)
            !! returned 2d array. It must be unallocated.

        integer astat, n_aos, i

        !nothing done when array_2d is already allocated
        if (allocated(array_2d)) return

        n_aos = size(arr_of_struct)

        allocate(array_2d(3,n_aos), stat=astat)

        if ( astat == 0 ) then
            do i = 1, n_aos
                array_2d(1,i) = arr_of_struct(i)%x
                array_2d(2,i) = arr_of_struct(i)%y
                array_2d(3,i) = arr_of_struct(i)%z
            end do
        else
            return    
        end if
        
    end subroutine

end module struct_array_m