module test_class_mod
    implicit none
    type myclass
        integer :: i = 0
        character(32) :: char = "I am a parent"
        integer,allocatable :: vals_dynamic(:)
        integer,allocatable :: vals_dynamic_internal(:)
        integer,dimension(3) :: vals_static
        contains
        procedure init
        procedure speak
    end type

    type,extends(myclass) :: subclass
        !nothing
    end type

    contains

    subroutine init(self, N, fillval)
        class(myclass),intent(inout) :: self
        integer,intent(in) :: N
        integer,intent(in) :: fillval

        allocate(self%vals_dynamic_internal(N), source = fillval)

    end subroutine

    subroutine speak(self)
        class(myclass),intent(in) :: self

        print*, self%char
        if (allocated(self%vals_dynamic)) then
            print"('dynamic: ', *(i0,:,','))", self%vals_dynamic
        else
            print *, "No array allocated"
        endif
        if (allocated(self%vals_dynamic_internal)) then
            print"('dynamic: ', *(i0,:,','))", self%vals_dynamic_internal
        else
            print *, "No array allocated"
        endif
        print"(*(i0,:,','))", self%vals_static
        
        print "('ID = ', i0)", self%i

    end subroutine

end module test_class_mod

module test_mod
    use test_class_mod
    implicit none

    type holder
        class(myclass),pointer :: ptr_ => null()
    end type


    type system
        integer :: nclass = 10
        type(holder),allocatable :: class_holder(:)
        integer,private :: nhold_ = 0
        contains
        procedure :: init => init_
        procedure add_class
        procedure hasnext
    end type

    contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_(this)
    class(system),intent(inout) :: this
    
    allocate(this%class_holder(this%nclass))

end subroutine

logical function hasnext(this, idx) result(exists)
    class(system),intent(in) :: this
    integer,intent(in) :: idx

    exists = associated(this%class_holder(idx)%ptr_)

end function

subroutine add_class(this, a_class)
    class(system),intent(inout) :: this
    class(myclass),intent(in),target :: a_class

    if ( this%nhold_+1 > this%nclass) then
        error stop "holder error"
    end if
    this%class_holder(this%nhold_+1)%ptr_ => a_class

    this%nhold_ = this%nhold_ + 1
end subroutine

subroutine delete_class(this)
    class(system),intent(inout) :: this

    integer n

    do n = 1, this%nhold_
        this%class_holder(n)%ptr_ => null()
    end do

end subroutine

    
end module test_mod

program main
    use test_class_mod
    use test_mod
    implicit none

    type(system) sys
    type(myclass) c1
    type(subclass) c2
    integer :: i = 1

    call sys%init()

    call c1%init(3, 1)
    c1%vals_dynamic = [1,2,3]

    call c2%init(3, 2)
    c2%char = "I am a child"
    c2%i = 10
    c2%vals_dynamic = [1,2,3,4,5]

    call sys%add_class(c1)
    call sys%add_class(c2)

    block
        type(subclass) c3

        call c3%init(3,3)
        c3%char = "I am a child in block"
        c3%i = -10
        c3%vals_dynamic = [-1,-2,-3,-4,-5]
        call sys%add_class(c3)
    end block

    do while (sys%hasnext(i))
        call sys%class_holder(i)%ptr_%speak()
        i = i + 1
    end do
        
end program main