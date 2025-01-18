module id_list_m
    implicit none
    integer,parameter,private :: INITIAL_SIZE = 30

    type, public :: id_list_t
        !! a simple integer list object
        !! enabling only four methods: add, push_back, resize and getter.
        !! 
        !! Referenced from: https://qiita.com/osada-yum/items/e2e5157602b927980080
        !!
        integer,private :: size_
            !! current index of the data
        integer,private :: capa_
            !! capacity of the data
        integer,allocatable,private :: data_(:)
            !! internal data
        contains
        procedure push_back_single
        procedure push_back_arr
        generic :: push_back => push_back_single, push_back_arr
        procedure get_at
        procedure get_data
        procedure resize
        procedure current_size
        procedure current_capacity
        final     dispose
    end type

    ! override default constructer
    interface id_list_t
        module procedure create_id_list_t
    end interface
    
contains
    
type(id_list_t) function create_id_list_t(n) result(obj)
    !! create a list object
    integer,intent(in),optional :: n
        !! initial size of the data

    integer capa_

    if ( present(n) ) then
        capa_ = n
    else
        capa_ = INITIAL_SIZE
    end if

    call obj%resize(capa_)

end function

subroutine dispose(this)
    !! destructer of id_list_t
    type(id_list_t),intent(inout) :: this

    if (allocated(this%data_)) deallocate(this%data_)
    this%size_ = 0
    this%capa_ = 0
end subroutine

pure subroutine push_back_single(this, val)
    !! add a value to the list
    class(id_list_t),intent(inout) :: this
    integer,intent(in) :: val

    if ( this%size_ == this%capa_ ) then
        call this%resize(this%size_*2)
    end if

    this%size_ = this%size_ + 1
    this%data_(this%size_) = val

end subroutine

pure subroutine push_back_arr(this, val)
    !! add values to the list
    class(id_list_t),intent(inout) :: this
    integer,intent(in) :: val(:)

    integer next_

    next_ = this%size_ + size(val) ! next size

    if ( next_ >= this%capa_ ) then
        call this%resize(next_*2)
    end if

    this%data_(this%size_+1:next_) = val(:)
    this%size_ = next_

end subroutine

impure integer function get_at(this, i) result(val)
    !! get a value in the list at the given index `i`
    class(id_list_t),intent(in) :: this
    integer,intent(in) :: i

    if ( i > this%size_ ) then
        error stop "id_list_t::index error: index out of range"
    end if

    val = this%data_(i)

end function


pure function get_data(this) result(data_)
    !! get internal data array as a deep copy.
    class(id_list_t),intent(in) :: this
    integer :: data_(this%size_)

    data_ = this%data_(1:this%size_)

end function


pure subroutine resize(this, n)
    !! resize a list
    class(id_list_t),intent(inout) :: this
    integer,intent(in) :: n
        !! new size of the data

    integer,allocatable :: tmp_(:)
    integer max_

    if ( .not. allocated(this%data_) ) then
        allocate(this%data_(n))
        this%capa_ = n
        this%size_ = 0
    else
        max_ = min(this%size_, n)
        allocate(tmp_(n))
        tmp_(1:max_) = this%data_(1:max_)
        call move_alloc(tmp_, this%data_)
        this%capa_ = n
        this%size_ = max_
    end if

end subroutine

pure integer function current_size(this)
    class(id_list_t),intent(in) :: this

    current_size = this%size_

end function

pure integer function current_capacity(this)
    class(id_list_t),intent(in) :: this

    current_capacity = this%capa_

end function

end module id_list_m