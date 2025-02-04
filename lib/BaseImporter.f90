module base_importer_m
    use kind_parameters_m
    use unstructured_mesh_m, only: ugrid_struct_t
    implicit none
    private

    integer,parameter :: IMPORT_ASCII  = 0
    integer,parameter :: IMPORT_BINARY = 1
    type,abstract:: ugrid_importer_t
        integer,private :: unit_
        logical,private :: shift_index_
        logical,private :: ascii_ = .true. 
        logical,private :: stream_ = .true.
        contains
        procedure,non_overridable :: open_stream_on
        procedure,non_overridable :: open_ascii_on
        procedure,non_overridable :: open_stream_off
        procedure,non_overridable :: open_ascii_off
        procedure,non_overridable :: open_file
        procedure,non_overridable :: get_current_unit
        ! procedure,non_overridable :: get_filename
        procedure(read_file_imp_),deferred :: read_file
        procedure,non_overridable :: close
    end type

    public ugrid_importer_t, &
           IMPORT_ASCII, &
           IMPORT_BINARY !, &
        !    FIELD_NAME_LEN

contains

! ~~~~~~~~~~~~~~~~~ importer ~~~~~~~~~~~~~~~~~~~~ !

subroutine open_stream_on(this)
    !! enable stream mode
    class(ugrid_importer_t),intent(inout) :: this
    this%stream_ = .true.
end subroutine

subroutine open_stream_off(this)
    !! disable stream mode
    class(ugrid_importer_t),intent(inout) :: this
    this%stream_ = .false.
end subroutine

subroutine open_ascii_on(this)
    !! enable ascii mode
    class(ugrid_importer_t),intent(inout) :: this
    this%ascii_ = .true.
end subroutine

subroutine open_ascii_off(this)
    !! disable ascii mode
    class(ugrid_importer_t),intent(inout) :: this
    this%ascii_ = .false.
end subroutine

subroutine open_file(this, filename)
    !! open mesh file
    !! raise error if file not found
    class(ugrid_importer_t),intent(inout) :: this
    character(*),intent(in) :: filename
        !! filename
    
    logical exist_
    character(:),allocatable :: form_, access_

    inquire(file=filename, exist=exist_)
    if ( .not. exist_ ) then
        error stop "BaseImporter_m/ugrid_importer_t::error:: file """//trim(filename)//""" not found"
    end if
    
    if ( this%ascii_ ) then
        form_ = "formatted"
    else
        form_ = "unformatted"
    end if

    if ( this%stream_ ) then
        access_ = "stream"
    else
        access_ = "sequential"
    end if

    open(newunit=this%unit_, file=filename, action="read", form=form_, access=access_, status="old")

    ! this%filename = filename
end subroutine


subroutine read_file_imp_(this, ugrid, shift_index)
    class(ugrid_importer_t),intent(inout) :: this
    type(ugrid_struct_t),intent(out) :: ugrid
    logical,intent(in) :: shift_index
        !! shift 0-start indices to 1-start ones

end subroutine

pure integer function get_current_unit(this)
    !! get current unit number related to the mesh file
    class(ugrid_importer_t),intent(in) :: this

    get_current_unit = this%unit_

end function

! character(FIELD_NAME_LEN) function get_filename(this)
!     class(ugrid_importer_t),intent(in) :: this

!     get_filename = this%filename

! end function

subroutine close(this)
    !! close the current mesh file
    class(ugrid_importer_t),intent(inout) :: this

    close(this%unit_)
    this%unit_ = 0

end subroutine

end module base_importer_m