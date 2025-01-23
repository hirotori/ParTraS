module field_updater_m
    use kind_parameters_m
    use flow_field_m
    use base_importer_m
    implicit none
    private

    !NOTE: クラスではなくprotectedなモジュール変数の構造体として運用するべきか. 

    type field_updater_t
        !! helper class for updating flow field
        private
        class(ugrid_importer_t),allocatable :: importer_
            !! importer
        type(cell_type_t) :: ct_
            !! cell type definitions
        type(face_vertex_def_t) :: fv_def_
            !! face-vertex connectivity definition
        real(DP) :: dt_flow_
            !! time stepping size in flow field
        real(DP) :: dt_part_
            !! time stepping size in particle simulation

        integer :: interval_
            !! time interval at which the flow field was dumped
            !! if interval = -1, no update is conducted
        character(:),allocatable :: basename_
        character(:),allocatable :: ext_
        integer:: pad_
        logical field_only_
        logical verts_only_

        logical is_assigned_

        integer :: counter_ = 1
        contains
        procedure,non_overridable,public :: assigned
        procedure,non_overridable,public :: construct_field_updater
        procedure,non_overridable,public :: disable_updater
        procedure,non_overridable,public :: update_field
        procedure,non_overridable,public :: should_update_field
        procedure,non_overridable,public :: get_current_step
        procedure,non_overridable,public :: get_current_filename
    end type

    integer,parameter,private :: NO_UPDATE = -1
        !! updater executes no update

    type(field_updater_t) mv_field_updater

    public mv_field_updater

contains

logical function assigned(this)
    class(field_updater_t),intent(in) :: this

    assigned = this%is_assigned_

end function

subroutine construct_field_updater(this, importer, cell_type_def, face_vert_def, dt_f, dt_p, &
                                   basename, pad, ext, field_only, verts_only, interval)
    !! create field_updater object
    class(field_updater_t),intent(inout) :: this
    class(ugrid_importer_t),intent(in) :: importer
        !! importer used in this updater
    type(cell_type_t),intent(in) :: cell_type_def
        !! cell type definitions for importer
    type(face_vertex_def_t),intent(in) :: face_vert_def
        !! face-vertex connectivity definitions for importer
    real(DP),intent(in) :: dt_f
        !! time stepping size for fluid flow
    real(DP),intent(in) :: dt_p
        !! time stepping size for particle simulator
    character(*),intent(in) :: basename
        !! filename for flow field without extention
    integer,intent(in) :: pad
        !! Number of zero-fillings in file name
    character(*),intent(in) :: ext
        !! extension of filename (for example: .vtk, .fph)
    logical,intent(in) :: field_only
        !! update cell velocity only
    logical,intent(in) :: verts_only
        !! update cell vertices only. This flag is valid only if field_only is false.
    integer,intent(in) :: interval
        !! Output interval of flow field file

    allocate(this%importer_, source=importer)
    this%ct_ = cell_type_def
    this%fv_def_ = face_vert_def

    this%dt_flow_ = dt_f
    this%dt_part_ = dt_p

    if ( dt_f < dt_p) then
        print "(A)", &
        "field_updater_t/construct_field_updater::WARNING:: the case dt_f < dt_p not supported. "//new_line("C")//&
        "In this case, updater may not work correctly."
    end if

    this%basename_ = basename
    this%pad_ = pad
    this%ext_ = ext

    this%field_only_ = field_only
    this%verts_only_ = verts_only
    this%interval_ = interval

    this%counter_ = 1

    this%is_assigned_ = .true.
end subroutine

subroutine disable_updater(this)
    !! Forcibly disable updates
    !! this process can be called before initializing this object
    class(field_updater_t),intent(inout) :: this
    this%interval_ = NO_UPDATE
    this%is_assigned_ = .true.

    if ( allocated(this%importer_) ) deallocate(this%importer_)
    if ( allocated(this%basename_))  deallocate(this%basename_)
    if ( allocated(this%ext_))  deallocate(this%ext_)
end subroutine

subroutine update_field(this, ncyc)
    !! update flow field. 
    class(field_updater_t),intent(inout) :: this
    integer,intent(in) :: ncyc

    type(ugrid_struct_t) ugrid_
    character(128) fname_

    if ( this%should_update_field(ncyc) ) then

        fname_ = get_filename(this%basename_, this%counter_*this%interval_, this%pad_, this%ext_)
        print "('field_updater_t/update_field::notice:: flow field is updated by ""',A, '""')", fname_

        call this%importer_%open_file(trim(fname_))

        call this%importer_%read_file(ugrid_, .true.)
        
        call update_flow_field(ugrid_, this%ct_, this%fv_def_, this%field_only_, this%verts_only_)
    
        call this%importer_%close()

        this%counter_ = this%counter_ + 1

    end if

end subroutine

pure integer function get_current_step(this)
    !! get current timestep of flow field
    class(field_updater_t),intent(in) :: this

    get_current_step = this%interval_*this%interval_

end function

pure function get_current_filename(this) result(fname)
    !! get current timestep of flow field
    class(field_updater_t),intent(in) :: this

    character(128) fname

    fname = get_filename(this%basename_, this%get_current_step(), this%pad_, this%ext_)

end function

logical function should_update_field(this, ncyc)
    !! notify whether the update timing is appropriate at the current time step `ncyc`.
    class(field_updater_t),intent(inout) :: this
    integer,intent(in) :: ncyc
        !! current cycle in a particle simulation

    if ( this%interval_ == NO_UPDATE ) then
        should_update_field = .false.
        return
    end if

    !NOTE: dt_part > dt_flow となるケースは想定していないため, 永遠に更新されない. 
    !NOTE: dt_part/dt_flowが割り切れない場合, 現在の時刻 nw*dt_part に最も近い時刻の流れ場で更新される. 
    !NOTE: 例えば dt_part = 0.03, dt_flow = 0.1である場合, ncyc=4で更新される. 
    if ( this%counter_ == int(this%dt_part_*ncyc/this%dt_flow_)) then
        should_update_field = .true.
    end if

end function
! end of definition for updater
! ========================================
! helper functions
!
pure function get_filename(basename, number, pad, ext) result(filename)
    !! get file name
    character(*),intent(in) :: basename
    integer,intent(in) :: number
    integer,intent(in) :: pad
    character(*),intent(in) :: ext

    character(128) filename
    character(:),allocatable :: pad_char_, fmt_

    if ( pad == 0 ) then
        fmt_ = "('"//basename//"',i0,'"//ext//"')"
        write(filename, fmt_) number
    else if (pad > 0) then
        allocate(character(get_digits_of(pad)):: pad_char_)
        write(pad_char_, "(i0)") pad
        fmt_ = "('"//basename//"',I"//pad_char_//"."//pad_char_//",'"//ext//"')"    
        write(filename, fmt_) number
    else
        error stop "update_m/get_filename::error:: negative pad number"
    end if

end function

pure integer function get_digits_of(num)
    integer,intent(in) :: num
    get_digits_of = int(log10(dble(num)))+1
end function

end module field_updater_m