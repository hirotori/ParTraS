module base_importer_m
    use kind_parameters_m
    implicit none
    private

    integer,parameter :: IMPORT_ASCII  = 0
    integer,parameter :: IMPORT_BINARY = 1


    type cell_type_t
        !! cell type definition. 
        integer :: cell_type_tetra
        integer :: cell_type_wedge
        integer :: cell_type_hexa
        integer :: cell_type_pyram
    end type

    type face_vertex_def_t
        !! face-vertex connectivity definition. 
        !! 1st dim is [n, v0, v1, v2, ..., vn] where n is the number of vertex on a face 
        !! 2nd dim of an array is face id. 
        !! NOTE: vn = -1 means the vertex vn is not found in the face
        integer :: face_def_tetra(3+1,4)
            !! definition for tetrahedron
        integer :: face_def_wedge(4+1,5)
            !! definition for wedge (prism)
        integer :: face_def_hexa(4+1,6)
            !! definition for hexahedron
        integer :: face_def_pyram(4+1,5)
            !! definition for pyramid
    end type

    type ugrid_struct_t
        !! basical unstructured grid data. 
        !! This is used for importing data from external file or define data manually
        integer(IP) ncell
        integer(IP) nvert
        real(DP),allocatable :: verts(:,:)
        integer(IP),allocatable :: conns(:)
        integer(IP),allocatable :: offsets(:)
        integer(IP),allocatable :: cell2verts(:,:)
        integer(IP),allocatable :: cell_types(:)
        real(DP),allocatable :: cell_velocity(:,:)

        ! CFDデータのサポート. 
        ! 特定のCFDデータは面ーセル接続関係をファイルに書き込む場合がある. 
        ! それらのデータを読み, 直接`flow_field_t`に渡せるようにメンバを追加. 
        ! ファイルにこれらのデータがない場合は割り付けてはならない.
        ! `flow_field_t`はこれらが割り付けられていない場合自分で構築する.  
        integer(IP),allocatable :: face2cells(:,:)
        integer(IP),allocatable :: face2verts(:,:)

        !幾何量
        real(DP),allocatable :: cell_centers(:,:)
        real(DP),allocatable :: face_centers(:,:)
        real(DP),allocatable :: face_normals(:,:)
    end type

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

    public ugrid_struct_t, ugrid_importer_t, &
           cell_type_t, face_vertex_def_t, &
           delete_ugrid, &
           IMPORT_ASCII, &
           IMPORT_BINARY !, &
        !    FIELD_NAME_LEN

contains

subroutine delete_ugrid(this)
    type(ugrid_struct_t),intent(inout) :: this

    if (allocated(this%verts))        deallocate(this%verts)
    if (allocated(this%conns))        deallocate(this%conns)
    if (allocated(this%offsets))      deallocate(this%offsets)
    if (allocated(this%cell_types))   deallocate(this%cell_types)
    if (allocated(this%cell_velocity))deallocate(this%cell_velocity)

    ! this%filename = FILENAME_NOT_ASSIGNED
    if ( allocated(this%face2cells) ) deallocate(this%face2cells) 
    if ( allocated(this%face2verts) ) deallocate(this%face2verts) 
    if ( allocated(this%face_centers) ) deallocate(this%face_centers) 
    if ( allocated(this%cell_centers) ) deallocate(this%cell_centers) 
end subroutine

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