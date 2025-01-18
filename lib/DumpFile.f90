module dump_file_m
    use particle_data_m, only: mv_pdata, export_pdata_ascii, export_pdata_binary, import_pdata_ascii, import_pdata_binary
    implicit none
    private

    public :: write_out_backup, load_from_backup, writeout_vtk

contains
    
pure function get_filename(basename, timestep, ext) result(fname)
    character(*),intent(in) :: basename
    integer,intent(in) :: timestep 
    character(*),intent(in) :: ext
    character(128) fname
    character(:),allocatable :: fmt_

    fmt_ = "('"//basename//"_',i0,'"//ext//"')"
    write(fname,fmt_) timestep

end function

subroutine write_out_backup(tag, timestep, ascii)
    !! write out data into file.
    !! file must be opend and closed in this process.
    !! Existing file is replaced.
    character(*),intent(in) :: tag
    integer,intent(in) :: timestep
    logical,intent(in) :: ascii

    integer unit_
    character(128) filename
    character(:),allocatable :: form_

    filename = get_filename(tag, timestep, ".pdata")

    if ( ascii ) then
        form_ = "formatted"
    else
        form_ = "unformatted"
    end if

    open(newunit=unit_, file=filename, status="replace", action="write", &
    form=form_, access="stream", position="asis")

    if ( ascii ) then
        call export_pdata_ascii(unit_)
    else
        call export_pdata_binary(unit_)
    end if

    close(unit_)

end subroutine

subroutine load_from_backup(filename, ascii)
    !! load backup and initialize module variable of particle data `mv_pdata`
    character(*),intent(in) :: filename
    logical,intent(in) :: ascii

    integer unit
    character(:),allocatable :: form_

    if ( ascii ) then
        form_ = "formatted"
    else
        form_ = "unformatted"
    end if
    
    open(newunit=unit, file=filename, form=form_, access="stream", status="old", action="read")

    if ( ascii ) then
        call import_pdata_ascii(unit)
    else
        call import_pdata_binary(unit)
    end if

    close(unit)

end subroutine


subroutine writeout_vtk(tag, timestep, ascii)
    !! write out particle data into file in legacy vtk format `.vtk`.
    !! file must be opend and closed in this process.
    !! Existing file is replaced.
    !! this file is used only for visualization by Paraview. 
    character(*),intent(in) :: tag
    integer,intent(in) :: timestep
    logical,intent(in) :: ascii

    integer unit_
    character(128) filename
    character(:),allocatable :: form_

    integer,allocatable :: offsets_(:)
    integer,allocatable :: conns_(:)

    filename = get_filename(tag, timestep, ".vtk")

    if ( ascii ) then
        form_ = "formatted"
    else
        form_ = "unformatted"
    end if

    open(newunit=unit_, file=filename, status="replace", action="write", &
    form=form_, access="stream", position="asis")

    if ( ascii ) then
        call write_vtk_ascii_(unit_)
    else
        ! call write_vtk_binary_(unit_)
    end if

    close(unit_)

end subroutine

subroutine write_vtk_ascii_(unit)
    integer,intent(in) :: unit

    character(:),allocatable :: points_form_, cells_form_, cell_types_form_, cell_data_form_
    integer,parameter :: ndigit = 8
        ! assume number of particles does not exceed 10^8
    character(ndigit) :: pts_, noff_, nconn_

    integer,allocatable,dimension(:) :: conns, offsets, cell_types, cell_data_scalar
    double precision,allocatable,dimension(:,:) :: xyz

    write(unit, "(A)") "# vtk DataFile Version 5.1"
    write(unit, "(A)") "vtk output"
    write(unit, "(A)") "ASCII"
    write(unit, "(A)") "DATASET UNSTRUCTURED_GRID"

    allocate(conns(mv_pdata%N_part))
    allocate(offsets(mv_pdata%N_part+1))
    allocate(xyz(3,mv_pdata%N_part))
    allocate(cell_types(mv_pdata%N_part), source=1)
    block
        integer i

        offsets(1) = 0
        do i = 1, mv_pdata%N_part
            conns(i) = i-1
            offsets(i+1) = i
            xyz(:,i) = mv_pdata%particles(i)%pos
        end do

    end block
    
    write(pts_, "(i0)") mv_pdata%N_part
    points_form_ = "POINTS "//trim(adjustl(pts_))//" double"

    ! POINTS
    write(unit, "(A)") points_form_
    write(unit, "(9(f12.5,1x))") xyz
    write(unit, *)
    
    ! CELLS
    write(noff_, "(i0)") mv_pdata%N_part+1
    write(nconn_, "(i0)") mv_pdata%N_part
    cells_form_ = "CELLS "//trim(adjustl(noff_))//" "//trim(adjustl(nconn_))
    write(unit, "(A)") cells_form_
    ! --- offsets
    write(unit, "(A)") "OFFSETS vtktypeint64"
    write(unit, "(9(i0,1x))") offsets
    ! --- connectivities
    write(unit, "(A)") "CONNECTIVITY vtktypeint64"
    write(unit, "(9(i0,1x))") conns
    write(unit,*)

    ! CELL_TYPES
    cell_types_form_ = "CELL_TYPES "//trim(adjustl(pts_)) ! ncells = npts
    write(unit, "(A)") cell_types_form_
    write(unit, "(i0)") cell_types
    write(unit,*)

    ! FIELD_DATA
    allocate(cell_data_scalar(mv_pdata%N_part))
    cell_data_form_ = "CELL_DATA "//trim(adjustl(pts_))
    write(unit, "(A)") cell_data_form_

    ! -- particle status
    write(unit, "(A)") "SCALARS status int"
    write(unit, "(A)") "LOOKUP_TABLE default"
    block 
        integer i
        do i = 1, mv_pdata%N_part
            cell_data_scalar(i) = mv_pdata%particles(i)%state
        end do
    end block
    write(unit, "(9(i0,1x))") cell_data_scalar
    write(unit,*)

    ! -- other fields such as radius, ref_cell can be appended. 

end subroutine

end module dump_file_m