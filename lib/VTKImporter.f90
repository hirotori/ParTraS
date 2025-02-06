module vtk_importer_m
    use kind_parameters_m
    use unstructured_mesh_m
    use base_importer_m
    use id_list_m
    implicit none
    private    
    type,extends(ugrid_importer_t) :: vtk_importer_t

        contains
        procedure :: read_file => import_vtk_ugrid_legacy
    end type

    public vtk_importer_t

contains


subroutine import_vtk_ugrid_legacy(this, ugrid, shift_index)
    class(vtk_importer_t),intent(inout) :: this
    type(ugrid_struct_t),intent(out) :: ugrid
    logical,intent(in) :: shift_index
    !! shift 0-start indices to 1-start ones

    character(12) form_

    inquire(unit=this%get_current_unit(), formatted=form_)

    select case (trim(form_))
    case ("YES")
        call read_ascii_(ugrid, this%get_current_unit())
    case ("NO")
        !call read_binary_
    case default
        error stop "vtk_importer_m/vtk_importer_t::error:: unknown form"
    end select

    if ( shift_index ) then
        ugrid%conns = ugrid%conns + 1
        ugrid%offsets = ugrid%offsets + 1
    end if

    ! ugrid%filename = this%get_filename()
    
end subroutine

subroutine read_ascii_(this, unit)
    type(ugrid_struct_t),intent(inout) :: this
    integer,intent(in) :: unit
    character(128) buff_

    integer istat
    integer :: noff, nconn
    integer file_maj_ver, file_min_ver
    integer,allocatable :: buffer_cells(:)
    type(id_list_t) conns_
    
    ! header
    read(unit,"(A)") buff_
    read(buff_(24:24), *) file_maj_ver
    read(buff_(26:), *) file_min_ver
    read(unit,"(A)")
    read(unit,"(A)")
    read(unit,"(A)") buff_
    if ( .not. find_string_("DATASET UNSTRUCTURED_GRID", trim(buff_))) then
        error stop "vtk_importer_m/read_ascii_::error:: DATASET WRONG"        
    end if

    do  
        read(unit,"(A)", iostat=istat) buff_
        if ( is_iostat_end(istat) ) then
            exit
        end if

        if ( find_string_("POINTS", buff_) ) then
            read(buff_(7:),*) this%nvert
            allocate(this%verts(3,this%nvert))
            read(unit,*) this%verts
        end if

        if ( find_string_("CELLS", buff_) ) then
            if ( file_maj_ver > 4 ) then
                read(buff_(6:),*) noff, nconn
                allocate(this%offsets(noff))
                allocate(this%conns(nconn))
                do  
                    read(unit,"(A)") buff_
                    if ( find_string_("OFFSETS", buff_) ) then                        
                        ! offsets
                        read(unit,*) this%offsets
                    end if
                    if ( find_string_("CONNECTIVITY", buff_) ) then                        
                        ! connectivity
                        read(unit,*) this%conns
                        exit
                    end if
                end do
            else 
                ! under ver. ~ 2
                read(buff_(6:),*) this%ncell, nconn
                allocate(buffer_cells(nconn))
                allocate(this%offsets(this%ncell+1))
                conns_ = id_list_t(n=this%ncell)
                read(unit,*) buffer_cells
                block
                    integer i, ist, n
                    ist = 1
                    this%offsets(1) = 0
                    do i = 2, this%ncell+1
                        n = buffer_cells(ist)
                        this%offsets(i) = this%offsets(i-1) + n
                        call conns_%push_back(buffer_cells(ist+1:ist+n))
                        ist = ist + n + 1
                    end do
                end block
                allocate(this%conns, source=conns_%get_data())
            end if
        end if

        if ( find_string_("CELL_TYPES", buff_) ) then
            read(buff_(11:),*) this%ncell
            allocate(this%cell_types(this%ncell))
            read(unit,*) this%cell_types
        end if

        if ( find_string_("CELL_DATA", buff_) ) then
            do
                read(unit, "(A)") buff_
                if ( find_string_("VECTORS", buff_) ) then
                    allocate(this%cell_velocity(3,this%ncell))
                    read(unit,*) this%cell_velocity
                    exit                    
                end if
            end do
        end if
    end do
end subroutine

logical function find_string_(string, line) result (found)
    character(*),intent(in) :: string, line
    found = index(line, string) > 0
end function

! subroutine read_binary_(unit)
!     integer,intent(in) :: unit

! end subroutine

end module vtk_importer_m