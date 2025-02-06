module init_m
    use iso_c_binding
    use util_m
    use kind_parameters_m
    use vtk_importer_m
    use afdet_importer_m
    use flow_field_m
    implicit none
    
contains

subroutine init_field(ugrid)
    !! initialize a flow field and delete ugrid.
    type(ugrid_struct_t),intent(inout) :: ugrid

    call construct_flow_field(ugrid)
    call delete_ugrid(ugrid)

end subroutine

subroutine init_field_vtk(c_filename, ascii) bind(c, name="init_field_vtk")
    character(c_char),dimension(*),intent(in) :: c_filename
    logical(c_bool),intent(in) :: ascii

    character(:),allocatable :: filename
    type(vtk_importer_t) vtk_importer_
    type(ugrid_struct_t) ugrid_

    filename = fstring(c_filename)

    if ( ascii ) then
        call vtk_importer_%open_ascii_on() 
    else
        call vtk_importer_%open_ascii_off()
    endif

    call vtk_importer_%open_stream_on()

    call vtk_importer_%open_file(filename)
    call vtk_importer_%read_file(ugrid_, .true.)
    call vtk_importer_%close()

    call init_field(ugrid_)

end subroutine

! subroutine init_field_scflow(filename)
!     character(*),intent(in) :: filename


! end subroutine


subroutine init_field_afdet(c_filename) bind(c, name="init_field_afdet")
    character(1, kind=c_char),dimension(*),intent(in) :: c_filename

    character(:),allocatable :: filename
    type(afdet_importer_t) importer_
    type(ugrid_struct_t) ugrid_

    filename = fstring(c_filename)
    
    call importer_%open_stream_on()
    call importer_%open_ascii_off()
    call importer_%open_file(filename)
    call importer_%read_file(ugrid_, .false.)
    call importer_%close()
    
    call init_field(ugrid_)

end subroutine

end module