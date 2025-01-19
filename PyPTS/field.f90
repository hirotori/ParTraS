module init_m
    use iso_c_binding
    use util_m
    use kind_parameters_m
    use base_importer_m
    use vtk_importer_m
    use flow_field_m
    use particle_data_m
    implicit none
    
contains

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

    call construct_flow_field(ugrid_, CELL_TYPE_VTK, FACE_VERT_DEF_VTK)

    call delete_ugrid(ugrid_)

end subroutine

! subroutine init_field_scflow(filename)
!     character(*),intent(in) :: filename


! end subroutine


end module