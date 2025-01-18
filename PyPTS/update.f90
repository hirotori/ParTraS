module update_m
    use iso_c_binding
    use util_m
    use base_importer_m
    use vtk_importer_m
    use field_updater_m
    implicit none
    
    public :: update_field_vtk
contains

subroutine update_field_vtk(dt_f, dt_p, basename, pad, field_only, interval, ascii) bind(c, name="update_field_vtk")
    real(c_double),intent(in) :: dt_f
    real(c_double),intent(in) :: dt_p
    character(1, kind=c_char),dimension(*),intent(in) :: basename
    integer(c_int),intent(in) :: pad
    logical(c_bool),intent(in) :: field_only
    integer(c_int),intent(in) :: interval
    logical(c_bool),intent(in) :: ascii

    type(vtk_importer_t) vtk_importer

    call vtk_importer%open_stream_on()

    if ( ascii ) then
        call vtk_importer%open_ascii_on()
    else
        call vtk_importer%open_ascii_off()
    end if

    call mv_field_updater%construct_field_updater(vtk_importer, CELL_TYPE_VTK, FACE_VERT_DEF_VTK, &
    dt_f, dt_p, fstring(basename), pad, ".vtk", logical(field_only), interval)

end subroutine

subroutine no_update_field() bind(c, name="no_update_field")

    call mv_field_updater%disable_updater()

end subroutine

end module update_m