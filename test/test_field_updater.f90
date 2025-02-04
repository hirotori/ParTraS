program test_field_updater
    use kind_parameters_m
    use unstructured_mesh_m
    use flow_field_m
    use vtk_importer_m
    use field_updater_m
    implicit none
    integer ncyc
    real(DP) dt_f, dt_p
    
    type(vtk_importer_t) vtk_importer, importer
    type(ugrid_struct_t) vtk_ugrid
    
    call vtk_importer%open_ascii_on()
    call vtk_importer%open_stream_on()
    call vtk_importer%open_file("flow_0.vtk")
    call vtk_importer%read_file(vtk_ugrid, shift_index=.true.)
    call vtk_importer%close()

    call construct_flow_field(vtk_ugrid)

    call importer%open_ascii_on()
    call importer%open_stream_on()

    dt_f = 0.01
    dt_p = 0.01
    call mv_field_updater%construct_field_updater(importer, & 
                                                  dt_f, dt_p, &
                                                  "flow_", 0, ".vtk", .true., .true., 1)

    do ncyc = 1, 9

        call mv_field_updater%update_field(ncyc)

        if ( all(int(mv_flow_field%velocity(3,:)) /= ncyc) ) then
            print*, ncyc
            error stop "error ncyc"
        end if

    end do

end program test_field_updater