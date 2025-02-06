program test_afdet_impoter
    use unstructured_mesh_m
    use afdet_importer_m
    use flow_field_m
    implicit none
    
    type(afdet_importer_t) importer
    type(ugrid_struct_t) grid

    call importer%open_stream_on()
    call importer%open_ascii_off()
    call importer%open_file("backup_data.bin")
    call importer%read_file(grid, .false.)
    call importer%close()

    if ( grid%ncell /= 32*16*16 ) then
        error stop "wrong cell num"
    end if

    !読み込めるかのテスト
    call construct_flow_field(grid)

end program 