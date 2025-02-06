program test
    use unstructured_mesh_m
    use vtk_importer_m
    implicit none
    type(vtk_importer_t) vtk_importer
    type(ugrid_struct_t) vtk_ugrid
    
    call vtk_importer%open_ascii_on()
    call vtk_importer%open_stream_on()
    call vtk_importer%open_file("sample_small2x2x2_old.vtk")
    call vtk_importer%read_file(vtk_ugrid, shift_index=.false.)
        
    if ( vtk_ugrid%ncell /= 8 ) then
        error stop "wrong number in ncell"
    end if
    
    if ( vtk_ugrid%nvert /= 27 ) then
        error stop "wrong number in nvert"
    end if

    if ( any(vtk_ugrid%offsets /= [0,8,16,24,32,40,48,56,64]) ) then
        error stop "wrong offsets"
    end if

    if ( any(vtk_ugrid%conns /= [ 0, 1, 4, 3, 9, 10, 13, 12, 1, &
                                  2, 5, 4, 10, 11, 14, 13, 3, 4, & 
                                  7, 6, 12, 13, 16, 15, 4,5, 8, & 
                                  7, 13, 14, 17, 16, 9 ,10 ,13 ,12 , &
                                 18, 19, 22, 21, 10, 11, 14, 13, 19, &
                                 20, 23, 22, 12, 13, 16, 15, 21, 22, &
                                 25, 24, 13, 14, 17, 16, 22, 23, 26, &
                                 25 ]) ) then
        error stop "wrong connectivity"        
    end if

    call delete_ugrid(vtk_ugrid)

end program test