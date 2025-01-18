program test
    use base_importer_m, only: ugrid_struct_t, delete_ugrid
    use vtk_importer_m
    use flow_field_m
    implicit none
    integer i
    real(8),allocatable :: v(:,:)
    type(vtk_importer_t) vtk_importer
    type(ugrid_struct_t) vtk_ugrid
    
    call vtk_importer%open_ascii_on()
    call vtk_importer%open_stream_on()
    call vtk_importer%open_file("sample_small2x2x2.vtk")
    call vtk_importer%read_file(vtk_ugrid, shift_index=.true.)
    call vtk_importer%close()

    allocate(v(3,vtk_ugrid%ncell), source=0.d0)
    call construct_flow_field(vtk_ugrid%ncell, vtk_ugrid%nvert, vtk_ugrid%conns, &
                              vtk_ugrid%offsets, vtk_ugrid%cell_types, vtk_ugrid%verts, v, &
                              CELL_TYPE_VTK, FACE_VERT_DEF_VTK)

    call delete_ugrid(vtk_ugrid)

    do i = 1, mv_flow_field%nface
        print*, mv_flow_field%face2verts(:,i)
        print*, mv_flow_field%face_centers(:,i)
    end do
    print*, "----"
    do i = 1, mv_flow_field%nface        
        print*, mv_flow_field%face_normals(:,i)
    end do
    print*, "----"
    do i = 1, mv_flow_field%ncell
        print*, mv_flow_field%cell_centers(:,i)
    end do
    print*, "----"
    do i = 1, mv_flow_field%nface
        print*, i, mv_flow_field%face2cells(:,i)
    end do
end program test