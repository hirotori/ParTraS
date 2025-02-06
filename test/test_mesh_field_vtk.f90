program test
    use vtk_importer_m
    use unstructured_mesh_m
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
    vtk_ugrid%cell_velocity = v
    call construct_flow_field(vtk_ugrid)

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