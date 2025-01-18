program test
    use base_importer_m, only: ugrid_struct_t, delete_ugrid
    use vtk_importer_m
    use flow_field_m
    implicit none
    integer,parameter :: nsample = 5
    integer :: locations(nsample) = [1, 12, 4, 7, 5]
    integer :: solutions(nsample)
    real(8) pts(3,nsample), rn_(3)
    integer i, cellId
    real(8),allocatable :: v(:,:)
    type(vtk_importer_t) vtk_importer
    type(ugrid_struct_t) vtk_ugrid

    block
        integer nseed
        integer,allocatable :: seeds(:)
        call random_seed(size=nseed)
        allocate(seeds(nseed), source=100)
        call random_seed(put=seeds)
    end block

    call vtk_importer%open_ascii_on()
    call vtk_importer%open_stream_on()
    call vtk_importer%open_file("sample_small.vtk")
    call vtk_importer%read_file(vtk_ugrid, shift_index=.true.)
    call vtk_importer%close()

    allocate(v(3,vtk_ugrid%ncell), source=0.d0)
    vtk_ugrid%cell_velocity = v
    call construct_flow_field(vtk_ugrid, CELL_TYPE_VTK, FACE_VERT_DEF_VTK)
    call delete_ugrid(vtk_ugrid)

    do i = 1, nsample
        call random_number(rn_)
        pts(:,i) = mv_flow_field%cell_centers(:,locations(i)) !+ rn_*0.1d0
    end do

    do i = 1, nsample
        print*, "test for", i
        cellId = 1
        call search_reference_cell(pts(:,i), pts(:,i), cellId, 10)
        solutions(i) = cellId
    end do
    
    print "(*(i0,1x))", locations
    print "(*(i0,1x))", solutions

    do i = 1, nsample
        if ( locations(i) /= solutions(i) ) then
            error stop
        end if
    end do
end program test