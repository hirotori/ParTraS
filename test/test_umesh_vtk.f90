program test
    use base_importer_m, only: ugrid_struct_t
    use vtk_importer_m
    use unstructured_mesh_m
    implicit none
    integer i, j, jf, sign_
    integer,allocatable :: face2cells(:,:), face2verts(:,:), boundary_faces(:)
    integer,allocatable :: cell_faces(:), cell_offsets(:)
    type(vtk_importer_t) vtk_importer
    type(ugrid_struct_t) vtk_ugrid

    call vtk_importer%open_ascii_on()
    call vtk_importer%open_stream_on()
    call vtk_importer%open_file("sample_small2x2x2_old.vtk")
    call vtk_importer%read_file(vtk_ugrid, shift_index=.true.)
    call vtk_importer%close()
    
    call construct_half_faces(vtk_ugrid%ncell, vtk_ugrid%nvert, vtk_ugrid%conns, vtk_ugrid%offsets, vtk_ugrid%cell_types, &
                              CELL_TYPE_VTK, FACE_VERT_DEF_VTK)
    print*, ""
    do i = 1, size(half_faces)
        print "('face = ', i0)", i
        print "('(owner, pair) = (', i0,1x,i0, ')')", half_faces(i)%owner, half_faces(i)%pair
        print "('verts = ', *(i0,:,',',1x), ']')", half_faces(i)%vertices(1:half_faces(i)%vert_count)
    end do

    call create_faces(face2cells,face2verts)

    do i = 1, size(face2cells,dim=2)
        print*, face2cells(:,i)
        print*, face2verts(:,i)
    end do

    call create_boundary_faces(face2cells, boundary_faces)

    print "(*(i0,1x))", boundary_faces(:)

    call create_cell_faces(vtk_ugrid%ncell, face2cells, cell_faces, cell_offsets)

    print "(*(i0,1x))", cell_offsets
    print "(*(i0,1x))", cell_faces

    do i = 1, vtk_ugrid%ncell
        print "('cell = ', i0)", i
    do j = cell_offsets(i), cell_offsets(i+1)-1
        jf = cell_faces(j)
        sign_ = (3 + sign(1, jf))/2
        print "('nbcell = ', i0)", face2cells(sign_, abs(jf))
    end do
    end do

end program test