program test
    use vtk_importer_m, only : CELL_TYPE_VTK, FACE_VERT_DEF_VTK
    use unstructured_mesh_m
    implicit none
    integer,parameter :: nvert = 12
    integer,parameter :: ncell = 2
    integer,parameter :: connectivity(16)  = [0,1,2,3,4,5,6,7, 1,8,9,2,5,10,11,6] + 1
    integer,parameter :: offset(3)         = [0, 8, 16] + 1
    integer,parameter :: cell_types(2)     = [12, 12]
    integer,allocatable :: face2cells(:,:), face2verts(:,:), boundary_faces(:), cell_offsets(:), cell_faces(:)
    integer i
    call construct_half_faces(ncell, nvert, connectivity, offset, cell_types, CELL_TYPE_VTK, FACE_VERT_DEF_VTK)
    print*, ""
    do i = 1, size(half_faces)
        print "('(owner, pair) = (', i0,1x,i0, ')')", half_faces(i)%owner-1, half_faces(i)%pair
        print "('verts = ', *(i0,:,',',1x), ']')", half_faces(i)%vertices(1:half_faces(i)%vert_count)-1
    end do

    call create_faces(face2cells,face2verts)

    do i = 1, size(face2cells,dim=2)
        print*, face2cells(:,i)
        print*, face2verts(:,i)
    end do

    call create_boundary_faces(face2cells, boundary_faces)

    print "(*(i0,1x))", boundary_faces(:)

    call create_cell_faces(ncell, face2cells, cell_faces, cell_offsets)

    print "(*(i0,1x))", cell_offsets
    print "(*(i0,1x))", cell_faces

end program test