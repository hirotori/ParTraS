program test
    use unstructured_mesh_m
    use flow_field_m
    implicit none
    integer,parameter :: nvert = 12
    integer,parameter :: ncell = 2
    integer,parameter :: connectivity(16)  = [0,1,2,3,4,5,6,7, 1,8,9,2,5,10,11,6] + 1
    integer,parameter :: offset(3)         = [0, 8, 16] + 1
    integer,parameter :: cell_types(2)     = [12, 12]
    real(8),parameter :: verts(3,nvert) = reshape([[-1.0, 0.0, 0.0], & !0
                                                   [ 0.0, 0.0, 0.0], & !1
                                                   [ 0.0, 1.0, 0.0], & !2
                                                   [-1.0, 1.0, 0.0], & !3
                                                   [-1.0, 0.0, 1.0], & !4
                                                   [ 0.0, 0.0, 1.0], & !5
                                                   [ 0.0, 1.0, 1.0], & !6
                                                   [-1.0, 1.0, 1.0], & !7
                                                   [ 1.0, 0.0, 0.0], & !8
                                                   [ 1.0, 1.0, 0.0], & !9
                                                   [ 1.0, 0.0, 1.0], & !10
                                                   [ 1.0, 1.0, 1.0]], shape=[3,nvert])
    real(8),parameter :: vv(3,ncell) = 0.d0
    integer i

    ! correct data
    integer,dimension(4,11) :: face2verts = reshape([[2,1,4,3, & !1
                                                         1,5,8,4, & !2
                                                         5,6,7,8, & !3
                                                         2,3,7,6, & !4 (= 8)
                                                         2,6,5,1, & !5
                                                         3,4,8,7, & !6
                                                         9,2,3,10, & !7
                                                         6,11,12,7, & !9
                                                         9,10,12,11, & !10
                                                         9,11,6,2, & !11, 12
                                                         10,3,7,12]], shape=[4,11])

    real(8),dimension(3,11) :: face_centers = reshape([[-0.5, 0.5, 0.0], &
                                                       [-1.0, 0.5, 0.5], &
                                                       [-0.5, 0.5, 1.0], &
                                                       [ 0.0, 0.5, 0.5], &
                                                       [-0.5, 0.0, 0.5], &
                                                       [-0.5, 1.0, 0.5], &
                                                       [ 0.5, 0.5, 0.0], &
                                                       [ 0.5, 0.5, 1.0], &
                                                       [ 1.0, 0.5, 0.5], &
                                                       [ 0.5, 0.0, 0.5], &
                                                       [ 0.5, 1.0, 0.5]], shape=[3,11])

    real(8),dimension(3,11) :: face_normals = reshape([[ 0.0, 0.0,-1.0], &
                                                       [-1.0, 0.0, 0.0], &
                                                       [ 0.0, 0.0, 1.0], &
                                                       [ 1.0, 0.0, 0.0], &
                                                       [ 0.0,-1.0, 0.0], &
                                                       [ 0.0, 1.0, 0.0], &
                                                       [ 0.0, 0.0,-1.0], &
                                                       [ 0.0, 0.0, 1.0], &
                                                       [ 1.0, 0.0, 0.0], &
                                                       [ 0.0,-1.0, 0.0], &
                                                       [ 0.0, 1.0, 0.0]], shape=[3,11])

    real(8),dimension(3,2) :: cell_centers = reshape([[-0.5, 0.5, 0.5],[0.5, 0.5, 0.5]], shape=[3,2])
    type(ugrid_struct_t) ugrid

    ugrid%ncell = ncell
    ugrid%nvert = nvert
    ugrid%conns = connectivity
    ugrid%offsets = offset
    ugrid%cell_types = cell_types
    ugrid%verts = verts
    ugrid%cell_velocity = vv
        
    call construct_flow_field(ugrid)
    
    if ( mv_flow_field%ncell /= 2 ) then
        error stop "ncell"
    end if

    if ( mv_flow_field%nface /= 11 ) then
        error stop "nface"
    end if

    if ( mv_flow_field%nvert /= 12 ) then
        error stop "nvert"
    end if

    do i = 1, mv_flow_field%nface
        if (any(mv_flow_field%face2verts(1:4,i) /= face2verts(:,i))) error stop "face2verts"
        if (any(abs(mv_flow_field%face_centers(:,i) - face_centers(:,i)) > epsilon(1.d0))) error stop "face_ceters"
    end do

    do i = 1, mv_flow_field%nface        
        if (any(abs(mv_flow_field%face_normals(:,i) - face_normals(:,i)) > epsilon(1.d0))) error stop "face_normals"
    end do

    do i = 1, mv_flow_field%ncell
       if ( any(abs(mv_flow_field%cell_centers(:,i) - cell_centers(:,i)) > epsilon(1.d0)) ) then
            error stop "cell_centers"
       end if
    end do

end program test