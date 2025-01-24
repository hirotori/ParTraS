module flow_field_m
    !! treat 3d unstructured flow field.
    use kind_parameters_m
    use base_importer_m, only : ugrid_struct_t, cell_type_t, face_vertex_def_t !, FIELD_NAME_LEN
    use unstructured_mesh_m
    use geometry_m
    implicit none
    private
    integer,parameter :: REFCELL_OUTBOUNDS = 0
    integer,parameter :: REFCELL_NOTFOUND  = -1

    type, public :: flow_field_t
        ! character(FIELD_NAME_LEN) :: filename
        integer(IP) :: ncell
        integer(IP) :: nface
        integer(IP) :: nvert

        ! geometry
        real(DP),allocatable :: verts(:,:)
        real(DP),allocatable :: face_centers(:,:)
        real(DP),allocatable :: face_normals(:,:)
        integer(IP),allocatable :: face2cells(:,:)
        integer(IP),allocatable :: face2verts(:,:)

        real(DP),allocatable :: cell_centers(:,:)
        integer(IP),allocatable :: cell2faces(:)
        integer(IP),allocatable :: cell_offsets(:)
        integer(IP),allocatable :: cell2verts(:,:)

        integer(IP),allocatable :: boundary_faces(:)

        ! field
        real(DP),allocatable :: velocity(:,:)

        logical :: is_assigned = .false.
    end type 

    type(flow_field_t),protected :: mv_flow_field
        !! flow field object (module variable)

    public mv_flow_field, &
           construct_flow_field, &
           update_flow_field, &
           search_reference_cell, &
           delete_mesh_field, &
           scale_field_variables, &
           REFCELL_NOTFOUND, REFCELL_OUTBOUNDS

contains

subroutine construct_flow_field(ugrid, cell_type_def, face_vert_def)
    type(ugrid_struct_t),intent(in) :: ugrid
        !! unstructured grid importer
    type(cell_type_t),intent(in) :: cell_type_def
        !! cell type definition
    type(face_vertex_def_t),intent(in) :: face_vert_def
        !! face-vertex connectivity definition

    ! これを呼び出す際はフィールド変数はすべて未割り付けと仮定. 
    ! すでに割り付けしているのにこれを呼び出した場合: 一旦delete処理が入るので普通に動く. 
    call update_flow_field(ugrid, cell_type_def, face_vert_def, .false., .false.)

end subroutine

subroutine update_members_(ugrid, field_only, verts_only)
    !! update members of flow_field_t. 
    type(ugrid_struct_t),intent(in) :: ugrid
        !! grid object
    logical,intent(in) :: field_only
        !! whether it updates field vairable (cell velocity) only or not
    logical,intent(in) :: verts_only
        !! whether it updates vertices only or not. Only valid if field_only == .false.

    ! 一旦すべて解放して再度割り付けるのは計算コストが高いと予想される. 
    ! 隣接関係は変わらないと仮定すれば, 流体の速度(cell_velocity)だけ変更すれば良いのでdeleteは不要. 
    if ( .not. field_only ) then 

        if ( verts_only ) then
            ! 格子点だけ更新する場合. 
            if ( size(ugrid%verts, dim=2) /= size(mv_flow_field%verts, dim=2) ) then
                error stop "flow_field_m/update_members_::ERROR::num of verts must be equal when updating verts only."
            end if
            mv_flow_field%verts(:,:) = ugrid%verts(:,:)
            return
        else
            ! 接続関係含めすべて. 

            ! 要素数が変わる可能性があるので, すべて解放する. 
            call delete_mesh_field()
            
            !geometry
            allocate(mv_flow_field%verts, source=ugrid%verts)

            !COMMENT: if分岐増えたせいで, とっちらかってきた. でもこれ以上は複雑にならないはず...
            ! 以下, ugridのメンバをコピーする形で構築. 割り付けられていない場合はこちらで構築. 

            if ( allocated(ugrid%cell2verts) ) then
                allocate(mv_flow_field%cell2verts, source=ugrid%cell2verts)
            else
                call create_cell_verts(ugrid%ncell, ugrid%conns, ugrid%offsets, mv_flow_field%cell2verts)
            end if

            ! 面-セル接続関係は両方が割り付けられている場合のみugridのデータを利用する. 
            if ( allocated(ugrid%face2cells) .and. allocated(ugrid%face2verts) ) then
                allocate(mv_flow_field%face2cells, source=ugrid%face2cells)
                allocate(mv_flow_field%face2verts, source=ugrid%face2verts)
            else
                call create_faces(mv_flow_field%face2cells, mv_flow_field%face2verts)
            end if

            call create_boundary_faces(mv_flow_field%face2cells, mv_flow_field%boundary_faces)
            call create_cell_faces(ugrid%ncell, mv_flow_field%face2cells, mv_flow_field%cell2faces, mv_flow_field%cell_offsets)
            
            mv_flow_field%ncell = ugrid%ncell
            mv_flow_field%nvert = ugrid%nvert
            mv_flow_field%nface = size(mv_flow_field%face2cells, dim=2)
            
            if ( allocated(ugrid%cell_centers) ) then
                allocate(mv_flow_field%cell_centers, source=ugrid%cell_centers)
            else
                allocate(mv_flow_field%cell_centers(3,mv_flow_field%ncell))
                call compute_cell_centers(mv_flow_field%cell2verts, mv_flow_field%verts, mv_flow_field%cell_centers)
            end if

            if ( allocated(ugrid%face_centers) .and. allocated(ugrid%face_normals) ) then
                allocate(mv_flow_field%face_centers, source=ugrid%face_centers)
                allocate(mv_flow_field%face_normals, source=ugrid%face_normals)                
            else
                allocate(mv_flow_field%face_centers(3,mv_flow_field%nface))
                allocate(mv_flow_field%face_normals(3,mv_flow_field%nface))
                call compute_face_centers_and_normals(mv_flow_field%face2verts, mv_flow_field%verts, mv_flow_field%face_centers, mv_flow_field%face_normals)
            end if
            
        end if

        !field variable
        allocate(mv_flow_field%velocity, source=ugrid%cell_velocity)

    else 
        ! すでに割り付け済みと仮定する. 
        mv_flow_field%velocity(:,:) = ugrid%cell_velocity(:,:)
    endif

end subroutine

subroutine update_flow_field(ugrid, cell_type_def, face_vert_def, field_only, verts_only)
    !! update flow field. 
    type(ugrid_struct_t),intent(in) :: ugrid
        !! unstructured grid importer
    type(cell_type_t),intent(in) :: cell_type_def
        !! cell type definition
    type(face_vertex_def_t),intent(in) :: face_vert_def
        !! face-vertex connectivity definition
    logical,intent(in) :: field_only
    logical,intent(in) :: verts_only


    !NOTE: メッシュがポリへドラルの場合, そもそもhalf_facesが構築できない. その場合half_faceの構築はskipするしかない. 
    !      したがって, ポリへドラルメッシュの場合面ーセル関係がugridで全て構築済みでなければならない. 
    if ( .not. field_only .and. .not. verts_only) call construct_half_faces(ugrid, cell_type_def, face_vert_def)

    call update_members_(ugrid, field_only, verts_only)

    call validate_flow_field()

    ! delete current data of half faces
    call delete_half_faces()
    
    mv_flow_field%is_assigned = .true.

    call print_field()

end subroutine

subroutine validate_flow_field()
    !! test if all members are allocated and have valid shapes
    integer ncellface

    if ( mv_flow_field%nvert <= 0 ) error stop "flow_field/validate::ERROR:: nvert <= 0" 
    if ( mv_flow_field%ncell <= 0 ) error stop "flow_field/validate::ERROR:: ncell <= 0" 
    if ( mv_flow_field%nface <= 0 ) error stop "flow_field/validate::ERROR:: nface <= 0" 

    if ( .not. validate_real_array(mv_flow_field%verts, [3, mv_flow_field%nvert])) then
        error stop "flow_field/validate::ERROR:: verts"
    endif

    if ( .not. validate_real_array(mv_flow_field%face_centers, [3, mv_flow_field%nface])) then
        error stop "flow_field/validate::ERROR:: face_centers"
    endif

    if ( .not. validate_real_array(mv_flow_field%face_normals, [3, mv_flow_field%nface])) then
        error stop "flow_field/validate::ERROR:: face_normals"
    endif

    if ( .not. validate_int_array(mv_flow_field%face2cells, [2, mv_flow_field%nface]) ) then
        error stop "flow_field/validate::ERROR:: face2cells"        
    end if

    if ( .not. validate_int_array(mv_flow_field%face2verts, [FACE_VERT_SIZE, mv_flow_field%nface]) ) then
        error stop "flow_field/validate::ERROR:: face2verts"        
    end if

    if ( .not. validate_real_array(mv_flow_field%cell_centers, [3, mv_flow_field%ncell])) then
        error stop "flow_field/validate::ERROR:: cell_centers"
    endif

    if ( .not. validate_int_array(mv_flow_field%cell2verts, [CELL_VERT_SIZE, mv_flow_field%ncell]) ) then
        error stop "flow_field/validate::ERROR:: cell2verts"        
    end if

    ncellface = mv_flow_field%cell_offsets(mv_flow_field%ncell+1)-1
    if ( .not. validate_int_array(mv_flow_field%cell2faces, [ncellface]) ) then
        error stop "flow_field/validate::ERROR:: cell2faces"                
    end if

    if ( .not. validate_int_array(mv_flow_field%cell_offsets, [mv_flow_field%ncell+1]) ) then
        error stop "flow_field/validate::ERROR:: cell_offsets"            
    end if

end subroutine

logical function validate_real_array(arr, arr_shape) result(valid)
    real(DP),allocatable,intent(in) :: arr(..)
    integer(IP),intent(in) :: arr_shape(rank(arr))

    valid = .true.
    if ( .not. allocated(arr) ) then
        valid = .false.
    else
        if ( all(shape(arr) /= arr_shape) ) then
            valid = .false.
        end if
    end if

end function

logical function validate_int_array(arr, arr_shape) result(valid)
    integer(IP),allocatable,intent(in) :: arr(..)
    integer(IP),intent(in) :: arr_shape(rank(arr))

    valid = .true.
    if ( .not. allocated(arr) ) then
        valid = .false.
    else
        if ( all(shape(arr) /= arr_shape) ) then
            valid = .false.
        end if
    end if

end function

subroutine print_field()

    print "('flow_field/print_field:: Status of current flow field: ')"
    print "(' number of cell = ',i0)", mv_flow_field%ncell
    print "(' number of vert = ',i0)", mv_flow_field%nvert
    print "(' number of face = ',i0)", size(mv_flow_field%face2cells, dim=2)

end subroutine

pure subroutine search_reference_cell(r0, ri, ref_cell, nrepeat)
    !! locating cell which contains a given particle i.
    real(DP),intent(in)    :: r0(3)
        !! particle coordinates (x,y,z) at the current step. 
    real(DP),intent(inout) :: ri(3)
        !! particle coordinates (x,y,z) at the next step. 
        !! if this exists out of boundary, then the position is moved back onto wall.  
    integer(IP),intent(in) :: nrepeat
        !! maximum number of iteration
    integer(IP),intent(inout) :: ref_cell
        !! reference cell. next reference cell is stored. 
        !! if new_ref_cell == -1, there are no cells which contain the particle

    integer np, nf, jf, jfa, jf_next, sign_
    real(DP) n_(3), rf_(3), proj_ij, proj_ijmax
    real(DP),parameter :: proj_ij_init = -100.d0
    real(DP),parameter :: epsilon_ = epsilon(1.d0)

    proj_ijmax = proj_ij_init
    do np = 1, nrepeat
        do nf = mv_flow_field%cell_offsets(ref_cell), mv_flow_field%cell_offsets(ref_cell+1)-1 
            jf = mv_flow_field%cell2faces(nf)
            jfa = abs(jf)
            n_(:) = mv_flow_field%face_normals(:,jfa)*sign(1, jf)
            rf_(:) = mv_flow_field%face_centers(:,jfa)

            proj_ij = sum(n_*(ri-rf_))
            if ( proj_ij > proj_ijmax ) then
                proj_ijmax = proj_ij
                jf_next = jf
            end if
        end do

        ! any proj_ij < proj_ijmax
        ! proj_ijmax < 0 yields proj_ij < 0 for all face j in cell i
        if ( proj_ijmax < epsilon_ ) then
            return
        else
            ! go to next cell
            sign_ = (3 + sign(1, jf_next))/2
            ref_cell = mv_flow_field%face2cells(sign_,abs(jf_next))
            proj_ijmax = proj_ij_init
        end if

        !
        if ( ref_cell == REFCELL_OUTBOUNDS ) then
            ! move particle back onto wall
            ! 粒子の現在位置r0と新しい位置riの軌道と面の交点を調べる．
            return
        end if

    end do

    !見つからない場合, 距離判定に切り替える. まだ実装してないので, NOTFOUNDを返す. 
    ref_cell = REFCELL_NOTFOUND
end subroutine

subroutine scale_field_variables(L_ref, U_ref)
    !! scale field variables by reference values
    real(DP),intent(in),optional :: L_ref
        !! reference length
    real(DP),intent(in),optional :: U_ref
        !! reference velocity

    if ( present(U_ref) ) then
        mv_flow_field%velocity(:,:) = mv_flow_field%velocity(:,:)/U_ref
    end if

end subroutine

subroutine delete_mesh_field()
    !! delete current data of flowfield

    if (allocated(mv_flow_field%face2cells)) deallocate(mv_flow_field%face2cells)
    if (allocated(mv_flow_field%face2verts)) deallocate(mv_flow_field%face2verts)
    if (allocated(mv_flow_field%face_centers)) deallocate(mv_flow_field%face_centers)
    if (allocated(mv_flow_field%face_normals)) deallocate(mv_flow_field%face_normals)
    if (allocated(mv_flow_field%cell2faces)) deallocate(mv_flow_field%cell2faces)
    if (allocated(mv_flow_field%cell2verts)) deallocate(mv_flow_field%cell2verts)
    if (allocated(mv_flow_field%cell_offsets)) deallocate(mv_flow_field%cell_offsets)
    if (allocated(mv_flow_field%cell_centers)) deallocate(mv_flow_field%cell_centers)
    if (allocated(mv_flow_field%boundary_faces)) deallocate(mv_flow_field%boundary_faces)
    if (allocated(mv_flow_field%velocity)) deallocate(mv_flow_field%velocity)
    if (allocated(mv_flow_field%verts)) deallocate(mv_flow_field%verts)

    mv_flow_field%ncell = 0
    mv_flow_field%nvert = 0
    mv_flow_field%nface = 0

    mv_flow_field%is_assigned = .false.

end subroutine

end module flow_field_m