module unstructured_mesh_m
    use id_list_m
    use base_importer_m
    implicit none
    private
    integer,parameter :: NULL_ID = 0
    integer,parameter :: CELL_OUTBOUND = 0
    integer,parameter :: CELL_VERT_SIZE = 10
    integer,parameter :: FACE_VERT_SIZE = 6
    
    type half_face_t
        integer :: pair  = NULL_ID
        integer :: owner = NULL_ID
        integer :: vertices(FACE_VERT_SIZE) = NULL_ID
        integer vert_count
    end type

    type(half_face_t),allocatable,protected :: half_faces(:)
        !! half-face array (module variable)

    interface construct_half_faces
        module procedure construct_half_faces_importer
        module procedure construct_half_faces_c2v
        module procedure construct_half_faces_conns
    end interface

    public half_faces, &
           construct_half_faces, &
           create_faces, &
           create_boundary_faces, &
           create_cell_faces, &
           create_cell_verts, &
           delete_half_faces

    contains

    subroutine construct_half_faces_importer(ugrid, cell_type_def, face_vert_def)
        !! construct half-face array. result is stored in a module variable `half_faces`.
        !! previous half_face data is deleted when this is called.

        type(ugrid_struct_t),intent(in) :: ugrid
            !! unstructured grid importer
        type(cell_type_t),intent(in) :: cell_type_def
            !! cell type definition
        type(face_vertex_def_t),intent(in) :: face_vert_def
            !! face-vertex connectivity definition

        call construct_half_faces_conns(ugrid%ncell, ugrid%nvert, ugrid%conns, ugrid%offsets, ugrid%cell_types, &
        cell_type_def, face_vert_def)

    end subroutine

    subroutine construct_half_faces_c2v(ncell, nvert, cell2verts, cell_types, cell_type_def, face_vert_def)
        !! construct half-face array. result is stored in a module variable `half_faces`.
        integer,intent(in) :: ncell
            !! number of cell
        integer,intent(in) :: nvert
            !! number of vertices
        integer,intent(in) :: cell2verts(:,:)
            !! cell offset array
        integer,intent(in) :: cell_types(:)
            !! cell type array
        type(cell_type_t),intent(in) :: cell_type_def
        type(face_vertex_def_t),intent(in) :: face_vert_def


        !convert cell2verts to connectivity and offset array
        !@NOTE むしろ connectivityとoffsetsからcell2verts作ったほうが効率的かも. 未実装. 

    end subroutine

    subroutine construct_half_faces_conns(ncell, nvert, connectivities, offsets, cell_types, cell_type_def, face_vert_def)
        !! construct half-face array. result is stored in a module variable `half_faces`.
        !! previous half_face data is deleted when this is called.
        integer,intent(in) :: ncell
            !! number of cell
        integer,intent(in) :: nvert
            !! number of vertices
        integer,intent(in) :: connectivities(:)
            !! cell connectivity array
        integer,intent(in) :: offsets(:)
            !! cell offset array
        integer,intent(in) :: cell_types(:)
            !! cell type array
        type(cell_type_t),intent(in) :: cell_type_def
            !! cell type definition
        type(face_vertex_def_t),intent(in) :: face_vert_def
            !! face-vertex connectivity definition

        integer,allocatable :: face_def_(:,:)
        integer,allocatable :: candidates(:)
        integer,dimension(CELL_VERT_SIZE) :: verts_ = NULL_ID
        integer,dimension(FACE_VERT_SIZE) :: ids_, ids_othres_
        integer ntet, nwed, nhex, npyr, nhface, nface, nvert_face, nvert_face_other, nvert_cell
        integer sumids, other
        integer ic, jf, jhf, kv, st_id

        type(id_list_t),allocatable :: vert2face(:)

        call delete_half_faces()

        ntet = count(cell_types == cell_type_def%cell_type_tetra)
        nwed = count(cell_types == cell_type_def%cell_type_wedge)
        nhex = count(cell_types == cell_type_def%cell_type_hexa)
        npyr = count(cell_types == cell_type_def%cell_type_pyram)
        nhface = ntet*4 + nwed*5 + nhex*6 + npyr*5

        allocate(half_faces(nhface))
        jhf = 1
        do ic = 1, ncell
            st_id = offsets(ic)
            nvert_cell = offsets(ic+1) - st_id
            verts_(:nvert_cell) = connectivities(st_id:st_id+nvert_cell-1)

            if (cell_types(ic) == cell_type_def%cell_type_tetra) then
                face_def_ = face_vert_def%face_def_tetra
                nface = 4
            else if (cell_types(ic) == cell_type_def%cell_type_wedge) then
                face_def_ = face_vert_def%face_def_wedge
                nface = 5
            else if (cell_types(ic) == cell_type_def%cell_type_hexa) then
                face_def_ = face_vert_def%face_def_hexa
                nface = 6
            else if (cell_types(ic) == cell_type_def%cell_type_pyram) then
                face_def_ = face_vert_def%face_def_pyram
                nface = 5
            else                
                error stop "error::construct_half_faces:: unknown cell type"
            end if

            do jf = 1, nface
                nvert_face = face_def_(1,jf)
                half_faces(jhf)%vertices(1:nvert_face) = verts_(face_def_(2:nvert_face+1,jf))
                half_faces(jhf)%owner = ic
                half_faces(jhf)%pair = NULL_ID
                half_faces(jhf)%vert_count = nvert_face
                ! if ( jhf == 75649 .or. jhf == 76761 ) then
                !     print*, verts_
                !     print*, jf
                !     print*, face_def_(2:nvert_face+1,jf)
                !     print*, half_faces(jhf)
                ! end if
                jhf = jhf + 1
            end do
        end do

        if ( jhf - 1 /= nhface ) then
            print "('unstructured_mesh_m/construct_half_faces::error:: jhf = ', i0, ', but nhface = ', i0)", jhf-1, nhface
            error stop "wrong count jhf"
        end if

        !!! search neighbor faces
        allocate(vert2face(nvert), source=id_list_t())
        do jf = 1, nhface
            do kv = 1, half_faces(jf)%vert_count
                call vert2face(half_faces(jf)%vertices(kv))%push_back(jf)
            end do
        end do

        ids_(:) = NULL_ID; ids_othres_(:) = NULL_ID
        do jhf = 1, nhface
            if ( half_faces(jhf)%pair /= NULL_ID ) cycle 

            nvert_face = half_faces(jhf)%vert_count
            ids_(:) = half_faces(jhf)%vertices(:)
            candidates = vert2face(ids_(1))%get_data()
            sumids = sum(ids_)

            ! if ( jhf == 75649 ) then
            !     print*, nvert_face
            !     print*, ids_
            !     print*, sumids
            ! end if
            do jf = 1, size(candidates)
                other = candidates(jf)

                if ( other == jhf ) cycle
                
                ids_othres_(:) = half_faces(other)%vertices(:)
                if (sumids /= sum(ids_othres_)) cycle
                ! if ( jhf == 75649 ) then
                !     print*, "----"
                !     print*, other
                !     print*, ids_othres_
                !     print*, sum(ids_othres_)
                !     print*, half_faces(other)%owner
                ! end if 

                !@NOTE: ids_の順番がids_others_の奇置換になっていないことを念頭に置いた処理.
                !       ふつうids_others_はids_と奇置換になっているはず. 
                !       昔のコードだとなっていないっぽい? (e.g. sax_flow.vtk)
                !       Tkiさんは奇置換前提でかいてる (is_same_face_)
                if ( is_same_set_(ids_, ids_othres_, nvert_face) ) then
                    half_faces(jhf)%pair = other
                    half_faces(other)%pair = jhf
                    exit
                endif                
            end do
        end do

    end subroutine

    pure logical function is_same_face_(idxs1, idxs2, n)
        !! 与えられた2つの節点集合が同一の面を指すかの判定. 
        !! idxs2はidxs1に対して反時計回りに並んでいると仮定されている. 
        integer, intent(in) :: idxs1(1:n), idxs2(1:n), n        
        integer i                
        is_same_face_ = .false.
        
        do i = 1, n
            if (sum(abs(idxs1 - cshift(idxs2(n:1:-1),i-1))) == 0) then  
                is_same_face_ = .true.
                exit
            end if
        end do
    end function

    pure logical function is_same_set_(idxs1, idxs2, n)
        !! 与えられた2つの集合 idxs1, idxs2が同一かを見る. 
        integer, intent(in) :: idxs1(1:n), idxs2(1:n), n        
        integer i1, i2
        integer count_
        is_same_set_ = .false.
        
        ! 考えられる方法: 
        ! 1) 昇順ソートした配列が同じなら同一. --> 要素数高々3~4なので, むしろ並び替えでコストか?
        ! 2) 総当たり. O(nm)
        ! ここでは簡単のため 2) を選択. 
        count_ = 0
        do i1 = 1, n
        do i2 = 1, n
            if ( idxs1(i1) == idxs2(i2) ) then
                count_ = count_ + 1                
            end if            
        end do            
        end do

        is_same_set_ = count_ == n

    end function

    subroutine create_faces(face2cells, face2verts)
        !! create face-to-cell and face-to-vertices connectivities.
        integer,allocatable,intent(out) :: face2cells(:,:)
            !! face-to-cell connectivities. shape=(2,nface)
            !! 1 in 1st dim: owner cell, 2 in 2nd dim: neighbor cell, where face2cell(1,face) < face2cell(2,face)
            !! face(2,face) == NULL_ID means the cell is not inner cell (i.e. this is a boundary face).  
        integer,allocatable,intent(out) :: face2verts(:,:)
            !! face-to-vertex connectivities.

        integer,allocatable :: f2c_(:,:)
        integer,allocatable :: f2v_(:,:)
        integer nhface, nvert
        integer jhf, pair
        integer nface

        nhface = size(half_faces)
        allocate(face2cells(             2,nhface), source=NULL_ID)
        allocate(face2verts(FACE_VERT_SIZE,nhface), source=NULL_ID)

        do jhf = 1, nhface
            face2verts(:,jhf) = half_faces(jhf)%vertices(:)
            face2cells(1,jhf) = half_faces(jhf)%owner
            pair = half_faces(jhf)%pair
            ! if (jhf == 75649) print*, half_faces(jhf)
            if ( pair /= NULL_ID ) then
                face2cells(2,jhf) = half_faces(pair)%owner
            else
                face2cells(2,jhf) = CELL_OUTBOUND
            end if

        end do

        ! use only faces where face2cells(1,face) < face2cells(2,face)
        nface = 0
        do jhf = 1, nhface
            if ( face2cells(1,jhf) < face2cells(2,jhf) .or. face2cells(2,jhf) == CELL_OUTBOUND) then
                nface = nface + 1
                face2cells(:,nface) = face2cells(:,jhf)
                face2verts(:,nface) = face2verts(:,jhf)
            end if            
        end do

        allocate(f2c_, source=face2cells(:,1:nface)); call move_alloc(f2c_, face2cells)
        allocate(f2v_, source=face2verts(:,1:nface)); call move_alloc(f2v_, face2verts)

    end subroutine

    subroutine create_boundary_faces(face2cells, boundary_faces)
        !! create an array pointing boundary face ids.
        integer,intent(in) :: face2cells(:,:)
        integer,allocatable,intent(out) :: boundary_faces(:)

        integer nface, jf
        type(id_list_t) :: boundary_faces_list_
        
        boundary_faces_list_ = id_list_t()
        nface = size(face2cells, dim=2)

        do jf = 1, nface
            if ( face2cells(2,jf) == CELL_OUTBOUND ) then
                call boundary_faces_list_%push_back(jf)
            end if            
        end do

        allocate(boundary_faces(boundary_faces_list_%current_size()), source=boundary_faces_list_%get_data())

    end subroutine

    subroutine create_cell_faces(ncell, face2cell, cell_faces, cell_offsets)
        !! create cell-to-face connectivities.
        ! cell_faces   = [f11,f12,f13,f14,f21,f22,f23,...] where fij denotes the j-th face composing the i-th cell.
        ! cell_offsets = [1, 5, ...]
        ! faces of the i-th cell is access via cell_faces and cell_offsets like below:
        ! do i = 1, ncell
        ! do ij = cell_offsets(i), cell_offsets(i+1)-1
        !  j = cell_faces(ij)
        !  ! if j < 0, face direction is inward the cell.
        !  ...
        ! enddo
        ! enddo 
        integer,intent(in) :: ncell
        integer,intent(in) :: face2cell(:,:)
        integer,allocatable,intent(out) :: cell_faces(:)
        integer,allocatable,intent(out) :: cell_offsets(:)

        integer ic, jf, ic1, ic2
        integer,allocatable :: offsets(:), counts(:)

        allocate(offsets(ncell+1))
        allocate(counts(0:ncell), source=0) !counts(0): ghost cells

        ! count number of faces for each cell
        do jf = 1, size(face2cell, dim=2)
            ic1 = face2cell(1,jf)
            ic2 = face2cell(2,jf)
            counts(ic1) = counts(ic1) + 1       
            counts(ic2) = counts(ic2) + 1       
        end do

        ! counts -> offsets
        offsets(1) = 1
        do ic = 2, ncell+1
            offsets(ic) = offsets(ic-1) + counts(ic-1)
        end do

        ! create cell_faces
        allocate(cell_faces(offsets(ncell+1)-1))
        counts(:) = 0
        do jf = 1, size(face2cell, dim=2)
            ic1 = face2cell(1,jf)
            ic2 = face2cell(2,jf)
            call sorting_ascending(offsets(ic1), counts(ic1), jf)
            counts(ic1) = counts(ic1) + 1
            if ( ic2 /= CELL_OUTBOUND ) then
                call sorting_ascending(offsets(ic2), counts(ic2), -jf)
                counts(ic2) = counts(ic2) + 1
            end if
        end do

        call move_alloc(offsets, cell_offsets)

        contains

        subroutine sorting_ascending(offsets_i, counts_i, jf_i)
            integer,intent(in) :: offsets_i, counts_i, jf_i
            
            integer i

            cell_faces(offsets_i+counts_i) = jf_i

            do i = offsets_i + counts_i, offsets_i + 1, -1
                if ( abs(cell_faces(i-1)) > abs(jf_i)) then
                    cell_faces(i) = cell_faces(i-1)
                    cell_faces(i-1) = jf_i
                else
                    exit
                end if
            end do

        end subroutine
    end subroutine

    subroutine create_cell_verts(ncell, connectivities, offsets, cell2verts)
        integer,intent(in) :: ncell
        integer,intent(in) :: connectivities(:)
        integer,intent(in) :: offsets(:)
        integer,allocatable,intent(out) :: cell2verts(:,:)

        integer ic, st_id, nvert_i

        allocate(cell2verts(CELL_VERT_SIZE,ncell), source=NULL_ID)

        do ic = 1, ncell
            st_id = offsets(ic)
            nvert_i = offsets(ic+1) - st_id
            cell2verts(1:nvert_i,ic) = connectivities(st_id:st_id+nvert_i-1)
        end do

    end subroutine


    subroutine assign_id_to_boundary(boundary_faces, face_centers, face_normals, proc, nzone, zone2boundary)
        interface
            integer function proc(face_centers, face_normals) result(id)
                real(8),intent(in) :: face_centers(:,:)
                real(8),intent(in) :: face_normals(:,:)

            end function
        end interface
        integer,intent(in) :: boundary_faces(:)
        real(8),intent(in) :: face_centers(:,:)
        real(8),intent(in) :: face_normals(:,:)
        integer,intent(in) :: nzone
        integer,intent(out) :: zone2boundary(:)


    end subroutine


    subroutine delete_half_faces()

        if ( allocated(half_faces) ) then
            deallocate(half_faces)
        end if

    end subroutine

end module unstructured_mesh_m