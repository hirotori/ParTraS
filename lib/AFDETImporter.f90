module afdet_importer_m
    use kind_parameters_m
    use base_importer_m
    use unstructured_mesh_m, only: ugrid_struct_t, VTK_HEXA, VTK_TETRA, VTK_PYRAMID, VTK_WEDGE
    implicit none
    
    type,extends(ugrid_importer_t) :: afdet_importer_t
        !! imports backup data written by afdet solver
        contains
        procedure :: read_file => read_backup_file
    end type

    interface read_arr_
        !! interface for read chunk data
        module procedure read_int_arr_r1
        module procedure read_int_arr_r2    
        module procedure read_real_arr_r2
    end interface

    integer,parameter,private :: AFDET_TETRA = 1
    integer,parameter,private :: AFDET_PYRAMID = 2
    integer,parameter,private :: AFDET_WEDGE = 3
    integer,parameter,private :: AFDET_HEXA = 4
    

    public :: afdet_importer_t

contains

subroutine read_backup_file(this, ugrid, shift_index)
    class(afdet_importer_t),intent(inout) :: this
    type(ugrid_struct_t),intent(out) :: ugrid
    logical,intent(in) :: shift_index
        !! shift 0-start indices to 1-start ones

    integer unit_

    real(DP) time_
    integer(IP) step_
    integer(IP) ver_
    integer(IP) face_
    real(DP) dummy_, dummy_3_(3)

    ! 以下3つは使わず, 読み込みのためだけに用意する. 
    ! 読み飛ばし処理書くよりも, フォーマット通り読むほうが楽なので. 
    integer(IP),allocatable :: boundary_offsets_(:)
    integer(IP),allocatable :: boundary_faces_(:)
    integer(IP),allocatable :: vertex_bounds_(:,:)

    !データはバイナリで書かれていること前提. 

    unit_ = this%get_current_unit()
    read(unit_) time_
    read(unit_) step_
    
    ! grid
    read(unit_) ver_
    if ( ver_ >= 1 ) then
        call read_arr_(unit_, ugrid%cell2verts)
        call read_arr_(unit_, ugrid%cell_types)
        call read_arr_(unit_, ugrid%face2cells)
        call read_arr_(unit_, ugrid%face2verts)   
        call read_arr_(unit_, ugrid%verts)
        call read_arr_(unit_, boundary_offsets_)
        call read_arr_(unit_, boundary_faces_)
    end if

    if (ver_ >= 2) then
        read(unit_) ugrid%nvert
        read(unit_) face_
        read(unit_) ugrid%ncell
        call read_arr_(unit_, vertex_bounds_)
    end if
    
    ! fluid
    read(unit_) ver_ !; print*, "ver = ", ver_
    read(unit_) dummy_
    read(unit_) dummy_
    read(unit_) dummy_
    read(unit_) dummy_
    read(unit_) dummy_
    read(unit_) dummy_
    read(unit_) dummy_3_
    call read_arr_(unit_, ugrid%cell_velocity)

    ! create offsets and connectivities
    ! cell2vertsから直にhal/faceを構築できればいいのだが, そうではないのでここで作成. 
    call create_conns_and_offsets_(ugrid%ncell, ugrid%cell2verts, ugrid%offsets, ugrid%conns)

    ! セルタイプを変換する
    where (ugrid%cell_types == AFDET_TETRA)
        ugrid%cell_types = VTK_TETRA
    elsewhere (ugrid%cell_types == AFDET_PYRAMID)
        ugrid%cell_types = VTK_PYRAMID
    elsewhere (ugrid%cell_types == AFDET_WEDGE)
        ugrid%cell_types = VTK_WEDGE
    elsewhere (ugrid%cell_types == AFDET_HEXA)
        ugrid%cell_types = VTK_HEXA
    end where

    ! ゴーストセルは除外する. 
    where (ugrid%face2cells > ugrid%ncell)
        ugrid%face2cells = 0
    end where
end subroutine

! helper
subroutine read_real_arr_r2(unit, arr)
    !! read data chunk.
    integer,intent(in) :: unit
    real(DP),allocatable :: arr(:,:)

    integer lower_(rank(arr))
    integer upper_(rank(arr))

    read(unit) lower_ !; print*, lower_
    read(unit) upper_ !; print*, upper_

    allocate(arr(lower_(1):upper_(1), lower_(2):upper_(2)))
    read(unit) arr

end subroutine

subroutine read_int_arr_r1(unit, arr)
    !! read data chunk.
    integer,intent(in) :: unit
    integer(IP),allocatable :: arr(:)

    integer lower_(rank(arr))
    integer upper_(rank(arr))

    read(unit) lower_ !; print*, lower_
    read(unit) upper_ !; print*, upper_

    allocate(arr(lower_(1):upper_(1)))
    read(unit) arr

end subroutine

subroutine read_int_arr_r2(unit, arr)
    !! read data chunk.
    integer,intent(in) :: unit
    integer(IP),allocatable :: arr(:,:)

    integer lower_(rank(arr))
    integer upper_(rank(arr))

    read(unit) lower_ !; print*, lower_
    read(unit) upper_ !; print*, upper_

    allocate(arr(lower_(1):upper_(1), lower_(2):upper_(2)))
    read(unit) arr

end subroutine

subroutine create_conns_and_offsets_(ncell, cell2verts, offsets, conns)
    !! cell2vertsからconnectivity と offsets を生成する. 
    integer,intent(in) :: ncell
    integer,intent(in) :: cell2verts(:,:)
    integer,allocatable,intent(inout) :: offsets(:)
    integer,allocatable,intent(inout) :: conns(:)

    integer i, nv

    allocate(offsets(ncell+1), source=0)
    offsets(1) = 1
    do i = 1, ncell
        nv = count(cell2verts(:,i) > 0)
        offsets(i+1) = offsets(i) + nv
    end do

    allocate(conns(offsets(ncell+1)-1))
    do i = 1, ncell
        nv = count(cell2verts(:,i) > 0)
        conns(offsets(i):offsets(i+1)-1) = cell2verts(1:nv,i)
    end do

end subroutine

end module afdet_importer_m