module init_m
    use iso_c_binding
    use util_m
    use kind_parameters_m
    use base_importer_m
    use vtk_importer_m
    use flow_field_m
    use particle_data_m
    implicit none
    
contains

subroutine init_field_vtk(c_filename, ascii) bind(c, name="init_field_vtk")
    character(c_char),dimension(*),intent(in) :: c_filename
    logical(c_bool),intent(in) :: ascii

    character(:),allocatable :: filename
    type(vtk_importer_t) vtk_importer_
    type(ugrid_struct_t) ugrid_

    filename = fstring(c_filename)

    if ( ascii ) then
        call vtk_importer_%open_ascii_on() 
    else
        call vtk_importer_%open_ascii_off()
    endif

    call vtk_importer_%open_stream_on()

    call vtk_importer_%open_file(filename)
    call vtk_importer_%read_file(ugrid_, .true.)
    call vtk_importer_%close()

    call construct_flow_field(ugrid_, CELL_TYPE_VTK, FACE_VERT_DEF_VTK)

    call delete_ugrid(ugrid_)

end subroutine

! subroutine init_field_scflow(filename)
!     character(*),intent(in) :: filename


! end subroutine


subroutine init_particle_data(n, particles) bind(c, name="init_particle_data")
    integer(c_int),intent(in) :: n
    type(particle_t),intent(in) :: particles(n)

    call construct_particle_data(n)

    mv_pdata%particles(:) = particles(:)

end subroutine

subroutine init_randomize_particle(center, width, seed, n, radius)
    real(c_double),dimension(3) :: center
    real(c_double),dimension(3) :: width
    integer(c_int),intent(in) :: seed
    integer(c_int),intent(in) :: n
    real(c_double),intent(in) :: radius

    integer i
    real(DP) dr_(3)
    type(particle_t) p_

    call construct_particle_data(n)

    call set_rand_seed()

    do i = 1, n
        call random_number(dr_)
        ! pos(i) = [center(i)-width(i)/2, center(i)-width(i)/2]
        p_%pos = center + width*(dr_ - 0.5d0)
        p_%vel = 0.d0
        p_%f   = 0.d0
        p_%radius = radius
        p_%ref_cell = 1
        p_%state = PARTICLE_ACTIVATE

        ! store
        mv_pdata%particles(i) = p_
    end do

    contains
    subroutine set_rand_seed()
        integer seed_size
        integer,allocatable :: seeds(:)

        ! set rand seed
        call random_seed(size=seed_size)
        allocate(seeds(seed_size), source=seed)
        call random_seed(put=seeds)

    end subroutine
end subroutine

end module