module particle_m
    use iso_c_binding
    use util_m
    use kind_parameters_m
    use particle_data_m
    implicit none
    
    type(particle_t),allocatable,private :: tmp_(:)
        !! temporary particle array used for passing data to python

contains
    
subroutine init_particle_data(n, particles) bind(c, name="init_particle_data")
    integer(c_int),intent(in) :: n
    type(particle_t),intent(in) :: particles(n)

    call construct_particle_data(n)

    mv_pdata%particles(:) = particles(:)

end subroutine


subroutine get_particle_data_size(n) bind(c, name="get_particle_data_size")
    integer(c_int),intent(out) :: n

    n = mv_pdata%N_part

end subroutine

subroutine get_particle_data(n, particles) bind(c, name="get_particle_data")
    integer(c_int),intent(in) :: n
    type(particle_t),intent(inout) :: particles(n)

    if ( n /= mv_pdata%N_part ) then
        error stop "particle_m/get_particle_data::ERROR:: number of particles doesn't match with current number"
    end if

    call get_particle_arrays(tmp_)

    particles(1:n) = tmp_(1:n)

    deallocate(tmp_)

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



end module particle_m