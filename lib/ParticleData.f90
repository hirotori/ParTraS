module particle_data_m
    !! treats particle data.
    !! Positions, velocity, reference cell ids and state flag of particles 
    !! are stored as members of the derived-types called "particle_data_t".
    !! These members are public, thus they are allowed to be changed and updated by user.
    !! 
    use iso_c_binding, only: c_int, c_double
    use kind_parameters_m
    implicit none
    private

    type,bind(c) :: particle_t
        !! presenting a spherical particle. 
        real(c_double) pos(3)
            !! position
        real(c_double) vel(3)
            !! velocity
        real(c_double) f(3)
            !! force (per unit mass)
        real(c_double) radius
            !! radius of particle
        integer(c_int) ref_cell
            !! reference cell
        integer(c_int) state
            !! particle state. 1:active, 0:inactive
        real(c_double) :: radius_min = -100.d0
            !! minimum radius of particle (used in virus droplet)
        real(c_double) :: dead_time = -100.d0
            !! time at which the particle is inactivated  (used in virus droplet)
    end type

    integer,parameter :: PARTICLE_ACTIVATE = 1
        !! particle is activated
    integer,parameter :: PARTICLE_INACTIVATE = 0
        !! particle is inactivated

    integer,parameter :: PSTATE_FLOATING = 1
        !! particle is floating in a field
    integer,parameter :: PSTATE_STUCK_TO_BOUNDARY = 0
        !! particle is stuck to boundary

    type particle_data_t
        integer(IP) :: N_part
            !! number of particles
        type(particle_t),allocatable :: particles(:)
            !! particle
    end type

    type(particle_data_t) :: mv_pdata
        !! particle data object (module variable)

    interface write(formatted)
        module procedure print_particle_ascii
    end interface


    public particle_t, particle_data_t, mv_pdata, PARTICLE_ACTIVATE, PARTICLE_INACTIVATE, &
           construct_particle_data, &
           delete_particle_data, &
           add_particle, &
           get_particle_arrays, &
           export_pdata_ascii, export_pdata_binary, &
           import_pdata_ascii, import_pdata_binary, &
           write(formatted)

contains

subroutine construct_particle_data(n)
    integer,intent(in) :: n
        !! number of particle

    call delete_particle_data()

    if ( n <= 0 ) error stop "particle_data_m/construct_particle_data::error:: invalid number of particles (< 0)" 

    mv_pdata%N_part = n
    allocate(mv_pdata%particles(n))

    print "('particle_data_m/construct_particle_data::INFO:: particle data constructed, n = ', i0)", mv_pdata%N_part

end subroutine

subroutine delete_particle_data()
    !! delete particle data
    if ( allocated(mv_pdata%particles) ) then
        deallocate(mv_pdata%particles)
        mv_pdata%N_part = 0
    end if    
end subroutine

subroutine add_particle(particles)
    !! add many new particles to the current particle data. 
    type(particle_t),allocatable :: particles(:)
        !! array of particle_t to be added to current particle data
        !! this variable will be deallocated after calling this
    integer n_
    type(particle_t),allocatable :: pdata_tmp_(:)

    if ( .not. allocated(mv_pdata%particles) ) then
        error stop "particle_data_m/add_particle::error:: pdata is not allocated."
    end if

    n_ = size(particles)
    allocate(pdata_tmp_(mv_pdata%N_part+n_))
    pdata_tmp_(1:mv_pdata%N_part) = mv_pdata%particles
    pdata_tmp_(mv_pdata%N_part+1:) = particles

    call move_alloc(from=pdata_tmp_, to=mv_pdata%particles)
    mv_pdata%N_part = mv_pdata%N_part + n_
    deallocate(particles)

end subroutine

subroutine get_particle_arrays(particles)
    !! get deep copy of `particle_data_t%particles`
    type(particle_t),allocatable :: particles(:)

    if ( allocated(particles) ) then
        error stop "particle_data_m/get_particle_arrays::ERROR:: array given in argument is already allocated"
    end if

    allocate(particles(mv_pdata%N_part))
    particles(:) = mv_pdata%particles(:) !deep copy

end subroutine

subroutine export_pdata_ascii(unit)
    !! export particle data
    integer,intent(in) :: unit
    integer i

    write(unit, "(A)") "number of particle"
    write(unit,"(i0)") mv_pdata%N_part
    write(unit, "(A)") "pos, vel, f, radius, ref_cell, state, radius_min, dead_time"
    do i = 1, mv_pdata%N_part
        write(unit, *) mv_pdata%particles(i)
    end do

end subroutine

subroutine export_pdata_binary(unit)
    !! export particle data
    integer,intent(in) :: unit
    integer i

    write(unit) mv_pdata%N_part
    ! do i = 1, pdata_%N_part
    !     write(unit) pdata_%particles(i)
    ! end do

    write(unit) mv_pdata%particles

end subroutine


subroutine import_pdata_ascii(unit)
    !! load particle data from file
    integer,intent(in) :: unit
    integer i, n

    read(unit, "(A)") 
    read(unit,   *  ) n
    read(unit, "(A)") 
    call construct_particle_data(n)
    ! do i = 1, pdata%N_part
    !     read(unit, *) pdata%particles(i)
    ! end do

    read(unit, *) mv_pdata%particles

end subroutine

subroutine import_pdata_binary(unit)
    !! load particle data from file
    integer,intent(in) :: unit
    integer i, n

    read(unit) n
    call construct_particle_data(n)
    ! do i = 1, pdata%N_part
    !     read(unit) pdata%particles(i)
    ! end do

    read(unit) mv_pdata%particles

end subroutine

subroutine print_particle_ascii(this, unit, iotype, arglist, iostatus, iomessage)
    type(particle_t),intent(in) :: this
    integer,intent(in) :: unit
    character(*),intent(in) :: iotype
    integer,intent(in) :: arglist(:)
    integer,intent(out) :: iostatus
    character(*),intent(inout) :: iomessage

    if ( iotype == "LISTDIRECTED" .or. size(arglist) < 2) then
        write(unit, fmt=*, iostat=iostatus, iomsg=iomessage) &
        this%pos, this%vel, this%f, this%radius, this%ref_cell, this%state, this%radius_min, this%dead_time
    else
        if ( iotype(3:) /= "particle_t" ) then
            error stop "particle_data_m/print_particle_ascii::ERROR:: type mismatch"
        else
            block
                character(2) width_tol, width_dec
                character(6) realspec
                character(:),allocatable :: fmt

                write(width_tol, "(I2)") arglist(1)
                write(width_dec, "(I2)") arglist(2)
                realspec = "F"//width_tol//"."//width_dec
                fmt = "(10("//realspec//",1x), 2(i0,1x), 2("//realspec//",1x))"
                write(unit, fmt=fmt, iostat=iostatus, iomsg=iomessage) &
                this%pos, this%vel, this%f, this%radius, this%ref_cell, this%state, this%radius_min, this%dead_time

            end block
        end if
    end if

end subroutine

end module