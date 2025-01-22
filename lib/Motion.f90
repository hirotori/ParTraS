module motion_m
    use kind_parameters_m
    use particle_data_m, only: particle_t, PARTICLE_ACTIVATE, PARTICLE_INACTIVATE, mv_pdata
    implicit none
    
    type,abstract:: motion_t
        real(DP),private :: dt_
        contains
        procedure construct_motion
        procedure,non_overridable :: get_dt
        procedure,non_overridable :: proceed_time_step
        procedure(integrate_one_step_impl),deferred :: integrate_one_step
        procedure(update_status_impl),deferred :: update_status
        procedure(compute_force_impl),deferred :: compute_force
    end type

contains

subroutine construct_motion(this, dt)
    class(motion_t),intent(inout) :: this
    real(DP),intent(in) :: dt

    this%dt_ = dt

end subroutine

pure real(DP) function get_dt(this) result(dt)
    class(motion_t),intent(in) :: this

    dt = this%dt_

end function

subroutine proceed_time_step(this, timestep)
    !! この処理はスレッドセーフではない. 
    !! @NOTE: pdataを引数にわたす形にすれば解決はするはず. 
    class(motion_t),intent(in) :: this
    integer,intent(in) :: timestep

    integer i
    type(particle_t) :: part_

    !$omp parallel do private(part_)
    do i = 1, mv_pdata%N_part

        ! load
        part_ = mv_pdata%particles(i)

        if ( part_%state == PARTICLE_INACTIVATE ) cycle 

        call this%update_status(part_, timestep)

        call this%integrate_one_step(part_)


        mv_pdata%particles(i) = part_

    end do

end subroutine

pure subroutine integrate_one_step_impl(this, part)
    !! integrate governing equation to get next velocity and position of a partile
    class(motion_t),intent(in) :: this
    type(particle_t),intent(inout) :: part

end subroutine

pure subroutine update_status_impl(this, part, timestep)
    !! update particle status
    !! this is called before `integrate_one_step`.
    !! this procedure must be thread-safe.
    class(motion_t),intent(in) :: this
    type(particle_t),intent(inout) :: part
    integer,intent(in) :: timestep

end subroutine

pure subroutine compute_force_impl(this, v, u, rad, f)
    !! implementation of `compute_force`.
    !! the unit of `f` is "force per mass".
    class(motion_t),intent(in) :: this
    real(DP),intent(in) :: v(3)
        !! particle velocity
    real(DP),intent(in) :: u(3)
        !! fluid velocity
    real(DP),intent(in) :: rad
        !! radius
    real(DP),intent(out) :: f(3)
        !! computed force

end subroutine

end module motion_m