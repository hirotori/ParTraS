module motion_m
    use kind_parameters_m
    use particle_data_m, only: particle_t, PARTICLE_ACTIVATE, PARTICLE_INACTIVATE, mv_pdata
    implicit none
    
    type,abstract:: motion_t
        real(DP),private :: dt_
            !! time stepping size [T]
        real(DP),private :: Re_ref
            !! Reynolds number
        real(DP),dimension(3),private :: body_force_ = 0.d0
            !! body force (per mass) exerted on a particle.
            !! unit is [LT^-2]
    
        contains
        procedure construct_motion
        procedure,non_overridable :: get_body_force
        procedure,non_overridable :: get_dt
        procedure,non_overridable :: get_Re_ref
        procedure,non_overridable :: proceed_time_step
        procedure(integrate_one_step_impl),deferred :: integrate_one_step
        procedure(update_status_impl),deferred :: update_status
        procedure(compute_force_impl),deferred :: compute_force
    end type

contains

subroutine construct_motion(this, dt, Re, body_force)
    !! construct motion object.

    class(motion_t),intent(inout) :: this
    real(DP),intent(in) :: dt
        !! time stepping size [L/U]
    real(DP),intent(in) :: Re
        !! Reynolds number
    real(DP),dimension(3),intent(in) :: body_force
        !! body force (per mass) exerted on a particle. The unit is [U^2/L]

    this%dt_ = dt
    this%Re_ref = Re
    this%body_force_ = body_force

end subroutine


pure function get_body_force(this) result(bf)
    !! body force (per mass) exerted on a particle
    class(motion_t),intent(in) :: this
    real(DP),dimension(3) :: bf

    bf = this%body_force_

end function

pure real(DP) function get_dt(this) result(dt)
    class(motion_t),intent(in) :: this

    dt = this%dt_

end function

pure real(DP) function get_Re_ref(this) result(Re_ref)
    !! reference Raynolds number
    class(motion_t),intent(in) :: this

    Re_ref = this%Re_ref

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