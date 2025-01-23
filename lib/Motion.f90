module motion_m
    use kind_parameters_m
    use particle_data_m, only: particle_t, PARTICLE_ACTIVATE, PARTICLE_INACTIVATE, mv_pdata
    implicit none
    
    type,abstract:: motion_t
        real(DP),private :: dt_
        real(DP),private :: L_ref
        real(DP),private :: U_ref
        real(DP),private :: Rho_ref
        real(DP),private :: Mu_ref
        real(DP),private :: Re_ref
    
        contains
        procedure construct_motion
        procedure,non_overridable :: get_dt
        procedure,non_overridable :: get_L_ref
        procedure,non_overridable :: get_U_ref
        procedure,non_overridable :: get_Rho_ref
        procedure,non_overridable :: get_Mu_ref
        procedure,non_overridable :: get_Re_ref
        procedure,non_overridable :: proceed_time_step
        procedure(integrate_one_step_impl),deferred :: integrate_one_step
        procedure(update_status_impl),deferred :: update_status
        procedure(compute_force_impl),deferred :: compute_force
    end type

contains

subroutine construct_motion(this, dt, L_ref, U_ref, Rho_ref, Mu_ref)
    !! construct motion object.
    !! equations of motions are non-dimensionalized by four characteristic variables:
    !! length, velocity, density, and viscosity.
    !! We assume density and viscosity of fluid to be reference density and viscosity, respectively.
    class(motion_t),intent(inout) :: this
    real(DP),intent(in) :: dt
        !! time stepping size [s]
    real(DP),intent(in) :: L_ref
        !! reference length [m]
    real(DP),intent(in) :: U_ref
        !! reference velocity [m/s]
    real(DP),intent(in) :: Rho_ref
        !! reference density [kg/m3]
    real(DP),intent(in) :: Mu_ref
        !! reference viscosity [Pa*s]

    this%dt_ = dt/(L_ref/U_ref)
    this%L_ref = L_ref
    this%U_ref = U_ref
    this%Rho_ref = Rho_ref
    this%Re_ref = Rho_ref*U_ref*L_ref/Mu_ref

end subroutine

pure real(DP) function get_dt(this) result(dt)
    class(motion_t),intent(in) :: this

    dt = this%dt_

end function

pure real(DP) function get_L_ref(this) result(L_ref)
    !! reference length
    class(motion_t),intent(in) :: this

    L_ref = this%L_ref

end function

pure real(DP) function get_U_ref(this) result(U_ref)
    !! reference velocity
    class(motion_t),intent(in) :: this

    U_ref = this%U_ref

end function

pure real(DP) function get_Rho_ref(this) result(Rho_ref)
    !! reference density
    class(motion_t),intent(in) :: this

    rho_ref = this%Rho_ref

end function

pure real(DP) function get_Mu_ref(this) result(mu_ref)
    !! reference viscosity
    class(motion_t),intent(in) :: this

    Mu_ref = this%Mu_ref

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