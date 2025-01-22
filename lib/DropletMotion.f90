module droplet_motion_m
    use kind_parameters_m
    use particle_data_m, only: particle_t, PARTICLE_ACTIVATE, PARTICLE_INACTIVATE, mv_pdata
    use flow_field_m, only: mv_flow_field, search_reference_cell, REFCELL_OUTBOUNDS
    use motion_m
    implicit none
    private

    type,extends(motion_t) :: droplet_motion_t
        real(DP),private :: rho_f
            !! density of fluid
        real(DP),private :: mu_f
            !! viscosity of fluid
        real(DP),private :: rho_p
            !! density of particle
        integer,private :: n_rk
            !! number of steps in Runge-Kutta

        real(DP),private :: nu_f
            !! kinetic viscosity of fluid
        real(DP),private :: coeff_f

        real(DP),private :: g_(3) = [real(DP):: 0.0, 0.0, -9.81]

        contains
        procedure construct_droplet_motion
        procedure :: integrate_one_step 
        procedure update_status
        procedure :: compute_force => compute_force_droplet
    end type
    
    public droplet_motion_t

contains
    
subroutine construct_droplet_motion(this, rho_f, mu_f, rho_p, dt, RK_order)
    class(droplet_motion_t),intent(inout) :: this
    real(DP),intent(in) :: rho_f
    real(DP),intent(in) :: mu_f
    real(DP),intent(in) :: rho_p
    real(DP),intent(in) :: dt
    integer,intent(in) :: RK_order

    this%rho_f = rho_f
    this%rho_p = rho_p
    this%nu_f = mu_f/rho_f
    this%coeff_f = 3./8.*(rho_f/rho_p)
    call this%construct_motion(dt)
    this%n_rk = RK_order

end subroutine

pure subroutine integrate_one_step(this, part)
    !! integrate governing equation to get next velocity and position of a partile
    !! integrated by k-step low-storage Runge-Kutta scheme
    class(droplet_motion_t),intent(in) :: this
    type(particle_t),intent(inout) :: part

    real(DP),dimension(3) :: r_, v_, f_
    real(DP) alpha_k, dt
    integer k

    v_ = part%vel
    r_ = part%pos
    f_ = part%f
    dt = this%get_dt()
    do k = 1, this%n_rk
        alpha_k = 1.0/(this%n_rk + 1 - k)
        r_(:) = part%pos(:) + alpha_k*v_(:)*dt
        v_(:) = part%vel(:) + alpha_k*f_(:)*dt

        call search_reference_cell(part%pos(:), r_, part%ref_cell, 10)
        if ( part%ref_cell == REFCELL_OUTBOUNDS ) then
            part%state = PARTICLE_INACTIVATE
            !TODO move particle onto wall
            v_ = 0.d0
            f_ = 0.d0
            exit
        endif
        call this%compute_force(v_, mv_flow_field%velocity(:,part%ref_cell), part%radius, f_)
        
    end do
    
    part%vel = v_
    part%pos = r_
    part%f   = f_

end subroutine

pure subroutine update_status(this, part, timestep)
    class(droplet_motion_t),intent(in) :: this
    type(particle_t),intent(inout) :: part
    integer,intent(in) :: timestep

    ! no update
    return

end subroutine

pure subroutine compute_force_droplet(this, v, u, rad, f)
    class(droplet_motion_t),intent(in) :: this
    real(DP),intent(in) :: v(3)
        !! particle velocity
    real(DP),intent(in) :: u(3)
        !! fluid velocity
    real(DP),intent(in) :: rad
        !! radius
    real(DP),intent(out) :: f(3)
        !! computed force

    real(8) Re_p_, dv(3), dv_norm

    dv = u - v
    dv_norm = sqrt(sum(dv*dv))
    Re_p_ = dv_norm*2*rad/this%nu_f
    f(:) = this%coeff_f*Cd(Re_p_)/rad*dv_norm*dv + this%g_

end subroutine

pure real(8) function Cd(Re)
    real(8),intent(in) :: Re

    real(8),parameter :: min_eps = epsilon(1.d0)

    Cd = max(0.1, 24.d0/(Re+min_eps)*(1+0.15*Re**0.687))

end function

end module droplet_motion_m