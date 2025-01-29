module droplet_motion_m
    use kind_parameters_m
    use particle_data_m, only: particle_t, PARTICLE_ACTIVATE, PARTICLE_INACTIVATE, mv_pdata
    use flow_field_m, only: mv_flow_field, search_reference_cell, REFCELL_OUTBOUNDS
    use motion_m
    implicit none
    private

    type,extends(motion_t) :: droplet_motion_t
        !! 飛沫の運動. Runge-Kutta法で解く. 
        !! 半径が小さい (<= 1d-6)とき, 安定に計算できる時間刻み幅がとても小さい (dt <= 1e-6). 
        integer,private :: n_rk
            !! number of steps in Runge-Kutta

        real(DP),private :: coeff_f

        contains
        procedure construct_droplet_motion
        procedure :: integrate_one_step 
        procedure update_status
        procedure :: compute_force => compute_force_droplet
        procedure Cd
        procedure get_coeff_f
        procedure,non_overridable :: get_RK_order
    end type
    
    public droplet_motion_t

contains
    
subroutine construct_droplet_motion(this, dt, Re, rho_f, rho_p, RK_order, gravity)
    !! construct an object describing the motion of droplets.

    class(droplet_motion_t),intent(inout) :: this
    real(DP),intent(in) :: dt
        !! time stepping size [L/U]
    real(DP),intent(in) :: Re
        !! Reynolds number of flow field
    real(DP),intent(in) :: rho_f
        !! density of fluid [ρ]
    real(DP),intent(in) :: rho_p
        !! density of particle [ρ]
    integer,intent(in) :: RK_order
        !! the order of Runge-Kutta scheme
    real(DP),dimension(3),intent(in) :: gravity
        !! gravity [U^2/L]

    this%coeff_f = 3./8.*(rho_f/rho_p)
    call this%construct_motion(dt, Re, gravity)
    this%n_rk = RK_order

end subroutine

pure integer(IP) function get_RK_order(this)
    !! order of Runge-Kutta scheme
    class(droplet_motion_t),intent(in) :: this

    get_RK_order = this%n_rk

end function

pure real(DP) function get_coeff_f(this)
    !! coefficient for drag force
    class(droplet_motion_t),intent(in) :: this

    get_coeff_f = this%coeff_f

end function


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
        if ( part%ref_cell <= REFCELL_OUTBOUNDS ) then
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
    Re_p_ = dv_norm*2*rad*this%get_Re_ref()
    f(:) = this%coeff_f*this%Cd(Re_p_)/rad*dv_norm*dv + this%get_body_force()

end subroutine

pure real(8) function Cd(this, Re)
    !! drag coefficient. 
    !! Cd = 0.1 is the limit value for Re -> infinity (eulerian flow).
    class(droplet_motion_t),intent(in) :: this
    real(8),intent(in) :: Re

    real(8),parameter :: min_eps = epsilon(1.d0)

    Cd = max(0.1, 24.d0/(Re+min_eps)*(1+0.15*Re**0.687))

end function

end module droplet_motion_m