module droplet_motion_legacy_m
    use kind_parameters_m
    use motion_m
    use droplet_motion_m
    use particle_data_m
    use flow_field_m
    implicit none
    
    type,extends(droplet_motion_t) :: motion_legacy_t
        !! 今まで使われてきた離散化手法. 速度差の2次の項を異なる時間段階で近似. 
        !! したがって速度は実質1次精度と考えられる. 
        !! 半径の小さい粒子 (~ 1d-6)でも安定に計算できる.
        contains
        procedure :: construct_droplet_motion_legacy
        procedure :: integrate_one_step
        procedure :: update_status
    end type

    public :: motion_legacy_t

contains

subroutine construct_droplet_motion_legacy(this, dt, Re, rho_f, rho_p, gravity)
    !! construct an object describing the motion of droplets.

    class(motion_legacy_t),intent(inout) :: this
    real(DP),intent(in) :: dt
        !! time stepping size [L/U]
    real(DP),intent(in) :: Re
        !! Reynolds number of flow field
    real(DP),intent(in) :: rho_f
        !! density of fluid [ρ]
    real(DP),intent(in) :: rho_p
        !! density of particle [ρ]
    real(DP),dimension(3),intent(in) :: gravity
        !! gravity [U^2/L]

    call this%construct_droplet_motion(dt, Re, rho_f, rho_p, 0, gravity)

end subroutine

pure subroutine integrate_one_step(this, part)
    !! integrate governing equation to get next velocity and position of a partile
    !! integrated by k-step low-storage Runge-Kutta scheme
    class(motion_legacy_t),intent(in) :: this
    type(particle_t),intent(inout) :: part

    real(DP),dimension(3) :: r_, v_, u_, f_
    real(DP) rad_, dt, dv(3), dv_norm, Re_p_, coeff_
    integer k

    v_ = part%vel
    r_ = part%pos
    f_ = part%f
    rad_ = part%radius
    dt = this%get_dt()
    
    call search_reference_cell(part%pos(:), r_, part%ref_cell, 10)
    if ( part%ref_cell <= REFCELL_OUTBOUNDS ) then
        part%state = PARTICLE_INACTIVATE
        !TODO move particle onto wall
        v_ = 0.d0
        f_ = 0.d0
        return
    endif

    u_ = mv_flow_field%velocity(:,part%ref_cell)
    dv = u_ - v_
    dv_norm = sqrt(sum(dv*dv))
    Re_p_ = dv_norm*2*rad_*this%get_Re_ref()
    coeff_ = this%get_coeff_f()*this%Cd(Re_p_)*dv_norm/rad_

    v_(:) = (part%vel(:) + dt * (this%get_body_force() + coeff_ * u_(:))) / (1 + coeff_ * dt)
    r_(:) = part%pos(:) + (part%vel(:) + v_)*0.5d0*dt
    
    part%vel = v_
    part%pos = r_

    ! 力の計算は不要. ほぼデバッグ用. 
    call this%compute_force(v_, mv_flow_field%velocity(:,part%ref_cell), part%radius, f_)
    part%f   = f_

end subroutine

pure subroutine update_status(this, part, timestep)
    class(motion_legacy_t),intent(in) :: this
    type(particle_t),intent(inout) :: part
    integer,intent(in) :: timestep

    ! no update
    return

end subroutine

end module droplet_motion_legacy_m