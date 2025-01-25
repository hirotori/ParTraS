module motions_m
    use iso_c_binding
    use droplet_motion_m
    !use virus_droplet_motion_m
    use simulator_m, only: mv_simulator
    implicit none
    
contains

subroutine assign_droplet_motion(dt, L_ref, U_ref, rho_f, mu_f, rho_p, RK_order, body_force) bind(c, name="assign_droplet_motion")
    real(c_double),intent(in) :: dt
    real(c_double),intent(in) :: L_ref
    real(c_double),intent(in) :: U_ref
    real(c_double),intent(in) :: rho_f
    real(c_double),intent(in) :: mu_f
    real(c_double),intent(in) :: rho_p
    integer(c_int),intent(in) :: RK_order
    real(c_double),intent(in) :: body_force(3)

    type(droplet_motion_t) motion_

    call motion_%construct_droplet_motion(dt, L_ref, U_ref, rho_f, mu_f, rho_p, RK_order, body_force)

    call mv_simulator%set_motion(motion_)


end subroutine


! subroutine assign_virus_droplet_motion(rho_f, mu_f, rho_p, dt, RK_order)
!     real(8),intent(in) :: rho_f
!     real(8),intent(in) :: mu_f
!     real(8),intent(in) :: rho_p
!     real(8),intent(in) :: dt
!     integer,intent(in) :: RK_order

!     type(droplet_motion_t) motion_

!     call motion_%construct_droplet_motion(rho_f, mu_f, rho_p, dt, RK_order)

!     call mv_simulator%set_motion(motion_)


! end subroutine

end module motions_m