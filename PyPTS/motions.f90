module motions_m
    use iso_c_binding
    use droplet_motion_m
    use droplet_motion_legacy_m
    !use virus_droplet_motion_m
    use simulator_m, only: mv_simulator
    implicit none
    
contains

subroutine delete_motion() bind(c, name="delete_motion")

    call mv_simulator%reset_motion()

end subroutine

subroutine assign_droplet_motion(dt, Re, rho_f, rho_p, RK_order, body_force) bind(c, name="assign_droplet_motion")
    real(c_double),intent(in) :: dt
    real(c_double),intent(in) :: Re
    real(c_double),intent(in) :: rho_f
    real(c_double),intent(in) :: rho_p
    integer(c_int),intent(in) :: RK_order
    real(c_double),intent(in) :: body_force(3)

    type(droplet_motion_t) motion_
    type(motion_legacy_t) motion_lg_

    if ( RK_order > 0 ) then        
        call motion_%construct_droplet_motion(dt, Re, rho_f, rho_p, RK_order, body_force)
        call mv_simulator%set_motion(motion_)
    else
        print"(A)", "motions/assign_droplet_motion::INFO:: use Legacy Scheme"
        call motion_lg_%construct_droplet_motion_legacy(dt, Re, rho_f, rho_p, body_force)
        call mv_simulator%set_motion(motion_lg_)
    end if


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