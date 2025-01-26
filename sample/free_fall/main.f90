program main
    use base_importer_m, only: ugrid_struct_t, delete_ugrid
    use vtk_importer_m
    use flow_field_m
    use particle_data_m
    use droplet_motion_m
    !$ use omp_lib
    implicit none

    integer i, n
    real(8) :: zmax
    real(8),parameter :: L = 1.0d0
    real(8),parameter :: U = 1.0d0
    real(8),parameter :: RHO = 1.0d0
    real(8),parameter :: MU = 1d-5
    real(8),parameter :: Re = RHO*U*L/MU
    
    integer,parameter :: n_rk = 4
    real(8),parameter :: rho_p = 1000.d0/RHO
    real(8),parameter :: rho_f = 1.d0/RHO
    real(8),parameter :: dt = 0.001/(L/U)
    real(8),parameter :: g(3) = [0.0, 0.0, -9.81]/(U*U/L)

    type(vtk_importer_t) vtk_importer
    type(ugrid_struct_t) vtk_ugrid
    type(droplet_motion_t) motion
    !$ integer nthread = 1

    !$ call omp_set_num_threads(nthread)

    call vtk_importer%open_stream_on()
    call vtk_importer%open_ascii_on()
    call vtk_importer%open_file("quiescent.vtk")
    call vtk_importer%read_file(vtk_ugrid, shift_index=.true.)
    call construct_flow_field(vtk_ugrid, CELL_TYPE_VTK, FACE_VERT_DEF_VTK)
    call delete_ugrid(vtk_ugrid)

    call construct_particle_data(200)
    
    ! initialize particle status
    zmax = maxval(mv_flow_field%cell_centers(3,:))
    do i = 1, mv_pdata%N_part
        mv_pdata%particles(i)%pos = [0.4d0, 0.5d0, zmax]
        mv_pdata%particles(i)%vel = 0.d0
        mv_pdata%particles(i)%state = PARTICLE_ACTIVATE
        mv_pdata%particles(i)%radius = 1d-5
        mv_pdata%particles(i)%ref_cell = 1

        call search_reference_cell(mv_pdata%particles(i)%pos, mv_pdata%particles(i)%pos, mv_pdata%particles(i)%ref_cell, 100)

        call compute_force(mv_pdata%particles(i)%vel, mv_flow_field%velocity(:,mv_pdata%particles(i)%ref_cell), &
                           mv_pdata%particles(i)%radius, mv_pdata%particles(i)%f)

    end do

    call motion%construct_droplet_motion(dt, Re, rho_f, rho_p, n_rk, g)
    ! Heun法 (2次ルンゲ・クッタ法) で離散化している.
    ! dt = 0.01とすると粒子が上に移動する. 力が+-振動し, 振幅が大きくなっていく. 
    ! dt = 0.001とすると計算はできた. 
    ! Cd = const.の場合, 終端速度が解析解と一致した. 
    ! スキームを low-storage Runge-Kutta に変更した. 
    print "('pos, vel, vel (exa.)')"
    do n = 1, 100000
        
        call motion%proceed_time_step(n)
        
        if (mod(n,10000) == 0) print "(*(g0.5,1x))", mv_pdata%particles(1)%pos(3), mv_pdata%particles(1)%vel(3)
        !sqrt(8*rho_p*rad_p*9.8/(3*rho_f*Cd))

    end do

    contains

    pure subroutine compute_force(v, u, radi, f)
        real(8),dimension(3),intent(in) :: v
        real(8),dimension(3),intent(in) :: u
        real(8),intent(in) :: radi
        real(8),dimension(3),intent(out) :: f

        real(8) Re_p_, dv(3), dv_norm

        dv = u - v
        dv_norm = sqrt(sum(dv*dv))
        Re_p_ = rho_f*dv_norm*2*radi/MU
        f(:) = 3.0/8.0*rho_f*Cd(Re_p_)/(rho_p*radi)*dv_norm*dv
        f(3) = f(3) - 9.81d0

    end subroutine

    pure real(8) function Cd(Re)
        real(8),intent(in) :: Re
        
        real(8),parameter :: min_eps = epsilon(1.d0)

        Cd = max(0.1, 24.d0/(Re+min_eps)*(1+0.15*Re**0.687))

    end function

end program main