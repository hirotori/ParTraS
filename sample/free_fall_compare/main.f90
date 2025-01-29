program main
    use base_importer_m, only: ugrid_struct_t, delete_ugrid
    use vtk_importer_m
    use flow_field_m
    use particle_data_m
    use droplet_motion_m
    use droplet_motion_legacy_m
    !$ use omp_lib
    implicit none

    integer i, n
    real(8),parameter :: tend = 0.02
    
    real(8) :: zmax
    real(8),parameter :: L = 1.0d0
    real(8),parameter :: U = 1.0d0
    real(8),parameter :: RHO = 1.0d0
    real(8),parameter :: MU = 1d-5
    real(8),parameter :: Re = RHO*U*L/MU
    
    integer,parameter :: n_rk = 4
    real(8),parameter :: rad_p = 9.d-7/L !1d-6 = 1μm
    real(8),parameter :: rho_p = 1000.d0/RHO
    real(8),parameter :: rho_f = 1.d0/RHO
    real(8),parameter :: dt = 0.00001/(L/U)
    real(8),parameter :: g(3) = [0.0, 0.0, -9.81]/(U*U/L)

    integer,parameter :: nend = int(tend/dt)
    integer,parameter :: nwrite = 1

    type(vtk_importer_t) vtk_importer
    type(ugrid_struct_t) vtk_ugrid
    type(droplet_motion_t) motion
    type(motion_legacy_t) motion_lg
    !$ integer nthread = 1

    !$ call omp_set_num_threads(nthread)

    call vtk_importer%open_stream_on()
    call vtk_importer%open_ascii_on()
    call vtk_importer%open_file("../free_fall/quiescent.vtk")
    call vtk_importer%read_file(vtk_ugrid, shift_index=.true.)
    call construct_flow_field(vtk_ugrid, CELL_TYPE_VTK, FACE_VERT_DEF_VTK)
    call delete_ugrid(vtk_ugrid)


    !legacy
    call motion_lg%construct_droplet_motion_legacy(dt, Re, rho_f, rho_p, g)

    call free_fall_test(nend, nwrite, motion_lg, "Legacy.dat")

    !RK-4
    call motion%construct_droplet_motion(dt, Re, rho_f, rho_p, n_rk, g)

    call free_fall_test(nend, nwrite, motion, "RK4.dat")

    ! 結果: z方向速度 vs 時間のプロット
    ! 1. Cd = 0.3 = const.のとき. 
    ! dt = 0.1で, RKは大雑把に厳密解に乗る. Legacyは速度の変化がおかしい (増加して収束する.)
    ! dt = 0.01以下で, どちらも厳密解に乗る. 
    ! radによらず, 正確に計算できる. 

    ! 2. Cd = f(Re)のとき. (rad = 1d-5)
    ! dt >= 0.01 だと RK4は安定に計算できないっぽい. Legacyは安定. 
    ! dt <= 0.001 ではどちらも同じ軌道. dt = 0.0001では実質問題なさそう. 
    ! Rk4の計算安定性は半径に依存していそう. 
    ! - rad = 1d-5のとき: dt = 0.001以下ならOK
    ! - rad = 1d-6のとき: dt = 0.0001でもNG. 
    ! RKは0.1μmの粒子に対して不安定. Legacyはどちらも安定. ただし, Legacyの場合
    ! rad <= 1d-6 (1μm)で, dt = 0.001だと, 速度が上昇してから一定値に収束する. 最初の値が v = -0.981*0.01. 重力の大きさ?
    ! dt = 1d-6だと, LegacyとRKともにうまくいく. ただ t = 1e-4くらいで速度が収束していることを考えると, 安定に計算できればそれで良さそう. 
    ! rad = 1d-7の場合は dt = 1e-8まで下げる必要がある. 

    ! RKの対処療法としては, 新しい速度と前の速度との差が, 前の速度と流体の速度との差を超えないようにリミタをかける. (Löhner et al.)
    
    contains

    subroutine free_fall_test(nend, nw, motion_obj, fname)
        integer,intent(in) :: nend, nw
        class(motion_t),intent(in) :: motion_obj
        character(*),intent(in) :: fname

        real(DP) re_, u_(3), v_(3), dv(3), dv_norm, cd_

        call construct_particle_data(1)
    
        ! initialize particle status
        zmax = maxval(mv_flow_field%cell_centers(3,:))
        do i = 1, mv_pdata%N_part
            mv_pdata%particles(i)%pos = [0.4d0, 0.5d0, zmax]
            mv_pdata%particles(i)%vel = 0.d0
            mv_pdata%particles(i)%state = PARTICLE_ACTIVATE
            mv_pdata%particles(i)%radius = rad_p
            mv_pdata%particles(i)%ref_cell = 1
            mv_pdata%particles(i)%f = 0.d0
    
            call search_reference_cell(mv_pdata%particles(i)%pos, mv_pdata%particles(i)%pos, mv_pdata%particles(i)%ref_cell, 100)
            call motion_obj%compute_force(mv_pdata%particles(i)%vel, mv_flow_field%velocity(:,mv_pdata%particles(i)%ref_cell), &
                                         mv_pdata%particles(i)%radius, mv_pdata%particles(i)%f)
        end do
    
        print *, "ref cell = ", mv_pdata%particles(1)%ref_cell
    

        open(unit=99, file=fname)
        print "('n, count, Cd, pos, vel')"
        do n = 1, nend
            
            call motion_obj%proceed_time_step(n)
            
            if (mv_pdata%particles(1)%state == PARTICLE_ACTIVATE) then
                u_ = mv_flow_field%velocity(:,mv_pdata%particles(1)%ref_cell)
                dv = u_ - mv_pdata%particles(1)%vel(:)
                dv_norm = sqrt(sum(dv*dv))
                re_ = dv_norm*2*rad_p*motion_obj%get_Re_ref()
                ! cd_ = drag_coeff(re_)
                cd_ = 24.d0/(re_+epsilon(1.d0))*(1+0.15*re_**0.687)
            endif
            
            if (mod(n,nw) == 0) print "(2(i0,1x), 5(g0.5,1x))", n, count(mv_pdata%particles(:)%state == PARTICLE_ACTIVATE), re_, cd_, &
                                                                 mv_pdata%particles(1)%pos(3), mv_pdata%particles(1)%vel(3), mv_pdata%particles(1)%f(3)
                                
            !sqrt(8*rho_p*rad_p*9.8/(3*rho_f*Cd))
            write(99, "(*(g0,1x))") n*dt, mv_pdata%particles(1)%pos(3), exact_pos(zmax, n, dt), mv_pdata%particles(1)%vel(3), mv_pdata%particles(1)%f(3)
            
        end do
    
        close(99)
        call delete_particle_data()

    end subroutine

    real(DP) function exact_vel(n, dt)
        !! exact solution for velocity with constant drag coefficient
        integer,intent(in) :: n
        real(DP),intent(in) :: dt

        real(DP) D
        real(DP) :: Cd = 0.1

        D = 3.*rho_f/8./rho_p/rad_p*Cd
        exact_vel = sqrt(-g(3)/D)*tanh(-sqrt(-D*g(3))*n*dt)

    end function

    real(DP) function exact_pos(z0, n, dt)
        !! exact solution for position with constant drag coefficient
        real(DP),intent(in) :: z0
        integer,intent(in) :: n
        real(DP),intent(in) :: dt

        real(DP) D, sq_Dg
        real(DP) :: Cd = 0.1
        
        D = 3.*rho_f/8./rho_p/rad_p*Cd
        sq_Dg = sqrt(-D*g(3))
        exact_pos = z0 - 1./D*(log(2.d0*(exp(-sq_Dg*n*dt) + exp(sq_Dg*n*dt))))

    end function

    real(8) function drag_coeff(Re_p_)
        real(8),intent(in) :: Re_p_

        real(8),parameter :: min_eps = epsilon(1.d0)

        ! drag_coeff = max(0.1, 24.d0/(Re+min_eps)*(1+0.15*Re**0.687))
        drag_coeff = 0.1!24.d0/(Re+min_eps)*(1+0.15*Re**0.687)
    end function


end program main