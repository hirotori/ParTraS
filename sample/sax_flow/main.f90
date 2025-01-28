program main
    use kind_parameters_m
    use base_importer_m, only: ugrid_struct_t, delete_ugrid
    use vtk_importer_m
    use flow_field_m
    use particle_data_m
    use droplet_motion_m
    use field_updater_m
    use simulator_m
    !$ use omp_lib

    implicit none
    integer nstep, nactive
    integer,parameter :: Npart = 200
    real(DP),dimension(3) :: r0
    type(droplet_motion_t) motion
    type(ugrid_struct_t) vtk_ugrid
    type(vtk_importer_t) vtk_importer
    logical :: write_ascii = .true.
    !$ integer :: nthread = 1

    !$ call omp_set_num_threads(nthread)

    call vtk_importer%open_file("sax_flow.vtk")
    call vtk_importer%read_file(vtk_ugrid, .true.)
    call construct_flow_field(vtk_ugrid, CELL_TYPE_VTK, FACE_VERT_DEF_VTK)
    call vtk_importer%close()
    call delete_ugrid(vtk_ugrid)

    call construct_particle_data(n=Npart)

    block
        integer i
        real(DP) dr(3), r_(3)
        integer seed, seed_size
        integer,allocatable :: seeds(:)

        ! set rand seed
        seed = 100
        call random_seed(size=seed_size)
        allocate(seeds(seed_size), source=seed)
        call random_seed(put=seeds)
        
        r0 = mv_flow_field%cell_centers(:,500) + [0.d0, 0.2d0, 0.d0]
        do i = 1, Npart
            call random_number(dr)
            dr = 0.02*dr - 0.01! [0,1) to [-0.01,0.01)
            r_ = r0 + dr
            mv_pdata%particles(i)%pos = r_
            mv_pdata%particles(i)%vel = 0.0
            mv_pdata%particles(i)%radius = 1d-5
            mv_pdata%particles(i)%state = PARTICLE_ACTIVATE
            mv_pdata%particles(i)%f   = 0.0
            mv_pdata%particles(i)%ref_cell = 1
            call search_reference_cell(r0, r_, mv_pdata%particles(i)%ref_cell, 100)
            if (i == 1) print*, mv_pdata%particles(i)%ref_cell
        end do
    end block

    block
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

        call motion%construct_droplet_motion(dt, Re, rho_f, rho_p, n_rk, g)
    end block

    call mv_field_updater%disable_updater()

    call mv_simulator%construct_simulator(10000, write_ascii, "./data/trajectory", 10000, .false., "./data/backup")
    call mv_simulator%set_motion(motion)
    call mv_simulator%run(1, 100000)

    ! run again
    ! call mv_simulator%set_dump_settings(1000, .true.)
    ! call mv_simulator%run(1001, 2000)

end program main