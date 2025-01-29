module simulator_m
    use kind_parameters_m
    use motion_m, only: motion_t
    use particle_data_m
    use base_importer_m
    use flow_field_m
    use dump_file_m, only: write_out_backup, writeout_vtk
    use field_updater_m
    use version_m
    !$ use omp_lib
    implicit none
    private
    
    type,abstract:: callback_t
        !! treats callback process.
        contains
        procedure(process_callback_imp),deferred :: process_callback
    end type

    type simulator_t
        integer,private :: ncall_
            !! time interval at which this is called
        integer,private :: nwrite_
            !! time interval at which the data is written
        integer,private :: nback_
            !! time interval at which the backup data is written
        class(callback_t),allocatable :: callback_
            !! registered callback process object
        class(motion_t),allocatable :: motion_
            !! registered particle motion object

        logical backup_ascii_
        logical write_ascii_
        character(:),allocatable :: basename_vtk_
        character(:),allocatable :: basename_back_
        contains
        procedure,non_overridable :: construct_simulator
        procedure,non_overridable :: set_dump_settings
        procedure,non_overridable :: set_writeout_settings
        procedure,non_overridable :: set_callback
        procedure,non_overridable :: reset_callback
        procedure,non_overridable :: set_motion
        procedure,non_overridable :: reset_motion
        procedure,non_overridable :: run
    end type

    type(simulator_t) mv_simulator
        !! simulator object (module variable)


    public :: callback_t, mv_simulator
contains


! ~~~~~~~~~~~~~~~~~~ callback  ~~~~~~~~~~~~~~~~~~~~~~
subroutine process_callback_imp(this, particle_data)
    class(callback_t),intent(inout) :: this
    type(particle_data_t),intent(inout) :: particle_data

end subroutine

! ~~~~~~~~~~~~~~~~~~ simulator ~~~~~~~~~~~~~~~~~~~~~~
subroutine construct_simulator(this, nwrite, write_ascii, path_write, nback, backup_ascii, path_back)
    !! construct simulator object
    !! particle data and flow field data are module variables, thus we needn't pass them to this object.
    class(simulator_t),intent(inout) :: this
    integer,intent(in) :: nwrite
    logical,intent(in) :: write_ascii
    character(*),intent(in) :: path_write
    integer,intent(in) :: nback
    logical,intent(in) :: backup_ascii
    character(*),intent(in) :: path_back

    call this%set_dump_settings(nback, backup_ascii, path_write)
    call this%set_writeout_settings(nwrite, write_ascii, path_back)

end subroutine

subroutine set_writeout_settings(this, nwrite, write_ascii, path_write)
    !! construct simulator object
    class(simulator_t),intent(inout) :: this
    integer,intent(in) :: nwrite
    logical,intent(in) :: write_ascii
    character(*),intent(in) :: path_write
    logical exists
    
    this%nwrite_ = nwrite
    this%write_ascii_ = write_ascii
    inquire(file=path_write, exist=exists)
    if ( .not. exists ) then
        error stop "simulator/set_writeout_settings::ERROR:: directory """//path_write//""" not found."
    end if

    this%basename_vtk_ = path_write//get_sep__()//"trajectory"

end subroutine

subroutine set_dump_settings(this, nback, backup_ascii, path_back)
    !! construct simulator object
    class(simulator_t),intent(inout) :: this
    integer,intent(in) :: nback
    logical,intent(in) :: backup_ascii
    character(*),intent(in) :: path_back
    logical exists

    this%nback_  = nback
    this%backup_ascii_ = backup_ascii
    inquire(file=path_back, exist=exists)
    if ( .not. exists ) then
        error stop "simulator/set_backup_settings::ERROR:: directory """//path_back//""" not found."
    end if
    this%basename_back_ = path_back//get_sep__()//"backup"

end subroutine

subroutine set_motion(this, motion)
    !! set motion object
    !! raise error if motion is already allocated
    class(simulator_t),intent(inout) :: this
    class(motion_t),intent(in) :: motion

    if ( .not. allocated(this%motion_) ) then
        allocate(this%motion_, source=motion)
    else
        error stop "simulator/set_motion::ERROR:: motion_t is already assigned"
    end if

end subroutine

subroutine reset_motion(this)
    !! deallocate motion object if it is allocated
    class(simulator_t),intent(inout) :: this

    if ( allocated(this%motion_) ) then
        deallocate(this%motion_)
    end if

end subroutine

subroutine set_callback(this, callback, ncall)
    class(simulator_t),intent(inout) :: this
    class(callback_t),intent(in) :: callback
    integer,intent(in) :: ncall

    if ( .not. allocated(this%callback_) ) then
        allocate(this%callback_, source=callback)
    else
        error stop "simulator/set_motion::ERROR:: callback_t is already assigned"
    end if
    this%ncall_ = ncall

end subroutine


subroutine reset_callback(this)
    class(simulator_t),intent(inout) :: this

    if ( allocated(this%callback_) ) then
        deallocate(this%callback_)
    end if
    this%ncall_ = 0

end subroutine

subroutine run(this, nstart, nend)
    class(simulator_t),intent(inout) :: this
    integer,intent(in) :: nstart
    integer,intent(in) :: nend

    integer ncyc, i, nactive
    !$ integer nthread, maxthread

    print "(A)", "============================================== "
    print "(A)", "          ParTraS - Version "//get_version()
    print "(A)", "       Particle Trajectory Simulator           "
    print "(A)", "============================================== "
    print "(A)", " " 
    print "(A)", "Developed by: Takahiro Ikeda, AFDET lab.         "
    print "(A)", "Description:                                     "
    print "(A)", "  ParTras simulates the advection and tracking   " 
    print "(A)", "  of droplets using the Lagrange particle method."
    print "(A)", "  It is designed for studying particle dynamics  "
    print "(A)", "  in various environments.                       "
    print "(A)", " "
    print "(A)", "Copyright (C) 2025-                              "
    print "(A)", "Takahiro Ikeda, and                              "
    print "(A)", "Advanced Fluid Dynamics and Energy Transfer Lab. "
    print "(A)", "(AFDET) All rights reserved.                     "
    print "(A)", "==============================================   "

    
    ! validation
    if ( mv_pdata%N_part <= 0 ) then
        error stop "simulator_m/simulator_t%run::ERROR:: no particles initialized"        
    end if

    if ( .not. mv_flow_field%is_assigned ) then
        error stop "simulator_m/simulator_t%run::ERROR:: no flow field initialized"        
    end if

    if ( .not. allocated(this%motion_) ) then
        error stop "simulator_m/simulator_t%run::ERROR:: no motion_t object allocated"
    end if

    if ( .not. mv_field_updater%assigned() ) then
        error stop "simulator_m/simulator_t%run::ERROR:: field updater is not assigned"
    end if

    print "(A)", "simulator_m/simulator_t%run:: Simulation started"

    if ( nstart <= 1 ) then

        print "(A)", "simulator_m/simulator_t%run:: Search reference cells ..."
        ! reference cell search
        i = 1
        call search_reference_cell(mv_pdata%particles(i)%pos, mv_pdata%particles(i)%pos, &
                                   mv_pdata%particles(i)%ref_cell, 100)
        ! assume that other particles are near from 1 
        if ( mv_pdata%particles(1)%ref_cell <= REFCELL_OUTBOUNDS ) then
            error stop "simulator_m/simulator_t%run::ERROR:: proper reference cell not found. The particle '1' may be out of field."
        end if
        mv_pdata%particles(2:)%ref_cell = mv_pdata%particles(1)%ref_cell
        do i = 2, mv_pdata%N_part
            call search_reference_cell(mv_pdata%particles(i)%pos, mv_pdata%particles(i)%pos, &
            mv_pdata%particles(i)%ref_cell, 100)
        end do
    else

        ! ACTIVATE状態でref_cell<=0はおかしい. 
        do i = 1, mv_pdata%N_part
            if ( mv_pdata%particles(i)%ref_cell <= REFCELL_NOTFOUND .and. mv_pdata%particles(i)%state <= PARTICLE_ACTIVATE) then
                print "('simulator_m/simulator_t%run::ERROR:: active particle with invalid ref_cell ', i0)", mv_pdata%particles(i)%ref_cell
                error stop "simulator_m/simulator_t%run::ERROR:: active particle must be located in a cell"
            endif  
        end do
    
    end if

    !$ print "('simulator_m/simulator_t%run::INFO::Simulation with OpenMP')"
    !$omp parallel
    !$omp single
    !$ nthread = omp_get_num_threads()
    !$ maxthread = omp_get_max_threads()
    !$omp end single
    !$omp end parallel
    !$ print "('simulator_m/simulator_t%run::INFO::Num of thread = ', i0, '/', i0)", nthread, maxthread

    do ncyc = nstart, nend

        nactive = count(mv_pdata%particles(:)%state == PARTICLE_ACTIVATE)

        if ( mod(ncyc, this%nwrite_) == 0 ) then
            !TODO: コンソール出力のタイミング制御と, 何を表示するかの策定. 
            print "(A)",    "-----------------------------------------------"
            print "(A,i0,A,f16.5,A)", "Timestep = ", ncyc, " t = (", ncyc*this%motion_%get_dt(), ") "
            !$ print "(A,i0)", "Num Thread = ", nthread
            print "(A)",    "-----------------------------------------------"                
            print "('There are ', i0, ' active particles')", nactive
            print "('ref_cell not found = ', i0, ' / ', i0)", count(mv_pdata%particles(:)%ref_cell == REFCELL_OUTBOUNDS), mv_pdata%N_part
        end if

        call mv_field_updater%update_field(ncyc)
        
        if ( allocated(this%callback_) ) then
            if ( mod(ncyc, this%ncall_) == 0 ) then
                call this%callback_%process_callback(mv_pdata)
            end if
        end if

        call this%motion_%proceed_time_step(ncyc)

        if ( mod(ncyc, this%nback_) == 0 ) then
            call write_out_backup(this%basename_back_, ncyc, this%backup_ascii_)
        end if

        if ( mod(ncyc, this%nwrite_) == 0 ) then
            call writeout_vtk(this%basename_vtk_, ncyc, this%write_ascii_)
        end if

        if ( nactive == 0 ) then
            print "(A)", "simulator_m/simulator_t%run::INFO:: There are no active particles."
            exit
        end if

    end do

    print "(A)", "simulator_m/simulator_t%run::INFO:: Simulation finished"

end subroutine

! private function
character(1) function get_sep__() result(sep)
    !! get file separator. 
    !! The file separator is assumed to be “/” for MacOS and Linux, and “\” for Windows. 
    !! If the detection fails, the separator is interpreted as an underscore “_”. 
    character(256) path_env_
    integer i, istat
    call get_environment_variable("PATH", path_env_, status=istat)
    if ( istat > 2 ) then
        error stop "Environmental Error:: error for environment variable `PATH`="""//trim(path_env_)//""" "
    end if
    do i = 1, 256
        sep = path_env_(i:i)
        if ( sep == "/" .or. sep == "\") then
            return
        end if
    end do

    ! 上で判定できない場合. 
    sep = "_"
end function
end module simulator_m