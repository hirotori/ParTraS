module simulation_m
    use iso_c_binding
    use util_m
    use simulator_m

    contains

    subroutine initialize_simulation(nwrite, write_ascii, write_path, nback, backup_ascii, back_path) &
        bind(c, name="initialize_simulation")
        integer(c_int),intent(in) :: nwrite
        logical(c_bool),intent(in) :: write_ascii
        character(1,kind=c_char),dimension(*),intent(in) :: write_path
        integer(c_int),intent(in) :: nback
        logical(c_bool),intent(in) :: backup_ascii
        character(1,kind=c_char),dimension(*),intent(in) :: back_path
        
        call mv_simulator%construct_simulator(nwrite, logical(write_ascii), fstring(write_path), &
                                              nback, logical(backup_ascii), fstring(back_path))

    end subroutine

    subroutine run_simulation(nstart, nend) bind(c, name="run_simulation")
        integer(c_int),intent(in) :: nstart
        integer(c_int),intent(in) :: nend

        call mv_simulator%run(nstart, nend)

    end subroutine

    subroutine set_writeout_settings(nwrite, write_ascii, write_path) bind(c, name="set_writeout_settings")
        !! construct simulator object
        integer(c_int),intent(in) :: nwrite
        logical(c_bool),intent(in) :: write_ascii
        character(1,kind=c_char),dimension(*),intent(in) :: write_path

        call mv_simulator%set_writeout_settings(nwrite, logical(write_ascii), fstring(write_path))
    
    end subroutine
    
    subroutine set_dump_settings(nback, backup_ascii, back_path) bind(c, name="set_dump_settings")
        !! construct simulator object
        integer(c_int),intent(in) :: nback
        logical(c_bool),intent(in) :: backup_ascii
        character(1,kind=c_char),dimension(*),intent(in) :: back_path

        call mv_simulator%set_dump_settings(nback, logical(backup_ascii), fstring(back_path))
    
    end subroutine

end module simulation_m