program test
    use kind_parameters_m
    use particle_data_m
    implicit none

    type(particle_data_t) pdata
    
    integer(IP) :: N

    N = 10
    pdata = particle_data(N)


end program test