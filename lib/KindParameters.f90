module kind_parameters_m
    !! kind parameter utility module
    use, intrinsic :: iso_fortran_env
    implicit none
    private

    integer,parameter,public :: I8 = int8
        !!parameter for 8-bit integer
    integer,parameter,public :: IP = int32
        !!parameter for 32-bit integer
    integer,parameter,public :: LIP = int64
        !!parameter for 64-bit integer
    integer,parameter,public :: DP = real64
        !!parameter for 64-bit real
    integer,parameter,public :: SP = real32
        !!parameter for 32-bit real

contains
    
end module kind_parameters_m