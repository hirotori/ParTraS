module version_m
    !! treat version of this project
    implicit none
    integer,parameter,private :: MAJOR_VERSION = 0
    integer,parameter,private :: MINOR_VERSION = 3
    integer,parameter,private :: PATCH_VERSION = 0

    public :: get_version

contains
    
    character(5) function get_version()
        !! get current version of this project
        write(get_version, "(3(i0,:,'.'))") MAJOR_VERSION, MINOR_VERSION, PATCH_VERSION

    end function

end module version_m