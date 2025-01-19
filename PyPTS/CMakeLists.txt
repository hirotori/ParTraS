add_library(pypts SHARED
    field.f90
    motions.f90
    particle.f90
    simulation.f90
    update.f90
    util.f90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    target_compile_options(pypts PRIVATE
        $<$<CONFIG:Debug>:-O0 -warn all -check all -traceback -fpp>
        $<$<CONFIG:Release>:-O3 -xHOST -finline -fpp>
        $<$<CONFIG:RelWithDevInfo>:-O0 -fpp -p>
    )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    target_compile_options(pypts PRIVATE
        $<$<CONFIG:Debug>: -O0 -Wall -fcheck=all -fbacktrace -cpp>
        $<$<CONFIG:Release>: -O3 -march=native -cpp>
        $<$<CONFIG:RelWithDevInfo>:-O0 -cpp -pg>
    )
endif()

if(OpenMP_Fortran_FOUND)
    target_compile_options(pypts PRIVATE 
        "${OpenMP_Fortran_FLAGS}"
    )

    target_link_libraries(pypts 
        OpenMP::OpenMP_Fortran 
    )

else()
    # OpenMPをサポートしない場合.

endif()

set_target_properties(pypts PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_NAME pypts
)

target_link_libraries(pypts PRIVATE partras)