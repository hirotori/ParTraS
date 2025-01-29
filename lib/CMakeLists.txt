add_library(partras STATIC
    AFDETImporter.f90
    BaseImporter.f90
    DropletMotion.f90
    DropletMotionLegacy.f90
    DumpFile.f90
    FieldUpdater.f90
    Geometry.f90
    Idlist.f90
    KindParameters.f90
    Motion.f90
    FlowField.f90
    ParticleData.f90
    Simulator.f90
    UnstructuredMesh.f90
    VTKImporter.f90
    Version.f90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    target_compile_options(partras PRIVATE
        $<$<CONFIG:Debug>:-O0 -warn all -check all -traceback -fpp>
        $<$<CONFIG:Release>:-O3 -xHOST -finline -fpp>
        $<$<CONFIG:RelWithDevInfo>:-O0 -fpp -p>
    )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    target_compile_options(partras PRIVATE
        $<$<CONFIG:Debug>: -O0 -Wall -fcheck=all -fbacktrace -cpp>
        $<$<CONFIG:Release>: -O3 -march=native -cpp>
        $<$<CONFIG:RelWithDevInfo>:-O0 -cpp -pg>
    )
endif()

if(OpenMP_Fortran_FOUND)
    target_compile_options(partras PRIVATE 
        "${OpenMP_Fortran_FLAGS}"
    )

    target_link_libraries(partras
        OpenMP::OpenMP_Fortran 
    )

else()
    # OpenMPをサポートしない場合.

endif()