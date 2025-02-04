function(ADDTEST name)
    ### テスト用プログラムのビルド設定
    add_executable(test_${name} test_${name}.f90)

    ### 参照するmodファイルのディレクトリとリンクするライブラリの設定
    target_link_libraries(test_${name} partras)
    set_target_properties(test_${name}
        PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
    )
    if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
        target_compile_options(test_${name} PRIVATE
            $<$<CONFIG:Debug>:-O0 -warn all -check all -traceback -fpp>
            $<$<CONFIG:Release>:-O3 -fPIC -finline -fpp>
            $<$<CONFIG:RelWithDevInfo>:-O2 -fpp -p>
        )
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
        target_compile_options(test_${name} PRIVATE
            $<$<CONFIG:Debug>: -O0 -Wall -fcheck=all -fbacktrace -cpp>
            $<$<CONFIG:Release>: -O3 -cpp>
            $<$<CONFIG:RelWithDevInfo>:-O2 -cpp -pg>
        )
    endif()

    ### テストの登録
    add_test(NAME ${name}
             COMMAND $<TARGET_FILE:test_${name}> ${CMAKE_CURRENT_BINARY_DIR}
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endfunction(ADDTEST)

# add_subdirectory(integrated)

ADDTEST(umesh)
ADDTEST(umesh_vtk)
ADDTEST(mesh_field)
ADDTEST(mesh_field_vtk)
ADDTEST(vtk_importer)
ADDTEST(afdet_importer)
ADDTEST(refcell_search)
# ADDTEST(field_updater)

if(OpenMP_Fortran_FOUND)
    ADDTEST(refcell_search_parallel)
endif()