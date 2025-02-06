cmake_minimum_required(VERSION 3.0)
enable_language(Fortran)

project(partras Fortran)

option(use_OpenMP "enable OpenMP directives" OFF)

# OpenMPのconfigure
if(use_OpenMP)
    find_package(OpenMP REQUIRED)
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mod)

add_subdirectory(lib)

# for Python interface
# add_subdirectory(PyPTS)

# test
enable_testing()
add_subdirectory(test)

# システムへのインストール
set(PROJ_DIR "${CMAKE_INSTALL_PREFIX}/partras")
install(DIRECTORY DESTINATION ${PROJ_DIR})
install(TARGETS partras ARCHIVE DESTINATION ${PROJ_DIR}/lib)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION ${PROJ_DIR})
# src/メイン関数はプロジェクトに含めない.