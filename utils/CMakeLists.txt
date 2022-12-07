# CMakeList for the Utils Module

message("Compile and Link Utils Library")

set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

set( UTILS_SUBMODULES "psfun_d_utils_mod.f90" "psfun_z_utils_mod.f90")
set( EXTERNAL_LIBRARIES "libraries/elliptic_integral.f90" "libraries/ogpf.f90")


add_library(psfun_utils psfun_utils_mod.f90 ${UTILS_SUBMODULES} ${EXTERNAL_LIBRARIES})
add_compile_definitions(${PSBFDEFINES})
target_include_directories(psfun_utils PUBLIC ${AMG4PSBLAS_MODULES} ${PSBLAS_MODULES} ${AUXILINCLUDE})
target_include_directories(psfun_utils PUBLIC ${AMG4PSBLAS_INCLUDE} ${PSBLAS_INCLUDE})
target_include_directories(psfun_utils PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
target_include_directories(psfun_utils PUBLIC ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
target_link_libraries(psfun_utils PUBLIC ${LINKED_LIBRARIES} ${LINK_AMG4PSBLAS} ${LINK_PSBLAS} ${AUXILLINK})

if(MPI_COMPILE_FLAGS)
  set_target_properties(psfun_utils PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(psfun_utils PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
