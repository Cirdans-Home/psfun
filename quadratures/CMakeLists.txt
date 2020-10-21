# CMakeList for the Quadrature Module

message("Compile and Link Quadrature Library")

set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

set( QUADRATURE_SUBMODULES  "psfun_d_quadrature_mod.f90"
                            "psfun_z_quadrature_mod.f90"
                            "psfun_z_quadrules_mod.f90"
                            "psfun_base_quadrature_mod.f90")

add_library(psfun_quadrature psfun_quadrature_mod.f90 ${QUADRATURE_SUBMODULES})
add_compile_definitions(${PSBFDEFINES})
target_include_directories(psfun_quadrature PUBLIC ${PSBLAS_MODULES} ${AUXILINCLUDE})
target_include_directories(psfun_quadrature PUBLIC ${PSBLAS_INCLUDE})
target_include_directories(psfun_krylov PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
target_include_directories(psfun_krylov PUBLIC ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
target_link_libraries(psfun_quadrature PUBLIC ${LINKED_LIBRARIES} ${LINK_PSBLAS} ${AUXILLINK} -lpsfun_utils)

if(MPI_COMPILE_FLAGS)
  set_target_properties(psfun_quadrature PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(psfun_quadrature PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
