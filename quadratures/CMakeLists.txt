# CMakeList for the Quadrature Module

message("Compile and Link Utils Library")

set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

add_library(psfun_quadrature psfun_d_quadrature_mod.f90)
add_compile_definitions(${PSBFDEFINES})
target_include_directories(psfun_quadrature PUBLIC ${PSBLAS_MODULES} ${AUXILINCLUDE})
target_include_directories(psfun_quadrature PUBLIC ${PSBLAS_INCLUDE})
target_link_libraries(psfun_quadrature PUBLIC ${LINKED_LIBRARIES} ${LINK_PSBLAS} ${AUXILLINK})

if(MPI_COMPILE_FLAGS)
  set_target_properties(psfun_quadrature PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(psfun_quadrature PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
