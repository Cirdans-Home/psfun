# CMakeList for the Krylov module

# We first have to compile the outher module by compiling against all the
# submodules in the various folders, therefore we use some variable with local
# scope containing the files in the subdirectories

set(SERIAL_DIR "${CMAKE_BINARY_DIR}/serial")
set(POLYNOMIAL_METHODS polynomial/psfun_d_arnoldi.f90 polynomial/psfun_d_lanczos.f90 )

# Now we can compile the library with all its dependencies
add_library(psfun_krylov psfun_d_krylov_mod.f90 ${POLYNOMIAL_METHODS})
add_compile_definitions(${PSBFDEFINES})
target_include_directories(psfun_krylov PUBLIC ${PSBLAS_MODULES})
target_include_directories(psfun_krylov PUBLIC ${PSBLAS_INCLUDE})
target_include_directories(psfun_krylov PUBLIC ${SERIAL_DIR})
target_link_libraries(psfun_krylov PUBLIC ${LINKED_LIBRARIES} ${LINK_PSBLAS} psfun_serial)

if(MPI_COMPILE_FLAGS)
  set_target_properties(psfun_krylov PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(psfun_krylov PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
