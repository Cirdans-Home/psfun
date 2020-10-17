# CMakeList for the Krylov module

message("Compile and Link Serial Library")

# We first have to compile the outher module by compiling against all the
# submodules in the various folders, therefore we use some variable with local
# scope containing the files in the subdirectories

set(SERIAL_SUBMODULES exponential/r8lib.f90
                      exponential/taylor.f90
                      exponential/dgpadm.f
                      exponential/dnchbv.f
                      exponential/dschbv.f
                      exponential/dgchbv.f)

# Now we can compile the library with all its dependencies
set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
add_library(psfun_serial psfun_d_serial_mod.f90 ${SERIAL_SUBMODULES})
add_compile_definitions(${PSBFDEFINES})
target_include_directories(psfun_serial PUBLIC ${PSBLAS_MODULES})
target_include_directories(psfun_serial PUBLIC ${PSBLAS_INCLUDE})
target_link_libraries(psfun_serial PUBLIC ${LINKED_LIBRARIES} ${LINK_PSBLAS})

if(MPI_COMPILE_FLAGS)
  set_target_properties(psfun_serial PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(psfun_serial PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
