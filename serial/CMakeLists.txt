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

# For the computation of the ϕ functions we use the code from
# Koikari, Souji. "Algorithm 894: On a block Schur--Parlett algorithm
# for ϕ-functions based on the sep-inverse estimate." ACM Transactions
# on Mathematical Software (TOMS) 36.2 (2009): 1-20.


# Now we can compile the library with all its dependencies
set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
add_library(psfun_serial psfun_d_serial_mod.F90 ${SERIAL_SUBMODULES})
add_compile_definitions(${PSBFDEFINES})
target_include_directories(psfun_serial PUBLIC  ${AMG4PSBLAS_MODULES} ${PSBLAS_MODULES} ${AUXILINCLUDE})
target_include_directories(psfun_serial PUBLIC  ${AMG4PSBLAS_INCLUDE} ${PSBLAS_INCLUDE})
target_link_libraries(psfun_serial PUBLIC ${LINKED_LIBRARIES} ${LINK_AMG4PSBLAS} ${LINK_PSBLAS} ${AUXILLINK})

if(MPI_COMPILE_FLAGS)
  set_target_properties(psfun_serial PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(psfun_serial PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
