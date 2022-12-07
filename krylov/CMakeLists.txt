# CMakeList for the Krylov module

message("Compile and Link Krylov Library")

# We first have to compile the outher module by compiling against all the
# submodules in the various folders, therefore we use some variable with local
# scope containing the files in the subdirectories
set(POLYNOMIAL_METHODS polynomial/psfun_d_arnoldi.f90 polynomial/psfun_d_lanczos.f90 )
set(RATIONAL_METHODS rational/psfun_d_sai.f90)
set(SHIFTED_METHODS shifted/psfun_dcg.F90)

# Now we can compile the library with all its dependencies
set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
add_library(psfun_krylov psfun_d_krylov_mod.F90 psfun_krylov_mod.f90
  ${POLYNOMIAL_METHODS} ${RATIONAL_METHODS} ${SHIFTED_METHODS})
add_compile_definitions(${PSBFDEFINES})
target_include_directories(psfun_krylov PUBLIC ${AMG4PSBLAS_MODULES} ${PSBLAS_MODULES} ${AUXILINCLUDE})
target_include_directories(psfun_krylov PUBLIC ${AMG4PSBLAS_INCLUDE} ${PSBLAS_INCLUDE})
target_include_directories(psfun_krylov PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
target_include_directories(psfun_krylov PUBLIC ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
target_link_libraries(psfun_krylov PUBLIC ${LINKED_LIBRARIES} ${LINK_PSBLAS}
  ${AUXILLINK} ${LINK_AMG4PSBLAS} ${LINKED_LIBRARIES_AMG4PSBLAS}
  psfun_serial psfun_utils)

if(MPI_COMPILE_FLAGS)
  set_target_properties(psfun_krylov PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(psfun_krylov PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
