# CMakeLists for the Examples/serial folder

message("Compile and Link Polynomial Krylov Example")


add_executable(arnolditest arnolditest.f90)
# We need to compile against PSBLAS
add_compile_definitions(${PSBFDEFINES})
target_include_directories(arnolditest PUBLIC ${PSBLAS_MODULES})
target_include_directories(arnolditest PUBLIC ${PSBLAS_INCLUDE})
# We need to compile against the SERIAL and KRYLOV part of the library
target_include_directories(serialtest PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
target_include_directories(serialtest PUBLIC ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
target_link_libraries(arnolditest PUBLIC ${LINKED_LIBRARIES} ${LINK_PSBLAS})
target_link_libraries(arnolditest PUBLIC ${LINKED_LIBRARIES} ${LINK_PSBLAS} psfun_serial psfun_krylov)

add_custom_command(
        TARGET arnolditest POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/examples/polynomial/arnoldi.inp
                ${CMAKE_CURRENT_BINARY_DIR}/arnoldi.inp)
add_custom_command(
        TARGET arnolditest POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/examples/polynomial/laplacian.mtx
                ${CMAKE_CURRENT_BINARY_DIR}/laplacian.mtx)
