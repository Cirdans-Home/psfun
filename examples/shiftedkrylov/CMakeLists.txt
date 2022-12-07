# CMakeLists for the Examples/serial folder

message("Compile and Link Shifted Krylov method Examples")


add_executable(shiftedtest shiftedtest.f90 getp.f90)
# We need to compile against PSBLAS
add_compile_definitions(${PSBFDEFINES})
target_include_directories(shiftedtest PUBLIC ${PSBLAS_MODULES})
target_include_directories(shiftedtest PUBLIC ${PSBLAS_INCLUDE})
# We need to compile against the SERIAL and KRYLOV part of the library
target_include_directories(shiftedtest PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
target_include_directories(shiftedtest PUBLIC ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
target_link_libraries(shiftedtest PUBLIC ${LINKED_LIBRARIES} ${LINK_AMG4PSBLAS}
        ${LINK_PSBLAS} ${LINKED_LIBRARIES_AMG4PSBLAS}
        psfun_serial psfun_krylov psfun_utils)

add_custom_command(
        TARGET shiftedtest POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/examples/shiftedkrylov/shifted.inp
                ${CMAKE_CURRENT_BINARY_DIR}/shifted.inp)
add_custom_command(
        TARGET shiftedtest POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/examples/polynomial/laplacian.mtx
                ${CMAKE_CURRENT_BINARY_DIR}/laplacian.mtx)
