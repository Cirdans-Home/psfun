# CMakeLists for the Examples/serial folder

message("Compile and Link Quadrature Example")

add_executable(quadraturetest quadraturetest.F90)

# We need to compile against PSBLAS
add_compile_definitions(${PSBFDEFINES})
target_include_directories(quadraturetest PUBLIC ${PSBLAS_MODULES} ${AUXILINCLUDE})
target_include_directories(quadraturetest PUBLIC ${PSBLAS_INCLUDE})
# We need to compile against the UTIL and SERIAL part of the library
target_include_directories(serialtest PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
target_include_directories(serialtest PUBLIC ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
target_link_libraries(quadraturetest PUBLIC ${LINKED_LIBRARIES}
        ${LINK_PSBLAS} psfun_serial psfun_krylov psfun_utils psfun_quadrature)

add_custom_command(
        TARGET quadraturetest POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/examples/quadrature/quadrature.inp
                ${CMAKE_CURRENT_BINARY_DIR}/quadrature.inp)
