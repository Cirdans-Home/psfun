# CMakeLists for the Examples/serial folder

message("Compile and Link Serial Example")

add_executable(serialtest serialtest.f90)
# We need to compile against PSBLAS
add_compile_definitions(${PSBFDEFINES})
target_include_directories(serialtest PUBLIC ${PSBLAS_MODULES})
target_include_directories(serialtest PUBLIC ${PSBLAS_INCLUDE})
# We need to compile against the SERIAL part of the library
target_include_directories(serialtest PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
target_include_directories(serialtest PUBLIC ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
target_link_libraries(serialtest PUBLIC ${LINKED_LIBRARIES} ${LINK_PSBLAS})
target_link_libraries(serialtest PUBLIC ${LINKED_LIBRARIES} ${LINK_PSBLAS} psfun_serial)

add_custom_command(
        TARGET serialtest POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/examples/serial/serial.inp
                ${CMAKE_CURRENT_BINARY_DIR}/serial.inp)
add_custom_command(
        TARGET serialtest POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/examples/serial/laplacian.mtx
                ${CMAKE_CURRENT_BINARY_DIR}/laplacian.mtx)
