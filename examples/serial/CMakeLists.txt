# CMakeLists for the Examples/serial folder

message("Compile and Link Serial Example")

# We need to compile against the SERIAL part of the library
set(SERIAL_DIR "${CMAKE_BINARY_DIR}/serial")

add_executable(serialtest serialtest.f90)
add_compile_definitions(${PSBFDEFINES})
target_include_directories(serialtest PUBLIC ${PSBLAS_MODULES} ${SERIAL_DIR})
target_include_directories(serialtest PUBLIC ${PSBLAS_INCLUDE} ${SERIAL_DIR})
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
                ${CMAKE_CURRENT_BINARY_DIR}/bcsstk06.mtx)
