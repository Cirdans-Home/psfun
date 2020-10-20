

message("Compile and Link Util Example")

add_executable(utiltest utiltest.f90)

# We need to compile against PSBLAS
add_compile_definitions(${PSBFDEFINES})
target_include_directories(utiltest PUBLIC ${PSBLAS_MODULES} ${AUXILINCLUDE})
target_include_directories(utiltest PUBLIC ${PSBLAS_INCLUDE})
# We need to compile against the UTILS part of the library
target_include_directories(utiltest PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
target_include_directories(utiltest PUBLIC ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
target_link_libraries(utiltest PUBLIC ${LINKED_LIBRARIES} ${LINK_PSBLAS} ${AUXILLINK} psfun_utils)
