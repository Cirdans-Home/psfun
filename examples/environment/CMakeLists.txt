message("Compile and Link Polynomial Krylov Example")


add_executable(psblastest psblastest.f90)
# We need to compile against PSBLAS
add_compile_definitions(${PSBFDEFINES})
target_include_directories(psblastest PUBLIC ${PSBLAS_MODULES})
target_include_directories(psblastest PUBLIC ${PSBLAS_INCLUDE})
# We need to compile against the SERIAL and KRYLOV part of the library
target_include_directories(psblastest PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
target_include_directories(psblastest PUBLIC ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
target_link_libraries(psblastest PUBLIC ${LINKED_LIBRARIES} ${LINK_PSBLAS})
target_link_libraries(psblastest PUBLIC ${LINKED_LIBRARIES} ${LINK_PSBLAS})
