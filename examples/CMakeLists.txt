# CMakeLists for the Examples folder

IF(BUILD_EXAMPLES)
  message("Compile and Link Examples")
  add_subdirectory(environment)
  add_subdirectory(utils)
  add_subdirectory(serial)
  add_subdirectory(polynomial)
  add_subdirectory(quadrature)
  add_subdirectory(shiftedkrylov)
ENDIF()
