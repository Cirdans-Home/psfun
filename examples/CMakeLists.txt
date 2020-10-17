# CMakeLists for the Examples folder

IF(BUILD_EXAMPLES)
  message("Compile and Link Examples")
  add_subdirectory(serial)
  add_subdirectory(polynomial)
ENDIF()
