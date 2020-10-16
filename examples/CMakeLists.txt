# CMakeLists for the Examples folder

IF(BUILD_EXAMPLES)
  message("Compile and Link Examples")
  add_subdirectory(serial)
ENDIF()
