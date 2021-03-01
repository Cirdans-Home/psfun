# We look for a Python3 environment
find_package(Python COMPONENTS Interpreter Development)
# We check if FORD is available to build the Documentation
execute_process(
      COMMAND pip show ford
      RESULT_VARIABLE EXIT_CODE
      OUTPUT_QUIET
)
IF(NOT ${EXIT_CODE} EQUAL 0)
  message(
          WARNING
          "The \"FORD\" Python3 package is not installed. Please install it using the following command: \"pip3 install ford\"."
  )
  SET(BUILD_DOCUMENTATION FALSE)
ELSE()
  # We have found FORD so we can build the documentation
  configure_file(
  "${CMAKE_SOURCE_DIR}/docs/forddocs.md.in"
  "${CMAKE_BINARY_DIR}/forddocs.md"
  @ONLY)

  add_custom_target(psfun_html_docs ALL
    ford ${CMAKE_BINARY_DIR}/forddocs.md
    )

ENDIF()
