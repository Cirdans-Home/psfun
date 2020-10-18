find_program(SPHINX_EXECUTABLE
	NAMES sphinx-build
	DOC "Sphinx Documentation Builder (sphinx-doc.org)"
	HINTS $ENV{SPHINX_DIR}
	PATH ${SPHINX_PATH}
	PATH_SUFFIXES bin
	# NO_DEFAULT_PATH
	# NO_CMAKE_ENVIRONMENT_PATH
	# NO_CMAKE_PATH
	# NO_SYSTEM_ENVIRONMENT_PATH
	# NO_CMAKE_SYSTEM_PATH
	# NO_CMAKE_FIND_ROOT_PATH
)
if (SPHINX_EXECUTABLE-NOTFOUND)
	message(STATUS "Failed to find sphinx-build at ${SPHINX_PATH}.")
	message(STATUS "Sphinx documentation targets will not be created.")
	return()
endif()

#Get Sphinx version
execute_process(COMMAND "${SPHINX_EXECUTABLE}" "--version" OUTPUT_VARIABLE
	_SPHINX_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
#Parse output
if(_SPHINX_VERSION)
	if(_SPHINX_VERSION MATCHES ".*sphinx-build[^0-9.]*([0-9.]+).*")
		string(REGEX REPLACE ".*sphinx-build ([0-9.]+).*" "\\1" SPHINX_VERSION "${_SPHINX_VERSION}")
	endif()
endif()

#Set "standard" find module return values
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Sphinx REQUIRED_VARS SPHINX_EXECUTABLE
	SPHINX_VERSION VERSION_VAR SPHINX_VERSION)

mark_as_advanced(SPHINX_EXECUTABLE)
