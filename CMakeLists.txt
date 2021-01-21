cmake_minimum_required(VERSION 3.15)
project(psfun VERSION 0.1 LANGUAGES Fortran C HOMEPAGE_URL https://github.com/Cirdans-Home/psfun)

# Directory with CMAKE modules, used to search for libraries
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# This is a parallel code! We need MPI and stuff:
find_package(MPI REQUIRED)
if(NOT MPI_Fortran_FOUND)
    message(FATAL_ERROR "Could not find Fortran MPI.
			Please set MPI_Fortran_COMPILER to point to the mpifort wrapper.")
endif()
IF(MPI_Fortran_COMPILER)
  # use MPI wrapper as the compiler
  SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ELSE()
  # add MPI_INCLUDE_PATH to include directories
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF()

if(NOT MPI_C_FOUND)
    message(FATAL_ERROR "Could not find C MPI.
			Please set MPI_C_COMPILER to point to the mpicc wrapper.")
endif()
IF(MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  SET(CMAKE_C_COMPILER ${MPI_C_COMPILER})
ELSE()
  # add MPI_INCLUDE_PATH to include directories
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF()


include_directories(${MPI_Fortran_INCLUDE_PATH})
add_compile_options(${MPI_Fortran_COMPILE_FLAGS})

# The user has to tell us where PSBLAS is installed : everyting is built on that
# if PSBLAS is not found all shall fail
set(PSBLAS_DIR /usr/lib CACHE PATH "/usr/lib")
find_package(PSBLAS REQUIRED)

#------------------------------------------------------------------------------#
# Options to be set by the user
#------------------------------------------------------------------------------#
option(BUILD_EXAMPLES "Build example programs" OFF)
option(EXTERNAL_PHI "Use Koikari codes for Phi Functiions" OFF)
option(BUILD_DOCUMENTATION "Build Sphinx Documentation" ON)
option(GNUPLOT "Enable GNUPLOT plots" OFF)

#------------------------------------------------------------------------------#
# Look for auxiliary libraries
#------------------------------------------------------------------------------#
IF(EXTERNAL_PHI)
	set(PHIFUNCTION_DIR /usr/lib CACHE PATH "/usr/lib")
	find_package(PHIFUNCTION)
ENDIF()
add_compile_definitions(${PHIDEFINES})
set(AUXILINCLUDE ${INCLUDEPHI})
set(AUXILLINK ${LINKPHI})

IF(GNUPLOT)
	find_package(Gnuplot)
ENDIF()
IF(GNUPLOT_FOUND)
	message("Gnuplot found")
	add_compile_definitions(WITHGNUPLOTFORTRAN)
ENDIF()

separate_arguments(AUXILINCLUDE SEPARATE_ARGS PROGRAM UNIX_COMMAND "${AUXILINCLUDE}")

#------------------------------------------------------------------------------#
# Set the folders that will contain modules, library, and include files in the
# building phase
#------------------------------------------------------------------------------#
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# ---------------------------------------------------------------------------- #
# Now we can start the compilation of all the bit and pieces (hope for the best)
# ---------------------------------------------------------------------------- #
add_subdirectory(utils)
add_subdirectory(serial)
add_subdirectory(krylov)
add_subdirectory(quadratures)
add_subdirectory(examples)

# --------------------------------------------------------------------------- #
# Sphinx autodoc
# --------------------------------------------------------------------------- #
IF(BUILD_DOCUMENTATION)
	add_subdirectory(docs)
ENDIF()
