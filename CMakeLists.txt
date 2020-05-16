cmake_minimum_required(VERSION 3.15)
project(psfun VERSION 1.0 LANGUAGES Fortran C HOMEPAGE_URL https://github.com/Cirdans-Home/psfun)

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


# ---------------------------------------------------------------------------- #
# Now we can start the compilation of all the bit and pieces (hope for the best)
# ---------------------------------------------------------------------------- #
add_subdirectory(krylov)

# --------------------------------------------------------------------------- #
# Sphinx autodoc
# --------------------------------------------------------------------------- #
add_subdirectory(docs)
