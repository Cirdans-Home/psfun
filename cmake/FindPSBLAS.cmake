find_library(PSBLAS_BASE_LIBRARY libpsb_base.a PATHS ${PSBLAS_DIR}/lib REQUIRED)
find_library(PSBLAS_CBIND_LIBRARY libpsb_cbind.a PATHS ${PSBLAS_DIR}/lib REQUIRED)
find_library(PSBLAS_UTIL_LIBRARY libpsb_util.a PATHS ${PSBLAS_DIR}/lib REQUIRED)
find_library(PSBLAS_KRYLOV_LIBRARY libpsb_krylov.a PATHS ${PSBLAS_DIR}/lib REQUIRED)
find_library(PSBLAS_PREC_LIBRARY libpsb_prec.a PATHS ${PSBLAS_DIR}/lib REQUIRED)

# We read the Make.inc.psblas file to get the compilation option in PSBLAS, when
# PSBLAS will be modified to be installed by means of CMake this part will
# become fully automatized and less REGEX-terryfing
if(PSBLAS_BASE_LIBRARY AND PSBLAS_CBIND_LIBRARY AND PSBLAS_UTIL_LIBRARY AND PSBLAS_KRYLOV_LIBRARY AND PSBLAS_PREC_LIBRARY)
	set(regex "BLAS=.*")
	file(STRINGS ${PSBLAS_DIR}/include/Make.inc.psblas LINK_BLAS REGEX "${regex}")
	set(regex "BLAS=")
	string(REGEX REPLACE "${regex}" "" LINK_BLAS "${LINK_BLAS}")

	set(regex "METIS_LIB=.*")
	file(STRINGS ${PSBLAS_DIR}/include/Make.inc.psblas LINK_METIS_LIB REGEX "${regex}")
	set(regex "METIS_LIB=")
	string(REGEX REPLACE "${regex}" "" LINK_METIS_LIB "${LINK_METIS_LIB}")

	set(regex "AMD_LIB=.*")
	file(STRINGS ${PSBLAS_DIR}/include/Make.inc.psblas LINK_AMD_LIB REGEX "${regex}")
	set(regex "AMD_LIB=")
	string(REGEX REPLACE "${regex}" "" LINK_AMD_LIB "${LINK_AMD_LIB}")

	set(regex "PSBFDEFINES=.*")
	file(STRINGS ${PSBLAS_DIR}/include/Make.inc.psblas PSBFDEFINES REGEX "${regex}")
	set(regex "PSBFDEFINES=")
	string(REGEX REPLACE "${regex}" "" PSBFDEFINES "${PSBFDEFINES}")
	set(regex "-D")
	string(REGEX REPLACE "${regex}" "" PSBFDEFINES "${PSBFDEFINES}")
	separate_arguments(PSBFDEFINES)

	set(regex "PSBCDEFINES=.*")
	file(STRINGS ${PSBLAS_DIR}/include/Make.inc.psblas PSBCDEFINES REGEX "${regex}")
	set(regex "PSBCDEFINES=")
	string(REGEX REPLACE "${regex}" "" PSBCDEFINES "${PSBCDEFINES}")
	set(regex "-D")
	string(REGEX REPLACE "${regex}" "" PSBCDEFINES "${PSBFDEFINES}")
	separate_arguments(PSBCDEFINES)

	set(LINK_PSBLAS -L${PSBLAS_DIR}/lib -lpsb_cbind -lpsb_prec -lpsb_krylov -lpsb_util -lpsb_base)
	set(LINKED_LIBRARIES ${LINK_BLAS} ${LINK_METIS_LIB} ${LINK_AMD_LIB})
	set(PSBLAS_INCLUDE ${PSBLAS_DIR}/include/)
	set(PSBLAS_MODULES ${PSBLAS_DIR}/modules/)
else()
	message("I did not find PSBLAS! Everything will fail")
endif()
#------------------------------------------------------------------------------#
