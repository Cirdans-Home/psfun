find_library(AMG4PSBLAS_PREC_LIBRARY libamg_prec.a PATHS ${AMG4PSBLAS_DIR}/lib REQUIRED)
find_library(AMG4PSBLAS_CBIND_LIBRARY libamg_cbind.a PATHS ${AMG4PSBLAS_DIR}/lib REQUIRED)

if(AMG4PSBLAS_PREC_LIBRARY AND AMG4PSBLAS_CBIND_LIBRARY)
	set(regex "MUMPSLIBS=.*")
	file(STRINGS ${AMG4PSBLAS_DIR}/include/Make.inc.amg4psblas LINK_MUMPS_LIB REGEX "${regex}")
	set(regex ".*#;MUMPSLIBS=")
	string(REGEX REPLACE "${regex}" "" LINK_MUMPS_LIB "${LINK_MUMPS_LIB}")

	set(regex "MUMPSFLAGS=.*")
	file(STRINGS ${AMG4PSBLAS_DIR}/include/Make.inc.amg4psblas FLAGS_MUMPS_LIB REGEX "${regex}")
	set(regex ".*#;MUMPSFLAGS=")
	string(REGEX REPLACE "${regex}" "" FLAGS_MUMPS_LIB "${FLAGS_MUMPS_LIB}")

	set(regex "SLULIBS=.*")
	file(STRINGS ${AMG4PSBLAS_DIR}/include/Make.inc.amg4psblas LINK_SLU_LIB REGEX "${regex}")
	set(regex ".*#;SLULIBS=")
	string(REGEX REPLACE "${regex}" "" LINK_SLU_LIB "${LINK_SLU_LIB}")

	set(regex "SLUFLAGS=.*")
	file(STRINGS ${AMG4PSBLAS_DIR}/include/Make.inc.amg4psblas FLAGS_SLU_LIB REGEX "${regex}")
	set(regex ".*#;SLUFLAGS=")
	string(REGEX REPLACE "${regex}" "" FLAGS_SLU_LIB "${FLAGS_SLU_LIB}")

	set(regex "SLUDISTLIBS=.*")
	file(STRINGS ${AMG4PSBLAS_DIR}/include/Make.inc.amg4psblas LINK_SLUDIST_LIB REGEX "${regex}")
	set(regex ".*#;SLUDISTLIBS=")
	string(REGEX REPLACE "${regex}" "" LINK_SLUDIST_LIB "${LINK_SLUDIST_LIB}")

	set(regex "SLUDISTFLAGS=.*")
	file(STRINGS ${AMG4PSBLAS_DIR}/include/Make.inc.amg4psblas FLAGS_SLUDIST_LIB REGEX "${regex}")
	set(regex ".*#;SLUDISTFLAGS=")
	string(REGEX REPLACE "${regex}" "" FLAGS_SLUDIST_LIB "${FLAGS_SLUDIST_LIB}")

	set(regex "UMFLIBS=.*")
	file(STRINGS ${AMG4PSBLAS_DIR}/include/Make.inc.amg4psblas LINK_UMF_LIB REGEX "${regex}")
	set(regex ".*#;UMFLIBS=")
	string(REGEX REPLACE "${regex}" "" LINK_UMF_LIB "${LINK_UMF_LIB}")

	set(regex "UMFFLAGS=.*")
	file(STRINGS ${AMG4PSBLAS_DIR}/include/Make.inc.amg4psblas FLAGS_UMF_LIB REGEX "${regex}")
	set(regex ".*#;UMFFLAGS=")
	string(REGEX REPLACE "${regex}" "" FLAGS_UMF_LIB "${FLAGS_UMF_LIB}")

	set(regex "EXTRALIBS=.*")
	file(STRINGS ${AMG4PSBLAS_DIR}/include/Make.inc.amg4psblas LINK_EXTRA_LIB REGEX "${regex}")
	set(regex "EXTRALIBS=")
	string(REGEX REPLACE "${regex}" "" LINK_EXTRA_LIB "${LINK_EXTRA_LIB}")

	set(AMGCDEFINES "${FLAGS_MUMPS_LIB} ${FLAGS_SLU_LIB} ${FLAGS_SLUDIST_LIB} ${FLAGS_UMF_LIB}")
	set(LINK_AMG4PSBLAS "-L${AMG4PSBLAS_DIR}/lib -lamg_cbind -lamg_prec")
	set(LINKED_LIBRARIES_AMG4PSBLAS ${LINK_MUMPS_LIB} ${LINK_SLU_LIB} ${LINK_SLUDIST_LIB} ${LINK_UMF_LIB} ${LINK_EXTRA_LIB} -lstdc++)
	set(AMG4PSBLAS_INCLUDE ${AMG4PSBLAS_DIR}/include/)
	set(AMG4PSBLAS_MODULES ${AMG4PSBLAS_DIR}/modules/)

else()
	message("I did not find AMG4PSBLAS! Everything will fail")
endif()
