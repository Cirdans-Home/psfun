find_library(PHIFUNCTION_LIBRARY1 libmtrcfgphi.a PATHS ${PHIFUNCTION_DIR}/lib)
find_library(PHIFUNCTION_LIBRARY2 libxlacon.a PATHS ${PHIFUNCTION_DIR}/lib)

if(PHIFUNCTION_LIBRARY1 AND PHIFUNCTION_LIBRARY2)
  SET(LINKPHI -L${PHIFUNCTION_DIR}/lib -lmtrcfgphi -lxlacon)
  SET(INCLUDEPHI ${PHIFUNCTION_DIR}/modules/)
  SET(PHIDEFINES WITHPHILIBRARY)
else()
  message("I did not find PHIFUNCTION! Everything will fail")
endif()