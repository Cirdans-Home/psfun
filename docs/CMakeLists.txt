# --------------------------------------------------------------------------- #
# Sphinx autodoc
# --------------------------------------------------------------------------- #
find_package(Sphinx REQUIRED)
if(NOT DEFINED SPHINX_THEME)
    set(SPHINX_THEME default)
endif()

if(NOT DEFINED SPHINX_THEME_DIR)
    set(SPHINX_THEME_DIR)
endif()
# Look for PdfLatex
find_package(LATEX COMPONENTS PDFLATEX)
set (LATEX_PDFLATEX_FOUND "${LATEX_PDFLATEX_FOUND}" PARENT_SCOPE)

# configured documentation tools and intermediate build results
set(BINARY_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_build")

# Sphinx cache with pickled ReST documents
set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_static")

# HTML output directory
set(SPHINX_HTML_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")

# Latex output directory
set(SPHINX_LATEX_DIR "${CMAKE_CURRENT_BINARY_DIR}/latex")

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in"
    "${BINARY_BUILD_DIR}/conf.py"
    @ONLY)

add_custom_target(psfun_html_docs ALL
    ${SPHINX_EXECUTABLE}
        -q -b html
        -c "${BINARY_BUILD_DIR}"
        -d "${SPHINX_CACHE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}"
        "${SPHINX_HTML_DIR}"
    COMMENT "Building HTML documentation with Sphinx")

if (LATEX_PDFLATEX_FOUND)
	add_custom_target(psfun_pdf_docs ALL
	    ${SPHINX_EXECUTABLE}
	        -c "${BINARY_BUILD_DIR}"
			-b latex
	        -d "${SPHINX_CACHE_DIR}"
	        "${CMAKE_CURRENT_SOURCE_DIR}"
	        "${SPHINX_LATEX_DIR}"
	    COMMENT "Building LATEX documentation with Sphinx")
	# We then compile the latex documentation
	add_custom_command(
					TARGET psfun_pdf_docs POST_BUILD
        	COMMAND make
        	WORKING_DIRECTORY ${SPHINX_LATEX_DIR}
	)
	add_custom_command(
	        TARGET psfun_pdf_docs POST_BUILD
	        COMMAND ${CMAKE_COMMAND} -E copy
	                ${SPHINX_LATEX_DIR}/psfun.pdf
	                ${CMAKE_SOURCE_DIR}/docs/psfun_guide.pdf)
else()
	message(WARNING "PDF generation tools are missing. PDF document generation target is not created.")
endif()
