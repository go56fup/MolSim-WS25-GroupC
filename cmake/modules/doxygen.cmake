option(BUILD_DOCUMENTATION "Enable generation of Doxygen documentation" ON)
if(NOT BUILD_DOCUMENTATION)
	return()
endif()

find_package(Doxygen QUIET)
if(NOT Doxygen_FOUND)
	message(WARNING "Doxygen not found, documentation will not be built.")
	return()
endif()

configure_file(Doxyfile.in Doxyfile)
doxygen_add_docs(doc_doxygen CONFIG_FILE "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile")
