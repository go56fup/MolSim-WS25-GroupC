option(BUILD_DOCUMENTATION "Enable generation of Doxygen documentation" ON)
if(NOT BUILD_DOCUMENTATION)
	return()
endif()

find_package(Doxygen QUIET)
if(NOT Doxygen_FOUND)
	message(WARNING "Doxygen not found, documentation will not be built.")
	return()
endif()

function(add_doxygen_support)
	configure_file(Doxyfile.in Doxyfile)
	doxygen_add_docs(doc_doxygen CONFIG_FILE "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile")
endfunction()

function(__build_variadic_list retval start end sep)
	set(result "")
	foreach(i RANGE "${start}" "${end}")
            if(NOT i EQUAL "${start}")
                string(APPEND result "${sep}")
            endif()
            string(APPEND result "\\${i}")
        endforeach()
		set("${retval}" "${result}" PARENT_SCOPE)
endfunction()

function(doxygen_alias name argc pattern)
	if(argc EQUAL 0)
		set(alias_arg_declaration "")
	else()
		set(alias_arg_declaration "{${argc}}")
	endif()
	set(alias "ALIASES += ${name}${alias_arg_declaration}=\"${pattern}\"")
	# message(WARNING "${alias}")
endfunction()


function(doxygen_generate_variadic)
    cmake_parse_arguments(ARG "" "NAME;NORMAL_ARGC;VARIADIC_UP_TO;PATTERN" "" ${ARGV})
    set(result "")
	math(EXPR args_begin "${ARG_NORMAL_ARGC} + 1")
	math(EXPR args_end "${ARG_VARIADIC_UP_TO} + ${ARG_NORMAL_ARGC}")
	
	foreach(i RANGE ${args_begin} ${args_end})
		__build_variadic_list(args "${args_begin}" "${i}" "\\,")
		string(REPLACE "**" "${args}" generated_alias "${ARG_PATTERN}")
		# Total arguments = non-variadic + iteration count
		# Iteration count = current index + 1
		# Current index = current no - start no
		math(EXPR total_argc "${ARG_NORMAL_ARGC} + ${i} - ${args_begin} + 1")
		doxygen_alias("${ARG_NAME}" "${total_argc}" "${generated_alias}")
	endforeach()
endfunction()

function(define_param_shorthands name argc)
	set(shorthands str val)
	# TODO: make build variadic list just start index, end index
	__build_variadic_list(args 1 "${argc}" ",")
	foreach(shorthand IN LISTS shorthands)
		doxygen_alias("${name}_${shorthand}" "${argc}" "${name}{${args},${shorthand}1,${shorthand}2}")
	endforeach()
endfunction()

function(define_whole_section_shorthand name argc argname_prefix)
	__build_variadic_list(args 1 "${argc}" ",")
	doxygen_alias("${name}_with_params" "${argc}" "${name}{${args}}")
endfunction()
