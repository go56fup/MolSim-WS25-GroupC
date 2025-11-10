option(BUILD_DOCUMENTATION "Enable generation of Doxygen documentation" ON)
if(NOT BUILD_DOCUMENTATION)
	return()
endif()

find_package(Doxygen QUIET)
if(NOT Doxygen_FOUND)
	message(WARNING "Doxygen not found, documentation will not be built.")
	return()
endif()

function(generate_documentation_target)
	get_property(CMAKE_GENERATED_ALIASES GLOBAL PROPERTY __generated_aliases)
	configure_file(Doxyfile.in Doxyfile)
	doxygen_add_docs(doc_doxygen CONFIG_FILE "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile")
endfunction()

function(
	__build_variadic_list
	retval
	start
	end
	sep
)
	set(result "")
	foreach(i RANGE "${start}" "${end}")
		if(NOT i EQUAL "${start}")
			string(APPEND result "${sep}")
		endif()
		string(APPEND result "\\${i}")
	endforeach()
	set(
		"${retval}"
		"${result}"
		PARENT_SCOPE
	)
endfunction()

function(doxygen_alias name argc pattern)
	if(argc EQUAL 0)
		set(alias_arg_declaration "")
	else()
		set(alias_arg_declaration "{${argc}}")
	endif()
	set(alias "ALIASES += ${name}${alias_arg_declaration}=\"${pattern}\"")
	get_property(prev GLOBAL PROPERTY __generated_aliases)
	set_property(GLOBAL PROPERTY __generated_aliases "${prev}\n${alias}")
endfunction()

function(doxygen_generate_variadic)
	cmake_parse_arguments(
		ARG
		""
		"NAME;NORMAL_ARGC;VARIADIC_UP_TO;PATTERN;SEPARATOR"
		""
		${ARGV}
	)
	if(NOT ARG_SEPARATOR)
		set(ARG_SEPARATOR "\\,")
	endif()
	set(result "")
	math(EXPR args_begin "${ARG_NORMAL_ARGC} + 1")
	math(EXPR args_end "${ARG_VARIADIC_UP_TO} + ${ARG_NORMAL_ARGC}")

	foreach(i RANGE ${args_begin} ${args_end})
		__build_variadic_list(args "${args_begin}" "${i}" "${ARG_SEPARATOR}")
		string(REPLACE "**" "${args}" generated_alias "${ARG_PATTERN}")
		# Total arguments = non-variadic + iteration count Iteration count = current index + 1
		# Current index = current no - start no
		math(EXPR total_argc "${ARG_NORMAL_ARGC} + ${i} - ${args_begin} + 1")
		doxygen_alias("${ARG_NAME}" "${total_argc}" "${generated_alias}")
	endforeach()
endfunction()

function(define_param_shorthands name argc)
	set(shorthands str val)
	if(argc LESS 1)
		set(args "")
	else()
		__build_variadic_list(args 1 "${argc}" ",")
		string(APPEND args ",")
	endif()
	foreach(shorthand IN LISTS shorthands)
		doxygen_alias(
			"${name}_${shorthand}" "${argc}" "@${name}{${args}${shorthand}1,${shorthand}2}"
		)
	endforeach()
endfunction()

function(define_whole_section_shorthand name arg_type)
	set(argc 2)
	__build_variadic_list(args 1 "${argc}" ",")
	doxygen_alias("${name}_with_params" "${argc}" "@${name}{${args}}^^@gtest_params_${arg_type}")
endfunction()
