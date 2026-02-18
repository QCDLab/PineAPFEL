cmake_minimum_required(VERSION 3.14)
project(pineapfel VERSION 0.1.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Dependencies
find_package(PkgConfig REQUIRED)
pkg_check_modules(PINEAPPL REQUIRED pineappl_capi)
pkg_check_modules(YAML_CPP REQUIRED yaml-cpp)

# APFEL++ via apfelxx-config
find_program(APFELXX_CONFIG apfelxx-config REQUIRED)
execute_process(COMMAND ${APFELXX_CONFIG} --cflags OUTPUT_VARIABLE APFELXX_CFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE RESULT_VARIABLE _apfel_cflags_rc)
execute_process(COMMAND ${APFELXX_CONFIG} --ldflags OUTPUT_VARIABLE APFELXX_LDFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE RESULT_VARIABLE _apfel_ldflags_rc)
if(NOT _apfel_cflags_rc EQUAL 0 OR NOT _apfel_ldflags_rc EQUAL 0)
    message(FATAL_ERROR "apfelxx-config failed (cflags rc=${_apfel_cflags_rc}, ldflags rc=${_apfel_ldflags_rc})")
endif()
message(STATUS "APFEL++ cflags: ${APFELXX_CFLAGS}")
message(STATUS "APFEL++ ldflags: ${APFELXX_LDFLAGS}")

# Parse APFEL++ flags
separate_arguments(APFELXX_CFLAGS_LIST UNIX_COMMAND "${APFELXX_CFLAGS}")
separate_arguments(APFELXX_LDFLAGS_LIST UNIX_COMMAND "${APFELXX_LDFLAGS}")

# Extract APFEL++ include dirs from -I flags
set(APFELXX_INCLUDE_DIRS "")
foreach(flag ${APFELXX_CFLAGS_LIST})
    if(flag MATCHES "^-I(.*)")
        list(APPEND APFELXX_INCLUDE_DIRS "${CMAKE_MATCH_1}")
    endif()
endforeach()

# Fallback: if cflags produced no -I flags, derive include dir from apfelxx-config location
if(NOT APFELXX_INCLUDE_DIRS)
    get_filename_component(_apfelxx_bindir "${APFELXX_CONFIG}" DIRECTORY)
    get_filename_component(_apfelxx_prefix "${_apfelxx_bindir}" DIRECTORY)
    set(APFELXX_INCLUDE_DIRS "${_apfelxx_prefix}/include")
endif()
message(STATUS "APFEL++ include dirs: ${APFELXX_INCLUDE_DIRS}")

# Library
add_library(pineapfel
    pineapfel/src/cards.cpp
    pineapfel/src/evolution.cpp
    pineapfel/src/fill.cpp
    pineapfel/src/grid_gen.cpp
)
target_include_directories(pineapfel PUBLIC
    pineapfel/include
    ${PINEAPPL_INCLUDE_DIRS}
    ${YAML_CPP_INCLUDE_DIRS}
    ${APFELXX_INCLUDE_DIRS}
)
target_link_directories(pineapfel PUBLIC
    ${PINEAPPL_LIBRARY_DIRS}
    ${YAML_CPP_LIBRARY_DIRS}
)
target_link_libraries(pineapfel PUBLIC
    ${PINEAPPL_LIBRARIES}
    ${YAML_CPP_LIBRARIES}
    ${APFELXX_LDFLAGS_LIST}
)
target_compile_options(pineapfel PRIVATE -Wall -Wextra -O3)

# CLI executables
add_executable(pineapfel-evolve pineapfel_cli/main.cpp)
target_link_libraries(pineapfel-evolve PRIVATE pineapfel)
target_compile_options(pineapfel-evolve PRIVATE -Wall -Wextra -O3)

add_executable(pineapfel-build pineapfel_cli/build.cpp)
target_link_libraries(pineapfel-build PRIVATE pineapfel)
target_compile_options(pineapfel-build PRIVATE -Wall -Wextra -O3)

# Install targets
install(TARGETS pineapfel pineapfel-evolve pineapfel-build
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)
install(DIRECTORY pineapfel/include/ DESTINATION include)
