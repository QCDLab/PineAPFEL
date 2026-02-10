cmake_minimum_required(VERSION 3.14)
project(pineapfel VERSION 0.1.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Dependencies
find_package(PkgConfig REQUIRED)
pkg_check_modules(PINEAPPL REQUIRED pineappl_capi)
pkg_check_modules(YAML_CPP REQUIRED yaml-cpp)

# APFEL++ via apfelxx-config
execute_process(COMMAND apfelxx-config --cflags OUTPUT_VARIABLE APFELXX_CFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND apfelxx-config --ldflags OUTPUT_VARIABLE APFELXX_LDFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)

# Parse APFEL++ flags
separate_arguments(APFELXX_CFLAGS_LIST UNIX_COMMAND "${APFELXX_CFLAGS}")
separate_arguments(APFELXX_LDFLAGS_LIST UNIX_COMMAND "${APFELXX_LDFLAGS}")

# Extract APFEL++ include dirs and libraries
set(APFELXX_INCLUDE_DIRS "")
set(APFELXX_LINK_FLAGS "")
foreach(flag ${APFELXX_CFLAGS_LIST})
    if(flag MATCHES "^-I(.*)")
        list(APPEND APFELXX_INCLUDE_DIRS "${CMAKE_MATCH_1}")
    endif()
endforeach()

# Library
add_library(pineapfel
    src/cards.cpp
    src/evolution.cpp
)
target_include_directories(pineapfel PUBLIC
    include
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

# CLI executable
add_executable(pineapfel-evolve cli/main.cpp)
target_link_libraries(pineapfel-evolve PRIVATE pineapfel)
target_compile_options(pineapfel-evolve PRIVATE -Wall -Wextra -O3)

# Install targets
install(TARGETS pineapfel pineapfel-evolve
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)
install(DIRECTORY include/pineapfel DESTINATION include)
