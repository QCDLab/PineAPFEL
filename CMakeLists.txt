cmake_minimum_required(VERSION 3.14)
project(pineapfel VERSION 0.1.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

# Dependencies
find_package(PkgConfig REQUIRED)
pkg_check_modules(PINEAPPL REQUIRED pineappl_capi)

find_package(yaml-cpp REQUIRED)
find_package(apfelxx REQUIRED)

# Library
add_library(pineapfel SHARED
    pineapfel/src/cards.cpp
    pineapfel/src/evolution.cpp
    pineapfel/src/fill.cpp
    pineapfel/src/grid_gen.cpp
    pineapfel/src/sidis_init.cpp
)

target_include_directories(pineapfel
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/pineapfel/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        ${PINEAPPL_INCLUDE_DIRS}
)

target_link_libraries(pineapfel
    PUBLIC
        yaml-cpp::yaml-cpp
        apfelxx::apfelxx
        ${PINEAPPL_LINK_LIBRARIES}
)

target_compile_options(pineapfel PRIVATE -Wall -Wextra -O3)

# CLI executables
add_executable(pineapfel-evolve pineapfel_cli/main.cpp)
target_link_libraries(pineapfel-evolve PRIVATE pineapfel)
target_compile_options(pineapfel-evolve PRIVATE -Wall -Wextra -O3)

add_executable(pineapfel-build pineapfel_cli/build.cpp)
target_link_libraries(pineapfel-build PRIVATE pineapfel)
target_compile_options(pineapfel-build PRIVATE -Wall -Wextra -O3)

# Tests
enable_testing()
add_executable(test_grid_vs_apfelxx tests/test_grid_vs_apfelxx.cpp tests/sidis_helper.cpp)
target_link_libraries(test_grid_vs_apfelxx PRIVATE pineapfel)
target_compile_options(test_grid_vs_apfelxx PRIVATE -Wall -Wextra -O3)
add_test(NAME grid_vs_apfelxx COMMAND test_grid_vs_apfelxx
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(test_f2nc_bsf tests/test_f2nc_lhapdf.cpp)
target_link_libraries(test_f2nc_bsf PRIVATE pineapfel)
target_compile_options(test_f2nc_bsf PRIVATE -Wall -Wextra -O3)
add_test(NAME f2nc_bsf COMMAND test_f2nc_bsf
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(test_sidis_g1_bsf tests/test_sidis_g1_bsf.cpp tests/sidis_helper.cpp)
target_link_libraries(test_sidis_g1_bsf PRIVATE pineapfel)
target_compile_options(test_sidis_g1_bsf PRIVATE -Wall -Wextra -O3)
add_test(NAME sidis_g1_bsf COMMAND test_sidis_g1_bsf
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

# Diagnostic test: SIDIS G1 convolution with LHAPDF vs BuildSidisG1 reference.
# Uses lhapdf-config to find include/link flags.
execute_process(COMMAND lhapdf-config --cflags
    OUTPUT_VARIABLE LHAPDF_CFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND lhapdf-config --libs
    OUTPUT_VARIABLE LHAPDF_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
separate_arguments(LHAPDF_CFLAGS_LIST UNIX_COMMAND "${LHAPDF_CFLAGS}")
separate_arguments(LHAPDF_LIBS_LIST   UNIX_COMMAND "${LHAPDF_LIBS}")

add_executable(test_sidis_g1_lhapdf tests/test_sidis_g1_lhapdf.cpp)
target_compile_options(test_sidis_g1_lhapdf PRIVATE -Wall -Wextra -O3 ${LHAPDF_CFLAGS_LIST})
target_link_libraries(test_sidis_g1_lhapdf PRIVATE pineapfel ${LHAPDF_LIBS_LIST})
# Not added to ctest: this is a manual diagnostic, not an automated pass/fail test.

# Install targets
install(TARGETS pineapfel pineapfel-evolve pineapfel-build
    EXPORT pineapfelTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(DIRECTORY pineapfel/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT pineapfelTargets
    FILE pineapfelTargets.cmake
    NAMESPACE pineapfel::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pineapfel
)

# Package config
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/pineapfelConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/pineapfelConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pineapfel
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/pineapfelConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/pineapfelConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/pineapfelConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pineapfel
)
