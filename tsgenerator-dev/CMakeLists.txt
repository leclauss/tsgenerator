cmake_minimum_required(VERSION 3.10.2)

# setup project
set(PROJECT_HOMEPAGE_URL
  "https://gitlab.com/r.moczalla/TSGenerator/tree/master")
set(TSG_SOVERSION 3)
set(TSG_VERSION ${TSG_SOVERSION}.0.0)
set(PROJECT_LIB_NAME tsgenerator)
project(${PROJECT_LIB_NAME}-dev
  VERSION ${TSG_VERSION}
  DESCRIPTION
"Synthetic time series generator for evaluation purposes of time series \
motif discovery algorithms."
  LANGUAGES CXX)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif ()


# check requirements
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if (WIN32)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W4")
else ()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
  set(CMAKE_CXX_FLAGS_DEBUG "--enable-checking -g -O0 -v -da -Q")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()

if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8)
  message(FATAL_ERROR "This project requires gcc version 8 or higher")
endif ()

enable_testing()

add_test(NAME "unit_tests" COMMAND unittests)


# build
add_library(shared SHARED
  src/motifsetcollection
  src/freepositions
  src/basets
  src/tpm
  src/tsm
  src/tsgenerator)

target_include_directories(shared PRIVATE
  api include src)

set_target_properties(shared PROPERTIES
  OUTPUT_NAME ${PROJECT_LIB_NAME}
  VERSION ${TSG_VERSION}
  SOVERSION ${TSG_SOVERSION})

add_library(static STATIC
  src/motifsetcollection
  src/freepositions
  src/basets
  src/tpm
  src/tsm
  src/tsgenerator)

set_target_properties(static PROPERTIES
  OUTPUT_NAME ${PROJECT_LIB_NAME}
  VERSION ${TSG_VERSION}
  SOVERSION ${TSG_SOVERSION})

target_include_directories(static PRIVATE
  api include src)

add_executable(unittests
  test/unittests.cpp)

target_include_directories(unittests PRIVATE
  test include src)

target_link_libraries(unittests PUBLIC
  ${STDCPPFS}
  shared)


# install
include(GNUInstallDirs)

configure_file(share/pkgconfig/tsgenerator.pc.in tsgenerator.pc @ONLY)

install(TARGETS shared
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})

install(TARGETS static
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})

install(FILES
  include/tsgenerator.hpp
  include/basets.hpp
  include/freepositions.hpp
  include/motifsetcollection.hpp
  include/tpm.hpp
  include/tsm.hpp
  include/tsgtypes.hpp
  include/tsgutils.hpp
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_LIB_NAME})

install(FILES ${CMAKE_BINARY_DIR}/tsgenerator.pc
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)


# package
set(CPACK_GENERATOR DEB)
set(CPACK_PACKAGE_FILE_NAME ${PROJECT_NAME}-${TSG_VERSION})
set(CPACK_PACKAGE_CONTACT "Rafael Moczalla")

include(CPack)


# uninstall
if (NOT TARGET uninstall)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P
    ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif ()

