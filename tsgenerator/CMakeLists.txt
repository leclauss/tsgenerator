cmake_minimum_required(VERSION 3.10.2)

# setup project
set(PROJECT_HOMEPAGE_URL
  "https://gitlab.com/r.moczalla/TSGenerator/tree/master")
set(TSG_SOVERSION 1)
set(TSG_VERSION ${TSG_SOVERSION}.0.0)
project(tsgenerator VERSION ${TSG_VERSION} DESCRIPTION "TSGenerator library"
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

option(BUILD_DOC "Build documentation" ON)
find_package(Doxygen)
if(DOXYGEN_FOUND)
    set(DOXYGEN_CONF ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxygen.conf)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    configure_file(${DOXYGEN_CONF} ${DOXYGEN_OUT} @ONLY)
    message("Documentation build started")

    add_custom_target(doc COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMENT "Doxygen started"
      VERBATIM)
else(DOXYGEN_FOUND)
  message("Doxygen missing for documentation build")
endif(DOXYGEN_FOUND)

enable_testing()

add_test(NAME "unit tests" COMMAND unittests)

# check the dependencies
if (WIN32)
  # build project
  set(ARCH win32)
  if (CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(ARCH win64)
  endif()
  if (CMAKE_GENERATOR_PLATFORM MATCHES .*32.*)
    set(ARCH win32)
  endif()
  if (CMAKE_GENERATOR_PLATFORM MATCHES .*64.*)
    set(ARCH win64)
  endif()

  file(GLOB FFTW3_LIB_NAME "../libs/${ARCH}/fftw3/libfftw3-*.lib")
  string(REGEX MATCH "libfftw3-.*.lib" FFTW3_LIB_NAME ${FFTW3_LIB_NAME})

  set(FFTW3_INCLUDEDIR ../libs/${ARCH}/fftw3)
else ()
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(FFTW3 REQUIRED fftw3)
endif ()


# build
find_library(FFTW3_LOCATION NAMES fftw3 ${FFTW3_LIB_NAME}
  HINTS ../libs/${ARCH}/fftw3)

if (FFTW3_LOCATION EQUAL "FFTW3_LOCATION-NOTFOUND")
  message(FATAL_ERROR "Cannot find FFTW3")
endif ()

add_library(FFTW3 STATIC IMPORTED)
set_target_properties(FFTW3 PROPERTIES
  IMPORTED_LOCATION ${FFTW3_LOCATION})

add_library(shared SHARED
  src/motifsetcollection
  src/freepositions
  src/basets
  src/scrimpplusplus
  src/tsgenerator)

target_include_directories(shared PRIVATE
  api include src
  PUBLIC
  ${FFTW3_INCLUDEDIR})

target_link_libraries(shared PUBLIC
  FFTW3)

target_compile_options(shared PUBLIC
  ${FFTW3_CFLAGS_OTHER})

set_target_properties(shared PROPERTIES
  OUTPUT_NAME ${PROJECT_NAME}
  VERSION ${TSG_VERSION}
  SOVERSION ${TSG_SOVERSION})

add_library(static STATIC
  src/motifsetcollection
  src/freepositions
  src/basets
  src/scrimpplusplus
  src/tsgenerator)

set_target_properties(static PROPERTIES
  OUTPUT_NAME ${PROJECT_NAME}
  VERSION ${TSG_VERSION}
  SOVERSION ${TSG_SOVERSION})

target_include_directories(static PRIVATE
  api include src
  PUBLIC
  ${FFTW3_INCLUDEDIR})

target_link_libraries(static PUBLIC
  FFTW3)

add_executable(unittests
  test/unittests.cpp)

target_include_directories(unittests PRIVATE
  test include src
  PUBLIC
  ${FFTW3_INCLUDEDIR})

target_compile_options(unittests PUBLIC
  ${FFTW3_CFLAGS_OTHER})

target_link_libraries(unittests PUBLIC
  ${STDCPPFS}
  shared
  FFTW3)


# install
include(GNUInstallDirs)

configure_file(tsgenerator.pc.in tsgenerator.pc @ONLY)

install(TARGETS shared
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(TARGETS static
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(FILES
  include/tsgenerator.hpp
  include/basets.hpp
  include/freepositions.hpp
  include/motifsetcollection.hpp
  include/scrimpplusplus.hpp
  DESTINATION
  ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})

install(FILES ${CMAKE_BINARY_DIR}/tsgenerator.pc DESTINATION
  ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)


# package
set(CPACK_GENERATOR DEB)
set(CPACK_PACKAGE_FILE_NAME ${PROJECT_NAME}-${TSG_VERSION})
set(CPACK_PACKAGE_CONTACT "Rafael Moczalla")
set(CPACK_DEBIAN_PACKAGE_DEPENDS libfftw3-dev)

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

