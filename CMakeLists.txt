cmake_minimum_required(VERSION 3.10.2)

set(TSG_SOVERSION 1)
set(TSG_VERSION ${TSG_SOVERSION}.0.0)
set(LIB_NAME tsgenerator)
project(TSGenerator VERSION ${TSG_VERSION} DESCRIPTION "TSGenerator library"
  LANGUAGES CXX)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif ()

set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "--enable-checking -g -O0 -v -da -Q")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

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

add_test(NAME "Unit_tests" COMMAND unittests)
add_test(NAME "Help_text" COMMAND ${PROJECT_NAME} -h)
add_test(NAME "Version_text" COMMAND ${PROJECT_NAME} -v)


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

  set(TSG_INCLUDEDIR tsgenerator/include)
  set(FFTW3_INCLUDEDIR libs/${ARCH}/fftw3)

  file(GLOB FFTW3_LIB_NAME "libs/${ARCH}/fftw3/libfftw3-*.lib")
  string(REGEX MATCH "libfftw3-.*.lib" FFTW3_LIB_NAME ${FFTW3_LIB_NAME})

  file(GLOB FFTW3_DLL_NAME "libs/${ARCH}/fftw3/libfftw3-*.dll")
else ()
  set(STDCPPFS stdc++fs)

  find_package(PkgConfig REQUIRED)
  pkg_check_modules(FFTW3 REQUIRED fftw3)
endif ()

# check dependencies
find_package(PkgConfig REQUIRED)

pkg_check_modules(TSG REQUIRED tsgenerator)
pkg_check_modules(FFTW3 REQUIRED fftw3)

# add the libraries
find_library(TSG_LOCATION NAMES tsgenerator
  HINTS libs/${ARCH}/tsgenerator)

add_library(TSG SHARED IMPORTED)
set_target_properties(TSG PROPERTIES
  IMPORTED_LOCATION ${TSG_LOCATION})

find_library(FFTW3_LOCATION NAMES fftw3 ${FFTW3_LIB_NAME}
  HINTS libs/${ARCH}/fftw3)

if (FFTW3_LOCATION EQUAL "FFTW3_LOCATION-NOTFOUND")
  message(FATAL "Cannot find FFTW3")
endif ()

add_library(FFTW3 STATIC IMPORTED)
set_target_properties(FFTW3 PROPERTIES
  IMPORTED_LOCATION ${FFTW3_LOCATION})

# build
add_executable(${PROJECT_NAME}
  src/main.cpp
  src/global.cpp
  src/outputgenerator.cpp)

target_include_directories(${PROJECT_NAME} PRIVATE
  include src
  PUBLIC
  ${TSG_INCLUDEDIR}
  ${FFTW3_INCLUDEDIR})

target_compile_options(${PROJECT_NAME} PUBLIC
  ${TSG_CFLAGS_OTHER}
  ${FFTW3_CFLAGS_OTHER})

target_link_libraries(${PROJECT_NAME} PUBLIC
  STDCPPFS
  TSG
  FFTW3)

add_executable(unittests
  src/unittests.cpp
  src/global.cpp
  src/outputgenerator.cpp)

target_include_directories(unittests PRIVATE
  include src
  PUBLIC
  ${TSG_INCLUDEDIR}
  ${FFTW3_INCLUDEDIR})

target_compile_options(unittests PUBLIC
  ${TSG_CFLAGS_OTHER}
  ${FFTW3_CFLAGS_OTHER})

target_link_libraries(unittests PUBLIC
  STDCPPFS
  TSG
  FFTW3)

if (WIN32)
  add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD COMMAND ${CMAKE_COMMAND}
    -E copy ${FFTW3_DLL_NAME} $<TARGET_FILE_DIR:${PROJECT_NAME}>)

  add_custom_command(TARGET unittests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
    copy ${FFTW3_DLL_NAME} $<TARGET_FILE_DIR:unittest>)
endif ()
