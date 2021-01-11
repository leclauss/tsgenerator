cmake_minimum_required(VERSION 3.10.2)

set(TSG_SOVERSION 2)
set(TSG_VERSION ${TSG_SOVERSION}.0.0)
set(LIB_NAME tsgenerator)
project(tsgenerator VERSION ${TSG_VERSION}
  DESCRIPTION
"Synthetic time series generator for evaluation purposes of time series \
motif discovery algorithms."
  LANGUAGES CXX)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif ()

set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME})

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
  message("Documentation build started")

  configure_file(share/doc/Doxygen.in Doxygen @ONLY)

  add_custom_target(doc
    COMMAND ${DOXYGEN_EXECUTABLE} Doxygen
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Doxygen started"
    VERBATIM)
else(DOXYGEN_FOUND)
  message("Doxygen missing for documentation build")
endif(DOXYGEN_FOUND)

enable_testing()

add_test(NAME "unit_tests" COMMAND unittests)


# check dependencies
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

  set(TSG_INCLUDEDIR ../${LIB_NAME}-dev/include)
else ()
  set(STDCPPFS stdc++fs)

  find_package(PkgConfig REQUIRED)
  pkg_check_modules(TSG REQUIRED ${LIB_NAME})
endif ()


# add the libraries
find_library(TSG_LOCATION NAMES ${LIB_NAME}
  HINTS
  ${TSG_LIBRARY_DIRS}
  ../libs/${ARCH}/${LIB_NAME})

add_library(TSG STATIC IMPORTED)
set_target_properties(TSG PROPERTIES
  IMPORTED_LOCATION ${TSG_LOCATION})


# build
configure_file(include/global.hpp.in global.hpp @ONLY)

add_executable(${PROJECT_NAME}
  src/main.cpp
  src/global.cpp
  src/outputgenerator.cpp)

target_include_directories(${PROJECT_NAME} PRIVATE
  include src
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}
  ${TSG_INCLUDEDIR})

target_compile_options(${PROJECT_NAME} PUBLIC
  ${TSG_CFLAGS_OTHER})

target_link_libraries(${PROJECT_NAME} PUBLIC
  ${STDCPPFS}
  TSG)

add_executable(unittests
  test/unittests.cpp
  src/global.cpp
  src/outputgenerator.cpp)

target_include_directories(unittests PRIVATE
  test include src
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}
  ${TSG_INCLUDEDIR})

target_compile_options(unittests PUBLIC
  ${TSG_CFLAGS_OTHER})

target_link_libraries(unittests PUBLIC
  ${STDCPPFS}
  TSG)


# install
include(GNUInstallDirs)

configure_file(share/man/tsgenerator.3.in tsgenerator.3 @ONLY)

install(TARGETS ${PROJECT_NAME}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})

install(FILES ${CMAKE_BINARY_DIR}/tsgenerator.3
  DESTINATION ${CMAKE_INSTALL_MANDIR}/man3)


# package
set(CPACK_GENERATOR DEB)
set(CPACK_PACKAGE_FILE_NAME ${PROJECT_NAME}-${PROJECT_VERSION})
set(CPACK_PACKAGE_CONTACT "Rafael Moczalla")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "tsgenerator-dev")

include(CPack)

