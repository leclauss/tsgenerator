cmake_minimum_required(VERSION 3.10.2)

set(TSG_SOVERSION 1)
set(TSG_VERSION ${TSG_SOVERSION}.0.0)
set(LIB_NAME tsgenerator)
project(tsgenerator-gui VERSION ${TSG_VERSION}
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

  set(TSG_INCLUDEDIR tsgenerator/include)
  set(FFTW3_INCLUDEDIR ../libs/${ARCH}/fftw3)

  file(GLOB FFTW3_LIB_NAME "../libs/${ARCH}/fftw3/libfftw3-*.lib")
  string(REGEX MATCH "../libfftw3-.*.lib" FFTW3_LIB_NAME ${FFTW3_LIB_NAME})

  file(GLOB FFTW3_DLL_NAME "../libs/${ARCH}/fftw3/libfftw3-*.dll")
else ()
  set(STDCPPFS stdc++fs)

  file(GLOB_RECURSE MAC_QT5 LIST_DIRECTORIES true "/usr/local/Cellar/qt/Qt5")

  find_package(PkgConfig REQUIRED)
  pkg_check_modules(TSG REQUIRED tsgenerator)
  pkg_check_modules(FFTW3 REQUIRED fftw3)
endif ()


# add the libraries
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)
find_package(Qt5 COMPONENTS Widgets REQUIRED
  HINTS ${MAC_QT5})

find_library(TSG_LOCATION NAMES tsgenerator
  HINTS ../libs/${ARCH}/tsgenerator)

add_library(TSG STATIC IMPORTED)
set_target_properties(TSG PROPERTIES
  IMPORTED_LOCATION ${TSG_LOCATION})

find_library(FFTW3_LOCATION NAMES fftw3 ${FFTW3_LIB_NAME}
  HINTS ../libs/${ARCH}/fftw3)

if (FFTW3_LOCATION EQUAL "FFTW3_LOCATION-NOTFOUND")
  message(FATAL_ERROR "Cannot find FFTW3")
endif ()

add_library(FFTW3 STATIC IMPORTED)
set_target_properties(FFTW3 PROPERTIES
  IMPORTED_LOCATION ${FFTW3_LOCATION})

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
  ${Qt5Widgets_INCLUDE_DIRS}
  ${TSG_INCLUDEDIR}
  ${FFTW3_INCLUDEDIR})

target_compile_options(${PROJECT_NAME} PUBLIC
  -fPIC ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}
  ${TSG_CFLAGS_OTHER}
  ${FFTW3_CFLAGS_OTHER})

target_link_libraries(${PROJECT_NAME} PUBLIC
  ${STDCPPFS}
  ${Qt5Widgets_LIBRARIES}
  TSG
  FFTW3)


# install
include(GNUInstallDirs)

install(TARGETS ${PROJECT_NAME}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})


# package
set(CPACK_GENERATOR DEB)
set(CPACK_PACKAGE_FILE_NAME ${PROJECT_NAME}-${PROJECT_VERSION})
set(CPACK_PACKAGE_CONTACT "Rafael Moczalla")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "qt5-default,libfftw3-dev,tsgenerator-dev")

include(CPack)

if (WIN32)
  add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD COMMAND ${CMAKE_COMMAND}
    -E copy ${FFTW3_DLL_NAME} $<TARGET_FILE_DIR:${PROJECT_NAME}>)

  add_custom_command(TARGET unittests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
    copy ${FFTW3_DLL_NAME} $<TARGET_FILE_DIR:unittests>)
endif ()
