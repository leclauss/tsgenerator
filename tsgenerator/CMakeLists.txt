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

  add_library(fftw3 STATIC IMPORTED)

  get_filename_component(FFTW3_LIB_PATH ../libs/${ARCH}/fftw3/libfftw3-3.lib
    ABSOLUTE)

  set_target_properties(fftw3 PROPERTIES IMPORTED_LOCATION
    ${FFTW3_LIB_PATH})

  add_library(shared SHARED
    src/motifsetcollection
    src/freepositions
    src/basets
    src/scrimpplusplus
    src/tsgenerator)

  set_target_properties(shared PROPERTIES
    OUTPUT_NAME ${PROJECT_NAME})

  target_include_directories(shared PRIVATE
    api include src
    PUBLIC
    ../libs/${ARCH}/fftw3)

  target_link_libraries(shared
    fftw3)

  add_library(static STATIC
    src/motifsetcollection
    src/freepositions
    src/basets
    src/scrimpplusplus
    src/tsgenerator)

  set_target_properties(static PROPERTIES
    OUTPUT_NAME ${PROJECT_NAME})

  target_include_directories(static PRIVATE
    api include src
    PUBLIC
    ../libs/${ARCH}/fftw3)

  target_link_libraries(static
    fftw3)
else ()
  # check the dependencies
  find_package(PkgConfig REQUIRED)

  pkg_check_modules(FFTW3 REQUIRED fftw3)


  # build
  add_library(${PROJECT_NAME} SHARED
    src/motifsetcollection
    src/freepositions
    src/basets
    src/scrimpplusplus
    src/tsgenerator)

  target_include_directories(${PROJECT_NAME} PRIVATE
    api include src
    PUBLIC
    ${FFTW3_INCLUDEDIR})

  target_link_libraries(${PROJECT_NAME}
    ${FFTW3_LIBRARIES})

  target_compile_options(${PROJECT_NAME} PUBLIC
    ${FFTW3_CFLAGS_OTHER})

  set_target_properties(${PROJECT_NAME} PROPERTIES
    VERSION ${TSG_VERSION}
    SOVERSION ${TSG_SOVERSION}
    PUBLIC_HEADER
    api/tsgenerator.hpp)


  # install
  include(GNUInstallDirs)

  configure_file(tsgenerator.pc.in tsgenerator.pc @ONLY)

  install(TARGETS ${PROJECT_NAME}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})

  install(FILES ${CMAKE_BINARY_DIR}/tsgenerator.pc DESTINATION
    ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)


  # package
  set(CPACK_GENERATOR DEB)
  set(CPACK_PACKAGING_INSTALL_PREFIX "/usr")
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
endif ()

