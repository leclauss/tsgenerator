include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

if (WIN32)
  set(ARCH win32)
  if (CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(ARCH win64)
  endif()

  add_library(fftw STATIC IMPORTED)
  set_target_properties(fftw PROPERTIES IMPORTED_LOCATION
    ${PROJECT_SOURCE_DIR}/libs/${ARCH}/fftw3/libfftw3-3.lib)

  include_directories(${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/libs/${ARCH}/fftw3.h)
  link_directories(${PROJECT_SOURCE_DIR}/src)

  add_executable(${PROJECT_NAME} main.cpp global.cpp motifsetcollection.cpp
    outputgenerator.cpp freepositions.cpp basets.cpp scrimpplusplus.cpp
    tsgenerator.cpp)
  target_link_libraries(${PROJECT_NAME} fftw)

  add_executable(unittests unittests.cpp global.cpp motifsetcollection.cpp
    outputgenerator.cpp freepositions.cpp basets.cpp scrimpplusplus.cpp
    tsgenerator.cpp)
  target_link_libraries(unittests fftw)
else ()
  MESSAGE(STATUS "Installing fftw3")

  ExternalProject_Add(project_fftw
    #GIT_REPOSITORY  https://github.com/FFTW/fftw3
    URL http://www.fftw.org/fftw-3.3.2.tar.gz
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/fftw
    CONFIGURE_COMMAND
  ${CMAKE_CURRENT_BINARY_DIR}/fftw/src/project_fftw/configure
  --prefix=${CMAKE_CURRENT_BINARY_DIR}/fftw/install
  --libdir=${CMAKE_CURRENT_BINARY_DIR}/fftw/install/lib
    INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/fftw/install
    )

  add_library(fftw STATIC IMPORTED)
  set(lib_fftw_name
  ${CMAKE_STATIC_LIBRARY_PREFIX}fftw3${CMAKE_STATIC_LIBRARY_SUFFIX})
  set_target_properties(fftw PROPERTIES IMPORTED_LOCATION
    ${CMAKE_CURRENT_BINARY_DIR}/fftw/install/lib/${lib_fftw_name})

  include_directories(${PROJECT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_BINARY_DIR}/fftw/install/include)
  link_directories(${PROJECT_SOURCE_DIR}/src)

  add_executable(${PROJECT_NAME} main.cpp global.cpp motifsetcollection.cpp
    outputgenerator.cpp freepositions.cpp basets.cpp scrimpplusplus.cpp
    tsgenerator.cpp fftw)

  add_dependencies(${PROJECT_NAME} project_fftw)
  target_link_libraries(${PROJECT_NAME} PUBLIC stdc++fs fftw)

  add_executable(unittests unittests.cpp global.cpp motifsetcollection.cpp
    outputgenerator.cpp freepositions.cpp basets.cpp scrimpplusplus.cpp
    tsgenerator.cpp fftw)

  add_dependencies(unittests project_fftw)
  target_link_libraries(unittests PUBLIC stdc++fs fftw)
endif ()
