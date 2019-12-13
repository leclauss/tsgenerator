MESSAGE(STATUS "Trying to install fftw...")

include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

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


include_directories(${PROJECT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_BINARY_DIR}/fftw/install/include)
link_directories(${PROJECT_SOURCE_DIR}/src)

add_executable(${PROJECT_NAME} main.cpp global.cpp motifsetcollection.cpp
  outputgenerator.cpp freepositions.cpp basets.cpp scrimpplusplus.cpp
  tsgenerator.cpp fftw)
add_dependencies(${PROJECT_NAME} project_fftw)

if (WIN32)
  target_link_libraries(${PROJECT_NAME} fftw)
else ()
  set_target_properties(fftw PROPERTIES IMPORTED_LOCATION
    ${CMAKE_CURRENT_BINARY_DIR}/fftw/install/lib/${lib_fftw_name})
  target_link_libraries(${PROJECT_NAME} PUBLIC stdc++fs fftw)
endif ()

add_executable(unittests unittests.cpp global.cpp motifsetcollection.cpp
  outputgenerator.cpp freepositions.cpp basets.cpp scrimpplusplus.cpp
  tsgenerator.cpp fftw)

add_dependencies(unittests project_fftw)
if (WIN32)
  target_link_libraries(unittests fftw)
else ()
  set_target_properties(fftw PROPERTIES IMPORTED_LOCATION
    ${CMAKE_CURRENT_BINARY_DIR}/fftw/install/lib/${lib_fftw_name})
  target_link_libraries(unittests PUBLIC stdc++fs fftw)
endif ()
