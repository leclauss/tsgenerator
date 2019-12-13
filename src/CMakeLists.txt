MESSAGE(STATUS "Trying to install fftw...")

include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

ExternalProject_Add(project_fftw
  #GIT_REPOSITORY  https://github.com/FFTW/fftw3
  URL http://www.fftw.org/fftw-3.3.2.tar.gz
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/fftw
  CONFIGURE_COMMAND
${CMAKE_CURRENT_BINARY_DIR}/fftw/src/project_fftw/configure
--prefix=${CMAKE_CURRENT_BINARY_DIR}/fftw/install
  INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/fftw/install
  )

add_library(fftw STATIC IMPORTED)
set(lib_fftw_name
${CMAKE_STATIC_LIBRARY_PREFIX}fftw3${CMAKE_STATIC_LIBRARY_SUFFIX})
if (EXISTS ${CMAKE_CURRENT_BINARY_DIR}/fftw/install/lib/${lib_fftw_name})
  set_target_properties(fftw PROPERTIES IMPORTED_LOCATION
	  ${CMAKE_CURRENT_BINARY_DIR}/fftw/install/lib/${lib_fftw_name})
elseif (EXISTS ${CMAKE_CURRENT_BINARY_DIR}/fftw/install/lib64/${lib_fftw_name})
  set_target_properties(fftw PROPERTIES IMPORTED_LOCATION
	  ${CMAKE_CURRENT_BINARY_DIR}/fftw/install/lib64/${lib_fftw_name})
else ()
	message(SEND_ERROR "cannot find fftw3 lib folder")
endif()

include_directories(${PROJECT_SOURCE_DIR}/include ${CMAKE_CURRENT_BINARY_DIR}/fftw/install/include)
link_directories(${PROJECT_SOURCE_DIR}/src)

add_executable(${PROJECT_NAME} main.cpp global.cpp motifsetcollection.cpp
  outputgenerator.cpp freepositions.cpp basets.cpp tsgenerator.cpp fftw)
if (WIN32)
  add_dependencies(${PROJECT_NAME} project_fftw)
  target_link_libraries(${PROJECT_NAME} fftw)
else ()
  add_dependencies(${PROJECT_NAME} project_fftw)
  target_link_libraries(${PROJECT_NAME} PUBLIC stdc++fs fftw)
endif ()

add_executable(unittests unittests.cpp global.cpp motifsetcollection.cpp
  outputgenerator.cpp freepositions.cpp basets.cpp tsgenerator.cpp fftw)
if (WIN32)
  add_dependencies(unittests project_fftw)
  target_link_libraries(unittests fftw)
else ()
  add_dependencies(unittests project_fftw)
  target_link_libraries(unittests PUBLIC stdc++fs fftw)
endif ()
