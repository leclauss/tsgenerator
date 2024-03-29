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

  set(TSG_INCLUDEDIR ../tsgenerator-dev/include)

  if (NOT Qt5_PATH)
    file(GLOB_RECURSE WIN_QT5 LIST_DIRECTORIES true "C:/Qt/5.*.*/Qt5")
    file(GLOB_RECURSE WIN_QT5_WIDGETS_DLL "C:/Qt/5.*.*/Qt5Widgets.dll")
    file(GLOB_RECURSE WIN_QT5_CORE_DLL "C:/Qt/5.*.*/Qt5Core.dll")
    file(GLOB_RECURSE WIN_QT5_GUI_DLL "C:/Qt/5.*.*/Qt5Gui.dll")
    file(GLOB_RECURSE WIN_QT5_CHARTS_DLL "C:/Qt/5.*.*/Qt5Charts.dll")
  else ()
    file(GLOB_RECURSE WIN_QT5 LIST_DIRECTORIES true "${Qt5_PATH}/5.*.*/Qt5")
    file(GLOB_RECURSE WIN_QT5_WIDGETS_DLL "${Qt5_PATH}/5.*.*/Qt5Widgets.dll")
    file(GLOB_RECURSE WIN_QT5_CORE_DLL "${Qt5_PATH}/5.*.*/Qt5Core.dll")
    file(GLOB_RECURSE WIN_QT5_GUI_DLL "${Qt5_PATH}/5.*.*/Qt5Gui.dll")
    file(GLOB_RECURSE WIN_QT5_CHARTS_DLL "${Qt5_PATH}/5.*.*/Qt5Charts.dll")
  endif ()
else ()
  set(STDCPPFS stdc++fs)

  if (NOT Qt5_PATH)
    if (APPLE)
      file(GLOB_RECURSE UNIX_QT5 LIST_DIRECTORIES true
        "/usr/local/Cellar/qt/Qt5")
    else ()
      file(GLOB_RECURSE UNIX_QT5 LIST_DIRECTORIES true
        "$ENV{HOME}/Qt5.*.*/5.*.*/gcc_*/Qt5")
    endif ()
  else ()
    file(GLOB_RECURSE UNIX_QT5 LIST_DIRECTORIES true "${Qt5_PATH}/Qt5")
  endif ()
  set(QT5_CFLAGS -fPIC)

  find_package(PkgConfig REQUIRED)
  pkg_check_modules(TSG REQUIRED ${LIB_NAME})
endif ()


# add the libraries
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)
find_package(Qt5 COMPONENTS Widgets Charts REQUIRED
  HINTS ${UNIX_QT5}
  ${WIN_QT5})

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
  src/tsggui.cpp
  src/global.cpp
  src/outputgenerator.cpp
  src/chartview.cpp
  src/chart.cpp
  src/canvas.cpp)

target_include_directories(${PROJECT_NAME} PRIVATE
  include src
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}
  ${Qt5Widgets_INCLUDE_DIRS}
  ${Qt5Charts_INCLUDE_DIRS}
  ${TSG_INCLUDEDIR})

target_compile_options(${PROJECT_NAME} PUBLIC
  ${Qt5_CFLAGS}
  ${TSG_CFLAGS_OTHER})

target_link_libraries(${PROJECT_NAME} PUBLIC
  ${STDCPPFS}
  Qt5::Widgets
  Qt5::Charts
  TSG)


# install
include(GNUInstallDirs)

install(TARGETS ${PROJECT_NAME}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})


# package
set(CPACK_GENERATOR DEB)
set(CPACK_PACKAGE_FILE_NAME ${PROJECT_NAME}-${PROJECT_VERSION})
set(CPACK_PACKAGE_CONTACT "Rafael Moczalla")
set(CPACK_DEBIAN_PACKAGE_DEPENDS
  "qt5-default,libqt5charts5-dev,tsgenerator-dev")

include(CPack)

if (WIN32)
  add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD COMMAND ${CMAKE_COMMAND}
    -E copy ${WIN_QT5_WIDGETS_DLL} $<TARGET_FILE_DIR:${PROJECT_NAME}>)
  add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD COMMAND ${CMAKE_COMMAND}
    -E copy ${WIN_QT5_CORE_DLL} $<TARGET_FILE_DIR:${PROJECT_NAME}>)
  add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD COMMAND ${CMAKE_COMMAND}
    -E copy ${WIN_QT5_GUI_DLL} $<TARGET_FILE_DIR:${PROJECT_NAME}>)
  add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD COMMAND ${CMAKE_COMMAND}
    -E copy ${WIN_QT5_CHARTS_DLL} $<TARGET_FILE_DIR:${PROJECT_NAME}>)
endif ()
