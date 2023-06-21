cmake_minimum_required(VERSION 3.16...3.22)

include(cmake/BuildType.cmake)

project(
  amstools
  VERSION 0.4
  DESCRIPTION "Collection of fast and efficient biological sequence tools"
  LANGUAGES CXX
)

## prevent in-source builds
#
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
  message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

## must include after the project call due to GNUInstallDirs requiring a
## language be enabled
#
include(GNUInstallDirs)

## check for zlib
#
find_package(ZLIB)

## check for cxxopts
#
find_package(
  cxxopts CONFIG
  HINTS $ENV{HOME} $ENV{HOME}/.local /usr/local /usr
)
if(NOT cxxopts_FOUND)
    message(STATUS "Fetching cxxopts library...")
    include(FetchContent)
    FetchContent_Declare(
      cxxopts
      GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
      GIT_TAG        v3.1.1
    )
    FetchContent_MakeAvailable(cxxopts)
endif()

enable_testing()

add_subdirectory(src)
