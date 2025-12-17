cmake_minimum_required (VERSION 3.13)

project (GFEM)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

include(cmake/StandardProjectSettings.cmake)

# Enables NeoPZ library
find_package(NeoPZ REQUIRED HINTS ${CMAKE_SOURCE_DIR}/../neopz_install/ ${CMAKE_SOURCE_DIR}/neopz_install/)
# Enables Eigen3 library
find_package(Eigen3 REQUIRED)

# Catch2 must be available before fetching autodiff
find_package(Catch2 3.0.0 QUIET)
if(NOT Catch2_FOUND)
    message(STATUS "Catch2 not found, fetching from GitHub...")
    include(FetchContent)
    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.4.0
    )
    FetchContent_MakeAvailable(Catch2)
    # Ensure Catch2 can be found by subdependencies
    list(APPEND CMAKE_PREFIX_PATH "${catch2_BINARY_DIR}")
endif()

# Include autodiff fetcher (after Catch2 is available)
include(FetchAutodiff)

# Mac OS X specific code
if (APPLE)
  add_definitions(-DMACOSX)
endif()

#TODO: Define this using modern CMake standards. @Gustavo, do you know how to do this?
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    option(GFEM_DEBUG "Debug verifications" OFF)
    if(ERRORESTIMATION_DEBUG)
        add_definitions(-DGFEM_DEBUG)
    endif()
endif()

# This option enables a lot of warnings and treat them as errors, to ensure
# good programming practices are used. Since its behaviour is extreme, it
# should be turned off by default.
option(ENABLE_ALL_WARNINGS "Enable a lot of compiler warnings and treat them as errors" OFF)
if (ENABLE_ALL_WARNINGS)
    add_library(project_warnings INTERFACE)
    include(cmake/CompilerWarnings.cmake)
    set_project_warnings(project_warnings)
    set(CMAKE_CXX_FLAGS "")
endif ()

find_package(nlohmann_json 3.2.0 REQUIRED)

#include_directories(/Users/philippedevloo/GitHub/GFemResearch/neopz_GFEM/build-cmake/_deps/eigen-src)
include_directories(${CMAKE_SOURCE_DIR}/NeoPZ_extensions)

# Creates GFEMLib with common extension classes
add_subdirectory(NeoPZ_extensions)

# Creates ErrorEstimationLib with common error estimation file
add_subdirectory(2D)
add_subdirectory(3D)

# Add UnitTests
add_subdirectory(UnitTests)

# Created executable targets


