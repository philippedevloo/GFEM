# Try to find a pybind11 installation.
# It provides:
#   pybind11_INCLUDE_DIRS  - include directories
#   pybind11::module       - imported interface target

# Look for a directory containing pybind11/pybind11.h
find_path(pybind11_INCLUDE_DIRS
    NAMES pybind11/pybind11.h
    HINTS
        ${CMAKE_CURRENT_LIST_DIR}/../external/pybind11/include
        $ENV{CONDA_PREFIX}/include
        $ENV{VIRTUAL_ENV}/include
    PATH_SUFFIXES
        include
        site-packages/pybind11/include     # pip install --user
        python*/site-packages/pybind11/include
)

if(NOT pybind11_INCLUDE_DIRS)
    message(FATAL_ERROR "pybind11 not found. Please install via pip install pybind11")
endif()

# Create an imported interface library
add_library(pybind11::module INTERFACE)

# Add include dirs to the target
target_include_directories(pybind11::module INTERFACE
    ${pybind11_INCLUDE_DIRS}
)

# Optionally detect Python version (for embedding Python)
find_package(Python3 COMPONENTS Interpreter Development QUIET)

if(Python3_FOUND)
    target_include_directories(pybind11::module INTERFACE
        ${Python3_INCLUDE_DIRS}
    )
    target_link_libraries(pybind11::module INTERFACE
        ${Python3_LIBRARIES}
    )
endif()

mark_as_advanced(pybind11_INCLUDE_DIRS)
