include(FetchContent)

# Declare autodiff
FetchContent_Declare(
    autodiff
    GIT_REPOSITORY https://github.com/autodiff/autodiff.git
    GIT_TAG        v1.1.0   # or any version tag you prefer
)
set(AUTODIFF_BUILD_PYTHON OFF CACHE BOOL "" FORCE)
set(AUTODIFF_BUILD_TESTS OFF CACHE BOOL "" FORCE)
#set(AUTODIFF_USE_EIGEN OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(autodiff)

# autodiff is header-only → but let's create an interface target
add_library(autodiff2 INTERFACE)
target_include_directories(autodiff2 INTERFACE
    ${autodiff_SOURCE_DIR}
)