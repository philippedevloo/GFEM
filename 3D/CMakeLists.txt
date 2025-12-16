# This project estimates the error of a flow around a NACA profile
set(GFem3Dsrc
    GFem3D.cpp
)

add_executable(GFem3D ${GFem3Dsrc})
target_link_libraries(GFem3D NeoPZ::pz GFEM_library Eigen3::Eigen)
# target_compile_definitions(GFem3D PUBLIC USING_EIGEN)
# copy a file from source directory to binary directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Reference.msh ${CMAKE_CURRENT_BINARY_DIR}/Reference.msh COPYONLY)
