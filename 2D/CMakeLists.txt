# This project estimates the error of a flow around a NACA profile
set(GFem2Dsrc
    GFem2D.cpp
)

add_executable(GFem2D ${GFem2Dsrc})
target_link_libraries(GFem2D NeoPZ::pz GFEM_library)
# copy a file from source directory to binary directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/quadmesh.msh ${CMAKE_CURRENT_BINARY_DIR}/quadmesh.msh COPYONLY)
# configure_file(${CMAKE_CURRENT_SOURCE_DIR}/SquareFrac.msh ${CMAKE_CURRENT_BINARY_DIR}/SquareFrac.msh COPYONLY)