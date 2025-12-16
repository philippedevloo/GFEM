# This project estimates the error of a flow around a NACA profile
set(GFem2Dsrc
    GFem2D.cpp
)

set(SBFem2Dsrc
    SBFem2D.cpp
)

add_executable(GFem2D ${GFem2Dsrc})
#target_link_libraries(GFem2D NeoPZ::pz GFEM_library Eigen3::Eigen)
target_link_libraries(GFem2D NeoPZ::pz GFEM_library)
# target_compile_definitions(GFem2D PUBLIC USING_EIGEN)

find_package(nlohmann_json 3.2.0 REQUIRED)

add_executable(SBFem2D ${SBFem2Dsrc})
target_link_libraries(SBFem2D NeoPZ::pz GFEM_library nlohmann_json::nlohmann_json)

add_executable(TestGFemEnrichment TestGFemEnrichment.cpp)
target_link_libraries(TestGFemEnrichment NeoPZ::pz GFEM_library nlohmann_json::nlohmann_json)
# copy a file from source directory to binary directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/quadmesh.msh ${CMAKE_CURRENT_BINARY_DIR}/quadmesh.msh COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/SquareFrac2.msh ${CMAKE_CURRENT_BINARY_DIR}/SquareFrac.msh COPYONLY)
configure_file(CrackDef.json CrackDef.json COPYONLY)

