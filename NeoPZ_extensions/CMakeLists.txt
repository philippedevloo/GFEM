# This is not the best. Please create (the extense) list of needed files
file(GLOB common_headers *.h)
file(GLOB common_sources *.cpp)

# find_package(nlohmann_json 3.2.0 REQUIRED)
# find_package(Eigen3 REQUIRED)

# Create the iMRS library statically!
add_library(GFEM_library STATIC ${common_headers} ${common_sources})
target_include_directories(GFEM_library PUBLIC Eigen3::Eigen)

# Links the other libraries
target_link_libraries(GFEM_library PUBLIC NeoPZ::pz PUBLIC
        nlohmann_json::nlohmann_json)
target_include_directories(GFEM_library PUBLIC ${PZ_INCLUDE_DIRS})