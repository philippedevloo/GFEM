# This is not the best. Please create (the extense) list of needed files
file(GLOB common_headers *.h)
file(GLOB common_sources *.cpp)

include (../cmake/EnableEigen.cmake)

# Create the iMRS library statically!
add_library(GFEM_library STATIC ${common_headers} ${common_sources})
enable_eigen(GFEM_library)
target_compile_definitions(GFEM_library PUBLIC USING_EIGEN)
message("GFEM lib EIGEN3_INCLUDE_DIR: ${EIGEN3_INCLUDE_DIR}")
target_include_directories(GFEM_library PUBLIC ${EIGEN3_INCLUDE_DIR})

# Links the other libraries
target_link_libraries(GFEM_library PUBLIC NeoPZ::pz )
