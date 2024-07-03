# This is not the best. Please create (the extense) list of needed files
file(GLOB common_headers *.h)
file(GLOB common_sources *.cpp)

# Create the iMRS library statically!
add_library(GFEM_library STATIC ${common_headers} ${common_sources})

# Links the other libraries
target_link_libraries(GFEM_library PUBLIC NeoPZ::pz)
