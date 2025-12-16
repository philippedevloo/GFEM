# This is not the best. Please create (the extense) list of needed files
file(GLOB common_headers *.h)
file(GLOB common_sources *.cpp)
list(REMOVE_ITEM common_sources "${CMAKE_CURRENT_SOURCE_DIR}/TestAutodiff.cpp")
#include (../cmake/EnableEigen.cmake)


find_package(nlohmann_json 3.2.0 REQUIRED)
find_package(Eigen3 REQUIRED)

# Create the iMRS library statically!
add_library(GFEM_library STATIC ${common_headers} ${common_sources})
#enable_eigen(GFEM_library)
#target_compile_definitions(GFEM_library PUBLIC USING_EIGEN)
#message("GFEM lib EIGEN3_INCLUDE_DIR: ${EIGEN3_INCLUDE_DIR}")
target_include_directories(GFEM_library PUBLIC Eigen3::Eigen)

# Links the other libraries
target_link_libraries(GFEM_library PUBLIC NeoPZ::pz PUBLIC
        nlohmann_json::nlohmann_json)
add_executable(test_autodiff TestAutodiff.cpp)
target_link_libraries(test_autodiff PUBLIC autodiff2 PUBLIC NeoPZ::pz PUBLIC Eigen3::Eigen)

# Unit test for TPZGFemEnrichment (Catch2)
add_executable(test_TPZGFemEnrichment test_TPZGFemEnrichment.cpp)
# Link with Catch2 main to provide the test runner entry point (Catch2 v3 uses an imported target
# that provides main). If your Catch2 installation provides Catch2::Catch2WithMain, prefer that.
if(TARGET Catch2::Catch2WithMain)
        target_link_libraries(test_TPZGFemEnrichment PUBLIC GFEM_library NeoPZ::pz Catch2::Catch2WithMain Eigen3::Eigen)
else()
        target_link_libraries(test_TPZGFemEnrichment PUBLIC GFEM_library NeoPZ::pz Catch2::Catch2 Eigen3::Eigen)
endif()
add_test(NAME test_TPZGFemEnrichment COMMAND test_TPZGFemEnrichment)
