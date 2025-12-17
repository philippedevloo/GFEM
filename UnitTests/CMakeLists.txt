enable_testing()
#in this file (and ONLY in this file) we will define CATCH_CONFIG_MAIN
add_library(test_library STATIC test_main.cpp)
if (APPLE)
  target_compile_definitions(test_library PUBLIC MACOSX)
endif()
if(PZ_LOG)
  target_link_libraries(test_library PUBLIC ${Log4cxx_LIBRARY})
endif()

target_link_libraries(test_library PUBLIC GFEM_library PUBLIC Catch2::Catch2)
if(APPLE)
   target_link_libraries(test_library PRIVATE ${ACCELERATE_LIB})
endif(APPLE)

include(${CMAKE_SOURCE_DIR}/cmake/add_unit_test.cmake)

add_unit_test(Test_TPZGFEMEnrichment test_TPZGFemEnrichment.cpp)
add_unit_test(Test_Autodiff TestAutodiff.cpp)

ENABLE_TESTING()
include(CTest)