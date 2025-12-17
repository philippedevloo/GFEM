function(enable_Catch2)
    #perhaps catch was already downloaded when running cmake
    if(NOT TARGET Catch2::Catch2)
        find_package(Catch2 3 QUIET)
        if(NOT TARGET Catch2::Catch2)
            Include(FetchContent)
            if(NOT Catch2_SOURCE_DIR)
                message(STATUS "Downloading Catch2")
            endif()
            FetchContent_Declare(
                Catch2
                GIT_REPOSITORY https://github.com/catchorg/Catch2.git
                GIT_TAG        v3.4.0)
            FetchContent_MakeAvailable(Catch2)
            list(APPEND CMAKE_PREFIX_PATH "${catch2_BINARY_DIR}")
        endif()
    endif()
    if(TARGET Catch2::Catch2)
        message(STATUS "Catch2 found at ${Catch2_SOURCE_DIR}")
    else()
        message(FATAL_ERROR "Could not satisfy dependency: Catch2")
    endif()
endfunction()