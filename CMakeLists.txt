cmake_minimum_required(VERSION 3.10)

# define project name, version
project(PSEMolDyn_GroupB)

# let ccmake and cmake-gui offer the default build type options
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

# set Release as the default build type if it is not yet set.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# collect all cpp files
file(GLOB_RECURSE MY_SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
        # header don't need to be included but this might be necessary for some IDEs
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
        )

# create make target
#add_executable(MolSim ${MY_SRC} src/ForceCalculation.h src/Gravitation.h src/Gravitation.cpp src/XMLReader)
add_executable(MolSim ${MY_SRC} src/forceCalculation/LennardJones.cpp src/forceCalculation/LennardJones.h src/forceCalculation/MixedLennardJones.cpp src/forceCalculation/MixedLennardJones.h src/container/LinkedCells.cpp src/container/LinkedCells.h src/container/LinkedCells.cpp src/container/LinkedCells.h src/XMLReader/MolSimImpl.h src/CheckpointFileReader.cpp src/CheckpointFileReader.h src/CheckpointFileWriter.cpp src/CheckpointFileWriter.h src/container/ParallelLinkedCells.cpp src/container/ParallelLinkedCells.h src/densityAndVelocityProfile/DensityAndVelocityToCSVWriter.h src/densityAndVelocityProfile/DensityAndVelocityProfile.h)

# set cxx standard. You may raise this if you want.
target_compile_features(MolSim
        PRIVATE
        cxx_std_14
        )

include(FetchContent)
FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG        v1.9.2
)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(spdlog googletest)

target_include_directories(MolSim
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/libs/libxsd
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        )

target_link_libraries(MolSim
        # stuff that is used in headers and source files
        PUBLIC
        xerces-c
        spdlog
        -fopenmp
        -frounding-math
        )

# activate all compiler warnings. Clean up your code :P
# depending on the compiler different flags are used
target_compile_options(MolSim
        PRIVATE
        $<$<CXX_COMPILER_ID:GNU>:-Wall>
        # disable some spam warnings for icpc...
        $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>
        -fopenmp
        -frounding-math
        )

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

include(doxygen)

enable_testing()

add_executable(
        main_test
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/main_test.cpp
)


target_include_directories(main_test
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/libs/libxsd
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/tests
        )

target_compile_features(main_test
        PRIVATE
        cxx_std_14
        )



target_link_libraries(
        main_test
        PUBLIC
        gtest_main
        PUBLIC
        xerces-c
        PUBLIC
        -fopenmp
)



include(GoogleTest)
gtest_discover_tests(main_test)