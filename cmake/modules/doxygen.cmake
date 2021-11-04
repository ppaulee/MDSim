# make doc_doxygen optional if someone does not have / like doxygen

# TODO: create CMake build option for the target.
option(BUILD_DOXY "Build documentation" OFF)
# TODO: Add a custom target for building the documentation.
if (BUILD_DOXY)
    find_package(Doxygen)
    if (DOXYGEN_FOUND)
        add_custom_target(doxyDoc COMMAND ${DOXYGEN_EXECUTABLE} "${CMAKE_BINARY_DIR}/../Doxyfile"
                WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}
                COMMENT "Building user's documentation into doxyDoc build dir..."
                )
    else (DOXYGEN_FOUND)
        message("Doxygen need to be installed to generate the doxygen documentation")
    endif (DOXYGEN_FOUND)
endif (BUILD_DOXY)