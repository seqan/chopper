# ---------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
# ---------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.18)

## CUSTOMISE

# Define the application name and version.
project (chopper
         LANGUAGES CXX
         VERSION 1.0.0
)

## BUILD

# Optimization
option (CHOPPER_NATIVE_BUILD "Optimize build for current architecture." ON)

# Make Release default build type
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE
         Release
         CACHE STRING "Choose the type of build, options are: Debug Release Coverage RelWithDebInfo MinSizeRel." FORCE
    )
endif ()

set (CHOPPER_SUBMODULES_DIR
     "${CMAKE_CURRENT_LIST_DIR}/lib"
     CACHE STRING "Directory containing submodules."
)

# Specify the directories where to store the built archives, libraries and executables
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set (CMAKE_INSTALL_BINDIR "bin")

# Messages
string (ASCII 27 Esc)
set (FontBold "${Esc}[1m")
set (FontReset "${Esc}[m")

# Dependency: SeqAn3.
find_package (SeqAn3 REQUIRED QUIET HINTS ${CHOPPER_SUBMODULES_DIR}/seqan3/build_system)

# Allow to include CMake scripts from seqan3.
list (APPEND CMAKE_MODULE_PATH "${SEQAN3_CLONE_DIR}/test/cmake/")
list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/test/cmake/")

# Define cmake configuration flags to configure and build external projects with the same flags as specified for
# this project.
set (SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_INSTALL_LIBDIR=lib")

# Use ccache.
set (CHOPPER_USE_CCACHE
     ON
     CACHE BOOL "Use ccache if available."
)
if (CHOPPER_USE_CCACHE)
    include (chopper_require_ccache)
    chopper_require_ccache ()
endif ()

# Dependency: Sharg.
find_package (Sharg REQUIRED QUIET HINTS ${CHOPPER_SUBMODULES_DIR}/sharg-parser/build_system)

# Dependency: HIBF.
include (FetchContent)
if (NOT TARGET seqan::hibf)
    FetchContent_Declare (chopper_hibf_project SOURCE_DIR ${CHOPPER_SUBMODULES_DIR}/hibf)

    set (HIBF_NATIVE_BUILD
         ${CHOPPER_NATIVE_BUILD}
         CACHE BOOL "Optimize build for current architecture."
    )
    set (INSTALL_HIBF
         OFF
         CACHE BOOL "Install HIBF targets."
    )

    FetchContent_MakeAvailable (chopper_hibf_project)
endif ()

# Add the application.
set (CHOPPER_INSTALL
     ON
     CACHE BOOL "Install Chopper."
)

add_subdirectory (src)
message (STATUS "${FontBold}You can run `make` to build the application.${FontReset}")

## DOCUMENTATION

set (CHOPPER_BUILD_DOC
     ON
     CACHE BOOL "Build chopper documentation."
)
if (CHOPPER_BUILD_DOC)
    add_subdirectory (doc EXCLUDE_FROM_ALL)
endif ()

## TEST

set (CHOPPER_BUILD_TEST
     ON
     CACHE BOOL "Build chopper tests."
)
if (CHOPPER_BUILD_TEST)
    enable_testing ()
    add_subdirectory (test EXCLUDE_FROM_ALL)
endif ()
