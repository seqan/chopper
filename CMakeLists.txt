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

# LTO
include (CheckIPOSupported)
check_ipo_supported (RESULT CHOPPER_HAS_LTO OUTPUT CHOPPER_HAS_LTO_OUTPUT)
if (CHOPPER_HAS_LTO)
    set (CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif ()

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

set (CPM_INDENT "CMake Package Manager CPM: ")
include (${chopper_SOURCE_DIR}/cmake/CPM.cmake)
CPMUsePackageLock (${chopper_SOURCE_DIR}/cmake/package-lock.cmake)

CPMGetPackage (use_ccache)
CPMGetPackage (hibf)
CPMGetPackage (sharg)
CPMGetPackage (seqan3)

# Allow to include CMake scripts from seqan3.
list (APPEND CMAKE_MODULE_PATH "${seqan3_SOURCE_DIR}/test/cmake/")
list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/test/cmake/")

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
