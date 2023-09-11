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
         LANGUAGES CXX C
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

# Check CXX and C compiler versions (xxHash uses C).
string (REGEX MATCH "^[0-9]+" CHOPPER_CXX_MAJOR_VERSION "${CMAKE_CXX_COMPILER_VERSION}")
string (REGEX MATCH "^[0-9]+" CHOPPER_C_MAJOR_VERSION "${CMAKE_C_COMPILER_VERSION}")
if (NOT "${CHOPPER_CXX_MAJOR_VERSION}" STREQUAL "${CHOPPER_C_MAJOR_VERSION}")
    message (FATAL_ERROR "CXX and C compiler major versions differ (${CHOPPER_CXX_MAJOR_VERSION} and\
                          ${CHOPPER_C_MAJOR_VERSION})! This will likely result in linker errors. Please set both\
                          -DCMAKE_CXX_COMPILER and -DCMAKE_C_COMPILER to use the same major version, and/or set the\
                          environment variables CXX and CC accordingly."
    )
endif ()

set (CHOPPER_SUBMODULES_DIR
     "${CMAKE_CURRENT_LIST_DIR}/lib"
     CACHE STRING "Directory containing submodules."
)

# Flag checks
include (CheckCXXCompilerFlag)

check_cxx_compiler_flag ("-fopenmp" CHOPPER_HAS_OPENMP)
if (CHOPPER_HAS_OPENMP)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif ()

check_cxx_compiler_flag ("-fopenmp-simd" CHOPPER_HAS_OPENMP_SIMD)
if (CHOPPER_HAS_OPENMP_SIMD)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp-simd -DSIMDE_ENABLE_OPENMP")
endif ()

check_cxx_compiler_flag ("-fopenmp-simd" CHOPPER_SUPPRESS_GCC4_ABI)
if (CHOPPER_SUPPRESS_GCC4_ABI)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-psabi")
endif ()

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

# Dependency: xxHash
if (NOT TARGET xxHash::xxhash)
    FetchContent_Declare (chopper_xxhash_project SOURCE_DIR ${CHOPPER_SUBMODULES_DIR}/xxHash SOURCE_SUBDIR
                                                 "cmake_unofficial"
    )

    set (XXHASH_BUILD_ENABLE_INLINE_API ON)
    set (XXHASH_BUILD_XXHSUM OFF)
    set (XXHASH_BUNDLED_MODE ON)
    set (BUILD_SHARED_LIBS OFF)

    FetchContent_MakeAvailable (chopper_xxhash_project)
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
