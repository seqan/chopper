cmake_minimum_required (VERSION 3.11)

## CUSTOMISE

# Define the application name and version.
project (chopper VERSION 1.0.0)

## BUILD

# Make Release default build type
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING
         "Choose the type of build, options are: Debug Release Coverage RelWithDebInfo MinSizeRel."
         FORCE)
endif ()

set (CHOPPER_SUBMODULES_DIR "${CMAKE_CURRENT_LIST_DIR}/lib" CACHE STRING "Directory containing submodules.")

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
set (CMAKE_INSTALL_LIBDIR "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}" CACHE STRING "")
set (CMAKE_INSTALL_BINDIR "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}" CACHE STRING "")

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

# Dependency: xxHash
set (CHOPPER_XXHASH_PATH "${PROJECT_BINARY_DIR}/lib/libxxhash.a")

if (EXISTS "${CHOPPER_XXHASH_PATH}")
    add_custom_target (xxhash_project)
else ()
    set (CHOPPER_XXHASH_ARGS "${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS}")
    list (APPEND CHOPPER_XXHASH_ARGS "-DXXHASH_BUILD_XXHSUM=OFF")
    list (APPEND CHOPPER_XXHASH_ARGS "-DXXHASH_BUNDLED_MODE=OFF")
    list (APPEND CHOPPER_XXHASH_ARGS "-DBUILD_SHARED_LIBS=OFF")
    list (APPEND CHOPPER_XXHASH_ARGS "--no-warn-unused-cli")

    include (ExternalProject)
    ExternalProject_Add (
        xxhash_project
        PREFIX xxhash_project
        GIT_REPOSITORY https://github.com/Cyan4973/xxHash.git
        GIT_TAG "f9155bd4c57e2270a4ffbb176485e5d713de1c9b"
        CMAKE_ARGS "${CHOPPER_XXHASH_ARGS}"
        SOURCE_SUBDIR cmake_unofficial
        BUILD_BYPRODUCTS "${CHOPPER_XXHASH_PATH}"
        UPDATE_DISCONNECTED 1
    )
endif ()

add_library (xxhash STATIC IMPORTED)
add_dependencies (xxhash xxhash_project)
set_target_properties (xxhash PROPERTIES IMPORTED_LOCATION "${CHOPPER_XXHASH_PATH}")
target_include_directories (xxhash INTERFACE "${PROJECT_BINARY_DIR}/include")

# Use ccache.
set (CHOPPER_USE_CCACHE ON CACHE BOOL "Use ccache if available.")
if (CHOPPER_USE_CCACHE)
     include (seqan3_require_ccache)
     seqan3_require_ccache ()
endif ()

# Add the application.
add_subdirectory (src)
message (STATUS "${FontBold}You can run `make` to build the application.${FontReset}")

## DOCUMENTATION

add_subdirectory (doc EXCLUDE_FROM_ALL)

## TEST

enable_testing ()
add_subdirectory (test EXCLUDE_FROM_ALL)
