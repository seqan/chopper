# ---------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
# ---------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.18)

if (TARGET chopper::interface)
    return ()
endif ()

add_library (chopper_interface INTERFACE)
target_link_libraries (chopper_interface INTERFACE seqan::hibf sharg::sharg seqan3::seqan3)
target_include_directories (chopper_interface INTERFACE "${PROJECT_SOURCE_DIR}/include")
target_compile_options (chopper_interface INTERFACE "-pedantic" "-Wall" "-Wextra")
add_library (chopper::interface ALIAS chopper_interface)

add_library (chopper_shared STATIC configuration.cpp input_functor.cpp)
target_link_libraries (chopper_shared PUBLIC chopper_interface)
add_library (chopper::shared ALIAS chopper_shared)

add_library (chopper_lib STATIC chopper_layout.cpp set_up_parser.cpp)
target_link_libraries (chopper_lib PUBLIC chopper::layout chopper::sketch)
add_library (chopper::chopper ALIAS chopper_lib)

add_executable (chopper chopper.cpp)
target_link_libraries (chopper PUBLIC chopper::chopper)

add_subdirectory (layout)
add_subdirectory (sketch)
add_subdirectory (util)

if (CHOPPER_INSTALL)
    include (GNUInstallDirs)
    install (TARGETS chopper RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
endif ()
