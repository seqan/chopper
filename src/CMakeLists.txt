# ---------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
# ---------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.18)

if (NOT TARGET chopper_shared)
    add_library (chopper_shared STATIC configuration.cpp input_functor.cpp)

    target_link_libraries (chopper_shared PUBLIC seqan3::seqan3)
    target_link_libraries (chopper_shared PUBLIC sharg::sharg)
    target_link_libraries (chopper_shared PUBLIC seqan::hibf)
    target_include_directories (chopper_shared PUBLIC ../include)
    target_compile_options (chopper_shared PUBLIC "-pedantic" "-Wall" "-Wextra")
endif ()

add_library (chopper_lib STATIC set_up_parser.cpp)
target_link_libraries (chopper_lib PUBLIC chopper_layout chopper_sketch)

add_executable (chopper chopper.cpp)
target_link_libraries (chopper PUBLIC chopper_lib)

add_subdirectory (layout)
add_subdirectory (sketch)
add_subdirectory (util)

if (CHOPPER_INSTALL)
    install (TARGETS chopper RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
endif ()
