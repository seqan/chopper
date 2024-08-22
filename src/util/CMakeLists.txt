# ---------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
# ---------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.18)

if (TARGET measure_hyperloglog)
    return ()
endif ()

add_executable (measure_hyperloglog EXCLUDE_FROM_ALL measure_hyperloglog.cpp)
target_link_libraries (measure_hyperloglog PUBLIC chopper::chopper)

add_subdirectory (display_layout)
