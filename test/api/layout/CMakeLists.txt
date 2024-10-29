# ---------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
# ---------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.18)

add_api_test (ibf_query_cost_test.cpp)
add_api_test (execute_layout_test.cpp)

add_api_test (execute_with_estimation_test.cpp)
target_use_datasources (execute_with_estimation_test FILES seq1.fa)
target_use_datasources (execute_with_estimation_test FILES seq2.fa)
target_use_datasources (execute_with_estimation_test FILES seq3.fa)
target_use_datasources (execute_with_estimation_test FILES small.fa)

add_api_test (hibf_statistics_test.cpp)
# Bogus LTO warning on gcc-12
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"
    AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 12
    AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 13
)
    target_link_options (hibf_statistics_test PRIVATE -Wno-stringop-overflow)
endif ()
add_api_test (user_bin_io_test.cpp)
add_api_test (input_test.cpp)
