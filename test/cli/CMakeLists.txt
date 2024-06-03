# ---------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
# ---------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.18)

add_cli_test (cli_chopper_basic_test.cpp)

add_cli_test (cli_chopper_pipeline_test.cpp)
target_use_datasources (cli_chopper_pipeline_test FILES seq1.fa)
target_use_datasources (cli_chopper_pipeline_test FILES small.fa)
target_use_datasources (cli_chopper_pipeline_test FILES small2.fa)
target_use_datasources (cli_chopper_pipeline_test FILES small.split)

add_cli_test (cli_timing_output_test.cpp)
target_use_datasources (cli_timing_output_test FILES small.fa)

add_cli_test (util_display_layout_test.cpp)
target_use_datasources (util_display_layout_test FILES seq1.fa)
target_use_datasources (util_display_layout_test FILES seq2.fa)
target_use_datasources (util_display_layout_test FILES seq3.fa)
target_use_datasources (util_display_layout_test FILES small.fa)
add_dependencies (util_display_layout_test display_layout)
