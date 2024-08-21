# ---------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
# ---------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.18)

add_api_test (check_filenames_test.cpp)
target_use_datasources (check_filenames_test FILES seq1.fa)
target_use_datasources (check_filenames_test FILES seq2.fa)
target_use_datasources (check_filenames_test FILES seq3.fa)
target_use_datasources (check_filenames_test FILES small.minimiser)

add_api_test (read_data_file_test.cpp)
target_use_datasources (read_data_file_test FILES seqinfo.tsv)
