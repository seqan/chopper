cmake_minimum_required (VERSION 3.8)

add_api_test (minimizer_msa_test.cpp)
target_use_datasources (minimizer_msa_test FILES small_graph.dot)
target_use_datasources (minimizer_msa_test FILES small.fa)

add_api_test (chopper_split_test.cpp)
target_use_datasources (chopper_split_test FILES small.fa)
target_use_datasources (chopper_split_test FILES small2.fa)

add_api_test (filename_batches_range_test.cpp)
target_use_datasources (filename_batches_range_test FILES high_level_ibf.binning)
target_use_datasources (filename_batches_range_test FILES low_level_ibfs.binning)
