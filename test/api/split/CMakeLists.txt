cmake_minimum_required (VERSION 3.8)

add_api_test (minimizer_msa_test.cpp)
target_use_datasources (minimizer_msa_test FILES small_graph.dot)
target_use_datasources (minimizer_msa_test FILES small.fa)

add_api_test (chopper_split_test.cpp)
target_use_datasources (chopper_split_test FILES small.fa)
target_use_datasources (chopper_split_test FILES small_traverse.out)
