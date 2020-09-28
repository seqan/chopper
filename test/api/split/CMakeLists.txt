cmake_minimum_required (VERSION 3.8)

add_api_test (minimizer_msa_test.cpp)
target_use_datasources (minimizer_msa_test FILES small_graph.dot)
target_use_datasources (minimizer_msa_test FILES small.fa)
