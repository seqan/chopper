cmake_minimum_required (VERSION 3.8)

add_api_test (search_test.cpp)
target_use_datasources (search_test FILES seq1.fa)
target_use_datasources (search_test FILES seq2.fa)
target_use_datasources (search_test FILES seq3.fa)

add_api_test (chopper_search_test.cpp)
