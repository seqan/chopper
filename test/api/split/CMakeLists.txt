cmake_minimum_required (VERSION 3.8)

add_api_test (chopper_split_test.cpp)
target_use_datasources (chopper_split_test FILES small.fa)
target_use_datasources (chopper_split_test FILES small2.fa)

add_api_test (filename_batches_range_test.cpp)

add_api_test (neighbour_joining_test.cpp)
add_api_test (map_distance_matrix_test.cpp)
