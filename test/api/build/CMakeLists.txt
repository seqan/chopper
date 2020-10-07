cmake_minimum_required (VERSION 3.8)

add_api_test (create_ibfs_from_data_file.cpp)
target_use_datasources (create_ibfs_from_data_file FILES small.fa)
target_use_datasources (create_ibfs_from_data_file FILES small2.fa)

add_api_test (parse_traversal_file_line_test.cpp)

add_api_test (chopper_build_test.cpp)
target_use_datasources (chopper_build_test FILES small.fa)
target_use_datasources (chopper_build_test FILES small2.fa)
