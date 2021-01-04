cmake_minimum_required (VERSION 3.8)

add_api_test (create_ibfs_test.cpp)
target_use_datasources (create_ibfs_test FILES small.fa)
target_use_datasources (create_ibfs_test FILES small2.fa)

add_api_test (parse_chopper_split_file_line_test.cpp)

add_api_test (read_chopper_split_file_test.cpp)

add_api_test (chopper_build_test.cpp)
target_use_datasources (chopper_build_test FILES small.fa)
target_use_datasources (chopper_build_test FILES small2.fa)

