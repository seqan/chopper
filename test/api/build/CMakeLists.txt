cmake_minimum_required (VERSION 3.8)

# add_api_test (create_ibfs_from_chopper_split_test.cpp)
# target_use_datasources (create_ibfs_from_chopper_split_test FILES small.fa)
# target_use_datasources (create_ibfs_from_chopper_split_test FILES small2.fa)

add_api_test (create_ibfs_from_chopper_pack_test.cpp)
target_use_datasources (create_ibfs_from_chopper_pack_test FILES seq1.fa)
target_use_datasources (create_ibfs_from_chopper_pack_test FILES seq2.fa)
target_use_datasources (create_ibfs_from_chopper_pack_test FILES seq3.fa)

# add_api_test (read_chopper_split_file_test.cpp)

add_api_test (read_chopper_pack_file_test.cpp)

add_api_test (chopper_build_test.cpp)
target_use_datasources (chopper_build_test FILES small.fa)
target_use_datasources (chopper_build_test FILES small2.fa)

