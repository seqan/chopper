cmake_minimum_required (VERSION 3.15)

add_api_test (check_filenames_test.cpp)

add_api_test (detail_apply_prefix_test.cpp)

add_api_test (read_data_file_test.cpp)
target_use_datasources (read_data_file_test FILES seqinfo.tsv)

add_api_test (count_kmers_test.cpp)
target_use_datasources (count_kmers_test FILES small.fa)
target_use_datasources (count_kmers_test FILES small.minimiser)

add_api_test (execute_test.cpp)
target_use_datasources (execute_test FILES small.fa)
target_use_datasources (read_data_file_test FILES seqinfo.tsv)
