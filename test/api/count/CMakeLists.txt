cmake_minimum_required (VERSION 3.8)

add_api_test(read_data_file_test.cpp)
target_use_datasources (read_data_file_test FILES seqinfo.tsv)

add_api_test(count_kmers_test.cpp)
target_use_datasources (count_kmers_test FILES small.fa)

add_api_test(chopper_count_test.cpp)
target_use_datasources (chopper_count_test FILES small.fa)
target_use_datasources (read_data_file_test FILES seqinfo.tsv)
