cmake_minimum_required (VERSION 3.18)

add_api_test (toolbox_test.cpp)
target_use_datasources (toolbox_test FILES small.fa)
target_use_datasources (toolbox_test FILES small2.fa)
target_use_datasources (toolbox_test FILES small.hll)
target_use_datasources (toolbox_test FILES small2.hll)

add_api_test (compute_sketches_test.cpp)
target_use_datasources (compute_sketches_test FILES small.fa)

add_api_test (check_filenames_test.cpp)
target_use_datasources (check_filenames_test FILES seq1.fa)
target_use_datasources (check_filenames_test FILES seq2.fa)
target_use_datasources (check_filenames_test FILES seq3.fa)
target_use_datasources (check_filenames_test FILES small.minimiser)

add_api_test (read_data_file_test.cpp)
target_use_datasources (read_data_file_test FILES seqinfo.tsv)
