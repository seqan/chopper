cmake_minimum_required (VERSION 3.8)

add_cli_test (cli_chopper_split_test.cpp)
target_use_datasources (cli_chopper_split_test FILES small.fa)
target_use_datasources (cli_chopper_split_test FILES small2.fa)

add_cli_test (cli_count_kmers_per_bin_test.cpp)
target_use_datasources (cli_count_kmers_per_bin_test FILES small.fa)
target_use_datasources (cli_count_kmers_per_bin_test FILES small2.fa)
add_dependencies (cli_count_kmers_per_bin_test count_kmers_per_bin)
