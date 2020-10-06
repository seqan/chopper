cmake_minimum_required (VERSION 3.8)

add_cli_test (cli_chopper_split_test.cpp)
target_use_datasources (cli_chopper_split_test FILES small.fa)
target_use_datasources (cli_chopper_split_test FILES small2.fa)
