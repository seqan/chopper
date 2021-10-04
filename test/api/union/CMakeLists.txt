cmake_minimum_required (VERSION 3.11)

add_api_test (hyperloglog_test.cpp)
target_use_datasources (hyperloglog_test FILES small.fa)
