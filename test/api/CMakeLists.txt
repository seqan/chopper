include(${CMAKE_SOURCE_DIR}/lib/seqan3/test/cmake/add_subdirectories.cmake)

add_api_test (detail_hibf_user_bins_test.cpp)
add_api_test (detail_parse_chopper_pack_header_line_test.cpp)
add_api_test (detail_parse_chopper_pack_line_test.cpp)

add_subdirectories ()
