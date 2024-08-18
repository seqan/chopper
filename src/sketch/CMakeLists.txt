cmake_minimum_required (VERSION 3.18)

if (NOT TARGET chopper_sketch)
    add_library (chopper_sketch STATIC check_filenames.cpp output.cpp read_data_file.cpp)

    target_link_libraries (chopper_sketch PUBLIC chopper_shared)
endif ()
