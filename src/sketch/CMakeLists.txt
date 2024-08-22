cmake_minimum_required (VERSION 3.18)

if (TARGET chopper::sketch)
    return ()
endif ()

add_library (chopper_sketch STATIC check_filenames.cpp output.cpp read_data_file.cpp)
target_link_libraries (chopper_sketch PUBLIC chopper::shared)
add_library (chopper::sketch ALIAS chopper_sketch)
