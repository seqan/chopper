cmake_minimum_required (VERSION 3.18)

if (TARGET chopper::layout)
    return ()
endif ()

add_library (chopper_layout STATIC
             determine_best_number_of_technical_bins.cpp
             determine_split_bins.cpp
             execute.cpp
             hibf_statistics.cpp
             ibf_query_cost.cpp
             input.cpp
             output.cpp
)
target_link_libraries (chopper_layout PUBLIC chopper::shared)
add_library (chopper::layout ALIAS chopper_layout)
