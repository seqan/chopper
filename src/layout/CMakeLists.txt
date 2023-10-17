cmake_minimum_required (VERSION 3.18)

if (NOT TARGET chopper_layout)
    add_library (chopper_layout STATIC determine_best_number_of_technical_bins.cpp execute.cpp hibf_statistics.cpp
                                       ibf_query_cost.cpp input.cpp output.cpp
    )

    target_link_libraries (chopper_layout PUBLIC chopper_shared)
endif ()
