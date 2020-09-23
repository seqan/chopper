cmake_minimum_required (VERSION 3.8)

add_executable ("${PROJECT_NAME}" chopper.cpp)
target_link_libraries ("${PROJECT_NAME}" PUBLIC seqan3::seqan3)
target_include_directories ("${PROJECT_NAME}" PUBLIC ../include)
target_compile_definitions ("${PROJECT_NAME}" PUBLIC ${SEQAN_DEFINITIONS})
target_compile_options ("${PROJECT_NAME}" PUBLIC ${SEQAN_CXX_FLAGS_LIST})
target_link_libraries ("${PROJECT_NAME}" PUBLIC "${SEQAN_LIBRARIES}")
target_include_directories ("${PROJECT_NAME}" PUBLIC "${SEQAN_INCLUDE_DIRS}")

