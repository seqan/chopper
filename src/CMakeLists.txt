cmake_minimum_required (VERSION 3.11)

add_library (chopper_interface INTERFACE)
target_link_libraries (chopper_interface INTERFACE seqan3::seqan3)
target_include_directories (chopper_interface INTERFACE ../include)
target_include_directories (chopper_interface INTERFACE ../lib/lemon/include)
target_include_directories (chopper_interface INTERFACE ../lib/robin-hood-hashing/src/include)
target_include_directories (chopper_interface INTERFACE ../lib/simde/simde)
target_include_directories (chopper_interface INTERFACE "${SEQAN_INCLUDE_DIRS}")
target_compile_definitions (chopper_interface INTERFACE ${SEQAN_DEFINITIONS})
target_compile_options (chopper_interface INTERFACE ${SEQAN_CXX_FLAGS_LIST})
target_link_libraries (chopper_interface INTERFACE "${SEQAN_LIBRARIES}")
target_link_libraries (chopper_interface INTERFACE xxhash)
target_compile_options (chopper_interface INTERFACE "-pedantic" "-Wall" "-Wextra")
add_dependencies (chopper_interface xxhash)

add_library (chopper_count_lib STATIC chopper_count.cpp)
target_link_libraries (chopper_count_lib "chopper_interface")

add_library (chopper_layout_lib STATIC chopper_layout.cpp)
target_link_libraries (chopper_layout_lib "chopper_interface")

add_library (chopper_lib INTERFACE)
target_link_libraries (chopper_lib INTERFACE "chopper_interface" "chopper_count_lib" "chopper_layout_lib")

add_executable (chopper chopper.cpp)
target_link_libraries (chopper "chopper_lib")

add_executable (measure_hyperloglog measure_hyperloglog.cpp)
target_link_libraries (measure_hyperloglog "chopper_interface")
