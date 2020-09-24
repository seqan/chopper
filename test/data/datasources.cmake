cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

# copies file to <build>/data/<filename>
declare_datasource (FILE small_graph.dot
                    URL ${CMAKE_SOURCE_DIR}/test/data/small_graph.dot
                    URL_HASH SHA256=17b4ea76944561f781221d6d30421b5074624a2c057cea1b924743bbfabb16a7)
declare_datasource (FILE small.fa
                    URL ${CMAKE_SOURCE_DIR}/test/data/small.fa
		    URL_HASH SHA256=215a646ccd2d156eeabdd792108fc0ba385e244a25218afc095658189d2a3c0d)
