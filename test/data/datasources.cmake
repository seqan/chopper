cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

declare_datasource (FILE small_graph.dot
                    URL ${CMAKE_SOURCE_DIR}/test/data/small_graph.dot
                    URL_HASH SHA256=17b4ea76944561f781221d6d30421b5074624a2c057cea1b924743bbfabb16a7)
declare_datasource (FILE small.fa
                    URL ${CMAKE_SOURCE_DIR}/test/data/small.fa
                    URL_HASH SHA256=215a646ccd2d156eeabdd792108fc0ba385e244a25218afc095658189d2a3c0d)

declare_datasource (FILE small_traverse.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/small_traverse.out
                    URL_HASH SHA256=9cbdd4cf4f78c35bdece7868b77c03406317255231224e88d75dc06192f8bde8)

declare_datasource (FILE only_filenames.tsv
                    URL ${CMAKE_SOURCE_DIR}/test/data/only_filenames.tsv
                    URL_HASH SHA256=581870b9c418bac44191ea967376c5b81d42a14d0d8aaf4e17dc700e146e714b)
declare_datasource (FILE filenames_and_counts.tsv
                    URL ${CMAKE_SOURCE_DIR}/test/data/filenames_and_counts.tsv
                    URL_HASH SHA256=4eedd49f5fa13fc981d8dab8c72a859c1814c3d283f030b0a81913f4916b611c)
declare_datasource (FILE filenames_counts_and_extra_information.tsv
                    URL ${CMAKE_SOURCE_DIR}/test/data/filenames_counts_and_extra_information.tsv
                    URL_HASH SHA256=2df6817a35aa73b75e1c63fdaa09dcb554f5bedddff59404ab48fdeb48e32649)

declare_datasource (FILE high_level_ibf.binning
                    URL ${CMAKE_SOURCE_DIR}/test/data/high_level_ibf.binning
                    URL_HASH SHA256=0d198d13f96cc7d4c6984b795ea5a39f271124278540584c68519a07a5cfffe7)
declare_datasource (FILE low_level_ibfs.binning
                    URL ${CMAKE_SOURCE_DIR}/test/data/low_level_ibfs.binning
                    URL_HASH SHA256=1bf3530d47ae044eef071c725e9e23ec2cec139222aad3ee3bade8ba001b46fb)
