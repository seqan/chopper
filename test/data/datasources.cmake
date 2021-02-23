cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

declare_datasource (FILE small_graph.dot
                    URL ${CMAKE_SOURCE_DIR}/test/data/small_graph.dot
                    URL_HASH SHA256=f573bb45abe677bda78d5f700bd6a9b5f2415f405953f5c085be491e2598f02e)

declare_datasource (FILE small.fa
                    URL ${CMAKE_SOURCE_DIR}/test/data/small.fa
                    URL_HASH SHA256=215a646ccd2d156eeabdd792108fc0ba385e244a25218afc095658189d2a3c0d)

declare_datasource (FILE seq1.fa
                    URL ${CMAKE_SOURCE_DIR}/test/data/seq1.fa
                    URL_HASH SHA256=f60b985d9f3be1d3dd6105d62f4321eeb617bc52f4fb0435b016bf1873aa2eac)

declare_datasource (FILE seq2.fa
                    URL ${CMAKE_SOURCE_DIR}/test/data/seq2.fa
                    URL_HASH SHA256=4153dfac5993b6ddeebf76fabeaaaad838b68f4abe0c3c03cc95b9196832af7a)

declare_datasource (FILE seq3.fa
                    URL ${CMAKE_SOURCE_DIR}/test/data/seq3.fa
                    URL_HASH SHA256=aa440a61539216096680e53a1d68b245e5a7cac67e7eca315996d51b863d9915)

declare_datasource (FILE small2.fa
                    URL ${CMAKE_SOURCE_DIR}/test/data/small2.fa
                    URL_HASH SHA256=3a444c41f71586d62514eea8faf3c699ac5f2cfce94fd3fc7c3e47d3ff3a8383)

declare_datasource (FILE only_filenames.tsv
                    URL ${CMAKE_SOURCE_DIR}/test/data/only_filenames.tsv
                    URL_HASH SHA256=581870b9c418bac44191ea967376c5b81d42a14d0d8aaf4e17dc700e146e714b)

declare_datasource (FILE filenames_and_counts.tsv
                    URL ${CMAKE_SOURCE_DIR}/test/data/filenames_and_counts.tsv
                    URL_HASH SHA256=4eedd49f5fa13fc981d8dab8c72a859c1814c3d283f030b0a81913f4916b611c)

declare_datasource (FILE filenames_counts_and_extra_information.tsv
                    URL ${CMAKE_SOURCE_DIR}/test/data/filenames_counts_and_extra_information.tsv
                    URL_HASH SHA256=2df6817a35aa73b75e1c63fdaa09dcb554f5bedddff59404ab48fdeb48e32649)

declare_datasource (FILE seqinfo.tsv
                    URL ${CMAKE_SOURCE_DIR}/test/data/seqinfo.tsv
                    URL_HASH SHA256=c6e28cc4ebf4902c41b1c237b1410665994e0fd0b8d8473c38a8107041d78172)

declare_datasource (FILE small.split
                    URL ${CMAKE_SOURCE_DIR}/test/data/small.split
                    URL_HASH SHA256=995ed2a9ff71dbf7b04a443b643d9268e370aa346cc5a5afa0b7176fbec7d9d0)
