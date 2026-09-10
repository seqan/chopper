# ---------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
# ---------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.18)

include (cmake/app_datasources.cmake)

declare_datasource (FILE small_graph.dot URL ${CMAKE_SOURCE_DIR}/test/data/small_graph.dot URL_HASH
                    SHA256=f573bb45abe677bda78d5f700bd6a9b5f2415f405953f5c085be491e2598f02e
)

declare_datasource (FILE small.fa URL ${CMAKE_SOURCE_DIR}/test/data/small.fa URL_HASH
                    SHA256=41d7eace7e07335afe6aa2c3b313eacae828e0f42a876fe97bf74a4cef4323f9
)

declare_datasource (FILE seq1.fa URL ${CMAKE_SOURCE_DIR}/test/data/seq1.fa URL_HASH
                    SHA256=752f6ff662e9c325886aab84db998c63fc89d7ff2f4504ccd686fd853e67a003
)

declare_datasource (FILE seq2.fa URL ${CMAKE_SOURCE_DIR}/test/data/seq2.fa URL_HASH
                    SHA256=c9c402bcb982a7f3ecb0ef8810b8d4ae3050553faa295a6c3e538713199564dd
)

declare_datasource (FILE seq3.fa URL ${CMAKE_SOURCE_DIR}/test/data/seq3.fa URL_HASH
                    SHA256=24a275e4666f815e46d1ccfb135bce5ea794f6630c8ffd9117eb0dcecce9d414
)

declare_datasource (FILE seqinfo.tsv URL ${CMAKE_SOURCE_DIR}/test/data/seqinfo.tsv URL_HASH
                    SHA256=c6e28cc4ebf4902c41b1c237b1410665994e0fd0b8d8473c38a8107041d78172
)

declare_datasource (FILE small.minimiser URL ${CMAKE_SOURCE_DIR}/test/data/small.minimiser URL_HASH
                    SHA256=6262de00ad97113320469ff952e657930d07151bda5174ee3599560eccb6f0e1
)
