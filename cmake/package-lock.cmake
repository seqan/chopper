# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

# CPM Package Lock
# This file should be committed to version control

# hibf
set (CHOPPER_HIBF_VERSION e61e1e88f018d70a33752a204bc44dd20282d554)
CPMDeclarePackage (hibf
                   NAME hibf
                   GIT_TAG ${CHOPPER_HIBF_VERSION} # main
                   GITHUB_REPOSITORY seqan/hibf
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_HIBF OFF"
)
# sharg
set (CHOPPER_SHARG_VERSION 215da66217957776691b7f8a572cf5a51d368b03)
CPMDeclarePackage (sharg
                   NAME sharg
                   GIT_TAG ${CHOPPER_SHARG_VERSION} # main
                   GITHUB_REPOSITORY seqan/sharg-parser
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_SHARG OFF" "INSTALL_TDL OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING" "BUILD_TESTING OFF"
)
# seqan3
set (CHOPPER_SEQAN3_VERSION c5a8eb5b7ac726a29706da69433d2faf3d0e9065)
CPMDeclarePackage (seqan3
                   NAME seqan3
                   GIT_TAG ${CHOPPER_SEQAN3_VERSION} # main
                   GITHUB_REPOSITORY seqan/seqan3
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_SEQAN3 OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)
# googletest
set (CHOPPER_GOOGLETEST_VERSION 1.15.2)
CPMDeclarePackage (googletest
                   NAME GTest
                   VERSION ${CHOPPER_GOOGLETEST_VERSION}
                   GITHUB_REPOSITORY google/googletest
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "BUILD_GMOCK OFF" "INSTALL_GTEST OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
                           "CMAKE_CXX_STANDARD 20"
)
# use_ccache
set (USE_CCACHE_VERSION d2a54ef555b6fc2d496a4c9506dbeb7cf899ce37)
CPMDeclarePackage (use_ccache
                   NAME use_ccache
                   GIT_TAG ${USE_CCACHE_VERSION} # main
                   GITHUB_REPOSITORY seqan/cmake-scripts
                   SOURCE_SUBDIR ccache
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
)
