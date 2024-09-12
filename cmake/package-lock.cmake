# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

# CPM Package Lock
# This file should be committed to version control

# hibf
set (CHOPPER_HIBF_VERSION 169cafe1fda30f2b0546cababd5740b282c38aa6)
CPMDeclarePackage (hibf
                   NAME hibf
                   GIT_TAG ${CHOPPER_HIBF_VERSION}
                   GITHUB_REPOSITORY seqan/hibf
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_HIBF OFF"
)
# sharg
set (CHOPPER_SHARG_VERSION 420193a9642ab34c87c21a8b9ba1dd3e23d2c83a)
CPMDeclarePackage (sharg
                   NAME sharg
                   GIT_TAG ${CHOPPER_SHARG_VERSION}
                   GITHUB_REPOSITORY seqan/sharg-parser
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_SHARG OFF" "INSTALL_TDL OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING" "BUILD_TESTING OFF"
)
# seqan3
set (CHOPPER_SEQAN3_VERSION 1fcf9a66a5db46785dc90bc2d1fc6b4acca94fd9)
CPMDeclarePackage (seqan3
                   NAME seqan3
                   GIT_TAG ${CHOPPER_SEQAN3_VERSION}
                   GITHUB_REPOSITORY seqan/seqan3
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_SEQAN3 OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)
# googletest
set (CHOPPER_GOOGLETEST_VERSION 1.15.2)
CPMDeclarePackage (googletest
                   NAME googletest
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
                   GIT_TAG ${USE_CCACHE_VERSION}
                   GITHUB_REPOSITORY seqan/cmake-scripts
                   SOURCE_SUBDIR ccache
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
)
