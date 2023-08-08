// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#pragma once

/*!\brief Workaround bogus memcpy errors in GCC 12. (Wrestrict and Wstringop-overflow)
 * \see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=105545
 */
#ifndef CHOPPER_WORKAROUND_GCC_BOGUS_MEMCPY
#    if defined(__GNUC__) && !defined(__llvm__) && !defined(__INTEL_COMPILER) && (__GNUC__ == 12)
#        define CHOPPER_WORKAROUND_GCC_BOGUS_MEMCPY 1
#    else
#        define CHOPPER_WORKAROUND_GCC_BOGUS_MEMCPY 0
#    endif
#endif
