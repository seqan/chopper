// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#pragma once

#if defined(__GNUC__) && !defined(__llvm__) && !defined(__INTEL_COMPILER)
#    define CHOPPER_COMPILER_IS_GCC 1
#else
#    define CHOPPER_COMPILER_IS_GCC 0
#endif

/*!\brief Workaround bogus memcpy errors in GCC 12. (Wrestrict and Wstringop-overflow)
 * \see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=105545
 */
#ifndef CHOPPER_WORKAROUND_GCC_BOGUS_MEMCPY
#    if CHOPPER_COMPILER_IS_GCC && (__GNUC__ == 12)
#        define CHOPPER_WORKAROUND_GCC_BOGUS_MEMCPY 1
#    else
#        define CHOPPER_WORKAROUND_GCC_BOGUS_MEMCPY 0
#    endif
#endif

/*!\brief Workaround bogus memmov errors in GCC 13. (Warray-bounds and Wstringop-overflow)
 */
#ifndef CHOPPER_WORKAROUND_GCC_BOGUS_MEMMOV
#    if CHOPPER_COMPILER_IS_GCC && (__GNUC__ == 13)
#        define CHOPPER_WORKAROUND_GCC_BOGUS_MEMMOV 1
#    else
#        define CHOPPER_WORKAROUND_GCC_BOGUS_MEMMOV 0
#    endif
#endif
