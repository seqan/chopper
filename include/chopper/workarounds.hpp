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
