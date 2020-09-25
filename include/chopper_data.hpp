#pragma once

#include <seqan/sequence.h>

#include "minimizer.hpp"

struct chopper_data
{
    typedef seqan::String<minimizer> TSequence;
    seqan::StringSet<TSequence, seqan::Owner<> > sequences;
    seqan::StringSet<seqan::String<char> > ids;
    seqan::String<size_t> lengths;
};
