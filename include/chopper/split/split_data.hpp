#pragma once

#include <seqan/sequence.h>

#include <chopper/split/minimizer.hpp>

struct split_data
{
    typedef seqan::String<minimizer> TSequence;
    seqan::StringSet<TSequence, seqan::Owner<> > sequences;
    seqan::StringSet<seqan::String<char> > ids;
    seqan::String<size_t> lengths;
    std::vector<std::string> files_of_origin;

    std::ofstream * outstream{nullptr};
};
