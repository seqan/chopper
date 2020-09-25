#pragma once

#include <seqan/seq_io.h>

#include "minimizer.hpp"
#include "chopper_config.hpp"

template <typename TNameSet>
bool load_minimizer_sequences(seqan::StringSet<seqan::String<minimizer>, seqan::Owner<>> & minimizer_sequences,
                              TNameSet& fastaIDs,
                              seqan::String<size_t> & original_sequence_lengths,
                              chopper_config const & config,
                              const char *fileName)
{
    seqan::SeqFileIn inFile;
    if (!open(inFile, fileName))
    {
        std::cerr << "Could not open " << fileName << " for reading!" << std::endl;
        return false;
    }

    std::cerr << ">>> Processing file " << fileName << std::endl;
    seqan::StringSet<seqan::String<seqan::Dna>, seqan::Owner<>> sequences;
    {
        seqan::StringSet<seqan::String<seqan::Iupac>, seqan::Owner<>> iupac_sequences;
        seqan::readRecords(fastaIDs, iupac_sequences, inFile);
        seqan::resize(sequences, seqan::length(iupac_sequences));

        for (size_t idx = 0; idx < seqan::length(iupac_sequences); ++idx)
        {
            for (size_t jdx = 0; jdx < seqan::length(iupac_sequences[idx]); ++jdx)
            {
                seqan::appendValue(sequences[idx], iupac_sequences[idx][jdx]);
            }
            seqan::appendValue(original_sequence_lengths, seqan::length(iupac_sequences[idx]));
        }
    }

    // compute minimizers per sequence and store the corresponding chain in minimizer_sequences
    auto from = seqan::length(minimizer_sequences);
    resize(minimizer_sequences, seqan::length(minimizer_sequences) + seqan::length(sequences));
    Minimizer mini;
    mini.resize(config.kmer_size, config.window_size);
    for (size_t idx = from; idx < seqan::length(minimizer_sequences); ++idx)
    {
        minimizer_sequences[idx] = mini.getMinimizer(sequences[idx - from]);
        // seqan3::debug_stream << seqan::length(minimizer_sequences[idx]) << std::endl;
    }

    return (seqan::length(fastaIDs) > 0u);
}
