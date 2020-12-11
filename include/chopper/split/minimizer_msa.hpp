#pragma once

#include <cstdlib>

#include <seqan/basic.h>
#include <seqan/graph_msa.h>
#include <seqan/modifier.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>

#include <seqan3/std/ranges>

#include <chopper/split/split_config.hpp>
#include <chopper/split/split_data.hpp>
#include <chopper/split/map_distance_matrix.hpp>
#include <chopper/split/distance_matrix_initialiser.hpp>
#include <chopper/split/minimizer.hpp>
#include <chopper/split/seqan2_msa_alignment.hpp>

inline auto minimizer_msa(split_data & data, batch_config const & config)
{
    // Alignment of the sequences
    typedef seqan::Graph<seqan::Alignment<seqan::StringSet<seqan::String<minimizer>, seqan::Dependent<> >, void, seqan::WithoutEdgeId> > TGraph;
    TGraph gAlign;

    distance_matrix_initialiser initialiser{};
    auto distance_matrix = initialiser.mash_distance(data, config);

    // MSA
    try
    {
        seqan2_msa_alignment(gAlign, data.sequences, distance_matrix, config);
    }
    catch (const std::bad_alloc & exception)
    {
        std::cerr << "Allocation for globalAlignment failed. Use smaller data or try a seeded alignment. \n"
                  << exception.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    return gAlign;
}
