#pragma once

#include <seqan/graph_msa.h>

#include <chopper/split/split_config.hpp>
#include <chopper/split/map_distance_matrix.hpp>
#include <chopper/split/neighbour_joining.hpp>

template<typename TString, typename TSpec, typename TSegmentMatches, typename TScoreValues>
inline void append_all_to_all_matches(seqan::StringSet<TString, seqan::Dependent<TSpec> > const & sequenceSet,
                                      TSegmentMatches& matches,
                                      TScoreValues& scores)
{
    typedef seqan::StringSet<TString, seqan::Dependent<TSpec> > TStringSet;
    typedef typename seqan::Size<TStringSet>::Type TSize;

    // create hash_map with positions where minimizers can be found
    std::unordered_map<uint64_t, std::vector<std::pair<size_t, size_t>>> hash_map{};
    for (size_t i = 0; i < seqan::length(sequenceSet); ++i)
        for (size_t p1 = 0; p1 < seqan::length(sequenceSet[i]); ++p1)
            hash_map[sequenceSet[i][p1].value].emplace_back(i, p1);

    // Traverse hash map and create fragments
    for (auto & [minimizer, id_pos_pairs] : hash_map)
    {
        (void)minimizer;
        for (size_t i = 0; i < id_pos_pairs.size(); ++i)
        {
            for (size_t j = i + 1; j < id_pos_pairs.size(); ++j)
            {
                if (id_pos_pairs[i].first != id_pos_pairs[j].first) // do not include matches within the same seq
                {
                    TSize const from = seqan::length(matches);

                    seqan::appendValue(matches, seqan::Fragment<>(id_pos_pairs[i].first, id_pos_pairs[i].second,
                                                                  id_pos_pairs[j].first, id_pos_pairs[j].second, 1));

                    TSize const to = seqan::length(matches);
                    _recordScores(scores, 10, from, to); // will be rescored anyway ?!?!?!?!
                }
            }
        }
    }
}

template <typename TStringSet, typename TCargo, typename TSpec, typename TStringSet1>
void seqan2_msa_alignment(seqan::Graph<seqan::Alignment<TStringSet, TCargo, TSpec> > & gAlign,
                          TStringSet1 const & sequenceSet,
                          map_distance_matrix & distance_matrix)
{
    typedef seqan::Score<int> TScore;
    typedef typename seqan::Value<TScore>::Type TScoreValue;
    typedef typename seqan::Size<TStringSet>::Type TSize;
    typedef seqan::Graph<seqan::Alignment<TStringSet, TSize> > TGraph;
    typedef seqan::String<seqan::Fragment<> > TFragmentString;
    typedef seqan::String<TScoreValue> TScoreValues;

    std::cout << "Start" << std::endl;
    double current_time = seqan::sysTime();

    // Initialize alignment object
    // -------------------------------------------------------------------------
    seqan::clear(gAlign);
    seqan::assignStringSet(gAlign, sequenceSet);
    TStringSet & seqSet = seqan::stringSet(gAlign); // make it into a dependent string set

    // Segment match generation
    // -------------------------------------------------------------------------
    current_time = seqan::sysTime();
    // Containers for segment matches and corresponding scores
    TFragmentString matches;
    TScoreValues scores;
    append_all_to_all_matches(seqSet, matches, scores); // all-to-all no segment matches
    std::cout << std::setw(30) << std::left << "Segment-match generation:" << std::setw(10) << std::right << seqan::sysTime() - current_time << " s" << std::endl;

    // Build Alignment Graph from matches
    // -------------------------------------------------------------------------
    current_time = seqan::sysTime();
    // Use these segment matches for the initial alignment graph
    TGraph g(seqSet);
    seqan::buildAlignmentGraph(matches, scores, g, seqan::FractionalScore());
    seqan::clear(matches);
    seqan::clear(scores);
    std::cout << std::setw(30) << std::left << "Build alignment graph:" << std::setw(10) << std::right << seqan::sysTime() - current_time << " s" << std::endl;
    std::cout << std::setw(30) << std::left << "Number of vertices:" << std::setw(10) << std::right << seqan::numVertices(g) << std::endl;

    // Build guide Tree
    // -------------------------------------------------------------------------
    current_time = seqan::sysTime();
    // Guide tree
    assert(!distance_matrix.empty()); // Check if we have a valid distance matrix
    auto guide_tree = neighbour_joining(std::move(distance_matrix));
    std::cout << std::setw(30) << std::left << "Build Guide Tree:" << std::setw(10) << std::right << seqan::sysTime() - current_time << " s" << std::endl;

    // Triplet extension
    // -------------------------------------------------------------------------
    current_time = seqan::sysTime();
    // Some alignment constants
    TSize nSeq = seqan::length(seqSet);
    bool isDeepAlignment = (nSeq > 50);  // threshold for what is a deep alignment
    TSize threshold = (isDeepAlignment) ? 30 : 10;  // experimentally proved relation

    if (nSeq < threshold)
        seqan::tripletLibraryExtension(g);
    else
        seqan::tripletLibraryExtension(g, guide_tree, threshold / 2);
    std::cout << std::setw(30) << std::left << "Triplet extension:" << std::setw(10) << std::right << seqan::sysTime() - current_time << "  " << std::endl;

    // Progressive Alignment
    // -------------------------------------------------------------------------
    current_time = seqan::sysTime();
    seqan::progressiveAlignment(g, guide_tree, gAlign);
    seqan::clear(guide_tree);
    seqan::clear(g);
    std::cout << std::setw(30) << std::left << "Progressive Alignment:" << std::setw(10) << std::right << seqan::sysTime() - current_time << " s" << std::endl;
}
