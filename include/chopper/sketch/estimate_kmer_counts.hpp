#pragma once

#include <vector>

#include <chopper/data_store.hpp>

namespace chopper::sketch
{

inline void estimate_kmer_counts(data_store & store)
{
    assert(!store.filenames.empty());
    assert(store.filenames.size() == store.all_sketches.size());

    store.kmer_counts.resize(store.filenames.size());

    for (size_t i = 0; i < store.all_sketches.size(); ++i)
        store.kmer_counts[i] = store.all_sketches[i].estimate();
}

} // namespace chopper::sketch
