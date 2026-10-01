// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides chopper::layout::phibf::partition_user_bins.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <cstddef>
#include <vector>

#include <chopper/configuration.hpp>

#include <hibf/sketch/hyperloglog.hpp>
#include <hibf/sketch/minhashes.hpp>

namespace chopper::layout::phibf
{

//!\brief The partitioning approaches of the partitioned HIBF. Selected by `configuration::partitioning_approach`.
enum partitioning_scheme
{
    blocked,       // 0
    sorted,        // 1
    folded,        // 2
    weighted_fold, // 3
    similarity,    // 4
    lsh,           // 5
    lsh_sim        // 6
};

/*!\brief Distributes all user bins onto `config.number_of_partitions` partitions of a partitioned HIBF.
 * \param[in]  config           The configuration. Uses `number_of_partitions`, `partitioning_approach`, and
 *                              `hibf_config.{sketch_bits, number_of_user_bins, tmax}`.
 * \param[in]  cardinalities    The cardinality of each user bin.
 * \param[in]  sketches         The HyperLogLog sketch of each user bin.
 * \param[in]  minHash_sketches The MinHash tables of each user bin. Only used by the `lsh` and `lsh_sim` approaches.
 * \param[out] partitions       Must have size `config.number_of_partitions` on entry. On return, `partitions[i]` holds
 *                              the user bin indices assigned to partition `i`. No partition is empty.
 * \throws std::invalid_argument If `config.partitioning_approach` is not a phibf::partitioning_scheme, or if there are
 *                               fewer user bins than partitions.
 * \throws std::logic_error If not all user bins have been assigned to a partition.
 * \throws std::runtime_error If a partition is empty.
 *
 * Each partition is laid out as a separate HIBF. In contrast to the fast layout's
 * chopper::layout::partition_user_bins, user bins are never split across partitions.
 */
void partition_user_bins(chopper::configuration const & config,
                         std::vector<size_t> const & cardinalities,
                         std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                         std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches,
                         std::vector<std::vector<size_t>> & partitions);

} // namespace chopper::layout::phibf
