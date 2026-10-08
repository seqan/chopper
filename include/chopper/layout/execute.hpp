// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#pragma once

#include <string>
#include <vector>

#include <chopper/configuration.hpp>

#include <hibf/sketch/hyperloglog.hpp>
#include <hibf/sketch/minhashes.hpp>

namespace chopper::layout
{

/*!\brief Computes the layout for the given user bins and writes it to `config.output_filename`.
 * \param[in,out] config           The configuration. `config.hibf_config` is validated and completed
 *                                 (`validate_and_set_defaults`), and the timers are updated.
 * \param[in]     filenames        The file names of each user bin. They are written to the layout file.
 * \param[in]     sketches         The HyperLogLog sketch of each user bin.
 * \param[in]     minHash_sketches The MinHash sketches of each user bin. Only used, and then required for every user
 *                                 bin, if `config.fast_layout` is set. May be empty otherwise.
 * \returns 0.
 * \throws std::invalid_argument If both `config.determine_best_tmax` and `config.fast_layout` are set, or if the
 *         computed layout is invalid (see seqan::hibf::layout::layout::validate).
 *
 * The layout is computed with
 * - determine_best_number_of_technical_bins if `config.determine_best_tmax` is set,
 * - fast_layout if `config.fast_layout` is set,
 * - the DP layout of the HIBF library (seqan::hibf::layout::compute_layout) otherwise.
 *
 * The layout is validated (seqan::hibf::layout::layout::validate) before it is written. Warnings and notes of the
 * validation are printed to `std::cerr`, notes only in debug builds of the HIBF library.
 *
 * Unless `config.determine_best_tmax` is set, `config.output_verbose_statistics` prints statistics of the layout to
 * `std::cout`.
 */
int execute(chopper::configuration & config,
            std::vector<std::vector<std::string>> const & filenames,
            std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
            std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches);

} // namespace chopper::layout
