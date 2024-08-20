// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#pragma once

#include <vector>

#include <cereal/types/vector.hpp>

#include <chopper/configuration.hpp>

#include <hibf/sketch/hyperloglog.hpp>
#include <hibf/sketch/minhashes.hpp>

namespace chopper::sketch
{

struct sketch_file
{
    chopper::configuration chopper_config{};
    std::vector<std::vector<std::string>> filenames{};
    std::vector<seqan::hibf::sketch::hyperloglog> hll_sketches{};
    std::vector<seqan::hibf::sketch::minhashes> minHash_sketches{};

private:
    friend class cereal::access;

    template <typename archive_t>
    void serialize(archive_t & archive)
    {
        uint32_t version{1};
        archive(CEREAL_NVP(version));

        archive(CEREAL_NVP(chopper_config));
        archive(CEREAL_NVP(filenames));
        archive(CEREAL_NVP(hll_sketches));
        archive(CEREAL_NVP(minHash_sketches));
    }
};

} // namespace chopper::sketch
