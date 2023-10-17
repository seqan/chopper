// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include <hibf/sketch/hyperloglog.hpp>

struct config
{
    std::filesystem::path input{};
    std::filesystem::path output{};
    uint8_t threads{1u};
};

void execute_general(config const & cfg);
void execute_sizes(config const & cfg);

void process_file(std::string const & filename,
                  std::vector<uint64_t> & current_kmers,
                  seqan::hibf::sketch::hyperloglog & sketch,
                  bool const fill_current_kmers,
                  uint8_t const kmer_size);