#pragma once

#include <cstdlib>
#include <string>

#include <chopper/build/build_data.hpp>
#include <chopper/detail_bin_prefixes.hpp>

void parse_chopper_pack_header_line(std::string const & line, build_data & data)
{
    if (line.substr(1, hibf_prefix.size()) == hibf_prefix)
    {
        assert(line.substr(hibf_prefix.size() + 2, 11) == "max_bin_id:");
        auto it = std::find(line.begin() + hibf_prefix.size() + 13, line.end(), '_'); // skip "MERGED"/"SPLIT"
        ++it; // skip "_"
        it = std::find(it, line.end(), '_'); // skip "BIN"
        ++it; // skip "_"
        data.hibf_max_bin = std::atoi(std::string(it, line.end()).c_str());
    }
    else if (line.substr(1, merged_bin_prefix.size()) == merged_bin_prefix)
    {
        std::string const hidx_str(line.begin() + 1 /*#*/ + merged_bin_prefix.size() + 1 /*_*/,
                                   std::find(line.begin() + merged_bin_prefix.size() + 2, line.end(), ' '));
        assert(line.substr(merged_bin_prefix.size() + hidx_str.size() + 3, 11) == "max_bin_id:");
        std::string const lidx_str = line.substr(merged_bin_prefix.size() + hidx_str.size() + 14,
                                                 line.size() - merged_bin_prefix.size() - hidx_str.size() - 14);

        size_t const hidx = std::atoi(hidx_str.c_str());
        size_t const lidx = std::atoi(lidx_str.c_str());

        data.merged_max_bin_map.emplace(hidx, lidx);
    }
}
