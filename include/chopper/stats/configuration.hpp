#pragma once

#include <filesystem>

#include <cereal/cereal.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

namespace chopper::stats
{

struct configuration
{
    std::filesystem::path layout_file;
};

} // namespace chopper::stats
