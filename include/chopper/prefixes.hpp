#pragma once

#include <string_view>

namespace chopper::prefix
{

constexpr std::string_view chopper{"chopper"};

constexpr std::string_view header{"#"};

constexpr std::string_view header_config{"#"};

constexpr std::string_view high_level{"HIGH_LEVEL_IBF"};

constexpr std::string_view first_header_line{"#HIGH_LEVEL_IBF"};
static_assert(first_header_line.starts_with(header));
static_assert(first_header_line.ends_with(high_level));

constexpr std::string_view merged_bin{"MERGED_BIN"};

constexpr std::string_view split_bin{"SPLIT_BIN"};

} // namespace chopper::prefix
