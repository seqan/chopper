#pragma once

#include <filesystem>
#include <string_view>

#include <chopper/prefixes.hpp>

namespace chopper::detail
{

/*!\brief Creates a filename and directory name from `prefix`.
 *
 * ### Example:
 * ```
 * foo     -> foo.count, foo_sketches
 * foo.txt -> foo.count, foo_sketches
 * path/   -> path/chopper.count, path/chopper_sketches
 * ```
 */
inline void apply_prefix(std::string_view const & prefix,
                         std::filesystem::path & filename,
                         std::filesystem::path & directory)
{
    // remove trailing slash if given
    std::filesystem::path safe_prefix{prefix};
    safe_prefix = (safe_prefix.has_filename()) ? safe_prefix.replace_extension("") : safe_prefix += prefix::chopper;

    filename += safe_prefix;
    filename += ".count";
    directory += safe_prefix;
    directory += "_sketches";
}

} // namespace chopper::detail
