# Chopper - partition your sequences [![build status][1]][2] [![codecov][3]][4]

[1]: https://github.com/seqan/chopper/actions/workflows/ci_linux.yml/badge.svg?branch=master
[2]: https://github.com/seqan/chopper/actions?query=branch%3Amaster
[3]: https://codecov.io/gh/seqan/chopper/branch/master/graph/badge.svg?token=SJVMYRUKW2
[4]: https://codecov.io/gh/seqan/chopper

## System requirements

* GCC Version >= 9
* CMake Version >= 3.11

## General setup

Set up the repository:

```
git clone --recurse-submodules https://github.com/seqan/chopper
```

Set up the build directory
```
mkdir chopper_build
cd chopper_build
cmake ../chopper
```

Build the test to check if everything works
```
make test
```

:warning: The following section is not yet adapted to the multi-level approach :warning:

## Chopper pack

This submodule uses a hierarchical DP algorithm to pack user bins into a given number technical bins,
optimizing the space consumption of a Hierarchical Binning Directory.

The app `chopper pack` needs an input file with filenames and weights as input
(see chopper count if you want to use kmer counts as weights).

The file has to be **tab separated** and looks like this:

```
file_path1    500
file_path2    1000
file_path3    10
...
```

It can also contain meta information like taxonomic ids which you can group by later on:

```
file_path1    500    taxID_1
file_path2    1000   taxID_2
file_path2    10     taxID_2
...
```

If you have the file you can use `chopper pack` like this:

```
./chopper pack -f fata.tsv --technical-bins 2 -o output_filename.txt
```

Given the example tsv file with the 3 lines above, it will create a file `output_filename.txt` which looks like this:

```
#MERGED_BIN_0 max_bin_id:0
#HIGH_LEVEL_IBF max_bin_id:SPLIT_BIN_1
#BIN_ID SEQ_IDS NUM_TECHNICAL_BINS  ESTIMATED_MAX_TB_SIZE
MERGED_BIN_0_0  file_path1  62  9
MERGED_BIN_0_62 file_path2  2   5
SPLIT_BIN_1 file_path2  1   1000
```

Every line that starts with `#` is a header line and is needed for `chopper split` and `chopper build`.

The output file has the following columns:

1. `BIN_ID`: An identification of the bin. **MERGED/SPLIT** indicates a bin that has been merged/split.
             Note that MERGED bins with the same number belong together.
             In the example `MERGED_BIN_0_0` and `MERGED_BIN_0_62` both belong to `MERGED_BIN_0`.
2. `SEQ_IDS`: The file paths that belong to this bin.
3. `NUM_TECHNICAL_BINS`: The number of technical bins.
                         For SPLIT bins, this indicates the number of technical bins in the High-Level IBF.
                         For MERGED bins, this indicates the number of technical bins in the Low-Level IBF
                         (MERGED bins have always exactly one bin in the High-Level IBF).
3. `ESTIMATED_MAX_TB_SIZE`: The estimated maximum bin size of this bin
                            (can be neglected and is used for internal calculations).

## To be continued...

More Examples and descriptions will follow soon.
