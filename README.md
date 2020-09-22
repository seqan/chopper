# SeqAn3 App Template

[![Build Status](https://travis-ci.com/seqan/app-template.svg?branch=master)](https://travis-ci.com/seqan/app-template)

This is a template for app developers with SeqAn3. 
You can easily clone this repository and modify the existing code to your needs. 
It provides the elementary set-up for all SeqAn3 applications.

The example application is a FastQ to FastA file format converter.
It demonstrates exemplarily the set-up of test cases, documentation, and build infrastructure.
Probably you want to name your app differently — simply replace `fastq_to_fasta` with your app name in the following.
Please note that the command line interface tests fail if you use an individual project name without adapting the
name in the test file.

Instructions:
1. clone this repository: `git clone --recurse-submodules https://github.com/seqan/app-template.git fastq_to_fasta`
2. edit the project name in the *project* command of `fastq_to_fasta/CMakeLists.txt`
3. create a build directory and visit it: `mkdir build && cd build`
4. run cmake: `cmake ../fastq_to_fasta`
5. build the application: `make`
6. optional: build and run the tests: `make test`
7. optional: build the api documentation: `make doc`
8. execute the app: `./fastq_to_fasta`
