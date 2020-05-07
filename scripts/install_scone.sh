#!/bin/bash
# This script must be executed from SCONE source root directory!
# It is indented to be used by CI Service

# This command prints command executed to terminal
# Also makes script terminate if any command ends with non-0 exit status
set -ex

# Create build folder (Remove if it exists)
# NOTE: MUst have space between " and ]]
# ./Build" ]] -> WORKS
# ./Build"]]  -> 'conditional binary operator expected'
[[ -d "./Build" ]] && rm -rf ./Build
mkdir ./Build

# Configuire with CMake
cmake -E chdir ./Build cmake ..

# Build SCONE
make -C ./Build
