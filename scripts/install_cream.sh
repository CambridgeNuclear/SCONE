#!/bin/bash
# This script must be executed from SCONE source root directory!
# It is indented to be used by CI Service

# This command prints command executed to terminal
# Also makes script terminate if any command ends with non-0 exit status
set -ex

# Change folder to cream
cd ./cream

# Install Cream
pip install -e .[test]

# Return to root directory
cd ./..
