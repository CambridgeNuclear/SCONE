#!/bin/bash

HELP="Create a SCONE nuclear data library file from raw ACE Neutron CE files.

make_ace_lib.sh <output> <MODE> <search-word> <path-to-directory>

Arguments:
  <output>      Library file to create.
  <MODE>        Choose tosearch for CE Neutron cards (ZZAAA.TTc id). or SAB cards (XXXXXX.TTt)
  <search-word> Keyword expected within the name of the ace files (e.g., ace, nc, t, ENDF, sssth)
  <path-to-directory> Path of the base directory where to search in

The script recursively searches inside <path-to-directory> for files with a name that includes the given
common <search-word>; in each file, it searches for the presence of an ID pattern (ZZAAA.TTc for CE; 
XXXXXX.TTt for SAB card) within the first line. By ACE definition this should match only in the first 
line of a header of an ACE card. For each match the script prints a line of the SCONE nuclear data library
file to <output> as:
  ZAID; LINE_NUMBER; FILE;
  ...
"

# Display help if number of arguments is wrong
if [ $# -ne 4 ]; then
  echo "${HELP}"
  exit 1
fi

# Pop first argument to be the library file
OUTNAME=$1
MODE=$2
WORD=$3
SEARCH_DIR=$4

# Check a path to directory is given correctly
if [ ! -d "$SEARCH_DIR" ]; then
  echo "Error: '$SEARCH_DIR' is not a directory"
  exit 1
fi

# If file already exists ask for confirmation and remove it
if [ -f $OUTNAME ]; then
  echo "Warning: File ${OUTNAME} exists and will be overriden."
  read -p "Continue [Y/n]? "  -r
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    exit 1
  fi
  rm $OUTNAME
fi

FILES=()
while IFS= read -r -d '' f; do
  FILES+=("$f")
done < <(
  find "$SEARCH_DIR" -type f -name "*${WORD}*" -print0
)

if [ "${#FILES[@]}" -eq 0 ]; then
  echo "No files with format '*${WORD}*' found in ${SEARCH_DIR}"
  exit 1
fi

# Process each ACE file with awk
for var in "${FILES[@]}"; 
do

  # We need full path to put in the file
  FULL_PATH=$(realpath "$var")
  echo "Processing file $var"

  if [[ $MODE =~ ^CE$ ]]; then
    awk -v FILE=$FULL_PATH \
    '/^[[:space:]]*[[:digit:]][[:digit:]][[:digit:]][[:digit:]][[:digit:]]?.[[:digit:]][[:digit:]]c/\
      {print $1 "; " NR "; " FILE "; "}' $var >> $OUTNAME

  elif [[ $MODE =~ ^SAB$ ]]; then
    awk -v FILE=$FULL_PATH \
    '/^[[:space:]]*([[:alnum:]]|_|\-)+.[[:digit:]][[:digit:]]t/\
      {print $1 "; " NR "; " FILE "; "}' $var >> $OUTNAME

  else
    echo "Unknown mode '${MODE}'. Must be CE or SAB"
    exit 1
  fi

done