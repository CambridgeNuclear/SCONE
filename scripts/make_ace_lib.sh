#!/bin/bash

HELP="Create a SCONE nuclear data library file from raw ACE Neutron CE files.

make_ace_lib.sh <output> <MODE> <ace-file>+

Arguments:
  <output>  Library file to create.
  <MODE>    Choose tosearch for CE Neutron cards (ZZAAA.TTc id). or SAB cards (XXXXXX.TTt)
  <ace-files>  List of files to search for the ID pattern

The script searches each of the <ace-file> for a presence of an ID pattern (ZZAAA.TTc for CE;
XXXXXX.TTt for SAB card) at a beginning of the line. By ACE definition this should match only
a first line of a header of an ACE card. For each match the script prints a line of the SCONE
nuclear data library file to <output> as:
  ZAID; LINE_NUMBER; FILE;
  ...
"

# Display help if number of arguments is wrong
if [ $# -le 2 ]; then
  echo "${HELP}"
  exit 1
fi

# Pop first argument to be the library file
OUTNAME=$1
MODE=$2
shift 2

# If file already exists ask for confirmation and remove it
if [ -f $OUTNAME ]; then
  echo "Warning: File ${OUTNAME} exists and will be overriden."
  read -p "Continue [Y/n]? "  -r
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    exit 1
  fi
  rm $OUTNAME
fi

# Process each ACE file with awk
for var in $@
do
  if [ ! -f $var ]; then
    echo -e "File '${var}' does not exist \U0001F61E"
    echo -e "Quitting \U0001F44B"
    exit 1
  fi

  # We need full path to put in the file
  FULL_PATH=$(realpath $var)

  echo "Processing file ${var}"
  if [[ $MODE =~ ^CE$ ]]; then
    awk -v FILE=$FULL_PATH \
    '/^[[:space:]]*[[:digit:]][[:digit:]][[:digit:]][[:digit:]][[:digit:]]?.[[:digit:]][[:digit:]]c/\
      {print $1 "; " NR "; " FILE "; "}' $var >> $OUTNAME

  elif [[ $MODE =~ ^SAB$ ]]; then
    awk -v FILE=$FULL_PATH \
    '/^[[:space:]]*[[:alnum:]]+.[[:digit:]][[:digit:]]t/\
      {print $1 "; " NR "; " FILE "; "}' $var >> $OUTNAME

  else
    echo "Unknown mode '${MODE}'. Must be CE or SAB"
    exit 1
  fi

done
