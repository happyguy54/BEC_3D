#!/bin/bash

# Check if the first argument is "-c2"
if [ "$1" == "-c2" ]; then
  shift  # Remove the first argument
  indices=$1  # The second argument is now the first
  shift  # Remove the indices argument
  rej_osl_args=$1  # The third argument is now the second
else
  echo "Invalid option: $1" 1>&2
  exit 1
fi

echo "Indices: $indices"  # Debug line
echo "rej_OSL arguments: $rej_osl_args"  # Debug line

# Split the indices string into an array
IFS=',' read -ra IDX <<< "$indices"

echo "IDX: ${IDX[@]}"  # Debug line

# Loop over the indices and call the commands
for i in "${IDX[@]}"; do
  echo "Running commands with index $i"
  #sourceroot
  echo "./do_fit.py -u --unf -n --c2-index $i $rej_osl_args"
  ./do_fit.py -u --unf -n --c2-index $i $rej_osl_args
  echo "./do_fit.py -u -k --c2-index $i $rej_osl_args"
  ./do_fit.py -u -k --c2-index $i $rej_osl_args
  echo "./do_fit.py -u --c2-index $i $rej_osl_args"
  ./do_fit.py -u --c2-index $i $rej_osl_args
done