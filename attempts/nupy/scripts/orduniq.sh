#!/bin/sh

# This script should be runned as './orduniq.sh output.txt > sort.txt'

EXECDIR="$(dirname $(realpath "$0"))"
FILEINPUT="$(realpath "$1")"

"$EXECDIR"/sort.py "$FILEINPUT" | uniq
