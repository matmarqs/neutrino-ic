#!/bin/sh

PROG="main"
PROGDIR="$(dirname "$(realpath $0)")"

make -s clean
cd "$PROGDIR"
make -s
"$PROGDIR"/"$PROG" # | "$PROGDIR"/genfig.py
make -s clean
