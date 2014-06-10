#!/bin/bash

#  2014_06_10
#  ryne beeson

#  purpose: move files around after leapfrog.py script is run
#  arguments:
#  $1 path to file(s) that need cleanup (includes .emtgopt and .qsub files)
#  $2 name of new directory to store files
#  $3 probe number (e.g. 1 or 2 or 3 etc...)

#  make a new directory to move files to
mkdir $2

#  mv file(s) into new directory
mv $1/*_$3_* $2
mv $1/EMTG_v8_results $2

