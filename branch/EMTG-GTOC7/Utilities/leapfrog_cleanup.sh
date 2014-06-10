#!/bin/bash

#  2014_06_10
#  ryne beeson

#  purpose: move files around after leapfrog.py script is run
#  arguments:
#  $1 path to file(s) that need cleanup
#  $2 name of new directory to store files
#  $3 probe number

#  make a new directory to move files to
mkdir $2

#  mv file(s) into new directory
mv $1/*_$3_* $2
mv $1/EMTG_v8_results $2

