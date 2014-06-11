#!/bin/bash

#  just print a list of numbers
#+ to a new file
for i in {2..16257}; do
    echo $i >> Full_AsteroidList.asteroidlist
done

#  exit shell script
exit
