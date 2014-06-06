#!/bin/bash

# 2014_05_28
# ryne beeson
#
# this bash script can take 2 arguments
# $1: is the directory path where .qsub and .emtgopt files live
# $2: the second argument, if set to 'test' will not run qsub()
#     but instead just print out what would have run

# change to directory with .qsub and .emtgopt files
cwd=$(pwd)
cd $1

# count the number of .qsub files
numqsub=$(ls -1 | grep '.qsub' | wc -l)
# count the number of .emtgopt files
numopt=$(ls -1 | grep '.emtgopt' | wc -l)

# print the directory where submission will happen
echo ""
echo "Jobs Directory is: $1"
echo "----------------------------------------"
# print out the number of .qsub and .emtgopt files found
echo "The number of .qsub    files is: $numqsub"
echo "The number of .emtgopt files is: $numopt" 
echo ""

# if number of files are incorrect, give a message and exit
# otherwise, set numfiles
if [ "$numqsub" != "$numopt" ]; then
    echo "number of .qsub and .emtgopt files are not the same in $1, exiting bash script"
    exit       
else
    numjobs=$numqsub
fi

# save list of files to run
list=$(ls -1 | grep '.qsub')
# log name
log=job_log.txt

# if-statement to check if actually queueing jobs or just testing
if [[ $2 == 'test' ]]; then
    echo "Testing......"
else
    echo "Submitting Jobs to Queue......"
fi

# for loop to submit jobs to queue and
# generate a list of jobs submitted
i=1
for m in $list; do 
    # print job to terminal for user to see progress
    echo $i $m
    # write job line to text file
    echo $i $m >> $log
    # update counter, i
    i=$(($i + 1))	
    # strip windows return characters
    $(tr -d '\015'<$m >'Temp.qsub')
    $(rm $m)
    $(mv 'Temp.qsub' $m)
    if [[ $2 != 'test' ]]; then
        # submit job using qsub()
        qsub $m
    fi
done    

# copy job submission file to original directory
mv $log $cwd
cd $cwd

# tell user where job submission file is
# and what it is called
echo ""
echo "A log file, $log, has been created here: $cwd" 

exit

