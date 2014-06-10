#!/bin/bash

#  2014_05_29
#  ryne beeson

#
#  INSTRUCTIONS:
#  this bash shell script only takes one argument.
#  the argument is the directory where all the 
#+ result subdirectories reside
#
#  e.g. $ ./emtg_result_cleanup.sh /My_EMTG_Result_Folder
#  

#  get current directory
cwd=$(pwd)

#  change to main results directory
cd $1

#  count the number of result directories
numdir=$(ls -1 | wc - l)

#  generate list of directories to step through
#+ and grab results from the tree_summary.LRTS
list=$(ls -1)

#  name of tree_summary.LRTS file
file='tree_summary.LRTS'
#  name of summary file
log='LRTS_summary.csv'
echo This is an LRTS Summary File >> $log
echo 'FileName, Initial Mass, Final Mass, Starting Epoch, Final Epoch, Epoch Diff, Total Score, Sequence (Asteroid - Epoch)' >> $log

#  using a for-loop to step into each directory,
#+ find the results information in the $file,
#+ record/store the info and step back out
#  counter
i=1
for m in $list; do
	#  skip over $log if it existed previously
	if [[ $m == $log ]]; then continue; fi
    #  step into a result directory
    cd $m
    #  get the best sequence
    tseq1=$(cat $file | grep -A 1 'Best')
    tseq2=$(cat $file | grep -A 2 'Best')
    tseq3=$(cat $file | grep -A 3 'Best')
    #  set the offset to start reading entries
    declare -i offset=-7
    #  count the number of entries minus 7 for
    #+ the string 'Best Sequence Found:' 
    #  also create a new list called 
    #+ seq1, which has unnecessary information stripped
    seq1='' 
    declare -i numasteroids=$offset
    for n in $tseq1; do
        numasteroids=$(($numasteroids + 1))
        if [[ $numasteroids -eq 1 ]]; then seq1="$n "; fi
        if [[ $numasteroids -gt 1 ]]; then seq1="$seq1 $n"; fi
    done
    #  strip information from tseq2 -> seq2
    seq2=''
    declare -i count=0
    for n in $tseq2; do
		epoch=${n%.*}
		epoch=${epoch%,*}
    	if [[ $count -eq $(($numasteroids + -$offset)) ]]; then
    		epoch=$(($epoch/86400))
    		seq2="$epoch" iepoch="$epoch"; 
    	fi
    	if [[ $count -gt $(($numasteroids + -$offset)) ]]; then
			epoch=$(($epoch/86400))
    		seq2="$seq2, $epoch"; 
    	fi
	    if [[ $count -eq $((2*$numasteroids + -$offset - 1)) ]]; then
	    	fepoch="$epoch"
	    	depoch=$(($fepoch - $iepoch))
	    fi
    	count=$(($count + 1))
    done
	#  strip information from tseq3 -> seq3
	seq3=''
    declare -i count=0
	for n in $tseq3; do
    	if [[ $count -eq $((2*$numasteroids + -$offset)) ]]; then seq3="$n " imass="$n"; fi
    	if [[ $count -gt $((2*$numasteroids + -$offset)) ]]; then seq3="$seq3 $n"; fi
		if [[ $count -eq $((3*$numasteroids + -$offset - 1)) ]]; then fmass="$n"; fi
    	count=$(($count + 1))
    done
    #  append the filename and total score
    if [[ $numasteroids -eq 1 ]]; then
        seq="$m, $imass, $fmass, $iepoch, $fepoch, $depoch, $numasteroids, "
    else
        seq="$m, $imass $fmass, $iepoch, $fepoch, $depoch, $numasteroids, "
    fi
    #  append the first asteroid then its epoch, the second asteroid then its epoch
    #+ etc.....
    declare -i ncount=1
    for n in $seq1; do
    	declare -i kcount=1
    	for k in $seq2; do
    		if [[ $ncount -eq $kcount ]] && [[ $ncount -ne $numasteroids ]]; then seq="$seq $n $k "; fi
    		if [[ $ncount -eq $kcount ]] && [[ $ncount -eq $numasteroids ]]; then seq="$seq $n, $k "; fi
    		kcount=$(($kcount + 1))
    	done
    	ncount=$(($ncount + 1))
    done
    #  step back a directory
	cd ..
	#  echo the best sequence
	echo $seq
	echo $seq >> $log
    #  increment the counter
    i=$(($i + 1))
done

#  exit shell script
exit

