#!/usr/bin/python

#  ryne beeson - rbeeson2@illinois.edu
#  2014_06_10

#  import statements
import sys
import csv

def main(argv=None):
	#  parse variable arguments in
	if argv is None:
		argv = sys.argv
		
	#  the first argument to main() should be the directory path that includes the 
	#+ .emtgopt and .qsub files
	main_dir = argv[1]
	#  the second argument to main() should be the stagger time for the new 
	#+ starting epoch
	epoch_stagger = float(argv[2])
	#  the third argument needs to be removed in a future version, but for now
	#+ it is used for simplicity. 
	#  the third argument is the number of nodes*cores to be used on BW
	num_cores = int(argv[3])
	#  a subdirectory called EMTG_v8_results should exist
	#  if .emtg_result_cleanup.sh was run on this directory
	#+ then it will have the LRTS_summary.csv file
	#  create a variable holding its path
	sub_dir = argv[1] + 'EMTG_v8_results'	 
	#  name of the result file
	result_filename = 'LRTS_summary.csv'
	#  path to the result file
	result_filepath = sub_dir + '/' + result_filename

	#  read in the result file
	data = []
	with open(result_filepath, 'r') as csvfile:
		input = csv.reader(csvfile, delimiter = ',')
		for row in input:
			data.append(row)
	#  cut off the file labeling in the first 2 rows
	data = data[2:]
	#  number of columns before sequence starts (asteroid, epoch, etc..)
	data_fixcols = 7
	#  convert data to appropriate numerical types (integers, floats)
	for i in range(len(data)):
		for j in range(1,len(data[i])):
			if j == 6: data[i][j] = int(data[i][j])
			elif j > 6 and j%2 == 1: data[i][j] = int(data[i][j])
			else: data[i][j] = float(data[i][j])
	#  close the result file
	csvfile.close()

	#  the number of missions run 
	#+ (note: not all LRTS missions complete, if len(asteroid_sequence) == 1 then 
	#+  it did not complete)
	number_of_missions = len(data)
	
	#  missions that were submitted
	
	#  for each mission run, read in its results from the result file
	#+ and then find its .emtgopt, .qsub, and .asteroidlist file
	#+ read those files in, copy to a new file, and make appropriate changes
	for mission in data:
		#  grab the number of asteroids visited in the sequence
		num_asteroids = mission[6]
		if num_asteroids == 1: continue		
		#  grab the mission name
		mission_name = mission[0]
		#  create a new mission name
		if '_1_LRTS' in mission_name: 
			new_mission_name = mission_name.replace('_1_LRTS', '_2_LRTS')
		elif '_2_LRTS' in mission_name:
			new_mission_name = mission_name.replace('_2_LRTS', '_3_LRTS')
		#  grab the starting asteroid location
		starting_asteroid = mission[7]
		#  segregate the asteroids and epochs into lists
		asteroid_sequence = []
		epoch_sequence = []
		for i in range(7, len(mission)):
			if i%2 == 1: asteroid_sequence.append(mission[i])
			if i%2 == 0: epoch_sequence.append(mission[i])
		
		''' .emtgopt file generation '''
		#  the options file path
		opt_file_path = main_dir + mission_name + '.emtgopt'
		new_opt_file_path = main_dir + new_mission_name + '.emtgopt'
		#  open the options file for the mission
		#  open a new options file for the next mission
		opt_file = open(opt_file_path, 'r')
		new_opt_file = 	open(new_opt_file_path, 'w')
		
		#  read each line of the option file, make appropriate changes
		#+ and write the new options file
		for line in opt_file:
			if 'lazy_race_tree_target_list_file' in line:
				line = 'lazy_race_tree_target_list_file ' + new_mission_name + '.asteroidlist\n'
			elif 'mission_name' in line:
				line = 'mission_name ' + new_mission_name + '\n'
			elif 'launch_window_open_date' in line: 
				current_epoch = float(line[24:])
				new_epoch = current_epoch + epoch_stagger
				line = 'launch_window_open_date ' + str(new_epoch) + '\n'
			new_opt_file.write(line)
		
		#  close the options files
		opt_file.close()
		new_opt_file.close()
	
		''' .qsub file generation '''
		#  the .qsub file path
		qsub_file_path = main_dir + mission_name + '.qsub'
		new_qsub_file_path = main_dir + new_mission_name + '.qsub'
		#  open the qsub file for the mission
		#  open a new qsub file for the next mission
		qsub_file = open(qsub_file_path, 'r')
		new_qsub_file = open(new_qsub_file_path, 'w')	
			
		#  read each line of the qsub file, make appropriate changes
		#+ and write the new qsub file
		for line in qsub_file:
			if '#PBS -N' in line:
				line = '#PBS -N ' + new_mission_name + '\n'
			elif 'aprun -n' in line:
				line = 'aprun -n ' + str(num_cores) + ' ../emtg ' + new_mission_name + '.emtgopt\n'
			new_qsub_file.write(line)
		
		#  close the qsub files
		qsub_file.close()
		new_qsub_file.close()
		
		''' .asteroidlist file generation '''
		#  the .asteroidlist file path
		asteroidlist_path = main_dir + mission_name + '.asteroidlist'
		new_asteroidlist_path = main_dir + new_mission_name + '.asteroidlist'
		#  open the asteroidlist file for the mission
		#  open a new asteroidlist file for the next mission
		asteroidlist = open(asteroidlist_path, 'r')
		new_asteroidlist = open(new_asteroidlist_path, 'w')	
		
		#  read each line of the .asteroid file, remove asteroids if they 
		#+ were already visited by the previous probe, and write a new 
		#+ .asteroid file
		for line in asteroidlist:
			if line.isspace(): continue 
			if int(line) in asteroid_sequence: continue
			else: new_asteroidlist.write(line)
			
		#  close the .asteroid files
		asteroidlist.close()
		new_asteroidlist.close()
		
	
		


		
if __name__ == "__main__":
	main()

