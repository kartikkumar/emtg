#tool for reading NSGA-II population and archive files
#for use with EMTG-NSGAII outer-loop by Vavrina and Englander
#Python interface by Jacob Englander begun 3-16-2014

import os

class NSGAII_outerloop_solution(object):
    #fields
    Xouter = [] #outer-loop decision vector
    Xinner = [] #inner-loop decision vector
    objective_values = [] #vector of objective function values
    power_system_size = [] #if applicable, power level of the array in kW
    thruster = []
    number_of_thrusters = []
    launch_vehicle = []
    launch_date = []
    description = [] # case name, can be parsed for data
    mission_sequence = []
    generation_found = [] #what generation was this solution found?
    timestamp = [] #at what time, in seconds from program start, was this solution found?
    
    #constructor
    def __init__(self, input_line, column_headers):
        self.parse_input_line(input_line, column_headers)

    #line parser
    def parse_input_line(self, input_line, column_headers):
        #strip off the newline character
        input_line = input_line.strip("\n")

        #strip the input line by commas
        input_cell = input_line.split(',')

        #declare arrays of launch vehicle and thruster names
        LV_names = ['AV401','AV411','AV421','AV431','AV501','AV511','AV521','AV531','AV541','AV551',
                    'F910','F911','AV551s48','F9H','D4H','SLSb1']
        thruster_names = ['NSTAR','XIPS25','BPT4000HIsp','BPT4000Hthrust','BPT4000XHIsp',
                          'NEXTHIspv9','VASIMRargon','VSIxenonhall','NEXTHIspv10','NEXTHthrustv10',
                          'BPT4000MALTO','NEXIS','H6MS','BHT20K','HiVHAc']

        for column_index in range(0, len(column_headers)):
            if column_headers[column_index] == 'Generation found':
                self.generation_found = int(input_cell[column_index])

            elif column_headers[column_index] == 'timestamp':
                self.timestamp = int(input_cell[column_index])
            
            elif column_headers[column_index] == 'Description':
                self.description = input_cell[column_index]

                #find the mission sequence descriptor
                left_parenthesis_index = self.description.find('(')
                self.mission_sequence = self.description[left_parenthesis_index+1:].strip(')')
                
                #reconstruct the full mission description from the case name
                descriptioncell = self.description.split('_')

                for descriptionitem in descriptioncell:
                    if descriptionitem.find('kW') > 0: #this entry encodes power system size
                        self.power_system_size = float(descriptionitem.strip('kW'))

                    if descriptionitem.find('nTh') > 0:
                        self.number_of_thrusters = descriptionitem.strip('nTh')

                    for LV_name in LV_names:
                        if descriptionitem == LV_name:
                            self.launch_vehicle = descriptionitem

                    for thruster_name in thruster_names:
                        if descriptionitem == thruster_name:
                            self.thruster = thruster_name

            elif column_headers[column_index] == 'BOL power at 1 AU (kW)' \
                or column_headers[column_index] == 'Launch epoch (MJD)' \
                or column_headers[column_index] == 'Flight time (days)' \
                or column_headers[column_index] == 'Thruster preference' \
		        or column_headers[column_index] == 'Number of thrusters' \
		        or column_headers[column_index] == 'Launch vehicle preference' \
		        or column_headers[column_index] == 'Delivered mass to final target' \
		        or column_headers[column_index] == 'Final journey mass increment (for maximizing sample return)' \
                or column_headers[column_index] == 'First journey departure C3 (km^2/s^2)' \
                or column_headers[column_index] == 'Final journey arrival C3 (km^2/s^2)' \
                or column_headers[column_index] == 'Total delta-v (km/s)':
                #this entry is an objective function value
                self.objective_values.append(float(input_cell[column_index]))

            elif column_headers[column_index].find('Gene ') > 0:
                self.Xouter.append(int(input_cell[column_index]))

            else: #this entry is a member of the inner-loop decision vector
                self.Xinner.append(float(input_cell[column_index]))

#top-level container of NSGAII_outerloop_solution objects
class NSGAII_outerloop_population(object):
    #fields
    solutions = [] #vector of NSGAII_outerloop_solution objects
    global_column_headers = []
    gene_column_headers = []
    objective_column_headers = []
    number_of_feasible_solutions = []
    
    #constructor
    def __init__(self, population_file_name):
        self.parse_population_file(population_file_name)

    def clear(self):
        self.solutions = [] #vector of NSGAII_outerloop_solution objects
        self.global_column_headers = []
        self.gene_column_headers = []
        self.objective_column_headers = []
        self.number_of_feasible_solutions = []

    #method to read a population file
    def parse_population_file(self, population_file_name):
        #Step 1: attempt to open a population file
        if os.path.isfile(population_file_name):
            inputfile = open(population_file_name, "r")
            self.success = 1
        else:
            print "Unable to open", population_file_name, "EMTG Error"
            self.success = 0
            return

        #Step 2: scan through the file
        linenumber = 0
        for line in inputfile:
            #strip off the newline character
            line = line.replace("\n","")
            linenumber = linenumber + 1

            #the fourth line of the population file contains the column headers
            if linenumber == 4:
               self.global_column_headers = line.split(',')
               for header in self.global_column_headers:
                if header== 'BOL power at 1 AU (kW)' \
                    or header== 'Launch epoch (MJD)' \
                    or header== 'Flight time (days)' \
                    or header== 'Thruster preference' \
                    or header== 'Number of thrusters' \
                    or header== 'Launch vehicle preference' \
                    or header== 'Delivered mass to final target' \
                    or header== 'Final journey mass increment (for maximizing sample return)' \
                    or header== 'First journey departure C3 (km^2/s^2)' \
                    or header== 'Final journey arrival C3 (km^2/s^2)' \
                    or header== 'Total delta-v (km/s)':
                    self.objective_column_headers.append(header)

            #the fifth line of the population file contains the gene names
            elif linenumber == 5:
                self.gene_column_headers = line.split(',')

            #every line after the fifth is a solution line
            elif linenumber > 5:
                self.solutions.append(NSGAII_outerloop_solution(line, self.global_column_headers))