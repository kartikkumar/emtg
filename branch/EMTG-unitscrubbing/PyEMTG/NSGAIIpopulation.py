#tool for reading NSGA-II population and archive files
#for use with EMTG-NSGAII outer-loop by Vavrina and Englander
#Python interface by Jacob Englander begun 3-16-2014

class NSGAII_outerloop_solution(object):
    #fields
    Xouter #outer-loop decision vector
    Xinner #inner-loop decision vector
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
            if column_headers[column_index].find('Gene'):
                self.Xouter.append(int(input_cell[column_index]))

            elif column_headers[column_index] == 'Generation found':
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
                    if descriptionitem.find('kW'): #this entry encodes power system size
                        self.power_system_size = double(descriptionitem.strip('kW'))

                    if descriptionitem.find('nTh'):
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
		        or column_headers[column_index] == 'Final journey mass increment (for maximizing sample return)':
                #this entry is an objective function value
                self.objective_values.append(float(input_cell[column_index]))

            else: #this entry is a member of the inner-loop decision vector
                self.Xinner.append(float(input_cell[column_index]))
