#tool for reading NSGA-II population and archive files
#for use with EMTG-NSGAII outer-loop by Vavrina and Englander
#Python interface by Jacob Englander begun 3-16-2014

import Universe

import os
import numpy as np
from scipy.integrate import ode
import matplotlib
import matplotlib.dates as dates
import matplotlib.ticker as ticker
import pylab
from mpl_toolkits.mplot3d import Axes3D
import copy
import wx
import datetime

class NSGAII_outerloop_solution(object):
    
    #constructor
    def __init__(self, input_line, column_headers):
        self.initialize()
        self.parse_input_line(input_line, column_headers)

    #clear function
    def initialize(self):
        self.Xouter = [] #outer-loop decision vector
        self.Xinner = [] #inner-loop decision vector
        self.objective_values = [] #vector of objective function values
        self.power_system_size = [] #if applicable, power level of the array in kW
        self.thruster = []
        self.number_of_thrusters = []
        self.launch_vehicle = []
        self.launch_date = []
        self.description = '' # case name, can be parsed for data
        self.mission_sequence = []
        self.generation_found = [] #what generation was this solution found?
        self.timestamp = [] #at what time, in seconds from program start, was this solution found?
        self.Legal_Solution = False

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
                          'BPT4000MALTO','NEXIS','H6MS','BHT20K','HiVHAc','13kWSTMDHallHisp','13kWSTMDHallHthrust',
                          'NEXT_TT11_Hisp','NEXT_TT11_Hthrust','NEXT_TT11_expanded',
                          '13kWSTMDHall_10_1_2014_Hisp','13kWSTMDHall_10_1_2014_Mthrust','13kWSTMDHall_10_1_2014_Hthrust']

        for column_index in range(0, len(column_headers)):
            if column_index < len(input_cell):
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
                        if descriptionitem.rfind('kW') == len(descriptionitem) - 2: #this entry encodes power system size
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
                    or column_headers[column_index] == 'First journey departure C3 (km^2/s^2)' \
                    or column_headers[column_index] == 'Final journey arrival C3 (km^2/s^2)' \
                    or column_headers[column_index] == 'Total delta-v (km/s)':
                    #this entry is an objective function value
                    self.objective_values.append(float(input_cell[column_index]))
                elif column_headers[column_index] == 'Thruster preference' \
                    or column_headers[column_index] == 'Number of thrusters' \
		            or column_headers[column_index] == 'Launch vehicle preference':
                    self.objective_values.append(int(input_cell[column_index]))
                elif column_headers[column_index] == 'Delivered mass to final target (kg)' \
                    or column_headers[column_index] == 'Final journey mass increment (for maximizing sample return)':
                    self.objective_values.append(-float(input_cell[column_index]))
                elif column_headers[column_index] == 'Point-group value':
                    self.objective_values.append(int(input_cell[column_index]))
                elif column_headers[column_index].find('Gene ') > 0:
                    self.Xouter.append(int(input_cell[column_index]))

                elif input_cell[column_index] != '': #this entry is a member of the inner-loop decision vector
                    self.Xinner.append(float(input_cell[column_index]))

    def plot_solution(self, PopulationAxes, PopulationFigure, ordered_list_of_objectives, colorbar, lowerbounds, upperbounds):
        if len(ordered_list_of_objectives) == 2: #2D
            self.point = PopulationAxes.scatter(self.objective_values[ordered_list_of_objectives[0]], self.objective_values[ordered_list_of_objectives[1]], s=50, c='b', marker='o', lw=0, picker=1)
        elif len(ordered_list_of_objectives) == 3: #3D
                self.point = PopulationAxes.scatter(self.objective_values[ordered_list_of_objectives[0]], self.objective_values[ordered_list_of_objectives[1]], self.objective_values[ordered_list_of_objectives[2]], s=50, c='b', marker='o', lw=0, picker=1)
        else: #4D
            if self.colorbar is None:
                self.point = PopulationAxes.scatter(self.objective_values[ordered_list_of_objectives[0]], self.objective_values[ordered_list_of_objectives[1]], self.objective_values[ordered_list_of_objectives[2]], s=50, c=self.objective_values[ordered_list_of_objectives[4]], marker='o', lw=0, picker=1)
                self.point.set_clim([lowerbounds[-1],upperbounds[-1]])
                colorbar = PopulationFigure.colorbar(self.point)
            else:
                self.point = PopulationAxes.scatter(self.objective_values[ordered_list_of_objectives[0]], self.objective_values[ordered_list_of_objectives[1]], self.objective_values[ordered_list_of_objectives[2]], s=50, c=self.objective_values[ordered_list_of_objectives[4]], marker='o', lw=0, picker=1)
                self.point.set_clim([lowerbounds[-1],upperbounds[-1]])

        self.picker = self.point.figure.canvas.mpl_connect('pick_event', self.onpick)

    def onpick(self, event):
        #description = []
        #for objective_index in ordered_list_of_objectives:
        #    description = description + self.objective_column_headers[objective_index] + ': ' + str(
        ind = event.ind[0]
        x, y, z = event.artist._offsets3d
        print self.description
        #print self.description, x[ind], y[ind], z[ind]
        #print ind


#top-level container of NSGAII_outerloop_solution objects
class NSGAII_outerloop_population(object):
        
    #constructor
    def __init__(self, population_file_name):
        self.clear()
        self.parse_population_file(population_file_name)

    def clear(self):
        self.solutions = [] #vector of NSGAII_outerloop_solution objects
        self.legal_solutions = []
        self.global_column_headers = []
        self.gene_column_headers = []
        self.objective_column_headers = []
        self.number_of_feasible_solutions = []
        self.points_array = []
        self.ordered_list_of_objectives = []

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
                if header == 'BOL power at 1 AU (kW)' \
                    or header == 'Launch epoch (MJD)' \
                    or header == 'Flight time (days)' \
                    or header == 'Thruster preference' \
                    or header == 'Number of thrusters' \
                    or header == 'Launch vehicle preference' \
                    or header == 'Delivered mass to final target (kg)' \
                    or header == 'Final journey mass increment (for maximizing sample return)' \
                    or header == 'First journey departure C3 (km^2/s^2)' \
                    or header == 'Final journey arrival C3 (km^2/s^2)' \
                    or header == 'Total delta-v (km/s)' \
                    or header == 'Inner-loop objective function' \
                    or header == 'Point-group value':
                    self.objective_column_headers.append(header)

            #the fifth line of the population file contains the gene names
            elif linenumber == 5:
                self.gene_column_headers = line.split(',')

            #every line after the fifth is a solution line
            elif linenumber > 5:
                tempSolution = NSGAII_outerloop_solution(line, self.global_column_headers)
                self.solutions.append(tempSolution)

    #method to plot the population
    #input is an ordered list of objectives, [x, y, z, color]. If there are two objectives, a monochrome 2D plot will be shown. If there are three objectives, a monochrome 3D plot will be shown.
    #if there are four, a colored 3D plot will be shown. If there are more than four there will be an error message.
    def plot_population(self, ordered_list_of_objectives, LowerBounds = None, UpperBounds = None, TimeUnit = 1, EpochUnit = 1, FontSize = 10, BaseMarkerSize = 20):
        self.ordered_list_of_objectives = ordered_list_of_objectives
        self.LowerBounds = LowerBounds
        self.UpperBounds = UpperBounds
        self.TimeUnit = TimeUnit
        self.EpochUnit = EpochUnit
        self.BaseMarkerSize = BaseMarkerSize

        #first check to see if the correct number of objective function indices were supplied
        if len(self.ordered_list_of_objectives) < 2 or len(self.ordered_list_of_objectives) > 5:
            print "NSGAII_outerloop_population::plot_population ERROR. You must specify between two and five objective functions to plot."
            return

        self.PopulationFigure = matplotlib.pyplot.figure()
        matplotlib.rcParams.update({'font.size': FontSize})
        self.PopulationFigure.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99)
        if len(ordered_list_of_objectives) == 2:
            self.PopulationAxes = self.PopulationFigure.add_subplot(111)
        else:
            self.PopulationAxes = self.PopulationFigure.add_subplot(111, projection='3d')

        #build up a list of objective values to be plotted
        self.objective_values_matrix = []
        self.legal_solutions = []
        for objective_index in range(0, len(self.ordered_list_of_objectives)):
            objective_values_vector = []
            for solution in self.solutions:
                if max(solution.objective_values) < 1.0e+99:
                    solution.Legal_Solution = True
                    if not (self.LowerBounds == None or self.UpperBounds == None):
                        #if bounds were supplied, check to see if the solution fits inside the bounds
                        for obj in range(0, len(self.ordered_list_of_objectives)):
                            if solution.objective_values[ordered_list_of_objectives[obj]] < self.LowerBounds[obj] or solution.objective_values[ordered_list_of_objectives[obj]] > self.UpperBounds[obj]:
                                solution.Legal_Solution = False
                                break
                    if solution.Legal_Solution:
                        self.legal_solutions.append(solution)

                        if self.objective_column_headers[self.ordered_list_of_objectives[objective_index]] == 'Flight time (days)' and self.TimeUnit == 0:
                            objective_values_vector.append(solution.objective_values[ordered_list_of_objectives[objective_index]] / 365.25)
                        elif self.objective_column_headers[self.ordered_list_of_objectives[objective_index]] == 'Launch epoch (MJD)' and self.EpochUnit == 0:
                            objective_values_vector.append(dates.date2num(datetime.datetime.fromtimestamp(wx.DateTimeFromJDN(solution.objective_values[ordered_list_of_objectives[objective_index]] + 2400000.5).GetTicks())))
                        else:
                            objective_values_vector.append(copy.deepcopy(float(solution.objective_values[ordered_list_of_objectives[objective_index]])))
            self.objective_values_matrix.append(np.array(objective_values_vector))

        #determine upper and lower bounds on each objective
        self.upperbounds = []
        self.lowerbounds = []
        for objective_index in range(0, len(self.ordered_list_of_objectives)):
            self.upperbounds.append(self.objective_values_matrix[objective_index].max())
            self.lowerbounds.append(self.objective_values_matrix[objective_index].min())

        #plot each solution
        self.plot_solution_points()

        #set the axes labels
        if self.objective_column_headers[self.ordered_list_of_objectives[0]] == 'Flight time (days)' and self.TimeUnit == 0:
            self.PopulationAxes.set_xlabel('Flight time (years)')
        elif self.objective_column_headers[self.ordered_list_of_objectives[0]] == 'Launch epoch (MJD)' and self.EpochUnit == 0:
            self.PopulationAxes.set_xlabel('Launch Epoch (TDB Gregorian)')
            self.PopulationAxes.w_xaxis.set_major_formatter(ticker.FuncFormatter(self.format_date))
        else:
            self.PopulationAxes.set_xlabel(self.objective_column_headers[self.ordered_list_of_objectives[0]])

        if self.objective_column_headers[self.ordered_list_of_objectives[1]] == 'Flight time (days)' and self.TimeUnit == 0:
            self.PopulationAxes.set_ylabel('Flight time (years)')
        elif self.objective_column_headers[self.ordered_list_of_objectives[1]] == 'Launch epoch (MJD)' and self.EpochUnit == 0:
            self.PopulationAxes.set_ylabel('Launch Epoch (TDB Gregorian)')
            self.PopulationAxes.w_yaxis.set_major_formatter(ticker.FuncFormatter(self.format_date))
        else:
            self.PopulationAxes.set_ylabel(self.objective_column_headers[self.ordered_list_of_objectives[1]])

        if len(ordered_list_of_objectives) > 2:
            if self.objective_column_headers[self.ordered_list_of_objectives[2]] == 'Flight time (days)' and self.TimeUnit == 0:
                self.PopulationAxes.set_zlabel('Flight time (years)')
            elif self.objective_column_headers[self.ordered_list_of_objectives[2]] == 'Launch epoch (MJD)' and self.EpochUnit == 0:
                self.PopulationAxes.set_zlabel('Launch Epoch (TDB Gregorian)')
                self.PopulationAxes.w_zaxis.set_major_formatter(ticker.FuncFormatter(self.format_date))
            else:
                self.PopulationAxes.set_zlabel(self.objective_column_headers[self.ordered_list_of_objectives[2]])
            self.PopulationAxes.autoscale_view(tight=True, scalex=True, scaley=True, scalez=True)
        else:
            self.PopulationAxes.autoscale_view(tight=True, scalex=True, scaley=True)
        self.PopulationAxes.grid(b=True)

        self.PopulationFigure.show()

    def plot_solution_points(self):
        self.colorbar = None
        self.solution_names = []
        X = []
        Y = []
        Z = []
        C = []
        S = []
        for solution in self.solutions:
            if solution.Legal_Solution:
                if self.objective_column_headers[self.ordered_list_of_objectives[0]] == 'Flight time (days)' and self.TimeUnit == 0:
                    X.append(solution.objective_values[self.ordered_list_of_objectives[0]] / 365.25)
                elif self.objective_column_headers[self.ordered_list_of_objectives[0]] == 'Launch epoch (MJD)' and self.EpochUnit == 0:
                    X.append(dates.date2num(datetime.datetime.fromtimestamp(wx.DateTimeFromJDN(solution.objective_values[0] + 2400000.5).GetTicks())))
                else:
                    X.append(solution.objective_values[self.ordered_list_of_objectives[0]])

                if self.objective_column_headers[self.ordered_list_of_objectives[1]] == 'Flight time (days)' and self.TimeUnit == 0:
                    Y.append(solution.objective_values[self.ordered_list_of_objectives[1]] / 365.25)
                elif self.objective_column_headers[self.ordered_list_of_objectives[1]] == 'Launch epoch (MJD)' and self.EpochUnit == 0:
                    Y.append(dates.date2num(datetime.datetime.fromtimestamp(wx.DateTimeFromJDN(solution.objective_values[1] + 2400000.5).GetTicks())))
                else:
                    Y.append(solution.objective_values[self.ordered_list_of_objectives[1]])

                if len(self.ordered_list_of_objectives) > 2:
                    if self.objective_column_headers[self.ordered_list_of_objectives[2]] == 'Flight time (days)' and self.TimeUnit == 0:
                        Z.append(solution.objective_values[self.ordered_list_of_objectives[2]] / 365.25)
                    elif self.objective_column_headers[self.ordered_list_of_objectives[2]] == 'Launch epoch (MJD)' and self.EpochUnit == 0:
                        Z.append(dates.date2num(datetime.datetime.fromtimestamp(wx.DateTimeFromJDN(solution.objective_values[2] + 2400000.5).GetTicks())))
                    else:
                        Z.append(solution.objective_values[self.ordered_list_of_objectives[2]])

                if len(self.ordered_list_of_objectives) > 3:
                    if self.objective_column_headers[self.ordered_list_of_objectives[3]] == 'Flight time (days)' and self.TimeUnit == 0:
                        C.append(solution.objective_values[self.ordered_list_of_objectives[3]] / 365.25)
                    elif self.objective_column_headers[self.ordered_list_of_objectives[3]] == 'Launch epoch (MJD)' and self.EpochUnit == 0:
                        C.append(dates.date2num(datetime.datetime.fromtimestamp(wx.DateTimeFromJDN(solution.objective_values[3] + 2400000.5).GetTicks())))
                    else:
                        C.append(solution.objective_values[self.ordered_list_of_objectives[3]])


                #if there is a fifth objective, size the markers to reflect it
                if len(self.ordered_list_of_objectives) > 4:
                    S.append(self.BaseMarkerSize * solution.objective_values[self.ordered_list_of_objectives[4]] / (self.upperbounds[4] - self.lowerbounds[4]))
                else:
                    S.append(self.BaseMarkerSize)

        solution.points = self.PopulationAxes.scatter(X, Y, Z, s=S, c=C, edgecolors='none', marker='o', picker=1)
        
        self.picker = self.PopulationFigure.canvas.mpl_connect('pick_event', self.onpick)
        if len(self.ordered_list_of_objectives) > 3:
            self.PopulationFigure.colorbar(solution.points, label=self.objective_column_headers[self.ordered_list_of_objectives[3]])
    
    
    def force_update(self, event):
        for solution in self.solutions:
            if solution.Legal_Solution:
                solution.points.changed()

    def onpick(self, event):
        ind = event.ind[0]
        if len(self.ordered_list_of_objectives) == 2: #2D
            print '2D picker not implemented'
        elif len(self.ordered_list_of_objectives) >= 3: #3D or 4D plot
            x, y, z = event.artist._offsets3d

            idx = np.where(self.objective_values_matrix[0] == x[ind])
            idy = np.where(self.objective_values_matrix[1] == y[ind])
            idz = np.where(self.objective_values_matrix[2] == z[ind])

            SolutionIndex = np.intersect1d(idx[0], np.intersect1d(idy[0], idz[0]))[0]
            ThisSolution = self.legal_solutions[SolutionIndex]

            print ThisSolution.description

            for objIndex in range(0, len(self.objective_column_headers)):
                if self.objective_column_headers[objIndex] == 'Flight time (days)' and self.TimeUnit == 0:
                    print 'Flight time (years)', ': ', ThisSolution.objective_values[objIndex] / 365.25
                elif self.objective_column_headers[objIndex] == 'Launch epoch (MJD)' and self.EpochUnit == 0:
                    dt = datetime.datetime.fromtimestamp(wx.DateTimeFromJDN(ThisSolution.objective_values[objIndex] + 2400000.5).GetTicks())
                    print 'Launch Epoch (TDB Gregorian):', dt.strftime('%m/%d/%Y')
                elif self.objective_column_headers[objIndex] == 'Thruster preference':
                    print 'Thruster type: ', ThisSolution.thruster
                elif self.objective_column_headers[objIndex] == 'Launch vehicle preference':
                    print 'Launch vehicle: ', ThisSolution.launch_vehicle
                else:
                    print self.objective_column_headers[objIndex], ': ', ThisSolution.objective_values[objIndex]

            print '---------------------------------------------------------------------------------------------'

    def format_date(self, x, pos=None):
        return dates.num2date(x).strftime('%Y-%m-%d')

    def generate_body_prevalence_report(self, Universe):
        #this method generates a list of tuples, (BodyName, NumberOfOccurrences)
        #pass in a Universe object, return a list of tuples
        
        prevalence_report = []
        for body in Universe.bodies:
            number_of_occurrences = 0
            for solution in self.solutions:
                #only book-keep feasible solutions

                if solution.objective_values[0] < 1.0e+99 and body.shortname in solution.description:
                    number_of_occurrences += 1
            prevalence_report.append((body.name, number_of_occurrences))

        return prevalence_report