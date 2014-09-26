#class for reading and parsing electric propulsion throttle tables

import os
import math

class ThrottleSetting(object):
    #constructor for parsing a string from a throttle file
    def __init__(self):
        self.clear()

    def initialize_from_string(self, throttlestring):
        self.clear()
        self.parse_throttle_line(throttlestring)

    #constructor when the propulsion information is known
    def initialize_from_input_data(self, Powerin, Thrust, Mdot, Isp, efficiency):
        self.clear()
        self.TL = -1
        self.Powerin = Powerin
        self.Thrust = Thrust
        self.Mdot = Mdot
        self.Isp = Isp
        self.efficiency = efficiency

    def clear(self):
        self.TL = [] 
        self.Powerin = [] #kW 
        self.Thrust = [] #mN      
        self.Mdot = [] #mg/s
        self.Isp = []#s
        self.efficiency = []

    #parse an input line
    def parse_throttle_line(self, throttlestring):
        throttlecell = throttlestring.split(',')
        self.TL = int(throttlecell[0].lstrip('TL'))
        self.Powerin = float(throttlecell[1])
        self.Thrust = float(throttlecell[2])
        self.Mdot = float(throttlecell[3])
        self.Isp = float(throttlecell[4])
        self.efficiency = float(throttlecell[5].strip('\n'))

    #function to compare this throttle setting to another
    def compare(self, otherThrottleSetting):
        diff_Powerin = self.Powerin - otherThrottleSetting.Powerin
        diff_Thrust = self.Thrust - otherThrottleSetting.Thrust
        diff_Mdot = self.Mdot - otherThrottleSetting.Mdot
        diff_Isp = self.Isp - otherThrottleSetting.Isp
        diff_efficiency = self.efficiency - otherThrottleSetting.Isp

        return diff_Powerin, diff_Thrust, diff_Mdot, diff_Isp, diff_efficiency


class ThrottleTable(object):
    def __init__(self, filename):
        self.clear()
        self.parse_throttle_file(filename)

    def clear(self):
        self.ThrottleSettings = []
        self.PPUefficiency = 1.0

    #function to parse a throttle file
    def parse_throttle_file(self, filename):
        self.clear()

        #read the throttle file
        if os.path.isfile(filename):
            inputfile = open(filename, "r")
        else:
            print "File ", inputfile, " does not exist!"
            return

        for line in inputfile:
            if line[0:3] == 'PPU':
                linecell = line.split(',')
                self.PPUefficiency = float(linecell[1])

            if line[0:2] == 'TL':
                NewThrottleSetting = ThrottleSetting()
                NewThrottleSetting.initialize_from_string(line)
                self.ThrottleSettings.append(NewThrottleSetting)

    #function to find the closest throttle setting to a reference performance point
    #this can be used in conjunction with an EMTG output file to find the throttle setting closest to what the optimizer is asking for at that time step
    def find_closest_throttle_setting(self, ReferenceThrottleSetting):
        #start by initializing the throttle level to zero
        TL = 0

        #we want the solution to the minimum-norm problem comparing the reference throttle setting to the throttle table
        #first we need to compute the normalized difference vector between the reference throttle setting and each throttle setting in the table
        normalized_difference_vectors = []
        normalized_distance_array = []

        for setting in self.ThrottleSettings:
            difference_vector = setting.compare(ReferenceThrottleSetting)
            normalized_difference_vectors.append([difference_vector[0] / ReferenceThrottleSetting.Powerin,
                                                  difference_vector[1] / ReferenceThrottleSetting.Thrust,
                                                  difference_vector[2] / ReferenceThrottleSetting.Mdot,
                                                  difference_vector[3] / ReferenceThrottleSetting.Isp])
            normalized_distance_array.append(math.sqrt(normalized_difference_vectors[-1][0]**2 +
                                                       normalized_difference_vectors[-1][1]**2 +
                                                       normalized_difference_vectors[-1][2]**2 +
                                                       normalized_difference_vectors[-1][3]**2))
        
        
        return self.ThrottleSettings[normalized_distance_array.index(min(normalized_distance_array))].TL