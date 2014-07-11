#classes for small body list
#for use in finding science flybys
#Jacob Englander 7-9-2014
#expected interface with Ryne Beeson's bubble search code

import math
import os
from kepler import kepler
from numpy import array

class SmallBody(object):
    def __init__(self, name, SPICEID, ReferenceEpoch, SMA, ECC, INC, RAAN, AOP, MA, Tholen, H, Diameter):
        self.set_values(name, SPICEID, ReferenceEpoch, SMA, ECC, INC, RAAN, AOP, MA, Tholen, H, Diameter)

    def set_values(self, name, SPICEID, ReferenceEpoch, SMA, ECC, INC, RAAN, AOP, MA, Tholen, H, Diameter):
        self.name = name
        self.SPICEID = SPICEID
        self.ReferenceEpoch = ReferenceEpoch
        self.SMA = SMA
        self.ECC = ECC
        self.INC = INC
        self.RAAN = RAAN
        self.AOP = AOP
        self.MA = MA
        self.Tholen = Tholen
        self.H = H
        self.Diameter = Diameter

        self.Aphelion = self.SMA * (1.0 + self.ECC)
        self.Perihelion = self.SMA * (1.0 - self.ECC)

    def locate_body(self, JD, mu, AU):
        #call Ryne's Kepler solver
        r, v = kepler(self.SMA * AU, self.ECC, self.INC, self.RAAN, self.AOP, self.MA, self.ReferenceEpoch, JD, mu)
       
        return r, v

class SmallBodyList(object):
    def __init__(self, filename, mu, LU):
        self.body_list = []
        self.LU = LU
        self.mu = mu
        self.read_small_body_list(filename)

    def read_small_body_list(self, filename):
        
        #read the small body list file
        if os.path.isfile(filename):
            inputfile = open(filename, "r")
        else:
            print "File ", inputfile, " does not exist!"
            return

        for line in inputfile:
            linecell = line.split(',')
            #change empty values to 0 so that float() and int() don't blow up
            for i in range(0, len(linecell)):
                if linecell[i].strip(' ') == '':
                    linecell[i] = '-1'
            if linecell[-1] == '\n':
                linecell.pop()
            while len(linecell) < 12:
                linecell.append('-1')
            linecell[-1] = linecell[-1] + '\n'

            #skip the first line
            if not linecell[0] == 'full_name':
                self.body_list.append(SmallBody(linecell[0].lstrip(), #name
                                                int(linecell[1]), #SPICE ID
                                                float(linecell[2]), #reference epoch
                                                float(linecell[3]), #SMA (AU)
                                                float(linecell[4]), #ECC
                                                float(linecell[5]) * math.pi/180.0, #INC
                                                float(linecell[6]) * math.pi/180.0, #RAAN
                                                float(linecell[7]) * math.pi/180.0, #AOP
                                                float(linecell[8]) * math.pi/180.0, #MA
                                                linecell[9],  #Tholen spectral type
                                                float(linecell[10]), #absolute magnitude
                                                float(linecell[11]))) #diameter (km)


    #function to find bodies that are near a reference state at a reference epoch
    def find_bubble_targets(self, ReferenceState, ReferenceJD, RelativePositionFilterMagnitude, RelativeVelocityFilterMagnitude, MaximumMagnitude):
        #initialize the list of bodies
        list_of_acceptable_bodies = []
        
        #find the reference sun distance
        ReferenceSunDistance = math.sqrt(ReferenceState[0]**2 + ReferenceState[1]**2 + ReferenceState[2]**2) / self.LU

        for body in self.body_list:
           #first filter is to discard all objects whose aphelion is smaller than the reference sun distance or whose perihelion is larger than the reference sun distance
           if body.Aphelion >= (ReferenceSunDistance - RelativePositionFilterMagnitude / self.LU) \
               and body.Perihelion <= (ReferenceSunDistance + RelativePositionFilterMagnitude / self.LU) \
               and body.H <= MaximumMagnitude:
               
               #second and third filters are on relative position and velocity, so we have to get the current state vector
               BodyR, BodyV = body.locate_body(ReferenceJD, self.mu, self.LU)

               #compute relative position
               RelativePosition = [BodyR[0] - ReferenceState[0], BodyR[1] - ReferenceState[1], BodyR[2] - ReferenceState[2]]
               body.RelativePositionMagnitude = math.sqrt(RelativePosition[0]**2 + RelativePosition[1]**2 + RelativePosition[2]**2)

               #compute relative velocity
               RelativeVelocity = [BodyV[0] - ReferenceState[3], BodyV[1] - ReferenceState[4], BodyV[2] - ReferenceState[5]]
               body.RelativeVelocityMagnitude = math.sqrt(RelativeVelocity[0]**2 + RelativeVelocity[1]**2 + RelativeVelocity[2]**2)

               if body.RelativePositionMagnitude < RelativePositionFilterMagnitude and body.RelativeVelocityMagnitude < RelativeVelocityFilterMagnitude:
                   list_of_acceptable_bodies.append(body)

        return list_of_acceptable_bodies