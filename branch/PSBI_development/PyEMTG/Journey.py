import MissionEvent
import EOM
import AstroFunctions

import math
import copy
import datetime
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import wx
import os
import warnings

import astropy.time

class Journey(object):
    def __init__(self):
        self.missionevents = []
        self.journey_name = "AJourney"
        self.central_body = "Sun"
        self.central_body_radius = 4.379e+6
        self.mu = 132712440017.99
        self.LU = 1.49597870691e+8
        self.TU = 5022642.890912973
        self.boundary_states = []
        self.flyby_periapse_states = []
        self.thruster_duty_cycle = 1.0

    def PlotJourney(self, JourneyAxes, PlotOptions):
        #plot each mission event
        for event in self.missionevents:
            event.PlotEvent(JourneyAxes, self.LU, self.TU, self.mu, PlotOptions)
        

    def PlotJourneyBoundaryOrbits(self, JourneyAxes):

        for boundarystate in self.boundary_states:
            BoundaryStateScaled = np.array(boundarystate) / self.LU
            BoundaryStateScaled[3:6] *= self.TU
            r = np.linalg.norm(BoundaryStateScaled[0:3])
            v = np.linalg.norm(BoundaryStateScaled[3:6])
            a = r / (2.0 - r*v*v)
            T = 2*math.pi*math.sqrt(a**3)

            StateIntegrateObject = ode(EOM.EOM_inertial_2body, jac=EOM.EOM_jacobian_intertial_2body).set_integrator('dop853', atol=1.0e-8, rtol=1.0e-8)
            StateIntegrateObject.set_initial_value(BoundaryStateScaled).set_f_params(1.0).set_jac_params(1.0)

            dt = T / 100
            StateHistory = []
            while StateIntegrateObject.successful() and StateIntegrateObject.t <= T * 1.01:
                StateIntegrateObject.integrate(StateIntegrateObject.t + dt)
                StateHistory.append(StateIntegrateObject.y * self.LU)

            X = []
            Y = []
            Z = []
            for StateLine in StateHistory:
                X.append(StateLine[0])
                Y.append(StateLine[1])
                Z.append(StateLine[2])

            JourneyAxes.plot(X, Y, Z, lw=2, c='0.75')

    def UpdateLabelPosition(self, Figure, Axes):
        for event in self.missionevents:
            event.UpdateLabelPosition(Figure, Axes)

    def PlotPhaseBoundariesOnDataPlot(self, DataAxes, PlotOptions, firstpass):
        #determine the Y limits of the plot
        Ybounds = DataAxes.get_ylim()

        #create vertical lines at important events
        date_string_vector = []
        boundarylegendflag = True
        burnlegendflag = True
        for event in self.missionevents:
            if event.EventType == 'upwr_flyby' or event.EventType == 'pwr_flyby' or event.EventType == 'LT_rndzvs' or event.EventType == 'rendezvous' or event.EventType == 'intercept' or event.EventType == 'insertion' or event.EventType == 'match-vinf' or event.EventType == 'launch' or event.EventType == 'departure' or event.EventType == "begin_spiral" or event.EventType == "end_spiral":
                event_epoch = datetime.datetime.strptime(event.GregorianDate,'%m/%d/%Y').date()
                if firstpass and boundarylegendflag:
                    DataAxes.plot([event_epoch]*2, Ybounds, c='k', marker='+', ls = '-.', lw=3, label='Phase boundary')
                    boundarylegendflag = False
                else:
                    DataAxes.plot([event_epoch]*2, Ybounds, c='k', marker='+', ls = '-.', lw=3)

            if event.EventType == 'chem_burn':
                event_epoch = datetime.datetime.strptime(event.GregorianDate,'%m/%d/%Y').date()
                if firstpass and burnlegendflag:
                    DataAxes.plot([event_epoch]*2, Ybounds, c='r', marker='+', ls = '-.', lw=3, label='Deep-Space Maneuver')
                    burnlegendflag = False
                else:
                    DataAxes.plot([event_epoch]*2, Ybounds, c='r', marker='+', ls = '-.', lw=3)


    def GenerateJourneyDataPlot(self, DataAxes, PlotOptions, firstpass):

        #generate a vector of dates
        date_string_vector = []
        for event in self.missionevents:
            if event.EventType != 'match_point' and event.EventType != 'upwr_flyby' and event.EventType != 'pwr_flyby':
                date_string_vector.append(event.GregorianDate)

        date_vector = [datetime.datetime.strptime(d,'%m/%d/%Y').date() for d in date_string_vector]

        #plot distance from central body
        if PlotOptions.PlotR:
            Rvector = []
            for event in self.missionevents:
                if event.EventType != 'match_point' and event.EventType != 'upwr_flyby' and event.EventType != 'pwr_flyby':
                    Rvector.append(math.sqrt(event.SpacecraftState[0]**2 + event.SpacecraftState[1]**2 + event.SpacecraftState[2]**2) / self.LU)
            if firstpass:
                DataAxes.plot(date_vector, Rvector, c='k', lw=2, label='Distance from central body (LU)')
            else:
                DataAxes.plot(date_vector, Rvector, c='k', lw=2)

        #plot velocity
        if PlotOptions.PlotV:
            Vvector = []
            for event in self.missionevents:
                if event.EventType != 'match_point' and event.EventType != 'upwr_flyby' and event.EventType != 'pwr_flyby':
                    Vvector.append(math.sqrt(event.SpacecraftState[3]**2 + event.SpacecraftState[4]**2 + event.SpacecraftState[5]**2) / self.LU * self.TU)
            if firstpass:
                DataAxes.plot(date_vector, Vvector, c='k', lw=2, ls='-.', label='Velocity magnitude (LU/TU)')
            else:
                DataAxes.plot(date_vector, Vvector, c='k', lw=2, ls='-.')

        #plot Thrust
        if PlotOptions.PlotThrust:
            Thrustvector = []
            for event in self.missionevents:
                if event.EventType != 'match_point' and event.EventType != 'upwr_flyby' and event.EventType != 'pwr_flyby':
                    Thrustvector.append(math.sqrt(event.Thrust[0]**2 + event.Thrust[1]**2 + event.Thrust[2]**2) * 10.0)
            if firstpass:
                DataAxes.plot(date_vector, Thrustvector, c='r', lw=2, ls='-', label='Applied thrust (0.1 N)')
            else:
                DataAxes.plot(date_vector, Thrustvector, c='r', lw=2, ls='-')

        #plot Isp
        if PlotOptions.PlotIsp:
            Ispvector = []
            for event in self.missionevents:
                if event.EventType != 'match_point' and event.EventType != 'upwr_flyby' and event.EventType != 'pwr_flyby':
                    Ispvector.append(event.Isp / 1000.0)
            if firstpass:
                DataAxes.plot(date_vector, Ispvector, c='c', lw=2, ls='-', label='Isp (1000 s)')
            else:
                DataAxes.plot(date_vector, Ispvector, c='c', lw=2, ls='-')

        #plot mass flow rate
        if PlotOptions.PlotMdot:
            Mdotvector = []
            for event in self.missionevents:
                if event.EventType != 'match_point' and event.EventType != 'upwr_flyby' and event.EventType != 'pwr_flyby':
                    Mdotvector.append(event.MassFlowRate * 1.0e6)
            if firstpass:
                DataAxes.plot(date_vector, Mdotvector, c='brown', lw=2, ls='-', label='Mass flow rate (mg/s)')
            else:
                DataAxes.plot(date_vector, Mdotvector, c='brown', lw=2, ls='-')

        #plot Efficiency
        if PlotOptions.PlotEfficiency:
            Efficiencyvector = []
            for event in self.missionevents:
                if event.EventType == 'SFthrust' or event.EventType == 'FBLTthrust' or event.EventType == "PSBIthrust":
                    Efficiencyvector.append(event.AvailableThrust * event.Isp * 9.80665 / (2000 * event.ActivePower))
                elif event.EventType != 'match_point' and event.EventType != 'upwr_flyby' and event.EventType != 'pwr_flyby':
                    Efficiencyvector.append(0.0)
            if firstpass:
                DataAxes.plot(date_vector, Efficiencyvector, c='DarkGreen', lw=2, ls='-', label='Propulsion system efficiency')
            else:
                DataAxes.plot(date_vector, Efficiencyvector, c='DarkGreen', lw=2, ls='-')

        #plot Throttle
        if PlotOptions.PlotThrottle:
            Throttlevector = []
            for event in self.missionevents:
                if event.EventType == 'SFthrust' or event.EventType == 'FBLTthrust' or event.EventType == "PSBIthrust":
                    Throttlevector.append(math.sqrt(event.Thrust[0]**2 + event.Thrust[1]**2 + event.Thrust[2]**2) / (event.AvailableThrust * self.thruster_duty_cycle))
                elif event.EventType != 'match_point' and event.EventType != 'upwr_flyby' and event.EventType != 'pwr_flyby':
                    Throttlevector.append(0.0)
            if firstpass:
                DataAxes.plot(date_vector, Throttlevector, c='r', lw=2, ls='--', label='Throttle')
            else:
                DataAxes.plot(date_vector, Throttlevector, c='r', lw=2, ls='--')

        #plot power
        if PlotOptions.PlotPower:
            Powervector = []
            for event in self.missionevents:
                if event.EventType != 'match_point' and event.EventType != 'upwr_flyby' and event.EventType != 'pwr_flyby':
                    if event.AvailablePower > 0.0:
                        Powervector.append(event.AvailablePower)
                    else:
                        Powervector.append(0.0)
            if firstpass:
                DataAxes.plot(date_vector, Powervector, c='Navy', lw=2, ls='-', label='Power produced by spacecraft (kW)')
            else:
                DataAxes.plot(date_vector, Powervector, c='Navy', lw=2, ls='-')

        #plot gamma
        if PlotOptions.PlotGamma:
            gammavector = []
            for event in self.missionevents:
                if event.EventType == 'SFthrust' or event.EventType == 'FBLTthrust' or event.EventType == "PSBIthrust":
                    gammavector.append(math.atan2(event.Thrust[1], event.Thrust[0]))
                elif event.EventType != 'match_point' and event.EventType != 'upwr_flyby' and event.EventType != 'pwr_flyby':
                    gammavector.append(0.0)
            if firstpass:
                DataAxes.plot(date_vector, gammavector, c='DarkGreen', lw=2, ls='--', label=r'$\gamma$ (radians)')
            else:
                DataAxes.plot(date_vector, gammavector, c='DarkGreen', lw=2, ls='--')

        #plot delta
        if PlotOptions.PlotDelta:
            deltavector = []
            for event in self.missionevents:
                if event.EventType == 'SFthrust' or event.EventType == 'FBLTthrust' or event.EventType == "PSBIthrust":
                    AppliedThrust = math.sqrt(event.Thrust[0]**2 + event.Thrust[1]**2 + event.Thrust[2]**2)
                    deltavector.append(math.asin(event.Thrust[2] / AppliedThrust))
                elif event.EventType != 'match_point' and event.EventType != 'upwr_flyby' and event.EventType != 'pwr_flyby':
                    deltavector.append(0.0)
            if firstpass:
                DataAxes.plot(date_vector, deltavector, c='LightGreen', lw=2, ls='--', label=r'$\delta$ (radians)')
            else:
                DataAxes.plot(date_vector, deltavector, c='LightGreen', lw=2, ls='--')


        #plot central body to thrust vector angle
        if PlotOptions.PlotCB_thrust_angle:
            CBthrustvector = []
            for event in self.missionevents:
                if event.EventType == 'SFthrust' or event.EventType == 'FBLTthrust' or event.EventType == "PSBIthrust":
                    r = math.sqrt(event.SpacecraftState[0]**2 + event.SpacecraftState[1]**2 + event.SpacecraftState[2]**2)
                    AppliedThrust = math.sqrt(event.Thrust[0]**2 + event.Thrust[1]**2 + event.Thrust[2]**2)
                    rdotT = event.SpacecraftState[0]*event.Thrust[0] + event.SpacecraftState[1]*event.Thrust[1] + event.SpacecraftState[2]*event.Thrust[2]
                    CBthrustvector.append( math.acos( rdotT / (r * AppliedThrust) ) )
                elif event.EventType != 'match_point' and event.EventType != 'upwr_flyby' and event.EventType != 'pwr_flyby':
                    CBthrustvector.append(0.0)
            if firstpass:
                DataAxes.plot(date_vector, CBthrustvector, c='Salmon', lw=2, ls='--', label='CB-thrust angle (radians)')
            else:
                DataAxes.plot(date_vector, CBthrustvector, c='Salmon', lw=2, ls='--')

        #plot mass
        if PlotOptions.PlotMass:
            mass = []
            for event in self.missionevents:
                if event.EventType != 'match_point' and event.EventType != 'upwr_flyby' and event.EventType != 'pwr_flyby':
                    mass.append(event.Mass * 1.0e-3)
            if firstpass:
                DataAxes.plot(date_vector, mass, c='DarkGrey', lw=2, ls='-', label='Mass (1000 kg)')
            else:
                DataAxes.plot(date_vector, mass, c='DarkGrey', lw=2, ls='-')

        #plot number of active thrusters
        if PlotOptions.PlotNumberOfEngines:
            numberofengines = []
            for event in self.missionevents:
                if event.EventType == 'SFthrust' or event.EventType == 'FBLTthrust' or event.EventType == "PSBIthrust":
                    numberofengines.append(event.Number_of_Active_Engines)
                elif event.EventType != 'match_point' and event.EventType != 'upwr_flyby' and event.EventType != 'pwr_flyby':
                    numberofengines.append(0)
            if firstpass:
                DataAxes.plot(date_vector, numberofengines, c='Orange', lw=2, ls='-', label='Number of active thrusters')
            else:
                DataAxes.plot(date_vector, numberofengines, c='Orange', lw=2, ls='-')

        #plot power actively used by the thrusters
        if PlotOptions.PlotActivePower:
            activepowervector = []
            for event in self.missionevents:
                if event.EventType == 'SFthrust' or event.EventType == 'FBLTthrust' or event.EventType == "PSBIthrust":
                    activepowervector.append(event.ActivePower)
                elif event.EventType != 'match_point' and event.EventType != 'upwr_flyby' and event.EventType != 'pwr_flyby':
                    activepowervector.append(0.0)

            if firstpass:
                DataAxes.plot(date_vector, activepowervector, c='Navy', lw=2, ls='--', label='Power used by the propulsion system (kW)')
            else:
                DataAxes.plot(date_vector, activepowervector, c='Navy', lw=2, ls='--')

        #plot waste heat from the propulsion system
        if PlotOptions.PlotWasteHeat:
            WasteHeatvector = []
            for event in self.missionevents:
                if event.EventType == 'SFthrust' or event.EventType == 'FBLTthrust' or event.EventType == "PSBIthrust":
                    WasteHeatvector.append( (1 - event.AvailableThrust * event.Isp * 9.80665 / (2000 * event.ActivePower)) * event.ActivePower )
                elif event.EventType != 'match_point' and event.EventType != 'upwr_flyby' and event.EventType != 'pwr_flyby':
                    WasteHeatvector.append(0.0)
            if firstpass:
                DataAxes.plot(date_vector, WasteHeatvector, c='Crimson', lw=2, ls='--', label='Waste heat from propulsion system (kW)')
            else:
                DataAxes.plot(date_vector, WasteHeatvector, c='Crimson', lw=2, ls='--')


    def OutputSTKEphemeris(self, MissionPanel):
        #Step 1: open the journey ephemeris file for writing
        dlg = wx.FileDialog(MissionPanel, "Save", MissionPanel.GetParent().dirname, self.journey_name, '.e', wx.SAVE | wx.FD_OVERWRITE_PROMPT)
        if dlg.ShowModal() == wx.ID_OK:
            saved = True
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
        else:
            saved = False
        
        dlg.Destroy()


        if saved:
            outputfile = open(os.path.join(self.dirname, self.filename), "w")
        else:
            return

        #Step 2: write the STK preamble information
        warnings.filterwarnings('ignore', category=UserWarning, append=True)
        epochstring = astropy.time.Time(self.missionevents[0].JulianDate, format='jd', scale='tdb', out_subfmt='date_hms').utc.iso
        warnings.resetwarnings()
        epocharray = epochstring.split(' ')
        datearray = epocharray[0].split('-')
        time = epocharray[1]
        months = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
        year = datearray[0]
        month = months[int(datearray[1]) - 1]
        day = datearray[2]
        printepochstring = day + ' ' + month + ' ' + year + ' ' + time

        outputfile.write('stk.v.9.0\n')
        outputfile.write('\n')
        outputfile.write('#Ephemeris file written by PyEMTG\n')
        outputfile.write('\n')
        outputfile.write('BEGIN Ephemeris\n')
        outputfile.write('\n')
        outputfile.write('ScenarioEpoch ' + printepochstring + '\n')
        outputfile.write('\n')
        outputfile.write('CentralBody ' + self.central_body + '\n')
        outputfile.write('\n') 
        outputfile.write('CoordinateSystem J2000\n')
        outputfile.write('\n')
        outputfile.write('InterpolationMethod Lagrange\n')
        outputfile.write('\n')
        outputfile.write('InterpolationSamplesM1 5\n')

        #Step 3: assemble a time-history of ephemeris data
        timehistory = []
        timehistory_fromzero = []
        statehistory = []
        rotatedstatehistory = []
        accelhistory = []
        rotatedaccelhistory = []

        for event in self.missionevents:
            event.OutputEphemerisData(self.LU, self.TU, self.mu, timehistory, statehistory, accelhistory)


        #Step 4: if applicable, rotate the state history to the J2000 Earth equatorial coordinate frame
        if self.central_body.lower() == 'sun':
            for line in statehistory:
                rotatedstatehistory.append(AstroFunctions.rotate_from_ecliptic_to_equatorial6(line))
            for line in accelhistory:
                rotatedaccelhistory.append(AstroFunctions.rotate_from_ecliptic_to_equatorial3(line))
        else:
            rotatedstatehistory = statehistory
            rotatedaccelhistory = accelhistory


        #Step 5: index the times from zero
        offset = copy.deepcopy(timehistory[0])
        for TimeEntry in timehistory:
            timehistory_fromzero.append(TimeEntry - offset)

        #Step 6: print the number of time steps
        outputfile.write('\n')
        outputfile.write('NumberOfEphemerisPoints ' + str(len(timehistory)) + '\n')
        outputfile.write('\n')
        outputfile.write('EphemerisTimePosVel\n')

        #Step 7: print the time and state history
        for index in range(0, len(timehistory)):
            outputfile.write('%1.14e'%(timehistory_fromzero[index]))
            for entry in rotatedstatehistory[index]:
                outputfile.write(' %1.14e'%entry)
            outputfile.write('\n')

        outputfile.write('END Ephemeris\n')

        outputfile.close()

        #Step 8: print an acceleration history file
        accelerationfilename = self.filename.split('.')[0] + '_acceleration.e'
        outputfile = open(os.path.join(self.dirname, accelerationfilename), "w")

        outputfile.write('stk.v.9.0\n')
        outputfile.write('\n')
        outputfile.write('#Acceleration file written by PyEMTG\n')
        outputfile.write('#position elements are actually acceleration components\n')
        outputfile.write('\n')
        outputfile.write('BEGIN Ephemeris\n')
        outputfile.write('\n')
        outputfile.write('ScenarioEpoch ' + printepochstring + '\n')
        outputfile.write('\n')
        outputfile.write('CentralBody ' + self.central_body + '\n')
        outputfile.write('\n') 
        outputfile.write('CoordinateSystem J2000\n')
        outputfile.write('\n')
        outputfile.write('InterpolationMethod Lagrange\n')
        outputfile.write('\n')
        outputfile.write('InterpolationSamplesM1 5\n')
        outputfile.write('\n')
        outputfile.write('NumberOfEphemerisPoints ' + str(len(timehistory)) + '\n')
        outputfile.write('\n')
        outputfile.write('EphemerisTimePos\n')

        for index in range(0, len(timehistory)):
            outputfile.write('%1.14e'%(timehistory_fromzero[index]))
            for entry in rotatedaccelhistory[index]:
                outputfile.write(' %1.14e'%entry)
            outputfile.write('\n')

        outputfile.write('END Ephemeris\n')

        outputfile.close()

