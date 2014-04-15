import EOM

import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
import copy


class MissionEvent(object):
    EventNumber = [] #index of event
    JulianDate = [] #Julian Ephemeris Date of event
    GregorianDate = [] #gregorian date (MM/DD/YYYY) of event
    EventType = [] #what class of event? Launch? Thrust? Coast? Flyby? etc
    Location = [] #where did the event take place
    TimestepLength = [] #how long is the time step in days?
    Altitude = [] #altitude for flybys
    BdotR = [] #for flybys
    BdotT = [] #for flybys
    RightAscension = [] #RA for maneuvers
    Declination = [] #DEC for maneuvers
    C3 = [] #C3 for arrivals and departures
    SpacecraftState = [] #6-vector
    DeltaVorThrustVectorControl = [] #3-vector
    Thrust = [] #3-vector
    DVmagorThrottle = [] #for impulses or thrust arcs
    AvailableThrust = [] #for thrust arcs
    Isp = [] #for all propulsive maneuvers
    AvailablePower = [] #for thrust arcs
    MassFlowRate = [] #kg/s
    Mass = [] #Mass after event occurs
    Number_of_Active_Engines = [] #number of thrusters firing at the center of this arc
    ActivePower = [] #how much power is currently being used by the thrust system

    def __init__(self, inputcell):
        self.parse_input_line(inputcell)

    def parse_input_line(self, inputcell):
        for i in range(0, len(inputcell)):
            inputcell[i] = inputcell[i].strip(' ')
            if inputcell[i] == '-' or inputcell[i] == '-\n':
                inputcell[i] = 0.0

        

        self.EventNumber = int(inputcell[0])
        self.JulianDate = float(inputcell[1])
        self.GregorianDate = inputcell[2]
        self.EventType = inputcell[3]
        self.Location = inputcell[4]
        self.TimestepLength = float(inputcell[5])
        self.Altitude = float(inputcell[6])
        self.BdotR = float(inputcell[7])
        self.BdotT = float(inputcell[8])
        self.RightAscension = float(inputcell[9])
        self.Declination = float(inputcell[10])
        self.C3 = float(inputcell[11])
        self.SpacecraftState = [float(inputcell[12]), float(inputcell[13]), float(inputcell[14]), float(inputcell[15]), float(inputcell[16]), float(inputcell[17])]
        self.DeltaVorThrustVectorControl = [float(inputcell[18]), float(inputcell[19]), float(inputcell[20])]
        self.Thrust = [float(inputcell[21]), float(inputcell[22]), float(inputcell[23])]
        self.DVmagorThrottle = float(inputcell[24])
        if inputcell[25] == "impulse":
            self.AvailableThrust = 0.0
        else:
            self.AvailableThrust = float(inputcell[25])
        if inputcell[26] == "LV-supplied":
            self.Isp = 0.0
        else:
            self.Isp = float(inputcell[26])
            
        self.AvailablePower = float(inputcell[27])
        self.MassFlowRate = float(inputcell[28])
        self.Mass = float(inputcell[29])

        if len(inputcell) >= 31:
            self.Number_of_Active_Engines = int(inputcell[30])
        
        if len(inputcell) >= 32:
            self.ActivePower = float(inputcell[31])


    def PlotEvent(self, GraphicsObject, LU, TU, mu, PlotOptions):
        #switch between different event types
        if self.EventType == "launch" or self.EventType == "departure":
            GraphicsObject.scatter(self.SpacecraftState[0], self.SpacecraftState[1], self.SpacecraftState[2], s=4, c='g', marker='^')
        elif self.EventType == "begin_spiral" or self.EventType == "end_spiral":
            GraphicsObject.scatter(self.SpacecraftState[0], self.SpacecraftState[1], self.SpacecraftState[2], s=4, c='orange', marker='^')
        elif self.EventType == "upwr_flyby":
            GraphicsObject.scatter(self.SpacecraftState[0], self.SpacecraftState[1], self.SpacecraftState[2], s=4, c='b', marker=r'$\circlearrowleft$')
        elif self.EventType == "pwr_flyby":
            GraphicsObject.scatter(self.SpacecraftState[0], self.SpacecraftState[1], self.SpacecraftState[2], s=4, c='r', marker=r'$\circlearrowleft$')
        elif self.EventType == "chem_burn":
            GraphicsObject.scatter(self.SpacecraftState[0], self.SpacecraftState[1], self.SpacecraftState[2], s=4, c='r', marker='o')
        elif self.EventType == "SFthrust" or self.EventType == "coast" or self.EventType == "force-coast":
            color = ''
            if self.EventType == "coast":
                linestyle = '--'
                color = 'k'
                weight = 1
            elif self.EventType == "SFthrust":
                linestyle = '-'
                color = 'k'
                weight = 2
                if PlotOptions.ShowThrustVectors:
                    ControlVector = np.array(self.Thrust) / self.AvailableThrust * LU * 0.1
                    GraphicsObject.plot(self.SpacecraftState[0] + np.array([0.0, ControlVector[0]]), self.SpacecraftState[1] + np.array([0.0, ControlVector[1]]), self.SpacecraftState[2] + np.array([0.0, ControlVector[2]]), c='r', lw=2)

            else:
                linestyle = '--'
                color = 'b'
                weight = 1

            GraphicsObject.scatter(self.SpacecraftState[0], self.SpacecraftState[1], self.SpacecraftState[2], s=1, c=color, marker='o')

            if PlotOptions.ShowPropagatedTrajectory:
                CenterPointState = np.zeros(6)
                CenterPointState[0:6] = np.array(self.SpacecraftState) / LU
                CenterPointState[3:6] *= TU

                ForwardIntegrateObject = ode(EOM.EOM_inertial_2body).set_integrator('dop853', atol=1.0e-8, rtol=1.0e-8)
                ForwardIntegrateObject.set_initial_value(CenterPointState).set_f_params(1.0)

                dt = self.TimestepLength * 86400 / TU / 10
                StateHistoryForward = []
                while ForwardIntegrateObject.successful() and ForwardIntegrateObject.t < self.TimestepLength * 86400 / TU / 2.0:
                    ForwardIntegrateObject.integrate(ForwardIntegrateObject.t + dt)
                    StateHistoryForward.append(ForwardIntegrateObject.y * LU)

                BackwardIntegrateObject = ode(EOM.EOM_inertial_2body).set_integrator('dop853', atol=1.0e-8, rtol=1.0e-8)
                BackwardIntegrateObject.set_initial_value(CenterPointState).set_f_params(1.0)

                dt = self.TimestepLength * 86400 / TU / 10
                StateHistoryBackward = []
                while BackwardIntegrateObject.successful() and BackwardIntegrateObject.t > -self.TimestepLength * 86400 / TU / 2.0:
                    BackwardIntegrateObject.integrate(BackwardIntegrateObject.t - dt)
                    StateHistoryBackward.append(BackwardIntegrateObject.y * LU)

                StateHistoryBackward.reverse()

                X = []
                Y = []
                Z = []
                for StateLine in StateHistoryBackward:
                    X.append(StateLine[0])
                    Y.append(StateLine[1])
                    Z.append(StateLine[2])
                for StateLine in StateHistoryForward:
                    X.append(StateLine[0])
                    Y.append(StateLine[1])
                    Z.append(StateLine[2])

                GraphicsObject.plot(X, Y, Z, lw=weight, c=color, ls=linestyle)

        elif self.EventType == "FBLTthrust" or self.EventType == "FBLTSthrust":
            GraphicsObject.scatter(self.SpacecraftState[0], self.SpacecraftState[1], self.SpacecraftState[2], s=2, c='k', marker='o')
            
            if PlotOptions.ShowThrustVectors:
                ControlVector = np.array(self.Thrust) / self.AvailableThrust * LU * 0.1
                GraphicsObject.plot(self.SpacecraftState[0] + np.array([0.0, ControlVector[0]]), self.SpacecraftState[1] + np.array([0.0, ControlVector[1]]), self.SpacecraftState[2] + np.array([0.0, ControlVector[2]]), c='r', lw=2)


            if PlotOptions.ShowPropagatedTrajectory:
                CenterPointState = np.zeros(7)
                CenterPointState[0:6] = np.array(self.SpacecraftState) / LU
                CenterPointState[3:6] *= TU
                CenterPointState[6] = 1.0
                ScaledThrust = np.array(self.Thrust) / self.Mass / LU / 1000* TU*TU
                ScaledMdot = self.MassFlowRate / self.Mass * TU

                ForwardIntegrateObject = ode(EOM.EOM_inertial_2bodyconstant_thrust).set_integrator('dop853', atol=1.0e-8, rtol=1.0e-8)
                ForwardIntegrateObject.set_initial_value(CenterPointState).set_f_params(ScaledThrust, ScaledMdot, 1.0)

                dt = self.TimestepLength * 86400 / TU / 10
                StateHistoryForward = []
                while ForwardIntegrateObject.successful() and ForwardIntegrateObject.t < self.TimestepLength * 86400 / TU / 2.0:
                    ForwardIntegrateObject.integrate(ForwardIntegrateObject.t + dt)
                    StateHistoryForward.append(ForwardIntegrateObject.y * LU)

                BackwardIntegrateObject = ode(EOM.EOM_inertial_2bodyconstant_thrust).set_integrator('dop853', atol=1.0e-8, rtol=1.0e-8)
                BackwardIntegrateObject.set_initial_value(CenterPointState).set_f_params(ScaledThrust, ScaledMdot, 1.0)

                dt = self.TimestepLength * 86400 / TU / 10
                StateHistoryBackward = []
                while BackwardIntegrateObject.successful() and BackwardIntegrateObject.t > -self.TimestepLength * 86400 / TU / 2.0:
                    BackwardIntegrateObject.integrate(BackwardIntegrateObject.t - dt)
                    StateHistoryBackward.append(BackwardIntegrateObject.y * LU)

                StateHistoryBackward.reverse()

                X = []
                Y = []
                Z = []
                for StateLine in StateHistoryBackward:
                    X.append(StateLine[0])
                    Y.append(StateLine[1])
                    Z.append(StateLine[2])
                for StateLine in StateHistoryForward:
                    X.append(StateLine[0])
                    Y.append(StateLine[1])
                    Z.append(StateLine[2])

                GraphicsObject.plot(X, Y, Z, lw=2, c='k')

        elif self.EventType == "insertion":
            GraphicsObject.scatter(self.SpacecraftState[0], self.SpacecraftState[1], self.SpacecraftState[2], s=4, c='r', marker='v')
        elif self.EventType == "LT_rndzvs":
            GraphicsObject.scatter(self.SpacecraftState[0], self.SpacecraftState[1], self.SpacecraftState[2], s=4, c='m', marker='o')
        elif self.EventType == "intercept":
            GraphicsObject.scatter(self.SpacecraftState[0], self.SpacecraftState[1], self.SpacecraftState[2], s=4, c='m', marker='d')
        elif self.EventType == "rendezvous":
            GraphicsObject.scatter(self.SpacecraftState[0], self.SpacecraftState[1], self.SpacecraftState[2], s=4, c='r', marker='d')
        elif self.EventType == "match-vinf":
            GraphicsObject.scatter(self.SpacecraftState[0], self.SpacecraftState[1], self.SpacecraftState[2], s=4, c='c', marker='s')
        elif self.EventType == "match_point":
            GraphicsObject.scatter(self.SpacecraftState[0], self.SpacecraftState[1], self.SpacecraftState[2], s=4, c='k', marker='s')

        if PlotOptions.ShowTextDescriptions and not (self.EventType == 'coast' or self.EventType == 'force-coast' or self.EventType == "SFthrust" or self.EventType == "FBLTthrust" or self.EventType == "FBLTSthrust" or self.EventType == "match_point"):
            self.LabelEvent(GraphicsObject, PlotOptions)

    def LabelEvent(self, GraphicsObject, PlotOptions):
        EventTypeFormatted = self.EventType
        if EventTypeFormatted == 'upwr_flyby':
            EventTypeFormatted = 'unpowered flyby'
        elif EventTypeFormatted == 'pwr_flyby':
            EventTypeFormatted = 'powered flyby'
        elif EventTypeFormatted == 'chem_burn':
            EventTypeFormatted = 'chemical burn'
        elif EventTypeFormatted == 'LT_rndzvs':
            EventTypeFormatted = 'LT rendezvous'
        elif EventTypeFormatted == 'begin_spiral':
            EventTypeFormatted = 'begin spiral'
        elif EventTypeFormatted == 'end_spiral':
            EventTypeFormatted = 'end spiral'

        description = EventTypeFormatted + '\n' + self.Location + '\n' + self.GregorianDate

        #for launches a C3 and DLA are needed
        if self.EventType == 'launch':
            description += '\nC3 = ' + "{0:.3f}".format(self.C3) + ' $km^2/s^2$'
            #add the LV to the description?
            description += '\nDLA = ' + "{0:.1f}".format(self.Declination) + '$^{\circ}$'

        #for non-launch departures only the C3 is needed
        if self.EventType == 'departure':
            description += '\nC3 = ' + "{0:.3f}".format(self.C3) + ' $km^2/s^2$'

        #for spirals output only the delta-v
        if self.EventType == 'begin_spiral' or self.EventType == 'end_spiral':
            description += '\n$\Delta v$ = ' + "{0:.3f}".format(self.DVmagorThrottle) + ' $km/s$'

        #for other events, output v-infinity and DLA
        if self.EventType == 'upwr_flyby' or self.EventType == 'pwr_flyby' or self.EventType == 'intercept' or self.EventType == 'interface' or self.EventType == 'insertion':
            description += '\n$v_\infty$ = ' + "{0:.3f}".format(math.sqrt(self.C3)) + ' $km/s$'
            description += '\nDEC = ' + "{0:.1f}".format(self.Declination) + '$^{\circ}$'

        #for flybys, altitude should be outputed
        if self.EventType == 'upwr_flyby' or self.EventType == 'pwr_flyby':
            description += '\naltitude = ' + "{0:.0f}".format(self.Altitude) + ' $km$'


        #for propulsive events, a deltaV is needed
        if self.EventType == 'departure' or self.EventType == 'pwr_flyby' or self.EventType == 'insertion' or self.EventType == 'chem_burn' or self.EventType == 'rendezvous':
                description += '\n$\Delta v$ = ' + "{0:.3f}".format(self.DVmagorThrottle) + ' $km/s$'

        #always append the spacecraft Mass
        description += '\nm = ' + "{0:.0f}".format(self.Mass) + ' $kg$'

        #draw the text
        #note, do not draw anything for chemical burns below 10 m/s
        if not (self.EventType == "chem_burn" and self.DVmagorThrottle < 0.001):
            x2D, y2D, _ = proj3d.proj_transform(self.SpacecraftState[0],self.SpacecraftState[1],self.SpacecraftState[2], GraphicsObject.get_proj())

            self.eventlabel = plt.annotate(description, xycoords = 'data', xy = (x2D, y2D), xytext = (20, 20), textcoords = 'offset points', ha = 'left', va = 'bottom',
                                           bbox = dict(boxstyle = 'round,pad=0.5', fc = 'white', alpha = 0.5), arrowprops = dict(arrowstyle = '->',
                                           connectionstyle = 'arc3,rad=0'), size=PlotOptions.FontSize)

            self.AnnotationHelper = self.eventlabel.draggable(use_blit=True)
            self.pcid = GraphicsObject.figure.canvas.mpl_connect('button_press_event', self.ClickAnnotation)
            self.rcid = GraphicsObject.figure.canvas.mpl_connect('button_release_event', self.ReleaseAnnotation)

    def UpdateLabelPosition(self, Figure, Axes):
        if not (self.EventType == 'coast' or self.EventType == 'force-coast' or self.EventType == "SFthrust" or self.EventType == "FBLTthrust" or self.EventType == "FBLTSthrust" or self.EventType == "match_point" or (self.EventType == "chem_burn" and self.DVmagorThrottle < 0.001)):
            x2, y2, _ = proj3d.proj_transform(self.SpacecraftState[0],self.SpacecraftState[1],self.SpacecraftState[2], Axes.get_proj())
            self.eventlabel.xy = x2,y2
            self.eventlabel.update_positions(Figure.canvas.renderer)

    def ClickAnnotation(self, event):
        if event.inaxes != self.eventlabel.axes: return

        contains, attrd = self.eventlabel.contains(event)
        if not contains: return
        self.eventlabel.axes.disable_mouse_rotation()

    def ReleaseAnnotation(self, event):
        if event.inaxes != self.eventlabel.axes: return

        contains, attrd = self.eventlabel.contains(event)
        if not contains: return
        self.eventlabel.axes.mouse_init()

    def OutputEphemerisData(self, LU, TU, mu, timehistory, statehistory, accelhistory):
        if self.EventType == 'coast' or self.EventType == 'force-coast' or self.EventType == "SFthrust" or self.EventType == "FBLTthrust" or self.EventType == "FBLTSthrust" or self.EventType == "match_point":

            #propagate the state forward and back
            CenterPointState = []
            CenterPointState = np.zeros(7)
            CenterPointState[0:6] = np.array(self.SpacecraftState) / LU
            CenterPointState[3:6] *= TU
            CenterPointState[6] = 1.0

            if self.EventType == 'coast' or self.EventType == 'force-coast' or self.EventType == "SFthrust":
                CenterPointState = np.zeros(6)
                CenterPointState[0:6] = np.array(self.SpacecraftState) / LU
                CenterPointState[3:6] *= TU

                ForwardIntegrateObject = ode(EOM.EOM_inertial_2body).set_integrator('dop853', atol=1.0e-8, rtol=1.0e-8)
                ForwardIntegrateObject.set_initial_value(CenterPointState).set_f_params(1.0)

                StateHistoryForward = []
                TimeHistoryForward = []
                while ForwardIntegrateObject.successful() and ForwardIntegrateObject.t < (self.TimestepLength/2.0 - 1) * 86400 / TU:
                    ForwardIntegrateObject.integrate(ForwardIntegrateObject.t + 86400 / TU)
                    StateHistoryForward.append(ForwardIntegrateObject.y * LU)
                    TimeHistoryForward.append(ForwardIntegrateObject.t * TU)

                BackwardIntegrateObject = ode(EOM.EOM_inertial_2body).set_integrator('dop853', atol=1.0e-8, rtol=1.0e-8)
                BackwardIntegrateObject.set_initial_value(CenterPointState).set_f_params(1.0)

                StateHistoryBackward = []
                TimeHistoryBackward = []
                while BackwardIntegrateObject.successful() and BackwardIntegrateObject.t > -(self.TimestepLength/2.0 - 1) * 86400 / TU:
                    BackwardIntegrateObject.integrate(BackwardIntegrateObject.t - 86400 / TU)
                    StateHistoryBackward.append(BackwardIntegrateObject.y * LU)
                    TimeHistoryBackward.append(BackwardIntegrateObject.t * TU)

                StateHistoryBackward.reverse()
                TimeHistoryBackward.reverse()

                for StateLine in StateHistoryBackward:
                    statehistory.append(np.array([StateLine[0], StateLine[1], StateLine[2], StateLine[3] / TU, StateLine[4] / TU, StateLine[5] / TU]) * 1000)
                    accelhistory.append(np.array([0.0, 0.0, 0.0]))
                for TimeLine in TimeHistoryBackward:
                    timehistory.append(TimeLine + self.JulianDate * 86400)

                statehistory.append(np.array([self.SpacecraftState[0], self.SpacecraftState[1], self.SpacecraftState[2], self.SpacecraftState[3], self.SpacecraftState[4], self.SpacecraftState[5]]) * 1000)
                accelhistory.append(np.array([0.0, 0.0, 0.0]))
                timehistory.append(self.JulianDate * 86400)

                for StateLine in StateHistoryForward:
                    statehistory.append(np.array([StateLine[0], StateLine[1], StateLine[2], StateLine[3] / TU, StateLine[4] / TU, StateLine[5] / TU]) * 1000)
                    accelhistory.append(np.array([0.0, 0.0, 0.0]))
                for TimeLine in TimeHistoryForward:
                    timehistory.append(TimeLine + self.JulianDate * 86400)
            elif self.EventType == 'FBLTthrust':
                CenterPointState = np.zeros(7)
                CenterPointState[0:6] = np.array(self.SpacecraftState) / LU
                CenterPointState[3:6] *= TU
                CenterPointState[6] = 1.0

                ScaledThrust = np.array(self.Thrust) / self.Mass / LU / 1000* TU*TU
                ScaledMdot = self.MassFlowRate / self.Mass * TU

                
                ForwardIntegrateObject = ode(EOM.EOM_inertial_2bodyconstant_thrust).set_integrator('dop853', atol=1.0e-8, rtol=1.0e-8)
                ForwardIntegrateObject.set_initial_value(CenterPointState).set_f_params(ScaledThrust, ScaledMdot, 1.0)

                StateHistoryForward = []
                TimeHistoryForward = []
                while ForwardIntegrateObject.successful() and ForwardIntegrateObject.t < (self.TimestepLength/2.0 - 1) * 86400 / TU:
                    ForwardIntegrateObject.integrate(ForwardIntegrateObject.t + 86400 / TU)
                    UnscaledState = copy.deepcopy(ForwardIntegrateObject.y)
                    UnscaledState[0:3] *= LU
                    UnscaledState[4:6] *= LU/TU
                    UnscaledState[6] *= self.Mass
                    StateHistoryForward.append(UnscaledState)
                    TimeHistoryForward.append(ForwardIntegrateObject.t * TU)

                BackwardIntegrateObject = ode(EOM.EOM_inertial_2bodyconstant_thrust).set_integrator('dop853', atol=1.0e-8, rtol=1.0e-8)
                BackwardIntegrateObject.set_initial_value(CenterPointState).set_f_params(ScaledThrust, ScaledMdot, 1.0)

                StateHistoryBackward = []
                TimeHistoryBackward = []
                while BackwardIntegrateObject.successful() and BackwardIntegrateObject.t > -(self.TimestepLength/2.0 - 1) * 86400 / TU:
                    BackwardIntegrateObject.integrate(BackwardIntegrateObject.t - 86400 / TU)
                    UnscaledState = copy.deepcopy(BackwardIntegrateObject.y)
                    UnscaledState[0:3] *= LU
                    UnscaledState[4:6] *= LU/TU
                    UnscaledState[6] *= self.Mass
                    StateHistoryBackward.append(UnscaledState)
                    TimeHistoryBackward.append(BackwardIntegrateObject.t * TU)

                StateHistoryBackward.reverse()
                TimeHistoryBackward.reverse()

                for StateLine in StateHistoryBackward:
                    statehistory.append(np.array([StateLine[0], StateLine[1], StateLine[2], StateLine[3], StateLine[4], StateLine[5]]) * 1000)
                    accelhistory.append(np.array([self.Thrust[0] / StateLine[6], self.Thrust[1] / StateLine[6], self.Thrust[2] / StateLine[6]]))
                for TimeLine in TimeHistoryBackward:
                    timehistory.append(TimeLine + self.JulianDate * 86400)

                statehistory.append(np.array([self.SpacecraftState[0], self.SpacecraftState[1], self.SpacecraftState[2], self.SpacecraftState[3], self.SpacecraftState[4], self.SpacecraftState[5]]) * 1000)
                accelhistory.append(np.array([self.Thrust[0] / self.Mass, self.Thrust[1] / self.Mass, self.Thrust[2] / self.Mass]))
                timehistory.append(self.JulianDate * 86400)

                for StateLine in StateHistoryForward:
                    statehistory.append(np.array([StateLine[0], StateLine[1], StateLine[2], StateLine[3], StateLine[4], StateLine[5]]) * 1000)
                    accelhistory.append(np.array([self.Thrust[0] / StateLine[6], self.Thrust[1] / StateLine[6], self.Thrust[2] / StateLine[6]]))
                for TimeLine in TimeHistoryForward:
                    timehistory.append(TimeLine + self.JulianDate * 86400)

        else:
            statehistory = np.array(self.SpacecraftState[0:6])
            accelhistory = np.array([0.0, 0.0, 0.0])
            timehistory = self.JulianDate * 86400

        return timehistory, statehistory, accelhistory