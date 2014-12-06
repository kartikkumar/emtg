class PlotOptions(object):
    def __init__(self):
        self.ShowBoundaryOrbits = True
        self.ShowPropagatedTrajectory = True
        self.ShowThrustVectors = True
        self.ShowTextDescriptions = True
        self.PlotR = False
        self.PlotV = False
        self.PlotThrust = False
        self.PlotIsp = False
        self.PlotMdot = False
        self.PlotEfficiency = False
        self.PlotThrottle = False
        self.PlotPower = False
        self.PlotGamma = False
        self.PlotDelta = False
        self.PlotCB_thrust_angle = False
        self.PlotMass = False
        self.PlotNumberOfEngines = False
        self.PlotActivePower = False
        self.PlotWasteHeat = False
        self.PlotCriticalEvents = False
        self.PlotEarthDistance = False
        self.PlotSunEarthSpacecraftAngle = False
        self.FontSize = 10

    def update_mission_panel(self, missionpanel):
        #trajectory plot options
        missionpanel.chkShowBoundaryOrbits.SetValue(self.ShowBoundaryOrbits)
        missionpanel.chkShowPropagatedTrajectory.SetValue(self.ShowPropagatedTrajectory)
        missionpanel.chkShowThrustVectors.SetValue(self.ShowThrustVectors)
        missionpanel.chkShowTextDescriptions.SetValue(self.ShowTextDescriptions)

        #data plot options
        missionpanel.chkPlotR.SetValue(self.PlotR)
        missionpanel.chkPlotV.SetValue(self.PlotV)
        missionpanel.chkPlotThrust.SetValue(self.PlotThrust)
        missionpanel.chkPlotIsp.SetValue(self.PlotIsp)
        missionpanel.chkPlotMdot.SetValue(self.PlotMdot)
        missionpanel.chkPlotEfficiency.SetValue(self.PlotEfficiency)
        missionpanel.chkPlotThrottle.SetValue(self.PlotThrottle)
        missionpanel.chkPlotPower.SetValue(self.PlotPower)
        missionpanel.chkPlotGamma.SetValue(self.PlotGamma)
        missionpanel.chkPlotDelta.SetValue(self.PlotDelta)
        missionpanel.chkPlotCB_thrust_angle.SetValue(self.PlotCB_thrust_angle)
        missionpanel.chkPlotMass.SetValue(self.PlotMass)
        missionpanel.chkPlotEarthDistance.SetValue(self.PlotEarthDistance)
        missionpanel.chkPlotSunEarthSpacecraftAngle.SetValue(self.PlotSunEarthSpacecraftAngle)

        #format options
        missionpanel.spnctrlFontSizeControl.SetValue(self.FontSize)