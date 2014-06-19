import Mission
import PlotOptions
import wx   

class MissionPanel(wx.Panel):
    #mission panel
    #contains post-processing options

    def __init__(self, parent, mission):
        wx.Panel.__init__(self, parent)

        self.mission = mission
        self.plotoptions = PlotOptions.PlotOptions()

        #Journey selection listbox
        #first we need to create an array of journey names
        self.journeynamelist = []
        for journey in self.mission.Journeys:
            self.journeynamelist.append(journey.journey_name)
        
        if len(self.journeynamelist) > 1:
            identicaljourneys = True
            for OtherJourney in self.mission.Journeys[1:len(self.mission.Journeys)]:
                if OtherJourney.central_body.lower() != self.mission.Journeys[0].central_body.lower():
                    identicaljourneys = False
            if identicaljourneys:
                self.journeynamelist.append("All journeys")

        self.mission.ActiveJourney = 0

        self.lblJourneyList = wx.StaticText(self, -1, "Choose a journey")
        self.JourneyListBox = wx.ListBox(self, -1, choices = self.journeynamelist, size=(300, -1))

        self.TrajectoryPlotBox = wx.StaticBox(self, -1, "Trajectory plot options")
        self.chkShowBoundaryOrbits = wx.CheckBox(self, -1, "Show boundary orbits", size=(300, -1))
        self.chkShowPropagatedTrajectory = wx.CheckBox(self, -1, "Show progagated trajectory", size=(300, -1))
        self.chkShowThrustVectors = wx.CheckBox(self, -1, "Show thust vectors", size=(300, -1))
        self.chkShowTextDescriptions = wx.CheckBox(self, -1, "Show text descriptions", size=(300, -1))

        self.btnPlotJourney = wx.Button(self, -1, "Plot trajectory", size = (300, -1))

        TrajectoryPlotBoxSizer = wx.StaticBoxSizer(self.TrajectoryPlotBox, wx.VERTICAL)
        TrajectoryPlotBoxSizer.AddMany([self.chkShowBoundaryOrbits, self.chkShowPropagatedTrajectory, self.chkShowThrustVectors, self.chkShowTextDescriptions, self.btnPlotJourney])
        font = self.GetFont()
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        self.TrajectoryPlotBox.SetFont(font)

        self.btnOutputSTKEphemeris = wx.Button(self, -1, "Output STK Ephemeris", (300, -1))
        
        LeftHandSizer = wx.BoxSizer(wx.VERTICAL)
        LeftHandSizer.AddMany([self.lblJourneyList, self.JourneyListBox, TrajectoryPlotBoxSizer, self.btnOutputSTKEphemeris])

        self.DataPlotBox = wx.StaticBox(self, -1, "Data Plots", size = (300, 300))
        DataPlotSizer = wx.StaticBoxSizer(self.DataPlotBox, wx.VERTICAL)

        self.chkPlotR = wx.CheckBox(self, -1, "distance from central body")
        self.chkPlotV = wx.CheckBox(self, -1, "velocity magnitude with respect to central body")
        self.chkPlotThrust = wx.CheckBox(self, -1, "applied thrust (N)")
        self.chkPlotIsp = wx.CheckBox(self, -1, "specific impulse (s)")
        self.chkPlotMdot = wx.CheckBox(self, -1, "mass flow rate (kg/s)")
        self.chkPlotEfficiency = wx.CheckBox(self, -1, "propulsion system efficiency")
        self.chkPlotThrottle = wx.CheckBox(self, -1, "throttle")
        self.chkPlotPower = wx.CheckBox(self, -1, "power produced by spacecraft (kW)")
        self.chkPlotGamma = wx.CheckBox(self, -1, "in-plane control angle (degrees)")
        self.chkPlotDelta = wx.CheckBox(self, -1, "out-of-plane control angle (degrees)")
        self.chkPlotCB_thrust_angle = wx.CheckBox(self, -1, "central body - thrust vector angle")
        self.chkPlotMass = wx.CheckBox(self, -1, "mass (kg)")
        self.chkPlotNumberOfEngines = wx.CheckBox(self, -1, "number of active thrusters")
        self.chkPlotActivePower = wx.CheckBox(self, -1, "power used by propulsion system (kW)")
        self.chkPlotWasteHeat = wx.CheckBox(self, -1, "propulsion system waste heat (kW)")
        self.chkPlotCriticalEvents = wx.CheckBox(self, -1, "mark critical events")


        DataPlotGrid = wx.GridSizer(5, 2, 5, 5)
        DataPlotGrid.AddMany([self.chkPlotR, self.chkPlotV,
                                self.chkPlotThrust, self.chkPlotIsp,
                                self.chkPlotMdot, self.chkPlotEfficiency,
                                self.chkPlotThrottle, self.chkPlotPower,
                                self.chkPlotGamma, self.chkPlotDelta,
                                self.chkPlotCB_thrust_angle, self.chkPlotMass,
                                self.chkPlotNumberOfEngines, self.chkPlotActivePower,
                                self.chkPlotWasteHeat, self.chkPlotCriticalEvents])

        self.btnGenerateDataPlot = wx.Button(self, -1, "Generate plot")

        DataPlotSizer.AddMany([DataPlotGrid, self.btnGenerateDataPlot])

        self.FormatBox = wx.StaticBox(self, -1, "Common plot options")
        self.FormatBox.SetFont(font)
        FormatBoxSizer = wx.StaticBoxSizer(self.FormatBox, wx.HORIZONTAL)
        self.lblFontSize = wx.StaticText(self, -1, "Font size")
        self.spnctrlFontSizeControl = wx.SpinCtrl(self, -1, min=1, max=100, initial=10, name="Font size")
        FormatBoxSizer.AddMany([self.lblFontSize, self.spnctrlFontSizeControl])

        RightHandSizer = wx.BoxSizer(wx.VERTICAL)
        RightHandSizer.AddMany([DataPlotSizer, FormatBoxSizer])

        mainbox = wx.BoxSizer(wx.HORIZONTAL)
        mainbox.AddMany([LeftHandSizer, RightHandSizer])

        self.SetSizer(mainbox)

        font = self.GetFont()
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        self.lblJourneyList.SetFont(font)
        self.DataPlotBox.SetFont(font)
        self.TrajectoryPlotBox.SetFont(font)

        self.JourneyListBox.SetSelection(self.mission.ActiveJourney)
        self.plotoptions.update_mission_panel(self)

        #trajectory plot bindings
        self.JourneyListBox.Bind(wx.EVT_LISTBOX, self.ClickJourneyListBox)
        self.btnPlotJourney.Bind(wx.EVT_BUTTON, self.ClickPlotMissionButton)
        self.chkShowBoundaryOrbits.Bind(wx.EVT_CHECKBOX, self.ChangeShowBoundaryOrbits)
        self.chkShowPropagatedTrajectory.Bind(wx.EVT_CHECKBOX, self.ChangeShowPropagatedTrajectory)
        self.chkShowThrustVectors.Bind(wx.EVT_CHECKBOX, self.ChangeShowThrustVectors)
        self.chkShowTextDescriptions.Bind(wx.EVT_CHECKBOX, self.ChangeShowTextDescriptions)

        #data plot bindings
        self.btnGenerateDataPlot.Bind(wx.EVT_BUTTON, self.ClickGenerateDataPlotButton)
        self.chkPlotR.Bind(wx.EVT_CHECKBOX, self.ChangePlotR)
        self.chkPlotV.Bind(wx.EVT_CHECKBOX, self.ChangePlotV)
        self.chkPlotThrust.Bind(wx.EVT_CHECKBOX, self.ChangePlotThrust)
        self.chkPlotIsp.Bind(wx.EVT_CHECKBOX, self.ChangePlotIsp)
        self.chkPlotMdot.Bind(wx.EVT_CHECKBOX, self.ChangePlotMdot)
        self.chkPlotEfficiency.Bind(wx.EVT_CHECKBOX, self.ChangePlotEfficiency)
        self.chkPlotThrottle.Bind(wx.EVT_CHECKBOX, self.ChangePlotThrottle)
        self.chkPlotPower.Bind(wx.EVT_CHECKBOX, self.ChangePlotPower)
        self.chkPlotGamma.Bind(wx.EVT_CHECKBOX, self.ChangePlotGamma)
        self.chkPlotDelta.Bind(wx.EVT_CHECKBOX, self.ChangePlotDelta)
        self.chkPlotCB_thrust_angle.Bind(wx.EVT_CHECKBOX, self.ChangePlotCB_thrust_angle)
        self.chkPlotMass.Bind(wx.EVT_CHECKBOX, self.ChangePlotMass)
        self.chkPlotNumberOfEngines.Bind(wx.EVT_CHECKBOX, self.ChangePlotNumberOfEngines)
        self.chkPlotActivePower.Bind(wx.EVT_CHECKBOX, self.ChangePlotActivePower)
        self.chkPlotWasteHeat.Bind(wx.EVT_CHECKBOX, self.ChangePlotWasteHeat)
        self.chkPlotCriticalEvents.Bind(wx.EVT_CHECKBOX, self.ChangePlotCriticalEvents)

        #format bindings
        self.spnctrlFontSizeControl.Bind(wx.EVT_SPINCTRL, self.ChangeFontSize)

        #ephemeris output bindings
        self.btnOutputSTKEphemeris.Bind(wx.EVT_BUTTON, self.ClickOutputSTKEphemeris)

    def ClickJourneyListBox(self, e):
        self.mission.ActiveJourney = self.JourneyListBox.GetSelection()

    def ClickPlotMissionButton(self, e):
        self.mission.PlotMission(self.plotoptions)

    def ChangeShowBoundaryOrbits(self, e):
        self.plotoptions.ShowBoundaryOrbits = self.chkShowBoundaryOrbits.GetValue()

    def ChangeShowPropagatedTrajectory(self, e):
        self.plotoptions.ShowPropagatedTrajectory = self.chkShowPropagatedTrajectory.GetValue()

    def ChangeShowThrustVectors(self, e):
        self.plotoptions.ShowThrustVectors = self.chkShowThrustVectors.GetValue()

    def ChangeShowTextDescriptions(self, e):
        self.plotoptions.ShowTextDescriptions = self.chkShowTextDescriptions.GetValue()

    def ClickGenerateDataPlotButton(self, e):
        #generate a custom plot if at least one checkbox is active
        if (self.plotoptions.PlotR or self.plotoptions.PlotV or self.plotoptions.PlotThrust or self.plotoptions.PlotIsp or self.plotoptions.PlotMdot
            or self.plotoptions.PlotEfficiency or self.plotoptions.PlotThrottle or self.plotoptions.PlotPower or self.plotoptions.PlotGamma or self.plotoptions.PlotDelta
            or self.chkPlotMass):
            self.mission.GenerateDataPlot(self.plotoptions)

    def ChangePlotR(self, e):
        self.plotoptions.PlotR = self.chkPlotR.GetValue()

    def ChangePlotV(self, e):
        self.plotoptions.PlotV = self.chkPlotV.GetValue()
    def ChangePlotThrust(self, e):
        self.plotoptions.PlotThrust = self.chkPlotThrust.GetValue()

    def ChangePlotIsp(self, e):
        self.plotoptions.PlotIsp = self.chkPlotIsp.GetValue()

    def ChangePlotMdot(self, e):
        self.plotoptions.PlotMdot = self.chkPlotMdot.GetValue()

    def ChangePlotEfficiency(self, e):
        self.plotoptions.PlotEfficiency = self.chkPlotEfficiency.GetValue()

    def ChangePlotThrottle(self, e):
        self.plotoptions.PlotThrottle = self.chkPlotThrottle.GetValue()

    def ChangePlotPower(self, e):
        self.plotoptions.PlotPower = self.chkPlotPower.GetValue()

    def ChangePlotGamma(self, e):
        self.plotoptions.PlotGamma = self.chkPlotGamma.GetValue()

    def ChangePlotDelta(self, e):
        self.plotoptions.PlotDelta = self.chkPlotDelta.GetValue()

    def ChangePlotCB_thrust_angle(self, e):
        self.plotoptions.PlotCB_thrust_angle = self.chkPlotCB_thrust_angle.GetValue()

    def ChangePlotMass(self, e):
        self.plotoptions.PlotMass = self.chkPlotMass.GetValue()

    def ChangePlotNumberOfEngines(self, e):
        self.plotoptions.PlotNumberOfEngines = self.chkPlotNumberOfEngines.GetValue()

    def ChangePlotActivePower(self, e):
        self.plotoptions.PlotActivePower = self.chkPlotActivePower.GetValue()

    def ChangePlotWasteHeat(self, e):
        self.plotoptions.PlotWasteHeat = self.chkPlotWasteHeat.GetValue()

    def ChangePlotCriticalEvents(self, e):
        self.plotoptions.PlotCriticalEvents = self.chkPlotCriticalEvents.GetValue()

    def ChangeFontSize(self, e):
        self.plotoptions.FontSize = self.spnctrlFontSizeControl.GetValue()

    def ClickOutputSTKEphemeris(self, e):
        self.mission.OutputSTKEphemeris(self)