import wx
import NSGAIIpopulation
import numpy
import copy

class NSGAIIPlotOptions:
    def __init__(self):
        self.UpperBounds = []
        self.LowerBounds = []
        self.TimeUnit = 0
        self.EpochUnit = 0

class NSGAIIpanel(wx.Panel):
    def __init__(self, parent, Population):
        wx.Panel.__init__(self, parent)

        self.NSGAIIpopulation = Population
        self.plotoptions = NSGAIIPlotOptions()
        self.Xobjective = 0
        self.Yobjective = 1
        self.Zobjective = 2
        self.Cobjective = 3

        self.plotoptions.LowerBounds = numpy.zeros(4)
        self.plotoptions.UpperBounds = 1.0e+50 * numpy.ones(4)

        #first we want an array of [Objective # label, combobox to select objective, lowerbound for objective display, upperbound for objective display]
        #we want one of these rows for every objective in the population file
        #each combobox after the second should have an option for "do not display"
        #put the array of objective selectors in a frame
        self.AxisOptionsBox = wx.StaticBox(self, -1, "Axis options")
        font = self.GetFont()
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        self.AxisOptionsBox.SetFont(font)
        AxisOptionsBoxSizer = wx.StaticBoxSizer(self.AxisOptionsBox, wx.VERTICAL)

        self.objective_selectors = []
        self.objective_upperbound_fields = []
        self.objective_lowerbound_fields = []
        self.objective_row_sizers = []
        
        #x axis
        xaxislabel = wx.StaticText(self, -1, "X axis", size=(100, -1))
        self.objective_selectors.append(wx.ComboBox(self, -1, choices = self.NSGAIIpopulation.objective_column_headers, style = wx.CB_READONLY, size=(200, -1)))
        self.objective_selectors[-1].SetSelection(self.Xobjective)
        self.objective_lowerbound_fields.append(wx.TextCtrl(self, -1, str(0.0)))
        self.objective_upperbound_fields.append(wx.TextCtrl(self, -1, str(1.0e+50)))
        self.objective_selectors[-1].Bind(wx.EVT_COMBOBOX, self.ChangeCObjective)
        self.objective_lowerbound_fields[-1].Bind(wx.EVT_KILL_FOCUS, self.ChangeXLowerBound)
        self.objective_upperbound_fields[-1].Bind(wx.EVT_KILL_FOCUS, self.ChangeXUpperBound)
        XaxisSizer = wx.BoxSizer(wx.HORIZONTAL)
        XaxisSizer.AddMany([xaxislabel, self.objective_selectors[-1], self.objective_lowerbound_fields[-1], self.objective_upperbound_fields[-1]])
        AxisOptionsBoxSizer.Add(XaxisSizer)

        #y axis
        yaxislabel = wx.StaticText(self, -1, "Y axis", size=(100, -1))
        self.objective_selectors.append(wx.ComboBox(self, -1, choices = self.NSGAIIpopulation.objective_column_headers, style = wx.CB_READONLY, size=(200, -1)))
        self.objective_selectors[-1].SetSelection(self.Yobjective)
        self.objective_lowerbound_fields.append(wx.TextCtrl(self, -1, str(0.0)))
        self.objective_upperbound_fields.append(wx.TextCtrl(self, -1, str(1.0e+50)))
        self.objective_selectors[-1].Bind(wx.EVT_COMBOBOX, self.ChangeCObjective)
        self.objective_lowerbound_fields[-1].Bind(wx.EVT_KILL_FOCUS, self.ChangeYLowerBound)
        self.objective_upperbound_fields[-1].Bind(wx.EVT_KILL_FOCUS, self.ChangeYUpperBound)
        YaxisSizer = wx.BoxSizer(wx.HORIZONTAL)
        YaxisSizer.AddMany([yaxislabel, self.objective_selectors[-1], self.objective_lowerbound_fields[-1], self.objective_upperbound_fields[-1]])
        AxisOptionsBoxSizer.Add(YaxisSizer)

        #z axis
        if len(self.NSGAIIpopulation.objective_column_headers) > 2:
            zaxislabel = wx.StaticText(self, -1, "Z axis", size=(100, -1))
            zaxischoices = copy.deepcopy(self.NSGAIIpopulation.objective_column_headers)
            zaxischoices.append('do not display')
            self.objective_selectors.append(wx.ComboBox(self, -1, choices = zaxischoices, style = wx.CB_READONLY, size=(200, -1)))
            self.objective_selectors[-1].SetSelection(self.Zobjective)
            self.objective_lowerbound_fields.append(wx.TextCtrl(self, -1, str(0.0)))
            self.objective_upperbound_fields.append(wx.TextCtrl(self, -1, str(1.0e+50)))
            self.objective_selectors[-1].Bind(wx.EVT_COMBOBOX, self.ChangeZObjective)
            self.objective_lowerbound_fields[-1].Bind(wx.EVT_KILL_FOCUS, self.ChangeZLowerBound)
            self.objective_upperbound_fields[-1].Bind(wx.EVT_KILL_FOCUS, self.ChangeZUpperBound)
            ZaxisSizer = wx.BoxSizer(wx.HORIZONTAL)
            ZaxisSizer.AddMany([zaxislabel, self.objective_selectors[-1], self.objective_lowerbound_fields[-1], self.objective_upperbound_fields[-1]])
            AxisOptionsBoxSizer.Add(ZaxisSizer)

        #color axis
        if len(self.NSGAIIpopulation.objective_column_headers) > 3:
            caxislabel = wx.StaticText(self, -1, "Color axis", size=(100, -1))
            caxischoices = copy.deepcopy(self.NSGAIIpopulation.objective_column_headers)
            caxischoices.append('do not display')
            self.objective_selectors.append(wx.ComboBox(self, -1, choices = caxischoices, style = wx.CB_READONLY, size=(200, -1)))
            self.objective_selectors[-1].SetSelection(self.Cobjective)
            self.objective_lowerbound_fields.append(wx.TextCtrl(self, -1, str(0.0)))
            self.objective_upperbound_fields.append(wx.TextCtrl(self, -1, str(1.0e+50)))
            self.objective_selectors[-1].Bind(wx.EVT_COMBOBOX, self.ChangeCObjective)
            self.objective_lowerbound_fields[-1].Bind(wx.EVT_KILL_FOCUS, self.ChangeCLowerBound)
            self.objective_upperbound_fields[-1].Bind(wx.EVT_KILL_FOCUS, self.ChangeCUpperBound)
            CaxisSizer = wx.BoxSizer(wx.HORIZONTAL)
            CaxisSizer.AddMany([caxislabel, self.objective_selectors[-1], self.objective_lowerbound_fields[-1], self.objective_upperbound_fields[-1]])
            AxisOptionsBoxSizer.Add(CaxisSizer)

        #next we want checkboxes for any other plot options
        self.lblTimeUnit = wx.StaticText(self, -1, "Display time unit", size=(200,-1))
        TimeUnitChoices = ['years','days']
        self.cmbTimeUnit = wx.ComboBox(self, -1, choices=TimeUnitChoices, style = wx.CB_READONLY)
        self.cmbTimeUnit.SetSelection(self.plotoptions.TimeUnit)
        self.cmbTimeUnit.Bind(wx.EVT_COMBOBOX, self.ChangeTimeUnit)
        TimeUnitSizer = wx.BoxSizer(wx.HORIZONTAL)
        TimeUnitSizer.AddMany([self.lblTimeUnit, self.cmbTimeUnit])

        self.lblEpochUnit = wx.StaticText(self, -1, "Display epoch unit", size=(200,-1))
        EpochUnitChoices = ['TDB Gregorian','TDB MJD']
        self.cmbEpochUnit = wx.ComboBox(self, -1, choices=EpochUnitChoices, style = wx.CB_READONLY)
        self.cmbEpochUnit.SetSelection(self.plotoptions.EpochUnit)
        self.cmbEpochUnit.Bind(wx.EVT_COMBOBOX, self.ChangeEpochUnit)
        EpochUnitSizer = wx.BoxSizer(wx.HORIZONTAL)
        EpochUnitSizer.AddMany([self.lblEpochUnit, self.cmbEpochUnit])


        self.PlotOptionsBox = wx.StaticBox(self, -1, "Plot Options", size = (300, 300))
        font = self.GetFont()
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        self.PlotOptionsBox.SetFont(font)
        PlotOptionsSizer = wx.StaticBoxSizer(self.PlotOptionsBox, wx.VERTICAL)
        PlotOptionsSizer.AddMany([TimeUnitSizer, EpochUnitSizer])
        
        #finally we want a button to make the plot
        self.btnPlotPopulation = wx.Button(self, -1, "Plot Population")
        self.btnPlotPopulation.Bind(wx.EVT_BUTTON, self.ClickPlotPopulation)

        #put everything in a big sizer

        mainsizer = wx.BoxSizer(wx.VERTICAL)
        mainsizer.AddMany([AxisOptionsBoxSizer, PlotOptionsSizer, self.btnPlotPopulation])

        self.SetSizer(mainsizer)

        
    #methods
    def ChangeXObjective(self, event):
        self.Xobjective = self.objective_selectors[0].GetSelection()

    def ChangeYObjective(self, event):
        self.Yobjective = self.objective_selectors[1].GetSelection()

    def ChangeZObjective(self, event):
        self.Zobjective = self.objective_selectors[2].GetSelection()

    def ChangeCObjective(self, event):
        self.Cobjective = self.objective_selectors[3].GetSelection()

    def ChangeXLowerBound(self, event):
        self.plotoptions.LowerBounds[0] = eval(self.objective_lowerbound_fields[0].GetValue)

    def ChangeXUpperBound(self, event):
        self.plotoptions.UpperBounds[0] = eval(self.objective_upperbound_fields[0].GetValue)

    def ChangeYLowerBound(self, event):
        self.plotoptions.LowerBounds[1] = eval(self.objective_lowerbound_fields[2].GetValue)

    def ChangeYUpperBound(self, event):
        self.plotoptions.UpperBounds[1] = eval(self.objective_upperbound_fields[2].GetValue)

    def ChangeZLowerBound(self, event):
        self.plotoptions.LowerBounds[2] = eval(self.objective_lowerbound_fields[2].GetValue)

    def ChangeZUpperBound(self, event):
        self.plotoptions.UpperBounds[2] = eval(self.objective_upperbound_fields[2].GetValue)

    def ChangeCLowerBound(self, event):
        self.plotoptions.LowerBounds[3] = eval(self.objective_lowerbound_fields[3].GetValue)

    def ChangeCUpperBound(self, event):
        self.plotoptions.UpperBounds[3] = eval(self.objective_upperbound_fields[3].GetValue)

    def ChangeTimeUnit(self, event):
        self.plotoptions.TimeUnit = self.cmbTimeUnit.GetSelection()

    def ChangeEpochUnit(self, event):
        self.plotoptions.EpochUnit = self.cmbEpochUnit.GetSelection()

    def ClickPlotPopulation(self, event):
        #first assemble the ordered list of objectives
        #note that if C is set but not Z, throw an error

        if self.Cobjective < len(self.NSGAIIpopulation.objective_column_headers) - 1 and self.Zobjective == len(self.NSGAIIpopulation.objective_column_headers):
            errordlg = wx.MessageDialog(self, "You cannot set the color axis without setting the Z axis first", "EMTG Error", wx.OK)
            errordlg.ShowModal()
            errordlg.Destroy()

        else:
            ordered_list_of_objectives = [self.Xobjective, self.Yobjective]

            if self.Zobjective < len(self.NSGAIIpopulation.objective_column_headers):
                ordered_list_of_objectives.append(self.Zobjective)

            if self.Cobjective < len(self.NSGAIIpopulation.objective_column_headers):
                ordered_list_of_objectives.append(self.Cobjective)

            #check for duplicate objectives. If present, throw an error
            s = set()
            if any(obj in s or s.add(obj) for obj in ordered_list_of_objectives):
                errordlg = wx.MessageDialog(self, "Objective axes must be unique", "EMTG Error", wx.OK)
                errordlg.ShowModal()
                errordlg.Destroy()
            
            else:
                self.NSGAIIpopulation.plot_population(ordered_list_of_objectives, LowerBounds = self.plotoptions.LowerBounds, UpperBounds = self.plotoptions.UpperBounds, TimeUnit = self.plotoptions.TimeUnit, EpochUnit = self.plotoptions.EpochUnit)
