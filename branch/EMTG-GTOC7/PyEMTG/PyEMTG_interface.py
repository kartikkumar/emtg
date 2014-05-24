import MissionOptions as MO
import JourneyOptions as JO
import Mission
import Universe
import BodyPicker
import Archive
import ArchivePanel
import NSGAIIpopulation
import NSGAIIpanel
import numpy as np
import OptionsNotebook
import UniverseNotebook
import MissionPanel
import webbrowser
import os
import subprocess
import platform
import wx
import wx.calendar
import copy

class PyEMTG_interface(wx.Frame):
    homedir = os.path.dirname(__file__)
    emtgpath = ''
    filename = ''
    dirname = ''
    mode = ''
    
    def __init__(self, *args, **kwargs):
        super(PyEMTG_interface, self).__init__(*args, **kwargs)

        self.read_PyEMTG_options_file()
        
        icon = wx.Icon("clemonaut.ico", wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)

        self.initialize_GUI()
        self.Maximize()

    def read_PyEMTG_options_file(self):
        input_file_name = os.path.join(self.homedir, "PyEMTG.options")
        if os.path.isfile(input_file_name):
            inputfile = open(input_file_name, "r")
        else:
            print "Unable to open", input_file_name, "EMTG Error"
            return

        for line in inputfile:
            linecell = line.split(' ')

            if linecell[0] == "EMTG_path":
                self.emtgpath = linecell[1]
            elif linecell[0] == "default_universe_path":
                self.default_universe_path = linecell[1]

        inputfile.close()
        
    def initialize_GUI(self):
        self.menubar = wx.MenuBar()
        self.fileMenu = wx.Menu()
        
        newMenu = wx.Menu()
        fnewmission = newMenu.Append(wx.ID_NEW, '&Mission\tCtrl+m')
        fnewuniverse = newMenu.Append(wx.ID_ANY, '&Universe\tCtrl+u')
        
        self.fileMenu.AppendMenu(wx.ID_ANY, '&New', newMenu)
        fopen = self.fileMenu.Append(wx.ID_OPEN, '&Open\tCtrl+o')
        fsave = self.fileMenu.Append(wx.ID_SAVE, '&Save\tCtrl+s')
        self.fileMenu.AppendSeparator()
        frun = self.fileMenu.Append(wx.ID_ANY, '&Run\tCtrl+r')
        self.fileMenu.AppendSeparator()
        fedit = self.fileMenu.Append(wx.ID_EDIT, 'Open file in &Editor\tCtrl+e')        
        self.fileMenu.AppendSeparator()
        fexit = self.fileMenu.Append(wx.ID_EXIT, 'E&xit\tCtrl+q')
        self.menubar.Append(self.fileMenu, '&File')
        self.SetMenuBar(self.menubar)
        
        self.Bind(wx.EVT_MENU, self.OnNewMission, fnewmission, id=wx.ID_NEW)
        self.Bind(wx.EVT_MENU, self.OnNewUniverse, fnewuniverse, id=wx.ID_NEW)

        self.Bind(wx.EVT_MENU, self.OnOpen, fopen, id=wx.ID_OPEN)
        self.Bind(wx.EVT_MENU, self.OnSave, fsave, id=wx.ID_SAVE)
        
        self.Bind(wx.EVT_MENU, self.OnRun, frun, id=wx.ID_ANY)

        self.Bind(wx.EVT_MENU, self.OnEdit, fedit, id=wx.ID_EDIT)
        
        self.Bind(wx.EVT_MENU, self.OnExit, fexit, id=wx.ID_EXIT)
        self.Bind(wx.EVT_CLOSE, self.OnExit)

        self.Bind(wx.EVT_SIZE, self.OnResize)
        self.Bind(wx.EVT_MAXIMIZE, self.OnResize)
        
        #Disable various menu options at program start
        self.fileMenu.Enable(wx.ID_SAVE, False)
        self.fileMenu.Enable(wx.ID_EDIT, False)

        #welcome message
        self.lblWelcome = wx.StaticText(self, -1, "Welcome to EMTG. Please select an option from the File menu to begin.", style = wx.ALIGN_CENTER)
        font = self.GetFont()
        font.SetPointSize(20)
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        self.lblWelcome.SetFont(font)
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(self.lblWelcome, 1, flag = wx.CENTER)
        self.SetSizer(sizer)

        #self.SetSize((800,600))
        self.SetTitle("EMTG Python Interface")
        
        self.Show()

    def OnResize(self, e):
        self.Layout()
        MySize = self.GetClientSize()
        if self.mode == "options":
            self.optionsnotebook.SetSize((MySize.x, MySize.y))
            self.optionsnotebook.Layout()

        elif self.mode == "mission":
            self.missionpanel.SetSize((MySize.x, MySize.y))
            self.missionpanel.Layout()

        elif self.mode == "universe":
            self.universenotebook.SetSize((MySize.x, MySize.y))
            self.universenotebook.Layout()

        elif self.mode == "archive":
            self.archivepanel.SetSize((MySize.x, MySize.y))
            self.archivepanel.Layout()

        elif self.mode == "NSGAII":
            self.NSGAIIpanel.SetSize((MySize.x, MySize.y))
            self.NSGAIIpanel.Layout()
        
    def OnNewMission(self, e):
        #If the GUI has not yet loaded anything then create a new mission. Otherwise as for permission first.
        if self.mode == "":
            self.dirname = self.homedir
            self.filename = "default.emtgopt"
            self.missionoptions = MO.MissionOptions(os.path.join(self.dirname, self.filename))
            
            if self.missionoptions.success == 1:
                self.mode = "options"
                self.fileMenu.Enable(wx.ID_SAVE, True)
                self.fileMenu.Enable(wx.ID_EDIT, True)
                self.lblWelcome.Show(False)
                self.InitializeMissionOptionsEditor()
        else:
            dlg = wx.MessageDialog(self,
                    "Do you really want to create a new mission? This will clear your current GUI settings.",
                    "Confirm New", wx.OK|wx.CANCEL|wx.ICON_QUESTION)
            result = dlg.ShowModal()
            dlg.Destroy()
            if result == wx.ID_OK:
                #destroy previous options memory block if you have one
                if self.mode == "options":
                    self.mode = ""
                    self.missionoptions = []
                    self.optionsnotebook.Destroy()

                elif self.mode == "mission":
                    self.mode = ""
                    self.mission = []
                    self.missionpanel.Destroy()

                elif self.mode == "universe":
                    self.mode = ""
                    self.universe = []
                    self.universenotebook.Destroy()

                elif self.mode == "archive":
                    self.mode = ""
                    self.archive = []
                    self.archivepanel.Destroy()

                elif self.mode == "NSGAII":
                    self.mode = ""
                    self.NSGAII = []
                    self.NSGAIIpanel.Destroy()

                #attempt to open a new options file
                self.dirname = self.homedir
                self.filename = "default.emtgopt"
                self.missionoptions = MO.MissionOptions(os.path.join(self.dirname, self.filename))
                
                if self.missionoptions.success == 1:
                    self.mode = "options"
                    self.fileMenu.Enable(wx.ID_SAVE, True)
                    self.fileMenu.Enable(wx.ID_EDIT, True)
                    self.lblWelcome.Show(False)
                    self.InitializeMissionOptionsEditor()
                    
    def OnNewUniverse(self, e):
        if self.mode == "":
            self.dirname = self.homedir
            self.filename = "default.emtg_universe"
            self.universe = Universe.Universe(os.path.join(self.dirname, self.filename))
            
            if self.universe.success == 1:
                self.mode = "universe"
                self.fileMenu.Enable(wx.ID_SAVE, True)
                self.fileMenu.Enable(wx.ID_EDIT, True)
                self.lblWelcome.Show(False)
                self.InitializeUniverseOptionsEditor()
        else:
            dlg = wx.MessageDialog(self,
                    "Do you really want to create a new universe? This will clear your current GUI settings.",
                    "Confirm New", wx.OK|wx.CANCEL|wx.ICON_QUESTION)
            result = dlg.ShowModal()
            dlg.Destroy()
            if result == wx.ID_OK:
                #destroy previous options memory block if you have one
                if self.mode == "options":
                    self.mode = ""
                    self.missionoptions = []
                    self.optionsnotebook.Destroy()

                elif self.mode == "mission":
                    self.mode = ""
                    self.mission = []
                    self.missionpanel.Destroy()

                elif self.mode == "universe":
                    self.mode = ""
                    self.universe = []
                    self.universenotebook.Destroy()

                elif self.mode == "archive":
                    self.mode = ""
                    self.archive = []
                    self.archivepanel.Destroy()

                elif self.mode == "NSGAII":
                    self.mode = ""
                    self.NSGAII = []
                    self.NSGAIIpanel.Destroy()

                #attempt to open a new universe file
                self.dirname = self.homedir
                self.filename = "default.emtg_universe"
                self.universe = Universe.Universe(os.path.join(self.dirname, self.filename))
                
                if self.universe.success == 1:
                    self.mode = "universe"
                    self.fileMenu.Enable(wx.ID_SAVE, True)
                    self.fileMenu.Enable(wx.ID_EDIT, True)
                    self.lblWelcome.Show(False)
                    self.InitializeUniverseOptionsEditor()
        
        
    def OnOpen(self, e):
        self.OpenFile(e)
        
    def OnSave(self, e):
        self.SetFocus()
        self.SaveFile(e)

    def OnRun(self, e):
        self.SetFocus()
        saved = self.SaveFile(e)

        

        if saved:
            if platform.system() == 'Windows':
                os.system('start ' + self.emtgpath.strip() + ' ' + os.path.join(self.dirname, self.filename))
            #DETACHED_PROCESS = 0x00000008

            #pid = subprocess.Popen([sys.executable, 'C://Users//Jacob//Documents//Projects//EMTG//bin//EMTG_v8 ' + os.path.join(self.dirname, self.filename)],
                       #creationflags=DETACHED_PROCESS).pid
        
    def OnEdit(self, e):
        webbrowser.open(os.path.join(self.dirname, self.filename))
        
    def OnExit(self, e):
        self.SetFocus()
        if self.mode == "options" or self.mode == "universe":
            dlg = wx.MessageDialog(self,
            "Save changes?",
            "Confirm Exit", wx.YES|wx.NO|wx.CANCEL|wx.ICON_QUESTION)
            result = dlg.ShowModal()
            dlg.Destroy()
            if result == wx.ID_NO:
                self.Destroy()
            elif result == wx.ID_YES:
                saved = self.SaveFile(e)
                if saved:
                    self.Destroy()
        
        else:
            dlg = wx.MessageDialog(self,
                "Do you really want to exit?",
                "Confirm Exit", wx.OK|wx.CANCEL|wx.ICON_QUESTION)
            result = dlg.ShowModal()
            dlg.Destroy()
            if result == wx.ID_OK:
                self.Destroy()
        
    
    def OpenFile(self, e):
        dlg = wx.FileDialog(self, "Open an EMTG file", self.dirname, "", "*.emtg*;*.NSGAII", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()

            fileparts = self.filename.split(".")
        

            #before we actually open the new file, we need to clear memory associated with whatever file we currently have open
            if self.mode == "options":
                self.missionoptions = []
                self.optionsnotebook.Destroy()
                
            elif self.mode == "mission":
                self.mission = []
                self.missionpanel.Destroy()

            elif self.mode == "universe":
                self.universe = []
                self.universenotebook.Destroy()

            elif self.mode == "archive":
                    self.archive = []
                    self.archivepanel.Destroy()

            elif self.mode == "NSGAII":
                    self.NSGAII = []
                    self.NSGAIIpanel.Destroy()

            self.mode = ""

            #next open the new file
            if fileparts[1] == "emtgopt":
                self.missionoptions = MO.MissionOptions(os.path.join(self.dirname, self.filename))
                if self.missionoptions.success == 1:
                    self.mode = "options"
                    self.lblWelcome.Show(False)
                    self.InitializeMissionOptionsEditor()
                    self.fileMenu.Enable(wx.ID_SAVE, True)
                    self.fileMenu.Enable(wx.ID_EDIT, True)

            elif fileparts[1] == "emtg":
                self.mission = Mission.Mission(os.path.join(self.dirname, self.filename))
                if self.mission.success == 1:
                    self.mode = "mission"
                    self.lblWelcome.Show(False)
                    self.missionpanel = MissionPanel.MissionPanel(self, self.mission)
                    self.missionpanel.SetSize(self.GetSize())
                    self.fileMenu.Enable(wx.ID_EDIT, True)
        
            elif fileparts[1] == "emtg_universe":
                self.universe = Universe.Universe(os.path.join(self.dirname, self.filename))
                if self.universe.success == 1:
                    self.mode = "universe"
                    self.lblWelcome.Show(False)
                    self.InitializeUniverseOptionsEditor()
                    self.fileMenu.Enable(wx.ID_SAVE, True)
                    self.fileMenu.Enable(wx.ID_EDIT, True)

            elif fileparts[1] == "emtgbatch":
                webbrowser.open(os.path.join(self.dirname, self.filename))

            elif fileparts[1] == "emtg_archive":
                self.archive = Archive.Archive(os.path.join(self.dirname, self.filename))
                if self.archive.success == 1:
                    self.mode = "archive"
                    self.lblWelcome.Show(False)
                    self.InitializeArchiveProcessor()
                    self.fileMenu.Enable(wx.ID_EDIT, True)
                    self.fileMenu.Enable(wx.ID_SAVE, False)

            elif fileparts[1] == "NSGAII":
                self.NSGAII = NSGAIIpopulation.NSGAII_outerloop_population(os.path.join(self.dirname, self.filename))
                if self.NSGAII.success == 1:
                    self.mode = "NSGAII"
                    self.lblWelcome.Show(False)
                    self.InitializeNSGAIIPanel()
                    self.fileMenu.Enable(wx.ID_EDIT, True)
                    self.fileMenu.Enable(wx.ID_SAVE, False)

            else:
                errordlg = wx.MessageDialog(self, "Unrecognized file type.", "EMTG Error", wx.OK)
                errordlg.ShowModal()
                errordlg.Destroy()

        dlg.Destroy()


    
    def SaveFile(self, e):
        extension = []
        if self.mode == "options":
            dialogtitle = self.missionoptions.mission_name
            extension = "*.emtgopt"
        elif self.mode == "universe":
            dialogtitle = self.universe.central_body_name
            extension = "*.emtg_universe"

        dlg = wx.FileDialog(self, "Save", self.dirname, dialogtitle, extension, wx.SAVE | wx.FD_OVERWRITE_PROMPT)
        if dlg.ShowModal() == wx.ID_OK:
            saved = True
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            
            if self.mode == "options":
                self.missionoptions.write_options_file(os.path.join(self.dirname, self.filename))
            if self.mode == "universe":
                self.universe.write_universe_file(os.path.join(self.dirname, self.filename))
        else:
            saved = False
        
        dlg.Destroy()

        return saved

    def InitializeUniverseOptionsEditor(self):
        #create and size a universe notebook object
        self.universenotebook = UniverseNotebook.UniverseNotebook(self, self.universe)
        self.universenotebook.SetSize(self.GetClientSize())

    def InitializeArchiveProcessor(self):
        #create and size an Archive panel object
        self.archivepanel = ArchivePanel.ArchivePanel(self, self.archive)
        self.archivepanel.SetSize(self.GetSize())

    def InitializeNSGAIIPanel(self):
        self.NSGAIIpanel = NSGAIIpanel.NSGAIIpanel(self, self.NSGAII)
        self.NSGAIIpanel.SetSize(self.GetSize())
        
    def InitializeMissionOptionsEditor(self):
        #create an options notebook object
        self.optionsnotebook = OptionsNotebook.OptionsBook(self)
        
        #update the latest information from the missionoptions object into the GUI
        self.missionoptions.update_all_panels(self.optionsnotebook)
        self.optionsnotebook.SetSize(self.GetClientSize())
        
        #bind the various GUI controls to events which modify the missionoptions object
        
        #global mission options
        self.optionsnotebook.tabGlobal.txtMissionName.Bind(wx.EVT_KILL_FOCUS, self.ChangeMissionName)
        self.optionsnotebook.tabGlobal.cmbMissionType.Bind(wx.EVT_COMBOBOX, self.ChangeMissionType)
        self.optionsnotebook.tabGlobal.cmbobjective_type.Bind(wx.EVT_COMBOBOX,self.Changeobjective_type)
        self.optionsnotebook.tabGlobal.chkinclude_initial_impulse_in_cost.Bind(wx.EVT_CHECKBOX,self.Changeinclude_initial_impulse_in_cost)
        self.optionsnotebook.tabGlobal.txtmax_phases_per_journey.Bind(wx.EVT_KILL_FOCUS,self.Changemax_phases_per_journey)
        self.optionsnotebook.tabGlobal.txtlaunch_window_open_date.Bind(wx.EVT_KILL_FOCUS,self.Changelaunch_window_open_date)
        self.optionsnotebook.tabGlobal.LaunchDateCalendar.Bind(wx.calendar.EVT_CALENDAR_SEL_CHANGED, self.ChangeLaunchDateCalendar)
        self.optionsnotebook.tabGlobal.txtnum_timesteps.Bind(wx.EVT_KILL_FOCUS,self.Changenum_timesteps)
        self.optionsnotebook.tabGlobal.cmbstep_size_distribution.Bind(wx.EVT_COMBOBOX,self.Changestep_size_distribution)
        self.optionsnotebook.tabGlobal.txtstep_size_stdv_or_scale.Bind(wx.EVT_KILL_FOCUS,self.Changestep_size_stdv_or_scale)
        self.optionsnotebook.tabGlobal.cmbcontrol_coordinate_system.Bind(wx.EVT_COMBOBOX, self.Changecontrol_coordinate_system)
        self.optionsnotebook.tabGlobal.txtDLA_bounds_lower.Bind(wx.EVT_KILL_FOCUS,self.ChangeDLA_bounds_lower)
        self.optionsnotebook.tabGlobal.txtDLA_bounds_upper.Bind(wx.EVT_KILL_FOCUS,self.ChangeDLA_bounds_upper)
        self.optionsnotebook.tabGlobal.chkglobal_timebounded.Bind(wx.EVT_CHECKBOX,self.Changeglobal_timebounded)
        self.optionsnotebook.tabGlobal.txttotal_flight_time_bounds_lower.Bind(wx.EVT_KILL_FOCUS,self.Changetotal_flight_time_bounds_lower)
        self.optionsnotebook.tabGlobal.txttotal_flight_time_bounds_upper.Bind(wx.EVT_KILL_FOCUS,self.Changetotal_flight_time_bounds_upper)
        self.optionsnotebook.tabGlobal.txtforced_post_launch_coast.Bind(wx.EVT_KILL_FOCUS,self.Changeforced_post_launch_coast)
        self.optionsnotebook.tabGlobal.txtforced_flyby_coast.Bind(wx.EVT_KILL_FOCUS,self.Changeforced_flyby_coast)
        self.optionsnotebook.tabGlobal.txtinitial_V_infinity_x.Bind(wx.EVT_KILL_FOCUS,self.Changeinitial_V_infinity_x)
        self.optionsnotebook.tabGlobal.txtinitial_V_infinity_y.Bind(wx.EVT_KILL_FOCUS,self.Changeinitial_V_infinity_y)
        self.optionsnotebook.tabGlobal.txtinitial_V_infinity_z.Bind(wx.EVT_KILL_FOCUS,self.Changeinitial_V_infinity_z)
        self.optionsnotebook.tabGlobal.txtminimum_dry_mass.Bind(wx.EVT_KILL_FOCUS,self.Changeminimum_dry_mass)
        self.optionsnotebook.tabGlobal.txtpost_mission_delta_v.Bind(wx.EVT_KILL_FOCUS,self.Changepost_mission_delta_v)

        
        #spacecraft options
        self.optionsnotebook.tabSpacecraft.txtmaximum_mass.Bind(wx.EVT_KILL_FOCUS, self.Changemaximum_mass)
        self.optionsnotebook.tabSpacecraft.chkallow_initial_mass_to_vary.Bind(wx.EVT_CHECKBOX, self.Changeallow_initial_mass_to_vary)
        self.optionsnotebook.tabSpacecraft.txtEP_dry_mass.Bind(wx.EVT_KILL_FOCUS, self.ChangeEP_dry_mass)
        self.optionsnotebook.tabSpacecraft.cmbLV_type.Bind(wx.EVT_COMBOBOX, self.ChangeLV_type)
        self.optionsnotebook.tabSpacecraft.txtIspDS.Bind(wx.EVT_KILL_FOCUS, self.ChangeIspDS)
        self.optionsnotebook.tabSpacecraft.txtIspChem.Bind(wx.EVT_KILL_FOCUS, self.ChangeIspChem)
        self.optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients0.Bind(wx.EVT_KILL_FOCUS, self.Changecustom_LV_coefficients0)
        self.optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients1.Bind(wx.EVT_KILL_FOCUS, self.Changecustom_LV_coefficients1)
        self.optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients2.Bind(wx.EVT_KILL_FOCUS, self.Changecustom_LV_coefficients2)
        self.optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients3.Bind(wx.EVT_KILL_FOCUS, self.Changecustom_LV_coefficients3)
        self.optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients4.Bind(wx.EVT_KILL_FOCUS, self.Changecustom_LV_coefficients4)
        self.optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients5.Bind(wx.EVT_KILL_FOCUS, self.Changecustom_LV_coefficients5)
        self.optionsnotebook.tabSpacecraft.txtcustom_LV_C3_bounds_lower.Bind(wx.EVT_KILL_FOCUS, self.Changecustom_LV_C3_bounds_lower)
        self.optionsnotebook.tabSpacecraft.txtcustom_LV_C3_bounds_upper.Bind(wx.EVT_KILL_FOCUS, self.Changecustom_LV_C3_bounds_upper)
        self.optionsnotebook.tabSpacecraft.txtparking_orbit_altitude.Bind(wx.EVT_KILL_FOCUS, self.Changeparking_orbit_altitude)
        self.optionsnotebook.tabSpacecraft.txtparking_orbit_inclination.Bind(wx.EVT_KILL_FOCUS, self.Changeparking_orbit_inclination)
        self.optionsnotebook.tabSpacecraft.txtpost_mission_Isp.Bind(wx.EVT_KILL_FOCUS, self.Changepost_mission_Isp)
        self.optionsnotebook.tabSpacecraft.chkenable_propellant_mass_constraint.Bind(wx.EVT_CHECKBOX, self.Changeenable_maximum_propellant_mass_constraint)
        self.optionsnotebook.tabSpacecraft.txtmaximum_propellant_mass.Bind(wx.EVT_KILL_FOCUS, self.Changemaximum_propellant_mass)
        self.optionsnotebook.tabSpacecraft.txtpropellant_margin.Bind(wx.EVT_KILL_FOCUS, self.Changepropellant_margin)
        self.optionsnotebook.tabSpacecraft.txtpower_margin.Bind(wx.EVT_KILL_FOCUS, self.Changepower_margin)
        self.optionsnotebook.tabSpacecraft.txtLV_margin.Bind(wx.EVT_KILL_FOCUS, self.ChangeLV_margin)
        self.optionsnotebook.tabSpacecraft.txtLV_adapter_mass.Bind(wx.EVT_KILL_FOCUS, self.Change_LV_adapter_mass)
        self.optionsnotebook.tabSpacecraft.cmbengine_type.Bind(wx.EVT_COMBOBOX, self.Changeengine_type)
        self.optionsnotebook.tabSpacecraft.txtnumber_of_engines.Bind(wx.EVT_KILL_FOCUS, self.Changenumber_of_engines)
        self.optionsnotebook.tabSpacecraft.cmbthrottle_logic_mode.Bind(wx.EVT_COMBOBOX, self.Changethrottle_logic_mode)
        self.optionsnotebook.tabSpacecraft.txtthrottle_sharpness.Bind(wx.EVT_KILL_FOCUS, self.Changethrottle_sharpness)
        self.optionsnotebook.tabSpacecraft.txtengine_duty_cycle.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_duty_cycle)
        self.optionsnotebook.tabSpacecraft.txtThrust.Bind(wx.EVT_KILL_FOCUS, self.ChangeThrust)
        self.optionsnotebook.tabSpacecraft.txtIspLT.Bind(wx.EVT_KILL_FOCUS, self.ChangeIspLT)
        self.optionsnotebook.tabSpacecraft.txtIspLT_minimum.Bind(wx.EVT_KILL_FOCUS, self.ChangeIspLT_minimum)
        self.optionsnotebook.tabSpacecraft.txtuser_defined_engine_efficiency.Bind(wx.EVT_KILL_FOCUS, self.Changeuser_defined_engine_efficiency)
        self.optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients0.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_thrust_coefficients0)
        self.optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients1.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_thrust_coefficients1)
        self.optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients2.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_thrust_coefficients2)
        self.optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients3.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_thrust_coefficients3)
        self.optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients4.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_thrust_coefficients4)
        self.optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients5.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_thrust_coefficients5)
        self.optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients6.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_thrust_coefficients6)
        self.optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients0.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_mass_flow_rate_coefficients0)
        self.optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients1.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_mass_flow_rate_coefficients1)
        self.optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients2.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_mass_flow_rate_coefficients2)
        self.optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients3.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_mass_flow_rate_coefficients3)
        self.optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients4.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_mass_flow_rate_coefficients4)
        self.optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients5.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_mass_flow_rate_coefficients5)
        self.optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients6.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_mass_flow_rate_coefficients6)
        self.optionsnotebook.tabSpacecraft.txtengine_input_power_bounds_lower.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_power_bounds_lower)
        self.optionsnotebook.tabSpacecraft.txtengine_input_power_bounds_upper.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_power_bounds_upper)
        self.optionsnotebook.tabSpacecraft.txtpower_at_1_AU.Bind(wx.EVT_KILL_FOCUS, self.Changepower_at_1_AU)
        self.optionsnotebook.tabSpacecraft.cmbpower_source_type.Bind(wx.EVT_COMBOBOX, self.Changepower_source_type)
        self.optionsnotebook.tabSpacecraft.txtsolar_power_gamma0.Bind(wx.EVT_KILL_FOCUS, self.Changesolar_power_gamma0)
        self.optionsnotebook.tabSpacecraft.txtsolar_power_gamma1.Bind(wx.EVT_KILL_FOCUS, self.Changesolar_power_gamma1)
        self.optionsnotebook.tabSpacecraft.txtsolar_power_gamma2.Bind(wx.EVT_KILL_FOCUS, self.Changesolar_power_gamma2)
        self.optionsnotebook.tabSpacecraft.txtsolar_power_gamma3.Bind(wx.EVT_KILL_FOCUS, self.Changesolar_power_gamma3)
        self.optionsnotebook.tabSpacecraft.txtsolar_power_gamma4.Bind(wx.EVT_KILL_FOCUS, self.Changesolar_power_gamma4)
        self.optionsnotebook.tabSpacecraft.cmbspacecraft_power_model_type.Bind(wx.EVT_COMBOBOX, self.Changespacecraft_power_model_type)
        self.optionsnotebook.tabSpacecraft.txtspacecraft_power_coefficients0.Bind(wx.EVT_KILL_FOCUS, self.Changespacecraft_power_coefficients0)
        self.optionsnotebook.tabSpacecraft.txtspacecraft_power_coefficients1.Bind(wx.EVT_KILL_FOCUS, self.Changespacecraft_power_coefficients1)
        self.optionsnotebook.tabSpacecraft.txtspacecraft_power_coefficients2.Bind(wx.EVT_KILL_FOCUS, self.Changespacecraft_power_coefficients2)
        self.optionsnotebook.tabSpacecraft.txtpower_decay_rate.Bind(wx.EVT_KILL_FOCUS, self.Changepower_decay_rate)

        
        #journey options
        self.optionsnotebook.tabJourney.JourneySelectBox.Bind(wx.EVT_LISTBOX, self.ChangeJourneySelectBoxChoice)
        self.optionsnotebook.tabJourney.btnAddNewJourney.Bind(wx.EVT_BUTTON, self.ClickAddNewJourney)
        self.optionsnotebook.tabJourney.btnDeleteJourney.Bind(wx.EVT_BUTTON, self.ClickDeleteJourney)
        self.optionsnotebook.tabJourney.btnMoveJourneyUp.Bind(wx.EVT_BUTTON, self.ClickMoveJourneyUp)
        self.optionsnotebook.tabJourney.btnMoveJourneyDown.Bind(wx.EVT_BUTTON, self.ClickMoveJourneyDown)

        self.optionsnotebook.tabJourney.txtjourney_names.Bind(wx.EVT_KILL_FOCUS,self.Changejourney_names)
        self.optionsnotebook.tabJourney.txtjourney_central_body.Bind(wx.EVT_KILL_FOCUS,self.Changejourney_central_body)
        self.optionsnotebook.tabJourney.txtdestination_list.Bind(wx.EVT_KILL_FOCUS,self.Changedestination_list)
        self.optionsnotebook.tabJourney.txtjourney_starting_mass_increment.Bind(wx.EVT_KILL_FOCUS,self.Changejourney_starting_mass_increment)
        self.optionsnotebook.tabJourney.chkjourney_variable_mass_increment.Bind(wx.EVT_CHECKBOX,self.Changejourney_variable_mass_increment)
        self.optionsnotebook.tabJourney.txtjourney_wait_time_bounds_lower.Bind(wx.EVT_KILL_FOCUS,self.Changejourney_wait_time_bounds_lower)
        self.optionsnotebook.tabJourney.txtjourney_wait_time_bounds_upper.Bind(wx.EVT_KILL_FOCUS,self.Changejourney_wait_time_bounds_upper)
        self.optionsnotebook.tabJourney.cmbjourney_timebounded.Bind(wx.EVT_COMBOBOX,self.Changejourney_timebounded)
        self.optionsnotebook.tabJourney.txtjourney_flight_time_bounds_lower.Bind(wx.EVT_KILL_FOCUS,self.Changejourney_flight_time_bounds_lower)
        self.optionsnotebook.tabJourney.txtjourney_flight_time_bounds_upper.Bind(wx.EVT_KILL_FOCUS,self.Changejourney_flight_time_bounds_upper)
        self.optionsnotebook.tabJourney.txtjourney_arrival_date_bounds_lower.Bind(wx.EVT_KILL_FOCUS,self.Changejourney_arrival_date_bounds_lower)
        self.optionsnotebook.tabJourney.ArrivalDateLowerCalendar.Bind(wx.calendar.EVT_CALENDAR_SEL_CHANGED, self.ChangeArrivalDateLowerCalendar)
        self.optionsnotebook.tabJourney.ArrivalDateUpperCalendar.Bind(wx.calendar.EVT_CALENDAR_SEL_CHANGED, self.ChangeArrivalDateUpperCalendar)
        self.optionsnotebook.tabJourney.txtjourney_arrival_date_bounds_upper.Bind(wx.EVT_KILL_FOCUS,self.Changejourney_arrival_date_bounds_upper)
        self.optionsnotebook.tabJourney.txtjourney_initial_impulse_bounds_lower.Bind(wx.EVT_KILL_FOCUS,self.Changejourney_initial_impulse_bounds_lower)
        self.optionsnotebook.tabJourney.txtjourney_initial_impulse_bounds_upper.Bind(wx.EVT_KILL_FOCUS,self.Changejourney_initial_impulse_bounds_upper)
        self.optionsnotebook.tabJourney.cmbjourney_departure_type.Bind(wx.EVT_COMBOBOX,self.Changejourney_departure_type)
        self.optionsnotebook.tabJourney.txtjourney_escape_spiral_starting_radius.Bind(wx.EVT_KILL_FOCUS, self.Changejourney_escape_spiral_radius)
        self.optionsnotebook.tabJourney.txtjourney_initial_velocity0.Bind(wx.EVT_KILL_FOCUS,self.Changejourney_initial_velocity0)
        self.optionsnotebook.tabJourney.txtjourney_initial_velocity1.Bind(wx.EVT_KILL_FOCUS,self.Changejourney_initial_velocity1)
        self.optionsnotebook.tabJourney.txtjourney_initial_velocity2.Bind(wx.EVT_KILL_FOCUS,self.Changejourney_initial_velocity2)
        self.optionsnotebook.tabJourney.cmbjourney_arrival_type.Bind(wx.EVT_COMBOBOX,self.Changejourney_arrival_type)
        self.optionsnotebook.tabJourney.txtjourney_capture_spiral_final_radius.Bind(wx.EVT_KILL_FOCUS, self.Changejourney_capture_spiral_radius)
        self.optionsnotebook.tabJourney.txtjourney_final_velocity0.Bind(wx.EVT_KILL_FOCUS,self.Changejourney_final_velocity0)
        self.optionsnotebook.tabJourney.txtjourney_final_velocity1.Bind(wx.EVT_KILL_FOCUS,self.Changejourney_final_velocity1)
        self.optionsnotebook.tabJourney.txtjourney_final_velocity2.Bind(wx.EVT_KILL_FOCUS,self.Changejourney_final_velocity2)
        self.optionsnotebook.tabJourney.chkjourney_arrival_declination_constraint_flag.Bind(wx.EVT_CHECKBOX, self.Changejourney_arrival_declination_constraint_flag)
        self.optionsnotebook.tabJourney.txtjourney_arrival_declination_bounds_lower.Bind(wx.EVT_KILL_FOCUS, self.Changejourney_arrival_declination_bounds_lower)
        self.optionsnotebook.tabJourney.txtjourney_arrival_declination_bounds_upper.Bind(wx.EVT_KILL_FOCUS, self.Changejourney_arrival_declination_bounds_upper)
        self.optionsnotebook.tabJourney.txtsequence.Bind(wx.EVT_KILL_FOCUS,self.Changesequence)
        self.optionsnotebook.tabJourney.txtjourney_perturbation_bodies.Bind(wx.EVT_KILL_FOCUS,self.Changejourney_perturbation_bodies)
        self.optionsnotebook.tabJourney.cmbjourney_departure_elements_type.Bind(wx.EVT_COMBOBOX,self.Changejourney_departure_elements_type)
        self.optionsnotebook.tabJourney.chkSMA_departure.Bind(wx.EVT_CHECKBOX,self.ChangevarySMA_departure)
        self.optionsnotebook.tabJourney.chkECC_departure.Bind(wx.EVT_CHECKBOX,self.ChangevaryECC_departure)
        self.optionsnotebook.tabJourney.chkINC_departure.Bind(wx.EVT_CHECKBOX,self.ChangevaryINC_departure)
        self.optionsnotebook.tabJourney.chkRAAN_departure.Bind(wx.EVT_CHECKBOX,self.ChangevaryRAAN_departure)
        self.optionsnotebook.tabJourney.chkAOP_departure.Bind(wx.EVT_CHECKBOX,self.ChangevaryAOP_departure)
        self.optionsnotebook.tabJourney.chkMA_departure.Bind(wx.EVT_CHECKBOX,self.ChangevaryMA_departure)
        self.optionsnotebook.tabJourney.txtSMA_departure.Bind(wx.EVT_KILL_FOCUS,self.ChangeSMA_departure)
        self.optionsnotebook.tabJourney.txtECC_departure.Bind(wx.EVT_KILL_FOCUS,self.ChangeECC_departure)
        self.optionsnotebook.tabJourney.txtINC_departure.Bind(wx.EVT_KILL_FOCUS,self.ChangeINC_departure)
        self.optionsnotebook.tabJourney.txtRAAN_departure.Bind(wx.EVT_KILL_FOCUS,self.ChangeRAAN_departure)
        self.optionsnotebook.tabJourney.txtAOP_departure.Bind(wx.EVT_KILL_FOCUS,self.ChangeAOP_departure)
        self.optionsnotebook.tabJourney.txtMA_departure.Bind(wx.EVT_KILL_FOCUS,self.ChangeMA_departure)
        self.optionsnotebook.tabJourney.txtSMA_departure0.Bind(wx.EVT_KILL_FOCUS,self.ChangeSMA_departure0)
        self.optionsnotebook.tabJourney.txtECC_departure0.Bind(wx.EVT_KILL_FOCUS,self.ChangeECC_departure0)
        self.optionsnotebook.tabJourney.txtINC_departure0.Bind(wx.EVT_KILL_FOCUS,self.ChangeINC_departure0)
        self.optionsnotebook.tabJourney.txtRAAN_departure0.Bind(wx.EVT_KILL_FOCUS,self.ChangeRAAN_departure0)
        self.optionsnotebook.tabJourney.txtAOP_departure0.Bind(wx.EVT_KILL_FOCUS,self.ChangeAOP_departure0)
        self.optionsnotebook.tabJourney.txtMA_departure0.Bind(wx.EVT_KILL_FOCUS,self.ChangeMA_departure0)
        self.optionsnotebook.tabJourney.txtSMA_departure1.Bind(wx.EVT_KILL_FOCUS,self.ChangeSMA_departure1)
        self.optionsnotebook.tabJourney.txtECC_departure1.Bind(wx.EVT_KILL_FOCUS,self.ChangeECC_departure1)
        self.optionsnotebook.tabJourney.txtINC_departure1.Bind(wx.EVT_KILL_FOCUS,self.ChangeINC_departure1)
        self.optionsnotebook.tabJourney.txtRAAN_departure1.Bind(wx.EVT_KILL_FOCUS,self.ChangeRAAN_departure1)
        self.optionsnotebook.tabJourney.txtAOP_departure1.Bind(wx.EVT_KILL_FOCUS,self.ChangeAOP_departure1)
        self.optionsnotebook.tabJourney.txtMA_departure1.Bind(wx.EVT_KILL_FOCUS,self.ChangeMA_departure1)
        self.optionsnotebook.tabJourney.cmbjourney_arrival_elements_type.Bind(wx.EVT_COMBOBOX,self.Changejourney_arrival_elements_type)
        self.optionsnotebook.tabJourney.chkSMA_arrival.Bind(wx.EVT_CHECKBOX,self.ChangevarySMA_arrival)
        self.optionsnotebook.tabJourney.chkECC_arrival.Bind(wx.EVT_CHECKBOX,self.ChangevaryECC_arrival)
        self.optionsnotebook.tabJourney.chkINC_arrival.Bind(wx.EVT_CHECKBOX,self.ChangevaryINC_arrival)
        self.optionsnotebook.tabJourney.chkRAAN_arrival.Bind(wx.EVT_CHECKBOX,self.ChangevaryRAAN_arrival)
        self.optionsnotebook.tabJourney.chkAOP_arrival.Bind(wx.EVT_CHECKBOX,self.ChangevaryAOP_arrival)
        self.optionsnotebook.tabJourney.chkMA_arrival.Bind(wx.EVT_CHECKBOX,self.ChangevaryMA_arrival)
        self.optionsnotebook.tabJourney.txtSMA_arrival.Bind(wx.EVT_KILL_FOCUS,self.ChangeSMA_arrival)
        self.optionsnotebook.tabJourney.txtECC_arrival.Bind(wx.EVT_KILL_FOCUS,self.ChangeECC_arrival)
        self.optionsnotebook.tabJourney.txtINC_arrival.Bind(wx.EVT_KILL_FOCUS,self.ChangeINC_arrival)
        self.optionsnotebook.tabJourney.txtRAAN_arrival.Bind(wx.EVT_KILL_FOCUS,self.ChangeRAAN_arrival)
        self.optionsnotebook.tabJourney.txtAOP_arrival.Bind(wx.EVT_KILL_FOCUS,self.ChangeAOP_arrival)
        self.optionsnotebook.tabJourney.txtMA_arrival.Bind(wx.EVT_KILL_FOCUS,self.ChangeMA_arrival)
        self.optionsnotebook.tabJourney.txtSMA_arrival0.Bind(wx.EVT_KILL_FOCUS,self.ChangeSMA_arrival0)
        self.optionsnotebook.tabJourney.txtECC_arrival0.Bind(wx.EVT_KILL_FOCUS,self.ChangeECC_arrival0)
        self.optionsnotebook.tabJourney.txtINC_arrival0.Bind(wx.EVT_KILL_FOCUS,self.ChangeINC_arrival0)
        self.optionsnotebook.tabJourney.txtRAAN_arrival0.Bind(wx.EVT_KILL_FOCUS,self.ChangeRAAN_arrival0)
        self.optionsnotebook.tabJourney.txtAOP_arrival0.Bind(wx.EVT_KILL_FOCUS,self.ChangeAOP_arrival0)
        self.optionsnotebook.tabJourney.txtMA_arrival0.Bind(wx.EVT_KILL_FOCUS,self.ChangeMA_arrival0)
        self.optionsnotebook.tabJourney.txtSMA_arrival1.Bind(wx.EVT_KILL_FOCUS,self.ChangeSMA_arrival1)
        self.optionsnotebook.tabJourney.txtECC_arrival1.Bind(wx.EVT_KILL_FOCUS,self.ChangeECC_arrival1)
        self.optionsnotebook.tabJourney.txtINC_arrival1.Bind(wx.EVT_KILL_FOCUS,self.ChangeINC_arrival1)
        self.optionsnotebook.tabJourney.txtRAAN_arrival1.Bind(wx.EVT_KILL_FOCUS,self.ChangeRAAN_arrival1)
        self.optionsnotebook.tabJourney.txtAOP_arrival1.Bind(wx.EVT_KILL_FOCUS,self.ChangeAOP_arrival1)
        self.optionsnotebook.tabJourney.txtMA_arrival1.Bind(wx.EVT_KILL_FOCUS,self.ChangeMA_arrival1)
        self.optionsnotebook.tabJourney.btndestination_list.Bind(wx.EVT_BUTTON,self.Clickdestination_list)
        self.optionsnotebook.tabJourney.btnjourney_central_body.Bind(wx.EVT_BUTTON,self.Clickjourney_central_body)
        self.optionsnotebook.tabJourney.btnsequence.Bind(wx.EVT_BUTTON,self.Clicksequence)
        self.optionsnotebook.tabJourney.btnjourney_perturbation_bodies.Bind(wx.EVT_BUTTON,self.Clickjourney_perturbation_bodies)

        
        #solver options
        self.optionsnotebook.tabSolver.cmbInnerLoopSolver.Bind(wx.EVT_COMBOBOX, self.ChangeInnerLoopSolver)
        self.optionsnotebook.tabSolver.cmbNLP_solver_type.Bind(wx.EVT_COMBOBOX, self.ChangeNLP_solver_type)
        self.optionsnotebook.tabSolver.cmbNLP_solver_mode.Bind(wx.EVT_COMBOBOX, self.ChangeNLP_solver_mode)
        self.optionsnotebook.tabSolver.chkquiet_NLP.Bind(wx.EVT_CHECKBOX, self.Changequiet_NLP)
        self.optionsnotebook.tabSolver.chkquiet_MBH.Bind(wx.EVT_CHECKBOX, self.Changequiet_MBH)
        self.optionsnotebook.tabSolver.chkMBH_two_step.Bind(wx.EVT_CHECKBOX, self.ChangeMBH_two_step)
        self.optionsnotebook.tabSolver.txtFD_stepsize.Bind(wx.EVT_KILL_FOCUS, self.ChangeFD_stepsize)
        self.optionsnotebook.tabSolver.txtFD_stepsize_coarse.Bind(wx.EVT_KILL_FOCUS, self.ChangeFD_stepsize_coarse)
        self.optionsnotebook.tabSolver.chkACE_feasible_point_finder.Bind(wx.EVT_CHECKBOX, self.ChangeACE_feasible_point_finder)
        self.optionsnotebook.tabSolver.txtMBH_max_not_improve.Bind(wx.EVT_KILL_FOCUS, self.ChangeMBH_max_not_improve)
        self.optionsnotebook.tabSolver.txtMBH_max_trials.Bind(wx.EVT_KILL_FOCUS, self.ChangeMBH_max_trials)
        self.optionsnotebook.tabSolver.txtMBH_max_run_time.Bind(wx.EVT_KILL_FOCUS, self.ChangeMBH_max_run_time)
        self.optionsnotebook.tabSolver.txtMBH_max_step_size.Bind(wx.EVT_KILL_FOCUS, self.ChangeMBH_max_step_size)
        self.optionsnotebook.tabSolver.cmbMBH_hop_distribution.Bind(wx.EVT_COMBOBOX, self.ChangeMBH_hop_distribution)
        self.optionsnotebook.tabSolver.txtMBH_Pareto_alpha.Bind(wx.EVT_KILL_FOCUS, self.ChangeMBH_Pareto_alpha)
        self.optionsnotebook.tabSolver.txtMBH_time_hop_probability.Bind(wx.EVT_KILL_FOCUS, self.ChangeMBH_time_hop_probability)
        self.optionsnotebook.tabSolver.txtsnopt_feasibility_tolerance.Bind(wx.EVT_KILL_FOCUS, self.Changesnopt_feasibility_tolerance)
        self.optionsnotebook.tabSolver.txtsnopt_major_iterations.Bind(wx.EVT_KILL_FOCUS, self.Changesnopt_major_iterations)
        self.optionsnotebook.tabSolver.txtsnopt_max_run_time.Bind(wx.EVT_KILL_FOCUS, self.Changesnopt_max_run_time)
        self.optionsnotebook.tabSolver.cmbderivative_type.Bind(wx.EVT_COMBOBOX, self.Changederivative_type)
        self.optionsnotebook.tabSolver.chkcheck_derivatives.Bind(wx.EVT_CHECKBOX, self.ChangeCheckDerivatives)
        self.optionsnotebook.tabSolver.chkseed_MBH.Bind(wx.EVT_CHECKBOX, self.ChangeSeedMBH)
        self.optionsnotebook.tabSolver.cmbinitial_guess_control_coordinate_system.Bind(wx.EVT_COMBOBOX, self.Changeinitial_guess_control_coordinate_system)
        self.optionsnotebook.tabSolver.chkinterpolate_initial_guess.Bind(wx.EVT_CHECKBOX, self.ChangeInterpolateInitialGuess)
        self.optionsnotebook.tabSolver.txtinitial_guess_num_timesteps.Bind(wx.EVT_KILL_FOCUS, self.Changeinitial_guess_num_timesteps)
        self.optionsnotebook.tabSolver.cmbinitial_guess_step_size_distribution.Bind(wx.EVT_COMBOBOX, self.Changeinitial_guess_step_size_distribution)
        self.optionsnotebook.tabSolver.txtinitial_guess_step_size_stdv_or_scale.Bind(wx.EVT_KILL_FOCUS, self.Changeinitial_guess_step_size_stdv_or_scale)
        self.optionsnotebook.tabSolver.cmbMBH_zero_control_initial_guess.Bind(wx.EVT_COMBOBOX, self.ChangeMBH_zero_control_initial_guess)
        
        self.optionsnotebook.tabSolver.cmbrun_outerloop.Bind(wx.EVT_COMBOBOX,self.Changerun_outerloop)
        self.optionsnotebook.tabSolver.txtouterloop_popsize.Bind(wx.EVT_KILL_FOCUS,self.Changeouterloop_popsize)
        self.optionsnotebook.tabSolver.txtouterloop_genmax.Bind(wx.EVT_KILL_FOCUS,self.Changeouterloop_genmax)
        self.optionsnotebook.tabSolver.txtouterloop_tournamentsize.Bind(wx.EVT_KILL_FOCUS,self.Changeouterloop_tournamentsize)
        self.optionsnotebook.tabSolver.txtouterloop_CR.Bind(wx.EVT_KILL_FOCUS,self.Changeouterloop_CR)
        self.optionsnotebook.tabSolver.txtouterloop_mu.Bind(wx.EVT_KILL_FOCUS,self.Changeouterloop_mu)
        self.optionsnotebook.tabSolver.txtouterloop_stallmax.Bind(wx.EVT_KILL_FOCUS,self.Changeouterloop_stallmax)
        self.optionsnotebook.tabSolver.txtouterloop_tolfit.Bind(wx.EVT_KILL_FOCUS,self.Changeouterloop_tolfit)
        self.optionsnotebook.tabSolver.txtouterloop_ntrials.Bind(wx.EVT_KILL_FOCUS,self.Changeouterloop_ntrials)
        self.optionsnotebook.tabSolver.txtouterloop_elitecount.Bind(wx.EVT_KILL_FOCUS,self.Changeouterloop_elitecount)
        self.optionsnotebook.tabSolver.txtouterloop_warmstart.Bind(wx.EVT_KILL_FOCUS,self.Changeouterloop_warmstart)

        self.optionsnotebook.tabSolver.chklazy_race_tree_allow_duplicates.Bind(wx.EVT_CHECKBOX, self.Changelazy_race_tree_allow_duplicates)
        self.optionsnotebook.tabSolver.txtlazy_race_tree_target_list_file.Bind(wx.EVT_KILL_FOCUS,self.Changelazy_race_tree_target_list_file)
        self.optionsnotebook.tabSolver.txtlazy_race_tree_start_location_ID.Bind(wx.EVT_KILL_FOCUS,self.Changelazy_race_tree_start_location_ID)
        self.optionsnotebook.tabSolver.btnlazy_race_tree_target_list_file.Bind(wx.EVT_BUTTON, self.Clicklazy_race_tree_target_list_file)
        self.optionsnotebook.tabSolver.txtlazy_race_tree_maximum_duration.Bind(wx.EVT_KILL_FOCUS,self.Changelazy_race_tree_maximum_duration)
        
        self.optionsnotebook.tabSolver.txttrialX.Bind(wx.EVT_KILL_FOCUS,self.ChangetrialX)
        self.optionsnotebook.tabSolver.btntrialX.Bind(wx.EVT_BUTTON, self.ClickTrialXButton)


        #physics options
        self.optionsnotebook.tabPhysics.cmbephemeris_source.Bind(wx.EVT_COMBOBOX,self.Changeephemeris_source)
        self.optionsnotebook.tabPhysics.txtSPICE_leap_seconds_kernel.Bind(wx.EVT_KILL_FOCUS,self.ChangeSPICE_leap_seconds_kernel)
        self.optionsnotebook.tabPhysics.txtSPICE_reference_frame_kernel.Bind(wx.EVT_KILL_FOCUS,self.ChangeSPICE_reference_frame_kernel)
        self.optionsnotebook.tabPhysics.txtuniverse_folder.Bind(wx.EVT_KILL_FOCUS,self.Changeuniverse_folder)
        self.optionsnotebook.tabPhysics.btnGetNewUniverseFolder.Bind(wx.EVT_BUTTON,self.GetNewUniverseFolder)
        self.optionsnotebook.tabPhysics.btnSetDefaultUniverse.Bind(wx.EVT_BUTTON,self.SetDefaultUniverse)
        self.optionsnotebook.tabPhysics.chkperturb_SRP.Bind(wx.EVT_CHECKBOX,self.Changeperturb_SRP)
        self.optionsnotebook.tabPhysics.chkperturb_thirdbody.Bind(wx.EVT_CHECKBOX,self.Changeperturb_thirdbody)
        self.optionsnotebook.tabPhysics.txtspacecraft_area.Bind(wx.EVT_KILL_FOCUS,self.Changespacecraft_area)
        self.optionsnotebook.tabPhysics.txtcoefficient_of_reflectivity.Bind(wx.EVT_KILL_FOCUS,self.Changecoefficient_of_reflectivity)
        self.optionsnotebook.tabPhysics.cmbspiral_model_type.Bind(wx.EVT_COMBOBOX, self.Changespiral_model_type)

        #output options
        self.optionsnotebook.tabOutput.chkcreate_GMAT_script.Bind(wx.EVT_CHECKBOX, self.Changecreate_GMAT_script)
        self.optionsnotebook.tabOutput.cmboutput_units.Bind(wx.EVT_COMBOBOX, self.Changeoutput_units)
        
        
    #event handlers for global mission options    
    def ChangeMissionName(self, e):
        self.missionoptions.mission_name = self.optionsnotebook.tabGlobal.txtMissionName.GetValue()
        
    def ChangeMissionType(self, e):
        self.missionoptions.mission_type = self.optionsnotebook.tabGlobal.cmbMissionType.GetSelection()
        self.missionoptions.update_all_panels(self.optionsnotebook)

    def Changeobjective_type(self, e):
        self.missionoptions.objective_type = self.optionsnotebook.tabGlobal.cmbobjective_type.GetSelection()
        self.missionoptions.update_all_panels(self.optionsnotebook)

    def Changeinclude_initial_impulse_in_cost(self, e):
        self.missionoptions.include_initial_impulse_in_cost = int(self.optionsnotebook.tabGlobal.chkinclude_initial_impulse_in_cost.GetValue())

    def Changemax_phases_per_journey(self, e):
        self.missionoptions.max_phases_per_journey = eval(self.optionsnotebook.tabGlobal.txtmax_phases_per_journey.GetValue())
        self.missionoptions.update_all_panels(self.optionsnotebook)

    def Changelaunch_window_open_date(self, e):
        self.missionoptions.launch_window_open_date = eval(self.optionsnotebook.tabGlobal.txtlaunch_window_open_date.GetValue())
        self.missionoptions.update_global_options_panel(self.optionsnotebook)

    def ChangeLaunchDateCalendar(self, e):
        CurrentLaunchDate = self.optionsnotebook.tabGlobal.LaunchDateCalendar.GetDate()
        CurrentLaunchDate = CurrentLaunchDate.FromUTC()
        self.missionoptions.launch_window_open_date = CurrentLaunchDate.GetMJD()
        self.missionoptions.update_global_options_panel(self.optionsnotebook)

    def Changenum_timesteps(self, e):
        self.missionoptions.num_timesteps = eval(self.optionsnotebook.tabGlobal.txtnum_timesteps.GetValue())
        self.missionoptions.update_global_options_panel(self.optionsnotebook)

    def Changestep_size_distribution(self, e):
        self.missionoptions.step_size_distribution = self.optionsnotebook.cmbstep_size_distribution.GetSelection()
        self.missionoptions.update_all_panels(self.optionsnotebook)

    def Changestep_size_stdv_or_scale(self, e):
        self.missionoptions.step_size_stdv_or_scale = eval(self.optionsnotebook.tabGlobal.txtstep_size_stdv_or_scale.GetValue())
        self.missionoptions.update_global_options_panel(self.optionsnotebook)

    def Changecontrol_coordinate_system(self, e):
        self.missionoptions.control_coordinate_system = self.optionsnotebook.tabGlobal.cmbcontrol_coordinate_system.GetSelection()
        self.missionoptions.update_global_options_panel(self.optionsnotebook)

    def ChangeDLA_bounds_lower(self, e):
        self.missionoptions.DLA_bounds[0] = eval(self.optionsnotebook.tabGlobal.txtDLA_bounds_lower.GetValue())
        self.missionoptions.update_global_options_panel(self.optionsnotebook)
        
    def ChangeDLA_bounds_upper(self, e):
        self.missionoptions.DLA_bounds[1] = eval(self.optionsnotebook.tabGlobal.txtDLA_bounds_upper.GetValue())
        self.missionoptions.update_global_options_panel(self.optionsnotebook)

    def Changeglobal_timebounded(self, e):
        self.missionoptions.global_timebounded = int(self.optionsnotebook.tabGlobal.chkglobal_timebounded.GetValue())
        self.missionoptions.update_all_panels(self.optionsnotebook)

    def Changetotal_flight_time_bounds_lower(self, e):
        self.missionoptions.total_flight_time_bounds[0] = eval(self.optionsnotebook.tabGlobal.txttotal_flight_time_bounds_lower.GetValue())
        self.missionoptions.update_global_options_panel(self.optionsnotebook)

    def Changetotal_flight_time_bounds_upper(self, e):
        self.missionoptions.total_flight_time_bounds[1] = eval(self.optionsnotebook.tabGlobal.txttotal_flight_time_bounds_upper.GetValue())
        self.missionoptions.update_global_options_panel(self.optionsnotebook)

    def Changeforced_post_launch_coast(self, e):
        self.missionoptions.forced_post_launch_coast = eval(self.optionsnotebook.tabGlobal.txtforced_post_launch_coast.GetValue())
        self.missionoptions.update_global_options_panel(self.optionsnotebook)

    def Changeforced_flyby_coast(self, e):
        self.missionoptions.forced_flyby_coast = eval(self.optionsnotebook.tabGlobal.txtforced_flyby_coast.GetValue())
        self.missionoptions.update_global_options_panel(self.optionsnotebook)

    def Changeinitial_V_infinity_x(self, e):
        self.missionoptions.initial_V_infinity[0] = eval(self.optionsnotebook.tabGlobal.txtinitial_V_infinity_x.GetValue())
        self.missionoptions.update_global_options_panel(self.optionsnotebook)

    def Changeinitial_V_infinity_y(self, e):
        self.missionoptions.initial_V_infinity[1] = eval(self.optionsnotebook.tabGlobal.txtinitial_V_infinity_y.GetValue())
        self.missionoptions.update_global_options_panel(self.optionsnotebook)

    def Changeinitial_V_infinity_z(self, e):
        self.missionoptions.initial_V_infinity[2] = eval(self.optionsnotebook.tabGlobal.txtinitial_V_infinity_z.GetValue())
        self.missionoptions.update_global_options_panel(self.optionsnotebook)

    def Changeminimum_dry_mass(self, e):
        self.missionoptions.minimum_dry_mass = eval(self.optionsnotebook.tabGlobal.txtminimum_dry_mass.GetValue())
        self.missionoptions.update_all_panels(self.optionsnotebook)

    def Changepost_mission_delta_v(self, e):
        self.missionoptions.post_mission_delta_v = eval(self.optionsnotebook.tabGlobal.txtpost_mission_delta_v.GetValue())
        self.missionoptions.update_global_options_panel(self.optionsnotebook)

    #event handlers for spacecraft options
    def Changemaximum_mass(self, e):
        self.missionoptions.maximum_mass = eval(self.optionsnotebook.tabSpacecraft.txtmaximum_mass.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changeallow_initial_mass_to_vary(self, e):
        self.missionoptions.allow_initial_mass_to_vary = int(self.optionsnotebook.tabSpacecraft.chkallow_initial_mass_to_vary.GetValue())
        self.missionoptions.update_all_panels(self.optionsnotebook)

    def ChangeEP_dry_mass(self, e):
        self.missionoptions.EP_dry_mass = eval(self.optionsnotebook.tabSpacecraft.txtEP_dry_mass.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def ChangeLV_type(self, e):
        self.missionoptions.LV_type = self.optionsnotebook.tabSpacecraft.cmbLV_type.GetSelection() - 2
        self.missionoptions.update_all_panels(self.optionsnotebook)

    def ChangeIspChem(self, e):
        self.missionoptions.IspChem = eval(self.optionsnotebook.tabSpacecraft.txtIspChem.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def ChangeIspDS(self, e):
        self.missionoptions.IspDS = eval(self.optionsnotebook.tabSpacecraft.txtIspDS.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changecustom_LV_coefficients0(self, e):
        self.missionoptions.custom_LV_coefficients[0] = eval(self.optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients0.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changecustom_LV_coefficients1(self, e):
        self.missionoptions.custom_LV_coefficients[1] = eval(self.optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients1.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changecustom_LV_coefficients2(self, e):
        self.missionoptions.custom_LV_coefficients[2] = eval(self.optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients2.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changecustom_LV_coefficients3(self, e):
        self.missionoptions.custom_LV_coefficients[3] = eval(self.optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients3.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changecustom_LV_coefficients4(self, e):
        self.missionoptions.custom_LV_coefficients[4] = eval(self.optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients4.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changecustom_LV_coefficients5(self, e):
        self.missionoptions.custom_LV_coefficients[5] = eval(self.optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients5.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changecustom_LV_C3_bounds_lower(self, e):
        self.missionoptions.custom_LV_C3_bounds[0] = eval(self.optionsnotebook.tabSpacecraft.txtcustom_LV_C3_bounds_lower.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changecustom_LV_C3_bounds_upper(self, e):
        self.missionoptions.custom_LV_C3_bounds[1] = eval(self.optionsnotebook.tabSpacecraft.txtcustom_LV_C3_bounds_upper.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changeparking_orbit_altitude(self, e):
        self.missionoptions.parking_orbit_altitude = eval(self.optionsnotebook.tabSpacecraft.txtparking_orbit_altitude.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changeparking_orbit_inclination(self, e):
        self.missionoptions.parking_orbit_inclination = eval(self.optionsnotebook.tabSpacecraft.txtparking_orbit_inclination.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changepost_mission_Isp(self, e):
        self.missionoptions.post_mission_Isp = eval(self.optionsnotebook.tabSpacecraft.txtpost_mission_Isp.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changeenable_maximum_propellant_mass_constraint(self, e):
        self.missionoptions.enable_maximum_propellant_mass_constraint = self.optionsnotebook.tabSpacecraft.chkenable_propellant_mass_constraint.GetValue()
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changemaximum_propellant_mass(self, e):
        self.missionoptions.maximum_propellant_mass = eval(self.optionsnotebook.tabSpacecraft.txtmaximum_propellant_mass.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changepropellant_margin(self, e):
        self.missionoptions.propellant_margin = eval(self.optionsnotebook.tabSpacecraft.txtpropellant_margin.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changepower_margin(self, e):
        self.missionoptions.power_margin = eval(self.optionsnotebook.tabSpacecraft.txtpower_margin.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def ChangeLV_margin(self, e):
        self.missionoptions.LV_margin = eval(self.optionsnotebook.tabSpacecraft.txtLV_margin.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Change_LV_adapter_mass(self, e):
        self.missionoptions.LV_adapter_mass = eval(self.optionsnotebook.tabSpacecraft.txtLV_adapter_mass.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changeengine_type(self, e):
        self.missionoptions.engine_type = self.optionsnotebook.tabSpacecraft.cmbengine_type.GetSelection()
        self.missionoptions.update_all_panels(self.optionsnotebook)

    def Changenumber_of_engines(self, e):
        self.missionoptions.number_of_engines = eval(self.optionsnotebook.tabSpacecraft.txtnumber_of_engines.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changethrottle_logic_mode(self, e):
        self.missionoptions.throttle_logic_mode = self.optionsnotebook.tabSpacecraft.cmbthrottle_logic_mode.GetSelection()
        self.missionoptions.update_all_panels(self.optionsnotebook)

    def Changethrottle_sharpness(self, e):
        self.missionoptions.throttle_sharpness = eval(self.optionsnotebook.tabSpacecraft.txtthrottle_sharpness.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changeengine_duty_cycle(self, e):
        self.missionoptions.engine_duty_cycle = eval(self.optionsnotebook.tabSpacecraft.txtengine_duty_cycle.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def ChangeThrust(self, e):
        self.missionoptions.Thrust = eval(self.optionsnotebook.tabSpacecraft.txtThrust.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def ChangeIspLT(self, e):
        self.missionoptions.IspLT = eval(self.optionsnotebook.tabSpacecraft.txtIspLT.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def ChangeIspLT_minimum(self, e):
        self.missionoptions.IspLT_minimum = eval(self.optionsnotebook.tabSpacecraft.txtIspLT_minimum.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changeuser_defined_engine_efficiency(self, e):
        self.missionoptions.user_defined_engine_efficiency = eval(self.optionsnotebook.tabSpacecraft.txtuser_defined_engine_efficiency.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changeengine_input_thrust_coefficients0(self, e):
        self.missionoptions.engine_input_thrust_coefficients[0] = eval(self.optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients0.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changeengine_input_thrust_coefficients1(self, e):
        self.missionoptions.engine_input_thrust_coefficients[1] = eval(self.optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients1.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changeengine_input_thrust_coefficients2(self, e):
        self.missionoptions.engine_input_thrust_coefficients[2] = eval(self.optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients2.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changeengine_input_thrust_coefficients3(self, e):
        self.missionoptions.engine_input_thrust_coefficients[3] = eval(self.optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients3.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changeengine_input_thrust_coefficients4(self, e):
        self.missionoptions.engine_input_thrust_coefficients[4] = eval(self.optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients4.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changeengine_input_thrust_coefficients5(self, e):
        self.missionoptions.engine_input_thrust_coefficients[5] = eval(self.optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients5.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changeengine_input_thrust_coefficients6(self, e):
        self.missionoptions.engine_input_thrust_coefficients[6] = eval(self.optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients6.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changeengine_input_mass_flow_rate_coefficients0(self, e):
        self.missionoptions.engine_input_mass_flow_rate_coefficients[0] = eval(self.optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients0.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changeengine_input_mass_flow_rate_coefficients1(self, e):
        self.missionoptions.engine_input_mass_flow_rate_coefficients[1] = eval(self.optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients1.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changeengine_input_mass_flow_rate_coefficients2(self, e):
        self.missionoptions.engine_input_mass_flow_rate_coefficients[2] = eval(self.optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients2.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changeengine_input_mass_flow_rate_coefficients3(self, e):
        self.missionoptions.engine_input_mass_flow_rate_coefficients[3] = eval(self.optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients3.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changeengine_input_mass_flow_rate_coefficients4(self, e):
        self.missionoptions.engine_input_mass_flow_rate_coefficients[4] = eval(self.optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients4.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changeengine_input_mass_flow_rate_coefficients5(self, e):
        self.missionoptions.engine_input_mass_flow_rate_coefficients[5] = eval(self.optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients5.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changeengine_input_mass_flow_rate_coefficients6(self, e):
        self.missionoptions.engine_input_mass_flow_rate_coefficients[6] = eval(self.optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients6.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changeengine_input_power_bounds_lower(self, e):
        self.missionoptions.engine_input_power_bounds[0] = eval(self.optionsnotebook.tabSpacecraft.txtengine_input_power_bounds_lower.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changeengine_input_power_bounds_upper(self, e):
        self.missionoptions.engine_input_power_bounds[1] = eval(self.optionsnotebook.tabSpacecraft.txtengine_input_power_bounds_upper.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changepower_at_1_AU(self, e):
        self.missionoptions.power_at_1_AU = eval(self.optionsnotebook.tabSpacecraft.txtpower_at_1_AU.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changepower_source_type(self, e):
        self.missionoptions.power_source_type = self.optionsnotebook.tabSpacecraft.cmbpower_source_type.GetSelection()
        self.missionoptions.update_all_panels(self.optionsnotebook)

    def Changesolar_power_gamma0(self, e):
        self.missionoptions.solar_power_gamma[0] = eval(self.optionsnotebook.tabSpacecraft.txtsolar_power_gamma0.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changesolar_power_gamma1(self, e):
        self.missionoptions.solar_power_gamma[1] = eval(self.optionsnotebook.tabSpacecraft.txtsolar_power_gamma1.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changesolar_power_gamma2(self, e):
        self.missionoptions.solar_power_gamma[2] = eval(self.optionsnotebook.tabSpacecraft.txtsolar_power_gamma2.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changesolar_power_gamma3(self, e):
        self.missionoptions.solar_power_gamma[3] = eval(self.optionsnotebook.tabSpacecraft.txtsolar_power_gamma3.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changesolar_power_gamma4(self, e):
        self.missionoptions.solar_power_gamma[4] = eval(self.optionsnotebook.tabSpacecraft.txtsolar_power_gamma4.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changespacecraft_power_model_type(self, e):
        self.missionoptions.spacecraft_power_model_type = self.optionsnotebook.tabSpacecraft.cmbspacecraft_power_model_type.GetSelection()
        self.missionoptions.update_all_panels(self.optionsnotebook)

    def Changespacecraft_power_coefficients0(self, e):
        self.missionoptions.spacecraft_power_coefficients[0] = eval(self.optionsnotebook.tabSpacecraft.txtspacecraft_power_coefficients0.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changespacecraft_power_coefficients1(self, e):
        self.missionoptions.spacecraft_power_coefficients[1] = eval(self.optionsnotebook.tabSpacecraft.txtspacecraft_power_coefficients1.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changespacecraft_power_coefficients2(self, e):
        self.missionoptions.spacecraft_power_coefficients[2] = eval(self.optionsnotebook.tabSpacecraft.txtspacecraft_power_coefficients2.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)

    def Changepower_decay_rate(self, e):
        self.missionoptions.power_decay_rate = eval(self.optionsnotebook.tabSpacecraft.txtpower_decay_rate.GetValue())
        self.missionoptions.update_spacecraft_options_panel(self.optionsnotebook)


    #event handlers for journey options
    def ChangeJourneySelectBoxChoice(self, e):
        #determine which journey is "active"
        self.missionoptions.ActiveJourney = self.optionsnotebook.tabJourney.JourneySelectBox.GetSelection()

        self.missionoptions.update_all_panels(self.optionsnotebook)

    
    def ClickAddNewJourney(self, e):
        #determine which journey is "active"
        self.missionoptions.ActiveJourney = self.optionsnotebook.tabJourney.JourneySelectBox.GetSelection()

        #add
        temp_JourneyOptions = JO.JourneyOptions(self.missionoptions.mission_type)
        temp_JourneyOptions.sequence = [[0]*self.missionoptions.max_phases_per_journey]*self.missionoptions.number_of_trial_sequences
        self.missionoptions.Journeys.append(temp_JourneyOptions)
        self.missionoptions.number_of_journeys += 1
        self.optionsnotebook.tabJourney.JourneySelectBox.SetSelection(-1)
        self.missionoptions.update_all_panels(self.optionsnotebook)

    def ClickDeleteJourney(self, e):
        #determine which journey is "active"
        self.missionoptions.ActiveJourney = self.optionsnotebook.tabJourney.JourneySelectBox.GetSelection()

        #delete
        self.missionoptions.Journeys.pop(self.missionoptions.ActiveJourney)
        self.missionoptions.number_of_journeys -= 1

        if self.missionoptions.ActiveJourney > self.missionoptions.number_of_journeys - 1:
            self.missionoptions.ActiveJourney -= 1

        self.optionsnotebook.tabJourney.JourneySelectBox.SetSelection(0)
        self.missionoptions.update_all_panels(self.optionsnotebook)

    def ClickMoveJourneyUp(self, e):
        #determine which journey is "active"
        self.missionoptions.ActiveJourney = self.optionsnotebook.tabJourney.JourneySelectBox.GetSelection()

        #move up
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney], self.missionoptions.Journeys[self.missionoptions.ActiveJourney-1] = self.missionoptions.Journeys[self.missionoptions.ActiveJourney-1], self.missionoptions.Journeys[self.missionoptions.ActiveJourney]
        self.missionoptions.ActiveJourney -= 1
        self.optionsnotebook.tabJourney.JourneySelectBox.SetSelection(self.missionoptions.ActiveJourney)
        self.missionoptions.update_all_panels(self.optionsnotebook)

    def ClickMoveJourneyDown(self, e):
        #determine which journey is "active"
        self.missionoptions.ActiveJourney = self.optionsnotebook.tabJourney.JourneySelectBox.GetSelection()

        #move down
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney], self.missionoptions.Journeys[self.missionoptions.ActiveJourney+1] = self.missionoptions.Journeys[self.missionoptions.ActiveJourney+1], self.missionoptions.Journeys[self.missionoptions.ActiveJourney]
        self.missionoptions.ActiveJourney += 1
        self.optionsnotebook.tabJourney.JourneySelectBox.SetSelection(self.missionoptions.ActiveJourney)
        self.missionoptions.update_all_panels(self.optionsnotebook)

    def Changejourney_names(self, e):
         namestring = self.optionsnotebook.tabJourney.txtjourney_names.GetValue().replace(' ', '_')
         self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_names = namestring
         self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def Changejourney_central_body(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_central_body = self.optionsnotebook.tabJourney.txtjourney_central_body.GetValue()

    def Changedestination_list(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].destination_list = eval(self.optionsnotebook.tabJourney.txtdestination_list.GetValue())
        self.missionoptions.update_all_panels(self.optionsnotebook)

    def Changejourney_starting_mass_increment(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_starting_mass_increment = eval(self.optionsnotebook.tabJourney.txtjourney_starting_mass_increment.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def Changejourney_variable_mass_increment(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_variable_mass_increment = int(self.optionsnotebook.tabJourney.chkjourney_variable_mass_increment.GetValue())

    def Changejourney_wait_time_bounds_lower(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_wait_time_bounds[0] = eval(self.optionsnotebook.tabJourney.txtjourney_wait_time_bounds_lower.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def Changejourney_wait_time_bounds_upper(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_wait_time_bounds[1] = eval(self.optionsnotebook.tabJourney.txtjourney_wait_time_bounds_upper.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def Changejourney_timebounded(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_timebounded = self.optionsnotebook.tabJourney.cmbjourney_timebounded.GetSelection()
        self.missionoptions.update_all_panels(self.optionsnotebook)

    def Changejourney_flight_time_bounds_lower(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_flight_time_bounds[0] = eval(self.optionsnotebook.tabJourney.txtjourney_flight_time_bounds_lower.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def Changejourney_flight_time_bounds_upper(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_flight_time_bounds[1] = eval(self.optionsnotebook.tabJourney.txtjourney_flight_time_bounds_upper.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def Changejourney_arrival_date_bounds_lower(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_date_bounds[0] = eval(self.optionsnotebook.tabJourney.txtjourney_arrival_date_bounds_lower.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeArrivalDateLowerCalendar(self, e):
        date = self.optionsnotebook.tabJourney.ArrivalDateLowerCalendar.GetDate()
        date = date.FromUTC()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_date_bounds[0] = date.GetMJD()
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeArrivalDateUpperCalendar(self, e):
        date = self.optionsnotebook.tabJourney.ArrivalDateUpperCalendar.GetDate()
        date = date.FromUTC()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_date_bounds[1] = date.GetMJD()
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def Changejourney_arrival_date_bounds_upper(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_date_bounds[1] = eval(self.optionsnotebook.tabJourney.txtjourney_arrival_date_bounds_upper.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def Changejourney_initial_impulse_bounds_lower(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_initial_impulse_bounds[0] = eval(self.optionsnotebook.tabJourney.txtjourney_initial_impulse_bounds_lower.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def Changejourney_initial_impulse_bounds_upper(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_initial_impulse_bounds[1] = eval(self.optionsnotebook.tabJourney.txtjourney_initial_impulse_bounds_upper.GetValue())
        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_initial_impulse_bounds[1] < 1.0e-8:
            self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_initial_impulse_bounds[1] = 1.0e-8

        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def Changejourney_departure_type(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_departure_type = self.optionsnotebook.tabJourney.cmbjourney_departure_type.GetSelection()
        self.missionoptions.update_all_panels(self.optionsnotebook)

    def Changejourney_escape_spiral_radius(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_escape_spiral_starting_radius = eval(self.optionsnotebook.tabJourney.txtjourney_escape_spiral_starting_radius.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def Changejourney_initial_velocity0(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_initial_velocity[0] = eval(self.optionsnotebook.tabJourney.txtjourney_initial_velocity0.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def Changejourney_initial_velocity1(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_initial_velocity[1] = eval(self.optionsnotebook.tabJourney.txtjourney_initial_velocity1.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def Changejourney_initial_velocity2(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_initial_velocity[2] = eval(self.optionsnotebook.tabJourney.txtjourney_initial_velocity2.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def Changejourney_arrival_type(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_type = self.optionsnotebook.tabJourney.cmbjourney_arrival_type.GetSelection()
        self.missionoptions.update_all_panels(self.optionsnotebook)

    def Changejourney_capture_spiral_radius(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_capture_spiral_final_radius = eval(self.optionsnotebook.tabJourney.txtjourney_capture_spiral_final_radius.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def Changejourney_final_velocity0(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_final_velocity[0] = eval(self.optionsnotebook.tabJourney.txtjourney_final_velocity0.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def Changejourney_final_velocity1(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_final_velocity[1] = eval(self.optionsnotebook.tabJourney.txtjourney_final_velocity1.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def Changejourney_final_velocity2(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_final_velocity[2] = eval(self.optionsnotebook.tabJourney.txtjourney_final_velocity2.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def Changejourney_arrival_declination_constraint_flag(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_declination_constraint_flag = int(self.optionsnotebook.tabJourney.chkjourney_arrival_declination_constraint_flag.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def Changejourney_arrival_declination_bounds_lower(self, e):
        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_declination_constraint_flag:
            self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_declination_bounds[0] = eval(self.optionsnotebook.tabJourney.txtjourney_arrival_declination_bounds_lower.GetValue())
            self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def Changejourney_arrival_declination_bounds_upper(self, e):
        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_declination_constraint_flag:
            self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_declination_bounds[1] = eval(self.optionsnotebook.tabJourney.txtjourney_arrival_declination_bounds_upper.GetValue())
            self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def Changesequence(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].sequence = eval(self.optionsnotebook.tabJourney.txtsequence.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def Changejourney_perturbation_bodies(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].sequence = eval(self.optionsnotebook.tabJourney.txtsequence.GetValue())

    def Changejourney_departure_elements_type(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_departure_elements_type = self.optionsnotebook.tabJourney.cmbjourney_departure_elements_type.GetSelection()
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangevarySMA_departure(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_departure_elements_vary_flag[0] = int(self.optionsnotebook.tabJourney.chkSMA_departure.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangevaryECC_departure(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_departure_elements_vary_flag[1] = int(self.optionsnotebook.tabJourney.chkECC_departure.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangevaryINC_departure(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_departure_elements_vary_flag[2] = int(self.optionsnotebook.tabJourney.chkINC_departure.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangevaryRAAN_departure(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_departure_elements_vary_flag[3] = int(self.optionsnotebook.tabJourney.chkRAAN_departure.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangevaryAOP_departure(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_departure_elements_vary_flag[4] = int(self.optionsnotebook.tabJourney.chkAOP_departure.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangevaryMA_departure(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_departure_elements_vary_flag[5] = int(self.optionsnotebook.tabJourney.chkMA_departure.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeSMA_departure(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_departure_elements[0] = eval(self.optionsnotebook.tabJourney.txtSMA_departure.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeECC_departure(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_departure_elements[1] = eval(self.optionsnotebook.tabJourney.txtECC_departure.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeINC_departure(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_departure_elements[2] = eval(self.optionsnotebook.tabJourney.txtINC_departure.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeRAAN_departure(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_departure_elements[3] = eval(self.optionsnotebook.tabJourney.txtRAAN_departure.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeAOP_departure(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_departure_elements[4] = eval(self.optionsnotebook.tabJourney.txtAOP_departure.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeMA_departure(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_departure_elements[5] = eval(self.optionsnotebook.tabJourney.txtMA_departure.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeSMA_departure0(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_departure_elements_bounds[0] = eval(self.optionsnotebook.tabJourney.txtSMA_departure0.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeECC_departure0(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_departure_elements_bounds[2] = eval(self.optionsnotebook.tabJourney.txtECC_departure0.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeINC_departure0(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_departure_elements_bounds[4] = eval(self.optionsnotebook.tabJourney.txtINC_departure0.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeRAAN_departure0(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_departure_elements_bounds[6] = eval(self.optionsnotebook.tabJourney.txtRAAN_departure0.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeAOP_departure0(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_departure_elements_bounds[8] = eval(self.optionsnotebook.tabJourney.txtAOP_departure0.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeMA_departure0(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_departure_elements_bounds[10] = eval(self.optionsnotebook.tabJourney.txtMA_departure0.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeSMA_departure1(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_departure_elements_bounds[1] = eval(self.optionsnotebook.tabJourney.txtSMA_departure1.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeECC_departure1(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_departure_elements_bounds[3] = eval(self.optionsnotebook.tabJourney.txtECC_departure1.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeINC_departure1(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_departure_elements_bounds[5] = eval(self.optionsnotebook.tabJourney.txtINC_departure1.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeRAAN_departure1(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_departure_elements_bounds[7] = eval(self.optionsnotebook.tabJourney.txtRAAN_departure1.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeAOP_departure1(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_departure_elements_bounds[9] = eval(self.optionsnotebook.tabJourney.txtAOP_departure1.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeMA_departure1(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_departure_elements_bounds[11] = eval(self.optionsnotebook.tabJourney.txtMA_departure1.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def Changejourney_arrival_elements_type(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_elements_type = self.optionsnotebook.tabJourney.cmbjourney_arrival_elements_type.GetSelection()
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangevarySMA_arrival(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_elements_vary_flag[0] = int(self.optionsnotebook.tabJourney.chkSMA_arrival.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangevaryECC_arrival(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_elements_vary_flag[1] = int(self.optionsnotebook.tabJourney.chkECC_arrival.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangevaryINC_arrival(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_elements_vary_flag[2] = int(self.optionsnotebook.tabJourney.chkINC_arrival.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangevaryRAAN_arrival(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_elements_vary_flag[3] = int(self.optionsnotebook.tabJourney.chkRAAN_arrival.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangevaryAOP_arrival(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_elements_vary_flag[4] = int(self.optionsnotebook.tabJourney.chkAOP_arrival.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangevaryMA_arrival(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_elements_vary_flag[5] = int(self.optionsnotebook.tabJourney.chkMA_arrival.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeSMA_arrival(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_elements[0] = eval(self.optionsnotebook.tabJourney.txtSMA_arrival.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeECC_arrival(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_elements[1] = eval(self.optionsnotebook.tabJourney.txtECC_arrival.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeINC_arrival(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_elements[2] = eval(self.optionsnotebook.tabJourney.txtINC_arrival.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeRAAN_arrival(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_elements[3] = eval(self.optionsnotebook.tabJourney.txtRAAN_arrival.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeAOP_arrival(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_elements[4] = eval(self.optionsnotebook.tabJourney.txtAOP_arrival.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeMA_arrival(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_elements[5] = eval(self.optionsnotebook.tabJourney.txtMA_arrival.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeSMA_arrival0(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_elements_bounds[0] = eval(self.optionsnotebook.tabJourney.txtSMA_arrival0.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeECC_arrival0(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_elements_bounds[2] = eval(self.optionsnotebook.tabJourney.txtECC_arrival0.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeINC_arrival0(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_elements_bounds[4] = eval(self.optionsnotebook.tabJourney.txtINC_arrival0.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeRAAN_arrival0(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_elements_bounds[6] = eval(self.optionsnotebook.tabJourney.txtRAAN_arrival0.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeAOP_arrival0(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_elements_bounds[8] = eval(self.optionsnotebook.tabJourney.txtAOP_arrival0.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeMA_arrival0(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_elements_bounds[10] = eval(self.optionsnotebook.tabJourney.txtMA_arrival0.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeSMA_arrival1(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_elements_bounds[1] = eval(self.optionsnotebook.tabJourney.txtSMA_arrival1.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeECC_arrival1(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_elements_bounds[3] = eval(self.optionsnotebook.tabJourney.txtECC_arrival1.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeINC_arrival1(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_elements_bounds[5] = eval(self.optionsnotebook.tabJourney.txtINC_arrival1.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeRAAN_arrival1(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_elements_bounds[7] = eval(self.optionsnotebook.tabJourney.txtRAAN_arrival1.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeAOP_arrival1(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_elements_bounds[9] = eval(self.optionsnotebook.tabJourney.txtAOP_arrival1.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def ChangeMA_arrival1(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_arrival_elements_bounds[11] = eval(self.optionsnotebook.tabJourney.txtMA_arrival1.GetValue())
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def Clickdestination_list(self, e):
        #call dialog to choose destination list
        self.universe = Universe.Universe(os.path.join(self.missionoptions.universe_folder, self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_central_body + ".emtg_universe"))
        dlg = BodyPicker.DestinationPicker(self, -1,
                                           self.universe,
                                           self.missionoptions.Journeys[self.missionoptions.ActiveJourney].destination_list[0],
                                           self.missionoptions.Journeys[self.missionoptions.ActiveJourney].destination_list[1])

        dlg.ShowModal()

        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].destination_list[0] = dlg.destination1
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].destination_list[1] = dlg.destination2

        dlg.Destroy()
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def Clickjourney_central_body(self, e):
        #call dialog to choose destination list
        dlg = wx.FileDialog(self, "Choose an emtg_universe file", self.missionoptions.universe_folder, self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_central_body+".emtg_universe", "*.emtg_universe", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
        else:
            filename = self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_central_body+".emtg_universe"

        fileparts = filename.split(".")
        dlg.Destroy()

        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_central_body = fileparts[0]

        self.missionoptions.update_journey_options_panel(self.optionsnotebook)

    def Clicksequence(self, e):
        #call dialog to choose sequence list
        self.universe = Universe.Universe(os.path.join(self.missionoptions.universe_folder, self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_central_body + ".emtg_universe"))
        dlg = BodyPicker.SequencePicker(self, -1,
                                        self.universe,
                                        self.missionoptions,
                                        self.missionoptions.ActiveJourney)

        dlg.ShowModal()

        dlg.Destroy()

        self.missionoptions.update_all_panels(self.optionsnotebook)

    def Clickjourney_perturbation_bodies(self, e):
        #call dialog to choose perturbation list
        self.universe = Universe.Universe(os.path.join(self.missionoptions.universe_folder, self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_central_body + ".emtg_universe"))

        dlg = wx.MultiChoiceDialog(self, "", "Choose perturbation bodies", choices=self.universe.perturbation_menu)

        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_perturbation_bodies != [0]:
            selections = []
            for i in self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_perturbation_bodies:
                selections.append(i-1)
            dlg.SetSelections(selections)

        if dlg.ShowModal() == wx.ID_OK:
            PertubationListIndices = dlg.GetSelections()
            if len(PertubationListIndices) > 0:
                self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_perturbation_bodies = []
                for i in PertubationListIndices:
                    self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_perturbation_bodies.append(self.universe.perturbation_indices[i]+1)
                self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_number_of_perturbation_bodies = len(PertubationListIndices)
            else:
                self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_perturbation_bodies = [0]
                self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_number_of_perturbation_bodies = 1

        dlg.Destroy()
        self.missionoptions.update_journey_options_panel(self.optionsnotebook)


    #event handlers for solver options
    def ChangeInnerLoopSolver(self, e):
        self.missionoptions.run_inner_loop = self.optionsnotebook.tabSolver.cmbInnerLoopSolver.GetSelection()
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)

    def ChangeNLP_solver_type(self, e):
        self.missionoptions.NLP_solver_type = self.optionsnotebook.tabSolver.cmbNLP_solver_type.GetSelection()
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)

    def ChangeNLP_solver_mode(self, e):
        self.missionoptions.NLP_solver_mode = self.optionsnotebook.tabSolver.cmbNLP_solver_mode.GetSelection()
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)

    def Changequiet_NLP(self, e):
        self.missionoptions.quiet_NLP = int(self.optionsnotebook.tabSolver.chkquiet_NLP.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)

    def Changequiet_MBH(self, e):
        self.missionoptions.quiet_basinhopping = int(self.optionsnotebook.tabSolver.chkquiet_MBH.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)

    def ChangeMBH_two_step(self, e):
        self.missionoptions.MBH_two_step = int(self.optionsnotebook.tabSolver.chkMBH_two_step.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)

    def ChangeFD_stepsize(self, e):
        self.missionoptions.FD_stepsize = eval(self.optionsnotebook.tabSolver.txtFD_stepsize.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)

    def ChangeFD_stepsize_coarse(self, e):
        self.missionoptions.FD_stepsize_coarse = eval(self.optionsnotebook.tabSolver.txtFD_stepsize_coarse.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)

    def ChangeACE_feasible_point_finder(self, e):
        self.missionoptions.ACE_feasible_point_finder = int(self.optionsnotebook.tabSolver.chkACE_feasible_point_finder.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)
        
    def ChangeMBH_max_not_improve(self, e):
        self.missionoptions.MBH_max_not_improve = eval(self.optionsnotebook.tabSolver.txtMBH_max_not_improve.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)

    def ChangeMBH_max_trials(self, e):
        self.missionoptions.MBH_max_trials = eval(self.optionsnotebook.tabSolver.txtMBH_max_trials.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)
        
    def ChangeMBH_max_run_time(self, e):
        self.missionoptions.MBH_max_run_time = eval(self.optionsnotebook.tabSolver.txtMBH_max_run_time.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)
        
    def ChangeMBH_max_step_size(self, e):
        self.missionoptions.MBH_max_step_size = eval(self.optionsnotebook.tabSolver.txtMBH_max_step_size.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)

    def ChangeMBH_hop_distribution(self, e):
        self.missionoptions.MBH_hop_distribution = self.optionsnotebook.tabSolver.cmbMBH_hop_distribution.GetSelection()
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)

    def ChangeMBH_Pareto_alpha(self, e):
        self.missionoptions.MBH_Pareto_alpha = self.optionsnotebook.tabSolver.txtMBH_Pareto_alpha.GetValue()
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)

    def ChangeMBH_time_hop_probability(self, e):
        self.missionoptions.MBH_time_hop_probability = eval(self.optionsnotebook.tabSolver.txtMBH_time_hop_probability.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)
      
    def Changesnopt_feasibility_tolerance(self, e):
        self.missionoptions.snopt_feasibility_tolerance = eval(self.optionsnotebook.tabSolver.txtsnopt_feasibility_tolerance.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)
                        
    def Changesnopt_major_iterations(self, e):
        self.missionoptions.snopt_major_iterations = eval(self.optionsnotebook.tabSolver.txtsnopt_major_iterations.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)
        
    def Changesnopt_max_run_time(self, e):
        self.missionoptions.snopt_max_run_time = eval(self.optionsnotebook.tabSolver.txtsnopt_max_run_time.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)
        
    def Changederivative_type(self, e):
        self.missionoptions.derivative_type = self.optionsnotebook.tabSolver.cmbderivative_type.GetSelection()
        
    def ChangeCheckDerivatives(self, e):
        self.missionoptions.check_derivatives = int(self.optionsnotebook.tabSolver.chkcheck_derivatives.GetValue())
        
    def ChangeSeedMBH(self, e):
        self.missionoptions.seed_MBH = int(self.optionsnotebook.tabSolver.chkseed_MBH.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)

    def Changeinitial_guess_control_coordinate_system(self, e):
        self.missionoptions.initial_guess_control_coordinate_system = self.optionsnotebook.tabSolver.cmbinitial_guess_control_coordinate_system.GetSelection()
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)
    
    def ChangeInterpolateInitialGuess(self, e):
        self.missionoptions.interpolate_initial_guess = int(self.optionsnotebook.tabSolver.chkinterpolate_initial_guess.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)
        
    def Changeinitial_guess_num_timesteps(self, e):
        self.missionoptions.initial_guess_num_timesteps = eval(self.optionsnotebook.tabSolver.txtinitial_guess_num_timesteps.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)
        
    def Changeinitial_guess_step_size_distribution(self, e):
        self.missionoptions.initial_guess_step_size_distribution = self.optionsnotebook.tabSolver.cmbinitial_guess_step_size_distribution.GetSelection()
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)
        
    def Changeinitial_guess_step_size_stdv_or_scale(self, e):
        self.missionoptions.initial_guess_step_size_stdv_or_scale = eval(self.optionsnotebook.tabSolver.txtinitial_guess_step_size_stdv_or_scale.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)
    
    def ChangeMBH_zero_control_initial_guess(self, e):
        self.missionoptions.MBH_zero_control_initial_guess = self.optionsnotebook.tabSolver.cmbMBH_zero_control_initial_guess.GetSelection()
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)
        
    def Changerun_outerloop(self, e):
        self.missionoptions.run_outerloop = int(self.optionsnotebook.tabSolver.cmbrun_outerloop.GetSelection())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)

    def Changelazy_race_tree_allow_duplicates(self, e):
        self.missionoptions.lazy_race_tree_allow_duplicates = int(self.optionsnotebook.tabSolver.chklazy_race_tree_allow_duplicates.GetValue)
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)

    def Changelazy_race_tree_target_list_file(self, e):
        self.missionoptions.lazy_race_tree_target_list_file = self.optionsnotebook.tabSolver.txtlazy_race_tree_target_list_file.GetValue
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)

    def Changelazy_race_tree_start_location_ID(self, e):
        self.missionoptions.lazy_race_tree_start_location_ID = int(self.optionsnotebook.tabSolver.txtlazy_race_tree_start_location_ID.GetValue)
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)

    def Changelazy_race_tree_maximum_duration(self, e):
        self.missionoptions.lazy_race_tree_maximum_duration = eval(self.optionsnotebook.tabSolver.txtlazy_race_tree_maximum_duration.GetValue)
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)

    def Clicklazy_race_tree_target_list_file(self, e):
        #call dialog to choose destination list
        dlg = wx.FileDialog(self, "Choose a lazy race-tree search input file", self.dirname, "", "*.txt", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
        else:
            filename = "none"

        dlg.Destroy()

        self.missionoptions.lazy_race_tree_target_list_file = filename

        self.missionoptions.update_solver_options_panel(self.optionsnotebook)
        
    def Changeouterloop_popsize(self, e):
        self.missionoptions.outerloop_popsize=eval(self.optionsnotebook.tabSolver.txtouterloop_popsize.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)
        
    def Changeouterloop_genmax(self, e):
        self.missionoptions.outerloop_genmax=eval(self.optionsnotebook.tabSolver.txtouterloop_genmax.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)
        
    def Changeouterloop_tournamentsize(self, e):
        self.missionoptions.outerloop_tournamentsize=eval(self.optionsnotebook.tabSolver.txtouterloop_tournamentsize.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)
        
    def Changeouterloop_CR(self, e):
        self.missionoptions.outerloop_CR=eval(self.optionsnotebook.tabSolver.txtouterloop_CR.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)
        
    def Changeouterloop_mu(self, e):
        self.missionoptions.outerloop_mu=eval(self.optionsnotebook.tabSolver.txtouterloop_mu.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)
        
    def Changeouterloop_stallmax(self, e):
        self.missionoptions.outerloop_stallmax=eval(self.optionsnotebook.tabSolver.txtouterloop_stallmax.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)
        
    def Changeouterloop_tolfit(self, e):
        self.missionoptions.outerloop_tolfit=eval(self.optionsnotebook.tabSolver.txtouterloop_tolfit.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)
        
    def Changeouterloop_ntrials(self, e):
        self.missionoptions.outerloop_ntrials=eval(self.optionsnotebook.tabSolver.txtouterloop_ntrials.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)
        
    def Changeouterloop_elitecount(self, e):
        self.missionoptions.outerloop_elitecount=eval(self.optionsnotebook.tabSolver.txtouterloop_elitecount.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)
        
    def Changeouterloop_warmstart(self, e):
        self.missionoptions.outerloop_warmstart=eval(self.optionsnotebook.tabSolver.txtouterloop_warmstart.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)
        
    def ChangetrialX(self, e):
        self.missionoptions.trialX = eval(self.optionsnotebook.tabSolver.txttrialX.GetValue())
        self.missionoptions.update_solver_options_panel(self.optionsnotebook)

    def ClickTrialXButton(self, e):
        dlg = wx.FileDialog(self, "Choose a mission file", self.dirname, "", "*.emtg", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()

            #now read the file into a mission object and extract the trialX entry
            tempmission = Mission.Mission(os.path.join(self.dirname, self.filename))
            self.missionoptions.trialX = [copy.deepcopy(tempmission.DecisionVector)]
            del tempmission

        self.missionoptions.update_solver_options_panel(self.optionsnotebook)

    #event handlers for physics options
    def Changeephemeris_source(self, e):
        self.missionoptions.ephemeris_source = self.optionsnotebook.tabPhysics.cmbephemeris_source.GetSelection()

    def ChangeSPICE_leap_seconds_kernel(self, e):
        self.missionoptions.SPICE_leap_seconds_kernel = self.optionsnotebook.tabPhysics.txtSPICE_leap_seconds_kernel.GetValue()

    def ChangeSPICE_reference_frame_kernel(self, e):
        self.missionoptions.SPICE_reference_frame_kernel = self.optionsnotebook.tabPhysics.txtSPICE_reference_frame_kernel.GetValue()

    def Changeuniverse_folder(self, e):
        self.missionoptions.universe_folder = self.optionsnotebook.tabPhysics.txtuniverse_folder.GetValue()

    def GetNewUniverseFolder(self, e):
        #file load dialog to get name of universe folder
        dlg = wx.DirDialog(self, "Choose a Universe folder", self.dirname)
        if dlg.ShowModal() == wx.ID_OK:
            self.missionoptions.universe_folder = dlg.GetPath()
            self.optionsnotebook.tabPhysics.txtuniverse_folder.SetValue(self.missionoptions.universe_folder)
        dlg.Destroy()

    def SetDefaultUniverse(self, e):
        self.missionoptions.universe_folder = self.default_universe_path
        self.optionsnotebook.tabPhysics.txtuniverse_folder.SetValue(self.default_universe_path)

    def Changeperturb_SRP(self, e):
        self.missionoptions.perturb_SRP = int(self.optionsnotebook.tabPhysics.chkperturb_SRP.GetValue())
        self.missionoptions.update_physics_options_panel(self.optionsnotebook)

    def Changeperturb_thirdbody(self, e):
        self.missionoptions.perturb_thirdbody = int(self.optionsnotebook.tabPhysics.chkperturb_thirdbody.GetValue())
        self.missionoptions.update_all_panels(self.optionsnotebook)

    def Changespacecraft_area(self, e):
        self.missionoptions.spacecraft_area = float(self.optionsnotebook.tabPhysics.txtspacecraft_area.GetValue())

    def Changecoefficient_of_reflectivity(self, e):
        self.missionoptions.coefficient_of_reflectivity = float(self.optionsnotebook.tabPhysics.txtcoefficient_of_reflectivity.GetValue())

    def Changespiral_model_type(self, e):
        self.missionoptions.spiral_model_type = self.optionsnotebook.tabPhysics.cmbspiral_model_type.GetSelection()

    #handlers for output options
    def Changecreate_GMAT_script(self, e):
        create_GMAT_script = int(self.optionsnotebook.tabOutput.chkcreate_GMAT_script.GetValue())

    def Changeoutput_units(self, e):
        output_units = int(self.optionsnotebook.tabOutput.cmboutput_units.GetSelection())