import wx
import wx.calendar
import wx.lib.scrolledpanel
import MissionOptions as MO


class GlobalOptionsPanel(wx.lib.scrolledpanel.ScrolledPanel):
    def __init__(self, parent):

        wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent)
        
        globaloptionsgrid = wx.FlexGridSizer(20,2,5,5)
        self.lblMissionName = wx.StaticText(self, -1, "Mission Name")
        self.txtMissionName = wx.TextCtrl(self, -1, "mission_name", size=(500,-1))

        self.lblMissionType = wx.StaticText(self, -1, "Mission Type")
        phasetypes = ['0: MGA','1: MGA-DSM','2: MGA-LT','3: FBLT','4: MGA-NDSM (experimental)']#,'6: DTLT']
                    #,'6: solver chooses (MGA, MGA-DSM)','7: solver chooses (MGA, MGA-LT)',
                    #'8: solver chooses (MGA-DSM, MGA-LT)','9: solver chooses (MGA, MGA-DSM, MGA-LT)']
        self.cmbMissionType = wx.ComboBox(self, -1, choices=phasetypes, style=wx.CB_READONLY)

        self.lblmaximum_number_of_lambert_revolutions = wx.StaticText(self, -1, "Maximum number of revolutions for solving Lambert's problem")
        self.txtmaximum_number_of_lambert_revolutions = wx.TextCtrl(self, -1, "maximum_number_of_lambert_revolutions")

        self.lblobjective_type = wx.StaticText(self, -1, "Include initial impulse in cost")
        objectivetypes = ['0: minimum deltaV','1: minimum time','2: maximum final mass','3: GTOC 1 asteroid deflection function',
                          '4: launch as late as possible in the window','5: launch as early as possible in the window',
                          '6: maximize orbit energy','7: minimize launch mass','8: arrive as early as possible',
                          '9: arrive as late as possible','10: minimum propellant (not the same as 2)','11: maximum dry/wet ratio',
                          '12: maximum arrival kinetic energy', '13: minimum BOL power']
        self.cmbobjective_type = wx.ComboBox(self, -1, choices=objectivetypes, style = wx.CB_READONLY)

        self.lblinclude_initial_impulse_in_cost = wx.StaticText(self, -1, "Include initial impulse in cost")
        self.chkinclude_initial_impulse_in_cost = wx.CheckBox(self, -1)

        self.lblmax_phases_per_journey = wx.StaticText(self, -1, "Maximum number of phases per journey")
        self.txtmax_phases_per_journey = wx.TextCtrl(self, -1, "max_phases_per_journey")
        
        self.lbllaunch_window_open_date = wx.StaticText(self, -1, "Launch window open date")
        self.txtlaunch_window_open_date = wx.TextCtrl(self, -1, "launch_window_open_date")
        self.LaunchDateCalendar = wx.calendar.CalendarCtrl(self, -1)
        calendarbox = wx.BoxSizer(wx.HORIZONTAL)
        calendarbox.AddMany([self.txtlaunch_window_open_date, self.LaunchDateCalendar])
        
        self.lblnum_timesteps = wx.StaticText(self, -1, "Number of time-steps")
        self.txtnum_timesteps = wx.TextCtrl(self, -1, "num_timesteps")

        self.lblstep_size_distribution = wx.StaticText(self, -1, "Step size distribution")
        distributionchoices = ["Uniform","Gaussian","Cauchy"]
        self.cmbstep_size_distribution = wx.ComboBox(self, -1, choices = distributionchoices, style=wx.CB_READONLY)

        self.lblstep_size_stdv_or_scale = wx.StaticText(self, -1, "Scale width/standard deviation")
        self.txtstep_size_stdv_or_scale = wx.TextCtrl(self, -1, "step_size_stdv_or_scale")

        self.lblcontrol_coordinate_system = wx.StaticText(self, -1, "Control coordinate system")
        control_coordinate_choices = ['Cartesian','Polar']
        self.cmbcontrol_coordinate_system = wx.ComboBox(self, -1, choices = control_coordinate_choices, style=wx.CB_READONLY)
                
        globaloptionsgrid.AddMany(  [self.lblMissionName, self.txtMissionName,
                                    self.lblMissionType, self.cmbMissionType,
                                    self.lblmaximum_number_of_lambert_revolutions, self.txtmaximum_number_of_lambert_revolutions,
                                    self.lblobjective_type, self.cmbobjective_type,
                                    self.lblinclude_initial_impulse_in_cost, self.chkinclude_initial_impulse_in_cost,
                                    self.lblmax_phases_per_journey, self.txtmax_phases_per_journey,
                                    self.lbllaunch_window_open_date, calendarbox,
                                    self.lblnum_timesteps, self.txtnum_timesteps,
                                    self.lblstep_size_distribution, self.cmbstep_size_distribution,
                                    self.lblstep_size_stdv_or_scale, self.txtstep_size_stdv_or_scale,
                                    self.lblcontrol_coordinate_system, self.cmbcontrol_coordinate_system])
        globaloptionsgrid.SetFlexibleDirection(wx.BOTH)

        #constraint fields
        constraintgrid = wx.FlexGridSizer(20, 2, 5, 5)

        self.lblDLA_bounds = wx.StaticText(self, -1, "DLA bounds (degrees)")
        self.txtDLA_bounds_lower = wx.TextCtrl(self, -1, "DLA_bounds[0]")
        self.txtDLA_bounds_upper = wx.TextCtrl(self, -1, "DLA_bounds[1]")
        DLAbox = wx.BoxSizer(wx.HORIZONTAL)
        DLAbox.AddMany([self.txtDLA_bounds_lower, self.txtDLA_bounds_upper])

        self.lblglobal_timebounded = wx.StaticText(self, -1, "Enable mission time bounds")
        self.chkglobal_timebounded = wx.CheckBox(self, -1)

        self.lbltotal_flight_time_bounds = wx.StaticText(self, -1, "Global flight time bounds")
        self.txttotal_flight_time_bounds_lower = wx.TextCtrl(self, -1, "total_flight_time_bounds[0]")
        self.txttotal_flight_time_bounds_upper = wx.TextCtrl(self, -1, "total_flight_time_bounds[1]")
        GlobalTimebox = wx.BoxSizer(wx.HORIZONTAL)
        GlobalTimebox.AddMany([self.txttotal_flight_time_bounds_lower, self.txttotal_flight_time_bounds_upper])

        self.lblforced_post_launch_coast = wx.StaticText(self, -1, "Forced post-launch coast duration (days)")
        self.txtforced_post_launch_coast = wx.TextCtrl(self, -1, "forced_post_launch_coast")

        self.lblforced_flyby_coast = wx.StaticText(self, -1, "Forced pre/post-flyby coast duration (days)")
        self.txtforced_flyby_coast = wx.TextCtrl(self, -1, "forced_post_launch_coast")
    
        self.lblinitial_V_infinity = wx.StaticText(self, -1, "Initial V-infinity in MJ2000 km/s")
        self.txtinitial_V_infinity_x = wx.TextCtrl(self, -1, "initial_V_infinity[0]")
        self.txtinitial_V_infinity_y = wx.TextCtrl(self, -1, "initial_V_infinity[1]")
        self.txtinitial_V_infinity_z = wx.TextCtrl(self, -1, "initial_V_infinity[2]")
        initial_V_infinity_box = wx.BoxSizer(wx.HORIZONTAL)
        initial_V_infinity_box.AddMany([self.txtinitial_V_infinity_x, self.txtinitial_V_infinity_y, self.txtinitial_V_infinity_z])

        self.lblminimum_dry_mass = wx.StaticText(self, -1, "Minimum dry mass (kg)")
        self.txtminimum_dry_mass = wx.TextCtrl(self, -1, "minimum_dry_mass")

        self.lblpost_mission_delta_v = wx.StaticText(self, -1, "Post-mission delta-v (km/s)")
        self.txtpost_mission_delta_v = wx.TextCtrl(self, -1, "post_mission_delta_v")

        constraintgrid.AddMany([self.lblDLA_bounds, DLAbox,
                                self.lblglobal_timebounded, self.chkglobal_timebounded,
                                self.lbltotal_flight_time_bounds, GlobalTimebox,
                                self.lblforced_post_launch_coast, self.txtforced_post_launch_coast,
                                self.lblforced_flyby_coast, self.txtforced_flyby_coast,
                                self.lblinitial_V_infinity, initial_V_infinity_box,
                                self.lblminimum_dry_mass, self.txtminimum_dry_mass,
                                self.lblpost_mission_delta_v, self.txtpost_mission_delta_v])

        vboxleft = wx.BoxSizer(wx.VERTICAL)
        vboxright = wx.BoxSizer(wx.VERTICAL)
        lblLeftTitle = wx.StaticText(self, -1, "Global mission options")
        lblRightTitle = wx.StaticText(self, -1, "Global mission constraints")
        vboxleft.Add(lblLeftTitle)
        vboxleft.Add(globaloptionsgrid)
        vboxright.Add(lblRightTitle)
        vboxright.Add(constraintgrid)
        
        font = self.GetFont()
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        lblLeftTitle.SetFont(font)
        lblRightTitle.SetFont(font)

        self.mainbox = wx.BoxSizer(wx.HORIZONTAL)
        
        self.mainbox.Add(vboxleft)
        self.mainbox.AddSpacer(20)
        self.mainbox.Add(vboxright)

        self.SetSizer(self.mainbox)
        self.SetupScrolling()


class SpacecraftOptionsPanel(wx.lib.scrolledpanel.ScrolledPanel):
    def __init__(self, parent):
        
        wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent)
        
        #spacecraft and launch vehicle fields
        spacecraftgrid = wx.FlexGridSizer(20,2,5,5)
        spacecraftgridtitle = wx.StaticText(self, -1, "Spacecraft and Launch Vehicle options")

        self.lblmaximum_mass = wx.StaticText(self, -1, "Maximum mass")
        self.txtmaximum_mass = wx.TextCtrl(self, -1, "maximum_mass")
        
        self.lblallow_initial_mass_to_vary = wx.StaticText(self, -1, "Allow initial mass to vary")
        self.chkallow_initial_mass_to_vary = wx.CheckBox(self, -1)

        self.lblEP_dry_mass = wx.StaticText(self, -1, "Propulsion stage dry mass")
        self.txtEP_dry_mass = wx.TextCtrl(self, -1, "EP_dry_mass")
        self.lblEP_dry_mass.Show(False)
        self.txtEP_dry_mass.Show(False)

        self.lblLV_type = wx.StaticText(self, -1, "Launch vehicle type")
        LV_choices = ['-2: custom launch vehicle','-1: burn with departure stage engine','0: fixed initial mass',
                      '1: Atlas V (401)  NLSII','2: Atlas V (411)  NLSII','3: Atlas V (421)  NLSII',
                      '4: Atlas V (431)  NLSII','5: Atlas V (501)  NLSII','6: Atlas V (511)  NLSII',
                      '7: Atlas V (521)  NLSII','8: Atlas V (531)  NLSII','9: Atlas V (541)  NLSII',
                      '10: Atlas V (551)  NLSII','11: Falcon 9 (v1.0) NLSII','12: Falcon 9 (v1.1)  NLSII',
                      '13: Atlas V (551) w/Star 48  NLSI','14: Falcon 9 Heavy  (notional)','15: Delta IV Heavy  NLSI',
                      '16: SLS Block 1  (notional)']
        self.cmbLV_type = wx.ComboBox(self, -1, choices=LV_choices, style=wx.CB_READONLY)

        self.lblIspDS = wx.StaticText(self, -1, "Departure stage Isp (s)")
        self.txtIspDS = wx.TextCtrl(self, -1, "IspDS")

        self.lblcustom_LV_coefficients = wx.StaticText(self, -1, "Custom launch vehicle coefficients")
        self.txtcustom_LV_coefficients0 = wx.TextCtrl(self, -1, "custom_LV_coefficients[0]")
        self.txtcustom_LV_coefficients1 = wx.TextCtrl(self, -1, "custom_LV_coefficients[1]")
        self.txtcustom_LV_coefficients2 = wx.TextCtrl(self, -1, "custom_LV_coefficients[2]")
        self.txtcustom_LV_coefficients3 = wx.TextCtrl(self, -1, "custom_LV_coefficients[3]")
        self.txtcustom_LV_coefficients4 = wx.TextCtrl(self, -1, "custom_LV_coefficients[4]")
        self.txtcustom_LV_coefficients5 = wx.TextCtrl(self, -1, "custom_LV_coefficients[5]")
        LV_coefficients_box = wx.BoxSizer(wx.HORIZONTAL)
        LV_coefficients_box.AddMany([self.txtcustom_LV_coefficients0, self.txtcustom_LV_coefficients1, self.txtcustom_LV_coefficients2,
                                     self.txtcustom_LV_coefficients3, self.txtcustom_LV_coefficients4, self.txtcustom_LV_coefficients5])

        self.lblcustom_LV_C3_bounds = wx.StaticText(self, -1, "Custom launch vehicle C3 bounds (km^2/s^2)")
        self.txtcustom_LV_C3_bounds_lower = wx.TextCtrl(self, -1, "custom_LV_C3_bounds[0]")
        self.txtcustom_LV_C3_bounds_upper = wx.TextCtrl(self, -1, "custom_LV_C3_bounds[1]")
        custom_LV_C3_bounds_box = wx.BoxSizer(wx.HORIZONTAL)
        custom_LV_C3_bounds_box.AddMany([self.txtcustom_LV_C3_bounds_lower, self.txtcustom_LV_C3_bounds_upper])

        self.lblLV_adapter_mass = wx.StaticText(self, -1, "Launch vehicle adapter mass (kg)")
        self.txtLV_adapter_mass = wx.TextCtrl(self, -1, "LV_margin")

        self.lblparking_orbit_altitude = wx.StaticText(self, -1, "Parking orbit altitude (km)")
        self.txtparking_orbit_altitude = wx.TextCtrl(self, -1, "parking_orbit_altitude")

        self.lblparking_orbit_inclination = wx.StaticText(self, -1, "Parking orbit inclination (degrees)")
        self.txtparking_orbit_inclination = wx.TextCtrl(self, -1, "parking_orbit_inclination")

        spacecraftgrid.AddMany([self.lblmaximum_mass, self.txtmaximum_mass,
                                self.lblallow_initial_mass_to_vary, self.chkallow_initial_mass_to_vary,
                                self.lblEP_dry_mass, self.txtEP_dry_mass,
                                self.lblLV_type, self.cmbLV_type,
                                self.lblLV_adapter_mass, self.txtLV_adapter_mass,
                                self.lblIspDS, self.txtIspDS,
                                self.lblcustom_LV_coefficients, LV_coefficients_box,
                                self.lblcustom_LV_C3_bounds, custom_LV_C3_bounds_box,
                                self.lblparking_orbit_altitude, self.txtparking_orbit_altitude,
                                self.lblparking_orbit_inclination, self.txtparking_orbit_inclination])

        spacecraftbox = wx.BoxSizer(wx.VERTICAL)
        spacecraftbox.AddMany([spacecraftgridtitle, spacecraftgrid])


        #terminal constraint/margining fields
        constraintsgrid = wx.FlexGridSizer(12,2,5,5)
        constraintsgridtitle = wx.StaticText(self, -1, "Margins and Constraints")

        self.lblpost_mission_Isp = wx.StaticText(self, -1, "Isp for post-mission delta-v (s)")
        self.txtpost_mission_Isp = wx.TextCtrl(self, -1, "post_mission_Isp")

        self.lblpropellant_margin = wx.StaticText(self, -1, "Propellant margin (fraction)")
        self.txtpropellant_margin = wx.TextCtrl(self, -1, "propellant_margin")

        self.lblpower_margin = wx.StaticText(self, -1, "Power margin (fraction)")
        self.txtpower_margin = wx.TextCtrl(self, -1, "power_margin")

        self.lblLV_margin = wx.StaticText(self, -1, "Launch vehicle margin (fraction)")
        self.txtLV_margin = wx.TextCtrl(self, -1, "LV_margin")

        self.lblenable_maximum_propellant_constraint = wx.StaticText(self, -1, "Enable maximum propellant constraint?")
        self.chkenable_propellant_mass_constraint = wx.CheckBox(self, -1)

        self.lblmaximum_propellant_mass = wx.StaticText(self, -1, "Maximum propellant mass (kg)")
        self.txtmaximum_propellant_mass = wx.TextCtrl(self, -1, "maximum_propellant_mass")

        constraintsgrid.AddMany([self.lblpropellant_margin, self.txtpropellant_margin,
                                 self.lblpower_margin, self.txtpower_margin,
                                 self.lblLV_margin, self.txtLV_margin,
                                 self.lblenable_maximum_propellant_constraint, self.chkenable_propellant_mass_constraint,
                                 self.lblmaximum_propellant_mass, self.txtmaximum_propellant_mass,
                                 self.lblpost_mission_Isp, self.txtpost_mission_Isp])

        constraintsbox = wx.BoxSizer(wx.VERTICAL)
        constraintsbox.AddMany([constraintsgridtitle, constraintsgrid])

        #propulsion
        propulsiongrid = wx.FlexGridSizer(26,2,5,5)
        propulsiongridtitle = wx.StaticText(self, -1, "Propulsion options")

        self.lblIspChem = wx.StaticText(self, -1, "Chemical Isp (s)")
        self.txtIspChem = wx.TextCtrl(self, -1, "IspChem")

        self.lblengine_type = wx.StaticText(self, -1, "Engine type")
        enginetypes = ['0: fixed thrust/Isp','1: constant Isp, efficiency, EMTG computes input power','2: choice of power model, constant efficiency, EMTG chooses Isp',
                       '3: choice of power model, constant efficiency and Isp','4: continuously-varying specific impulse','5: custom thrust and mass flow rate polynomial',
                       '6: NSTAR','7: XIPS-25','8: BPT-4000 High-Isp','9: BPT-4000 High-Thrust','10: BPT-4000 Ex-High-Isp','11: NEXT high-Isp v9',
                       '12: VASIMR (argon, using analytical model, not available in open-source)','13: Hall Thruster (Xenon, using analytical model, not available in open-source)','14: NEXT high-ISP v10',
                       '15: NEXT high-thrust v10','16: BPT-4000 MALTO','17: NEXIS Cardiff 8-15-201','18: H6MS Cardiff 8-15-2013','19: BHT20K Cardiff 8-16-2013','20: HiVHAC EM','21: 13 kW STMD Hall high-Isp (not available in open-source)','22: 13 kW STMD Hall high-thrust (not available in open-source)']
        self.cmbengine_type = wx.ComboBox(self, -1, choices = enginetypes, style=wx.CB_READONLY)

        self.lblnumber_of_engines = wx.StaticText(self, -1, "Number of thrusters")
        self.txtnumber_of_engines = wx.TextCtrl(self, -1, "number_of_engines")

        self.lblthrottle_logic_mode = wx.StaticText(self, -1, "Throttle logic mode")
        throttle_logic_types = ['maximum power use','maximum thrust','maximum Isp','maximum efficiency']
        self.cmbthrottle_logic_mode = wx.ComboBox(self, -1, choices = throttle_logic_types, style = wx.CB_READONLY)

        self.lblthrottle_sharpness = wx.StaticText(self, -1, "Throttle sharpness")
        self.txtthrottle_sharpness = wx.TextCtrl(self, -1, "throttle_sharpness")

        self.lblengine_duty_cycle = wx.StaticText(self, -1, "Thruster duty cycle")
        self.txtengine_duty_cycle = wx.TextCtrl(self, -1, "engine_duty_cycle")

        self.lblThrust = wx.StaticText(self, -1, "Electric thruster thrust (N)")
        self.txtThrust = wx.TextCtrl(self, -1, "Thrust")

        self.lblIspLT = wx.StaticText(self, -1, "Electric thruster Isp (s)")
        self.txtIspLT = wx.TextCtrl(self, -1, "IspLT")

        self.lblIspLT_minimum = wx.StaticText(self, -1, "Minimum Isp for VSI systems (s)")
        self.txtIspLT_minimum = wx.TextCtrl(self, -1, "IspLT_minimum")

        self.lbluser_defined_engine_efficiency = wx.StaticText(self, -1, "Thruster efficiency")
        self.txtuser_defined_engine_efficiency = wx.TextCtrl(self, -1, "user_defined_engine_efficiency")

        self.lblengine_input_thrust_coefficients = wx.StaticText(self, -1, "Custom thrust coefficients")
        self.txtengine_input_thrust_coefficients0 = wx.TextCtrl(self, -1, "engine_input_thrust_coefficients[0]")
        self.txtengine_input_thrust_coefficients1 = wx.TextCtrl(self, -1, "engine_input_thrust_coefficients[1]")
        self.txtengine_input_thrust_coefficients2 = wx.TextCtrl(self, -1, "engine_input_thrust_coefficients[2]")
        self.txtengine_input_thrust_coefficients3 = wx.TextCtrl(self, -1, "engine_input_thrust_coefficients[3]")
        self.txtengine_input_thrust_coefficients4 = wx.TextCtrl(self, -1, "engine_input_thrust_coefficients[4]")
        self.txtengine_input_thrust_coefficients5 = wx.TextCtrl(self, -1, "engine_input_thrust_coefficients[5]")
        self.txtengine_input_thrust_coefficients6 = wx.TextCtrl(self, -1, "engine_input_thrust_coefficients[6]")
        thrust_coefficients_box = wx.BoxSizer(wx.HORIZONTAL)
        thrust_coefficients_box.AddMany([self.txtengine_input_thrust_coefficients0, self.txtengine_input_thrust_coefficients1, self.txtengine_input_thrust_coefficients2,
                                         self.txtengine_input_thrust_coefficients3, self.txtengine_input_thrust_coefficients4, self.txtengine_input_thrust_coefficients5,
                                         self.txtengine_input_thrust_coefficients6])

        self.lblengine_input_mass_flow_rate_coefficients = wx.StaticText(self, -1, "Custom mass flow rate coefficients")
        self.txtengine_input_mass_flow_rate_coefficients0 = wx.TextCtrl(self, -1, "engine_input_mass_flow_rate_coefficients[0]")
        self.txtengine_input_mass_flow_rate_coefficients1 = wx.TextCtrl(self, -1, "engine_input_mass_flow_rate_coefficients[1]")
        self.txtengine_input_mass_flow_rate_coefficients2 = wx.TextCtrl(self, -1, "engine_input_mass_flow_rate_coefficients[2]")
        self.txtengine_input_mass_flow_rate_coefficients3 = wx.TextCtrl(self, -1, "engine_input_mass_flow_rate_coefficients[3]")
        self.txtengine_input_mass_flow_rate_coefficients4 = wx.TextCtrl(self, -1, "engine_input_mass_flow_rate_coefficients[4]")
        self.txtengine_input_mass_flow_rate_coefficients5 = wx.TextCtrl(self, -1, "engine_input_mass_flow_rate_coefficients[5]")
        self.txtengine_input_mass_flow_rate_coefficients6 = wx.TextCtrl(self, -1, "engine_input_mass_flow_rate_coefficients[6]")
        mass_flow_rate_coefficients_box = wx.BoxSizer(wx.HORIZONTAL)
        mass_flow_rate_coefficients_box.AddMany([self.txtengine_input_mass_flow_rate_coefficients0, self.txtengine_input_mass_flow_rate_coefficients1, self.txtengine_input_mass_flow_rate_coefficients2,
                                                 self.txtengine_input_mass_flow_rate_coefficients3, self.txtengine_input_mass_flow_rate_coefficients4, self.txtengine_input_mass_flow_rate_coefficients5,
                                                 self.txtengine_input_mass_flow_rate_coefficients6])

        self.lblengine_input_power_bounds = wx.StaticText(self, -1, "Thruster input power bounds")
        self.txtengine_input_power_bounds_lower = wx.TextCtrl(self, -1, "engine_input_power_bounds[0]")
        self.txtengine_input_power_bounds_upper = wx.TextCtrl(self, -1, "engine_input_power_bounds[1]")
        enginepowerbox = wx.BoxSizer(wx.HORIZONTAL)
        enginepowerbox.AddMany([self.txtengine_input_power_bounds_lower, self.txtengine_input_power_bounds_upper])

        propulsiongrid.AddMany([self.lblIspChem, self.txtIspChem,
                                self.lblengine_type, self.cmbengine_type,
                                self.lblnumber_of_engines, self.txtnumber_of_engines,
                                self.lblthrottle_logic_mode, self.cmbthrottle_logic_mode,
                                self.lblthrottle_sharpness, self.txtthrottle_sharpness,
                                self.lblengine_duty_cycle, self.txtengine_duty_cycle,
                                self.lblThrust, self.txtThrust,
                                self.lblIspLT, self.txtIspLT,
                                self.lblIspLT_minimum, self.txtIspLT_minimum,
                                self.lbluser_defined_engine_efficiency, self.txtuser_defined_engine_efficiency,
                                self.lblengine_input_thrust_coefficients, thrust_coefficients_box,
                                self.lblengine_input_mass_flow_rate_coefficients, mass_flow_rate_coefficients_box,
                                self.lblengine_input_power_bounds, enginepowerbox])

        propulsionbox = wx.BoxSizer(wx.VERTICAL)
        propulsionbox.AddMany([propulsiongridtitle, propulsiongrid])

        #power
        powergrid = wx.FlexGridSizer(20,2,5,5)
        self.powergridtitle = wx.StaticText(self, -1, "Power options")

        self.lblpower_at_1_AU = wx.StaticText(self, -1, "Power at 1 AU (kW)")
        self.txtpower_at_1_AU = wx.TextCtrl(self, -1, "power_at_1_AU")

        self.lblpower_source_type = wx.StaticText(self, -1, "Power source type")
        power_source_choices = ['0: solar','1: radioisotope']
        self.cmbpower_source_type = wx.ComboBox(self, -1, choices=power_source_choices, style=wx.CB_READONLY)

        self.lblsolar_power_gamma = wx.StaticText(self, -1, "Solar power coefficients")
        self.txtsolar_power_gamma0 = wx.TextCtrl(self, -1, "solar_power_gamma[0]")
        self.txtsolar_power_gamma1 = wx.TextCtrl(self, -1, "solar_power_gamma[1]")
        self.txtsolar_power_gamma2 = wx.TextCtrl(self, -1, "solar_power_gamma[2]")
        self.txtsolar_power_gamma3 = wx.TextCtrl(self, -1, "solar_power_gamma[3]")
        self.txtsolar_power_gamma4 = wx.TextCtrl(self, -1, "solar_power_gamma[4]")
        solarpowerbox = wx.BoxSizer(wx.HORIZONTAL)
        solarpowerbox.AddMany([self.txtsolar_power_gamma0, self.txtsolar_power_gamma1, self.txtsolar_power_gamma2, self.txtsolar_power_gamma3, self.txtsolar_power_gamma4])

        self.lblspacecraft_power_model_type = wx.StaticText(self, -1, "Spacecraft power model type")
        power_model_choices = ['0: P_sc = A + B/r + C/r^2','1: P_sc = A if P > A, A + B(C - P) otherwise']
        self.cmbspacecraft_power_model_type = wx.ComboBox(self, -1, choices=power_model_choices, style = wx.CB_READONLY)
        
        self.lblspacecraft_power_coefficients = wx.StaticText(self, -1, "Spacecraft power coefficients")
        self.txtspacecraft_power_coefficients0 = wx.TextCtrl(self, -1, "spacecraft_power_coefficients[0]")
        self.txtspacecraft_power_coefficients1 = wx.TextCtrl(self, -1, "spacecraft_power_coefficients[1]")
        self.txtspacecraft_power_coefficients2 = wx.TextCtrl(self, -1, "spacecraft_power_coefficients[2]")
        spacecraftpowerbox = wx.BoxSizer(wx.HORIZONTAL)
        spacecraftpowerbox.AddMany([self.txtspacecraft_power_coefficients0, self.txtspacecraft_power_coefficients1, self.txtspacecraft_power_coefficients2])

        self.lblpower_decay_rate = wx.StaticText(self, -1, "Power decay rate (fraction per year)")
        self.txtpower_decay_rate = wx.TextCtrl(self, -1, "power_decay_rate")

        powergrid.AddMany([self.lblpower_at_1_AU, self.txtpower_at_1_AU,
                           self.lblpower_source_type, self.cmbpower_source_type,
                           self.lblsolar_power_gamma, solarpowerbox,
                           self.lblspacecraft_power_model_type, self.cmbspacecraft_power_model_type,
                           self.lblspacecraft_power_coefficients, spacecraftpowerbox,
                           self.lblpower_decay_rate, self.txtpower_decay_rate])

        powerbox = wx.BoxSizer(wx.VERTICAL)
        powerbox.AddMany([self.powergridtitle, powergrid])

        #now tie everything together
        leftvertsizer = wx.BoxSizer(wx.VERTICAL)
        leftvertsizer.AddMany([spacecraftbox, propulsionbox, powerbox])
        
        self.mainbox = wx.BoxSizer(wx.HORIZONTAL)
        self.mainbox.AddMany([leftvertsizer, constraintsbox]) 

        font = self.GetFont()
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        spacecraftgridtitle.SetFont(font)
        self.powergridtitle.SetFont(font)
        constraintsgridtitle.SetFont(font)
        propulsiongridtitle.SetFont(font)

        self.SetSizer(self.mainbox)
        self.SetupScrolling()


class JourneyOptionsPanel(wx.lib.scrolledpanel.ScrolledPanel):
    def __init__(self, parent):
        
        wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent)

        self.JourneyList = []

        self.JourneySelectBox = wx.ListBox(self, -1, choices=self.JourneyList, size=(300,200), style=wx.LB_SINGLE)
        self.btnAddNewJourney = wx.Button(self, -1, "New Journey", size=(200,-1))
        self.btnDeleteJourney = wx.Button(self, -1, "Delete Journey", size=(200,-1))
        self.btnMoveJourneyUp = wx.Button(self, -1, "Move Journey Up", size=(200,-1))
        self.btnMoveJourneyDown = wx.Button(self, -1, "Move Journey Down", size=(200,-1))

        buttonstacksizer = wx.BoxSizer(wx.VERTICAL)
        buttonstacksizer.AddMany([self.btnAddNewJourney, self.btnDeleteJourney, self.btnMoveJourneyUp, self.btnMoveJourneyDown])

        JourneySelectionSizer = wx.BoxSizer(wx.HORIZONTAL)
        JourneySelectionSizer.Add(self.JourneySelectBox)
        JourneySelectionSizer.AddSpacer(5)
        JourneySelectionSizer.Add(buttonstacksizer)

        self.lbljourney_names = wx.StaticText(self, -1, "Journey name")
        self.txtjourney_names = wx.TextCtrl(self, -1, "journey_names", size=(300,-1))

        self.lbljourney_central_body = wx.StaticText(self, -1, "Central body")
        self.txtjourney_central_body = wx.TextCtrl(self, -1, "journey_central_body")
        self.btnjourney_central_body = wx.Button(self, -1, "...")
        journey_central_body_box = wx.BoxSizer(wx.HORIZONTAL)
        journey_central_body_box.Add(self.txtjourney_central_body)
        journey_central_body_box.AddSpacer(5)
        journey_central_body_box.Add(self.btnjourney_central_body)

        self.lbldestination_list = wx.StaticText(self, -1, "Destination list")
        self.txtdestination_list = wx.TextCtrl(self, -1, "destination_list")
        self.btndestination_list = wx.Button(self, -1, "...")
        destination_list_box = wx.BoxSizer(wx.HORIZONTAL)
        destination_list_box.Add(self.txtdestination_list)
        destination_list_box.AddSpacer(5)
        destination_list_box.Add(self.btndestination_list)

        self.lbljourney_starting_mass_increment = wx.StaticText(self, -1, "Starting mass increment (kg)")
        self.txtjourney_starting_mass_increment = wx.TextCtrl(self, -1, "journey_starting_mass_increment")

        self.lbljourney_variable_mass_increment = wx.StaticText(self, -1, "Variable mass increment")
        self.chkjourney_variable_mass_increment = wx.CheckBox(self, -1)

        self.lbljourney_wait_time_bounds = wx.StaticText(self, -1, "Wait time bounds")
        self.txtjourney_wait_time_bounds_lower = wx.TextCtrl(self, -1, "journey_wait_time_bounds[0]")
        self.txtjourney_wait_time_bounds_upper = wx.TextCtrl(self, -1, "journey_wait_time_bounds[1]")
        wait_time_sizer = wx.BoxSizer(wx.HORIZONTAL)
        wait_time_sizer.AddMany([self.txtjourney_wait_time_bounds_lower, self.txtjourney_wait_time_bounds_upper])

        self.lbljourney_timebounded = wx.StaticText(self, -1, "Journey time bounds")
        journey_time_bounds_choices = ['unbounded','bounded flight time','bounded arrival date','bounded aggregate flight time']
        self.cmbjourney_timebounded = wx.ComboBox(self, -1, choices=journey_time_bounds_choices, style=wx.CB_READONLY)

        self.lbljourney_flight_time_bounds = wx.StaticText(self, -1, "Journey flight time bounds")
        self.txtjourney_flight_time_bounds_lower = wx.TextCtrl(self, -1, "journey_flight_time_bounds[0]")
        self.txtjourney_flight_time_bounds_upper = wx.TextCtrl(self, -1, "journey_flight_time_bounds[1]")
        flight_time_sizer = wx.BoxSizer(wx.HORIZONTAL)
        flight_time_sizer.AddMany([self.txtjourney_flight_time_bounds_lower, self.txtjourney_flight_time_bounds_upper])

        self.lbljourney_arrival_date_bounds = wx.StaticText(self, -1, "Journey arrival date bounds")
        self.txtjourney_arrival_date_bounds_lower = wx.TextCtrl(self, -1, "journey_arrival_date_bounds[0]")
        self.txtjourney_arrival_date_bounds_upper = wx.TextCtrl(self, -1, "journey_arrival_date_bounds[1]")
        self.ArrivalDateLowerCalendar = wx.calendar.CalendarCtrl(self, -1)
        self.ArrivalDateUpperCalendar = wx.calendar.CalendarCtrl(self, -1)
        arrival_date_sizer = wx.BoxSizer(wx.HORIZONTAL)
        arrival_date_sizer.AddMany([self.txtjourney_arrival_date_bounds_lower, self.ArrivalDateLowerCalendar, self.txtjourney_arrival_date_bounds_upper, self.ArrivalDateUpperCalendar])

        self.lbljourney_initial_impulse_bounds = wx.StaticText(self, -1, "Journey initial impulse bounds")
        self.txtjourney_initial_impulse_bounds_lower = wx.TextCtrl(self, -1, "journey_initial_impulse_bounds[0]")
        self.txtjourney_initial_impulse_bounds_upper = wx.TextCtrl(self, -1, "journey_initial_impulse_bounds[1]")
        initial_impulse_sizer = wx.BoxSizer(wx.HORIZONTAL)
        initial_impulse_sizer.AddMany([self.txtjourney_initial_impulse_bounds_lower, self.txtjourney_initial_impulse_bounds_upper])

        self.lbljourney_departure_type = wx.StaticText(self, -1, "Journey departure type")
        journey_departure_type_choices = ['0: launch or direct insertion','1: depart from parking orbit','2: free direct departure',
                                        '3: flyby','4: flyby with fixed v-infinity-out','5: Spiral-out from circular orbit','6: zero-turn flyby (for small bodies)']
        self.cmbjourney_departure_type = wx.ComboBox(self, -1, choices=journey_departure_type_choices, style=wx.CB_READONLY)

        self.lbljourney_initial_velocity = wx.StaticText(self, -1, "Journey initial velocity vector")
        self.txtjourney_initial_velocity0 = wx.TextCtrl(self, -1, "journey_initial_velocity[0]")
        self.txtjourney_initial_velocity1 = wx.TextCtrl(self, -1, "journey_initial_velocity[1]")
        self.txtjourney_initial_velocity2 = wx.TextCtrl(self, -1, "journey_initial_velocity[2]")
        journey_initial_velocity_box = wx.BoxSizer(wx.HORIZONTAL)
        journey_initial_velocity_box.AddMany([self.txtjourney_initial_velocity0, self.txtjourney_initial_velocity1, self.txtjourney_initial_velocity2])

        self.lbljourney_escape_spiral_starting_radius = wx.StaticText(self, -1, "Orbital radius for beginning of escape spiral (km)")
        self.txtjourney_escape_spiral_starting_radius = wx.TextCtrl(self, -1, "journey_escape_spiral_starting_radius")

        self.lbljourney_arrival_type = wx.StaticText(self, -1, "Journey arrival type")
        journey_arrival_type_choices = ['0: insertion into parking orbit (use chemical Isp)','1: rendezvous (use chemical Isp)','2: intercept with bounded V_infinity',
                                        '3: low-thrust rendezvous (does not work if terminal phase is not low-thrust)','4: match final v-infinity vector',
                                        '5: match final v-infinity vector (low-thrust)','6: escape (E = 0)','7: capture spiral']
        self.cmbjourney_arrival_type = wx.ComboBox(self, -1, choices=journey_arrival_type_choices, style=wx.CB_READONLY)

        self.lbljourney_capture_spiral_final_radius = wx.StaticText(self, -1, "Orbital radius for end of capture spiral (km)")
        self.txtjourney_capture_spiral_final_radius = wx.TextCtrl(self, -1, "journey_capture_spiral_final_radius")

        self.lbljourney_final_velocity = wx.StaticText(self, -1, "Journey final velocity vector")
        self.txtjourney_final_velocity0 = wx.TextCtrl(self, -1, "journey_final_velocity[0]")
        self.txtjourney_final_velocity1 = wx.TextCtrl(self, -1, "journey_final_velocity[1]")
        self.txtjourney_final_velocity2 = wx.TextCtrl(self, -1, "journey_final_velocity[2]")
        journey_final_velocity_box = wx.BoxSizer(wx.HORIZONTAL)
        journey_final_velocity_box.AddMany([self.txtjourney_final_velocity0, self.txtjourney_final_velocity1, self.txtjourney_final_velocity2])

        self.lbljourney_arrival_declination_constraint_flag = wx.StaticText(self, -1, "Apply arrival declination constraint?")
        self.chkjourney_arrival_declination_constraint_flag = wx.CheckBox(self, -1)
        self.lbljourney_arrival_declination_bounds = wx.StaticText(self, -1, "Arrival Declination bounds")
        self.txtjourney_arrival_declination_bounds_lower = wx.TextCtrl(self, -1, "journey_arrival_declination_bounds[0]")
        self.txtjourney_arrival_declination_bounds_upper = wx.TextCtrl(self, -1, "journey_arrival_declination_bounds[1]")
        declination_bounds_box = wx.BoxSizer(wx.HORIZONTAL)
        declination_bounds_box.AddMany([self.txtjourney_arrival_declination_bounds_lower, self.txtjourney_arrival_declination_bounds_upper])

        self.lblsequence = wx.StaticText(self, -1, "Flyby sequence")
        self.txtsequence = wx.TextCtrl(self, -1, "sequence", size=(300,60), style=wx.TE_MULTILINE)
        self.btnsequence = wx.Button(self, -1, "...")
        sequence_box = wx.BoxSizer(wx.HORIZONTAL)
        sequence_box.Add(self.txtsequence)
        sequence_box.AddSpacer(5)
        sequence_box.Add(self.btnsequence)

        self.lbljourney_perturbation_bodies = wx.StaticText(self, -1, "Perturbation_bodies")
        self.txtjourney_perturbation_bodies = wx.TextCtrl(self, -1, "journey_perturbation_bodies", size=(300,-1))
        self.btnjourney_perturbation_bodies = wx.Button(self, -1, "...")
        journey_perturbation_bodies_box = wx.BoxSizer(wx.HORIZONTAL)
        journey_perturbation_bodies_box.Add(self.txtjourney_perturbation_bodies)
        journey_perturbation_bodies_box.AddSpacer(5)
        journey_perturbation_bodies_box.Add(self.btnjourney_perturbation_bodies)

        

        JourneyInformationGrid = wx.FlexGridSizer(40,2,5,5)
        JourneyInformationGrid.AddMany([self.lbljourney_names, self.txtjourney_names,
                                        self.lbljourney_central_body, journey_central_body_box,
                                        self.lbldestination_list, destination_list_box,
                                        self.lbljourney_starting_mass_increment, self.txtjourney_starting_mass_increment,
                                        self.lbljourney_variable_mass_increment, self.chkjourney_variable_mass_increment,
                                        self.lbljourney_wait_time_bounds, wait_time_sizer,
                                        self.lbljourney_timebounded, self.cmbjourney_timebounded,
                                        self.lbljourney_flight_time_bounds, flight_time_sizer,
                                        self.lbljourney_arrival_date_bounds, arrival_date_sizer,
                                        self.lbljourney_initial_impulse_bounds, initial_impulse_sizer,
                                        self.lbljourney_departure_type, self.cmbjourney_departure_type,
                                        self.lbljourney_escape_spiral_starting_radius, self.txtjourney_escape_spiral_starting_radius,
                                        self.lbljourney_initial_velocity, journey_initial_velocity_box,
                                        self.lbljourney_arrival_type, self.cmbjourney_arrival_type,
                                        self.lbljourney_capture_spiral_final_radius, self.txtjourney_capture_spiral_final_radius,
                                        self.lbljourney_final_velocity, journey_final_velocity_box,
                                        self.lbljourney_arrival_declination_constraint_flag, self.chkjourney_arrival_declination_constraint_flag,
                                        self.lbljourney_arrival_declination_bounds, declination_bounds_box,
                                        self.lblsequence, sequence_box,
                                        self.lbljourney_perturbation_bodies, journey_perturbation_bodies_box])

        JourneyInformationStacker = wx.BoxSizer(wx.VERTICAL)
        JourneyInformationStacker.AddMany([JourneySelectionSizer, JourneyInformationGrid])

        #custom departure elements
        self.boxjourney_departure_elements = wx.StaticBox(self, -1, "Journey departure elements")
        self.lbljourney_departure_elements_type= wx.StaticText(self, -1, "Journey departure elements type")
        journey_departure_elements_type_choices = ['0: inertial', '1: COE']
        self.cmbjourney_departure_elements_type = wx.ComboBox(self, -1, choices=journey_departure_elements_type_choices, style=wx.CB_READONLY)
        departure_elements_type_box = wx.BoxSizer(wx.HORIZONTAL)
        departure_elements_type_box.Add(self.lbljourney_departure_elements_type)
        departure_elements_type_box.AddSpacer(5)
        departure_elements_type_box.Add(self.cmbjourney_departure_elements_type)

        empty_departure_cell = wx.StaticText(self, -1, "")
        self.lblvarydepartureelements = wx.StaticText(self, -1, "Vary?")
        self.lbldepartureelementsvalue = wx.StaticText(self, -1, "Value")
        self.lbldepartureelementslower = wx.StaticText(self, -1, "Lower bound")
        self.lbldepartureelementsupper = wx.StaticText(self, -1, "Upper bound")
        self.lblSMA_departure = wx.StaticText(self, -1, "SMA")
        self.lblECC_departure = wx.StaticText(self, -1, "ECC")
        self.lblINC_departure = wx.StaticText(self, -1, "INC")
        self.lblRAAN_departure = wx.StaticText(self, -1, "RAAN")
        self.lblAOP_departure = wx.StaticText(self, -1, "AOP")
        self.lblMA_departure = wx.StaticText(self, -1, "MA")
        self.chkSMA_departure = wx.CheckBox(self, -1)
        self.chkECC_departure = wx.CheckBox(self, -1)
        self.chkINC_departure = wx.CheckBox(self, -1)
        self.chkRAAN_departure = wx.CheckBox(self, -1)
        self.chkAOP_departure = wx.CheckBox(self, -1)
        self.chkMA_departure = wx.CheckBox(self, -1)
        self.txtSMA_departure = wx.TextCtrl(self, -1, "SMA_val")
        self.txtECC_departure = wx.TextCtrl(self, -1, "ECC_val")
        self.txtINC_departure = wx.TextCtrl(self, -1, "INC_val")
        self.txtRAAN_departure = wx.TextCtrl(self, -1, "RAAN_val")
        self.txtAOP_departure = wx.TextCtrl(self, -1, "AOP_val")
        self.txtMA_departure = wx.TextCtrl(self, -1, "MA_val")
        self.txtSMA_departure0 = wx.TextCtrl(self, -1, "SMA_val0")
        self.txtECC_departure0 = wx.TextCtrl(self, -1, "ECC_val0")
        self.txtINC_departure0 = wx.TextCtrl(self, -1, "INC_val0")
        self.txtRAAN_departure0 = wx.TextCtrl(self, -1, "RAAN_val0")
        self.txtAOP_departure0 = wx.TextCtrl(self, -1, "AOP_val0")
        self.txtMA_departure0 = wx.TextCtrl(self, -1, "MA_val0")
        self.txtSMA_departure1 = wx.TextCtrl(self, -1, "SMA_val1")
        self.txtECC_departure1 = wx.TextCtrl(self, -1, "ECC_val1")
        self.txtINC_departure1 = wx.TextCtrl(self, -1, "INC_val1")
        self.txtRAAN_departure1 = wx.TextCtrl(self, -1, "RAAN_val1")
        self.txtAOP_departure1 = wx.TextCtrl(self, -1, "AOP_val1")
        self.txtMA_departure1 = wx.TextCtrl(self, -1, "MA_val1")
        DepartureElementsSizer = wx.FlexGridSizer(14,5,5,5)
        DepartureElementsSizer.AddMany([empty_departure_cell, self.lblvarydepartureelements, self.lbldepartureelementsvalue, self.lbldepartureelementslower, self.lbldepartureelementsupper, 
                                            self.lblSMA_departure, self.chkSMA_departure, self.txtSMA_departure, self.txtSMA_departure0, self.txtSMA_departure1,
                                            self.lblECC_departure, self.chkECC_departure, self.txtECC_departure, self.txtECC_departure0, self.txtECC_departure1,
                                            self.lblINC_departure, self.chkINC_departure, self.txtINC_departure, self.txtINC_departure0, self.txtINC_departure1,
                                            self.lblRAAN_departure, self.chkRAAN_departure, self.txtRAAN_departure, self.txtRAAN_departure0, self.txtRAAN_departure1,
                                            self.lblAOP_departure, self.chkAOP_departure, self.txtAOP_departure, self.txtAOP_departure0, self.txtAOP_departure1,
                                            self.lblMA_departure, self.chkMA_departure, self.txtMA_departure, self.txtMA_departure0, self.txtMA_departure1])
        self.DepartureElementsBox = wx.StaticBoxSizer(self.boxjourney_departure_elements, wx.VERTICAL)
        self.DepartureElementsBox.AddMany([departure_elements_type_box, DepartureElementsSizer])
        
        font = self.GetFont()
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        self.boxjourney_departure_elements.SetFont(font)


        #custom arrival elements
        self.boxjourney_arrival_elements = wx.StaticBox(self, -1, "Journey arrival elements")
        self.lbljourney_arrival_elements_type= wx.StaticText(self, -1, "Journey arrival elements type")
        journey_arrival_elements_type_choices = ['0: inertial', '1: COE']
        self.cmbjourney_arrival_elements_type = wx.ComboBox(self, -1, choices=journey_arrival_elements_type_choices, style=wx.CB_READONLY)
        arrival_elements_type_box = wx.BoxSizer(wx.HORIZONTAL)
        arrival_elements_type_box.Add(self.lbljourney_arrival_elements_type)
        arrival_elements_type_box.AddSpacer(5)
        arrival_elements_type_box.Add(self.cmbjourney_arrival_elements_type)

        empty_arrival_cell = wx.StaticText(self, -1, "")
        self.lblvaryarrivalelements = wx.StaticText(self, -1, "Vary?")
        self.lblarrivalelementsvalue = wx.StaticText(self, -1, "Value")
        self.lblarrivalelementslower = wx.StaticText(self, -1, "Lower bound")
        self.lblarrivalelementsupper = wx.StaticText(self, -1, "Upper bound")
        self.lblSMA_arrival = wx.StaticText(self, -1, "SMA")
        self.lblECC_arrival = wx.StaticText(self, -1, "ECC")
        self.lblINC_arrival = wx.StaticText(self, -1, "INC")
        self.lblRAAN_arrival = wx.StaticText(self, -1, "RAAN")
        self.lblAOP_arrival = wx.StaticText(self, -1, "AOP")
        self.lblMA_arrival = wx.StaticText(self, -1, "MA")
        self.chkSMA_arrival = wx.CheckBox(self, -1)
        self.chkECC_arrival = wx.CheckBox(self, -1)
        self.chkINC_arrival = wx.CheckBox(self, -1)
        self.chkRAAN_arrival = wx.CheckBox(self, -1)
        self.chkAOP_arrival = wx.CheckBox(self, -1)
        self.chkMA_arrival = wx.CheckBox(self, -1)
        self.txtSMA_arrival = wx.TextCtrl(self, -1, "SMA_val")
        self.txtECC_arrival = wx.TextCtrl(self, -1, "ECC_val")
        self.txtINC_arrival = wx.TextCtrl(self, -1, "INC_val")
        self.txtRAAN_arrival = wx.TextCtrl(self, -1, "RAAN_val")
        self.txtAOP_arrival = wx.TextCtrl(self, -1, "AOP_val")
        self.txtMA_arrival = wx.TextCtrl(self, -1, "MA_val")
        self.txtSMA_arrival0 = wx.TextCtrl(self, -1, "SMA_val0")
        self.txtECC_arrival0 = wx.TextCtrl(self, -1, "ECC_val0")
        self.txtINC_arrival0 = wx.TextCtrl(self, -1, "INC_val0")
        self.txtRAAN_arrival0 = wx.TextCtrl(self, -1, "RAAN_val0")
        self.txtAOP_arrival0 = wx.TextCtrl(self, -1, "AOP_val0")
        self.txtMA_arrival0 = wx.TextCtrl(self, -1, "MA_val0")
        self.txtSMA_arrival1 = wx.TextCtrl(self, -1, "SMA_val1")
        self.txtECC_arrival1 = wx.TextCtrl(self, -1, "ECC_val1")
        self.txtINC_arrival1 = wx.TextCtrl(self, -1, "INC_val1")
        self.txtRAAN_arrival1 = wx.TextCtrl(self, -1, "RAAN_val1")
        self.txtAOP_arrival1 = wx.TextCtrl(self, -1, "AOP_val1")
        self.txtMA_arrival1 = wx.TextCtrl(self, -1, "MA_val1")
        ArrivalElementsSizer = wx.FlexGridSizer(14,5,5,5)
        ArrivalElementsSizer.AddMany([empty_arrival_cell, self.lblvaryarrivalelements, self.lblarrivalelementsvalue, self.lblarrivalelementslower, self.lblarrivalelementsupper, 
                                            self.lblSMA_arrival, self.chkSMA_arrival, self.txtSMA_arrival, self.txtSMA_arrival0, self.txtSMA_arrival1,
                                            self.lblECC_arrival, self.chkECC_arrival, self.txtECC_arrival, self.txtECC_arrival0, self.txtECC_arrival1,
                                            self.lblINC_arrival, self.chkINC_arrival, self.txtINC_arrival, self.txtINC_arrival0, self.txtINC_arrival1,
                                            self.lblRAAN_arrival, self.chkRAAN_arrival, self.txtRAAN_arrival, self.txtRAAN_arrival0, self.txtRAAN_arrival1,
                                            self.lblAOP_arrival, self.chkAOP_arrival, self.txtAOP_arrival, self.txtAOP_arrival0, self.txtAOP_arrival1,
                                            self.lblMA_arrival, self.chkMA_arrival, self.txtMA_arrival, self.txtMA_arrival0, self.txtMA_arrival1])
        self.ArrivalElementsBox = wx.StaticBoxSizer(self.boxjourney_arrival_elements, wx.VERTICAL)
        self.ArrivalElementsBox.AddMany([arrival_elements_type_box, ArrivalElementsSizer])
        
        font = self.GetFont()
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        self.boxjourney_arrival_elements.SetFont(font)

        ElementsStacker = wx.BoxSizer(wx.VERTICAL)
        ElementsStacker.AddMany([self.DepartureElementsBox, self.ArrivalElementsBox])

        self.mainbox = wx.BoxSizer(wx.HORIZONTAL)
        self.mainbox.AddMany([JourneyInformationStacker, ElementsStacker])
        self.SetSizer(self.mainbox)
        self.SetupScrolling()
        
class SolverOptionsPanel(wx.lib.scrolledpanel.ScrolledPanel):
    def __init__(self, parent):

        wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent)
        
        innerloopgrid = wx.GridSizer(28,2,5,5)
        
        self.lblInnerLoopSolver = wx.StaticText(self, -1, "Inner-loop Solver Mode")
        innerloopsolvertypes = ['Evaluate trialX', 'Evaluate a batch of trialX vectors','Monotonic Basin Hopping',
                                'Adaptive Constrained Differential Evolution','SNOPT with initial guess']
        self.cmbInnerLoopSolver = wx.ComboBox(self, -1, choices = innerloopsolvertypes, style=wx.CB_READONLY)

        self.lblNLP_solver_type = wx.StaticText(self, -1, "NLP solver")
        NLP_solver_types = ['SNOPT','WORHP']
        self.cmbNLP_solver_type = wx.ComboBox(self, -1, choices = NLP_solver_types, style=wx.CB_READONLY)

        self.lblNLP_solver_mode = wx.StaticText(self, -1, "NLP solver mode")
        NLP_solver_modes = ['Feasible point','Optimize']
        self.cmbNLP_solver_mode = wx.ComboBox(self, -1, choices = NLP_solver_modes, style = wx.CB_READONLY)

        self.lblquiet_NLP = wx.StaticText(self, -1, "Quiet NLP solver?")
        self.chkquiet_NLP = wx.CheckBox(self, -1)

        self.lblquiet_MBH = wx.StaticText(self, -1, "Quiet MBH solver?")
        self.chkquiet_MBH = wx.CheckBox(self, -1)

        self.lblMBH_two_step = wx.StaticText(self, -1, "Two-step MBH?")
        self.chkMBH_two_step = wx.CheckBox(self, -1)

        self.lblFD_stepsize = wx.StaticText(self, -1, "Finite differencing step size")
        self.txtFD_stepsize = wx.TextCtrl(self, -1, "FD_stepsize")

        self.lblFD_stepsize_coarse = wx.StaticText(self, -1, "Finite differencing coarse step size")
        self.txtFD_stepsize_coarse = wx.TextCtrl(self, -1, "FD_stepsize_coarse")

        self.lblACE_feasible_point_finder = wx.StaticText(self, -1, "Enable ACE feasible point finder?")
        self.chkACE_feasible_point_finder = wx.CheckBox(self, -1)
        
        self.lblMBH_max_not_improve = wx.StaticText(self, -1, "MBH Impatience")
        self.txtMBH_max_not_improve = wx.TextCtrl(self, -1, "MBH_max_not_improve")
        
        self.lblMBH_max_trials = wx.StaticText(self, -1, "Maximum number of innerloop trials")
        self.txtMBH_max_trials = wx.TextCtrl(self, -1, "MBH_max_trials")
        
        self.lblMBH_max_run_time = wx.StaticText(self, -1, "Maximum run-time")
        self.txtMBH_max_run_time = wx.TextCtrl(self, -1, "MBH_max_run_time")
        
        self.lblMBH_max_step_size = wx.StaticText(self, -1, "MBH maximum perturbation size")
        self.txtMBH_max_step_size = wx.TextCtrl(self, -1, "MBH_max_step_size")

        self.lblMBH_hop_distribution = wx.StaticText(self, -1, "MBH hop probability distribution")
        hop_distribution_choices = ["Uniform","Cauchy","Pareto","Gaussian"]
        self.cmbMBH_hop_distribution = wx.ComboBox(self, -1, choices = hop_distribution_choices, style = wx.CB_READONLY)

        self.lblMBH_Pareto_alpha = wx.StaticText(self, -1, "MBH Pareto distribution alpha")
        self.txtMBH_Pareto_alpha = wx.TextCtrl(self, -1, "MBH_Pareto_alpha")

        self.lblMBH_time_hop_probability = wx.StaticText(self, -1, "Probability of MBH time hop")
        self.txtMBH_time_hop_probability = wx.TextCtrl(self, -1, "MBH_time_hop_probability")
        
        self.lblsnopt_feasibility_tolerance = wx.StaticText(self, -1, "Feasibility tolerance")
        self.txtsnopt_feasibility_tolerance = wx.TextCtrl(self, -1, "snopt_feasibility_tolerance")
        
        self.lblsnopt_major_iterations = wx.StaticText(self, -1, "SNOPT major iterations limit")
        self.txtsnopt_major_iterations = wx.TextCtrl(self, -1, "snopt_major_iterations")
        
        self.lblsnopt_max_run_time = wx.StaticText(self, -1, "SNOPT maximum run time")
        self.txtsnopt_max_run_time = wx.TextCtrl(self, -1, "snopt_max_run_time")
        
        self.lblderivative_type = wx.StaticText(self, -1, "Maximum number of innerloop trials")
        derivativechoices = ["Finite Differencing","Analytical flybys and objective function","Analytical all but time","All but current phase flight time derivatives","Fully analytical (experimental)"]
        self.cmbderivative_type = wx.ComboBox(self, -1, choices = derivativechoices, style = wx.CB_READONLY)
        
        self.lblcheck_derivatives = wx.StaticText(self, -1, "Check derivatives via finite differencing?")
        self.chkcheck_derivatives = wx.CheckBox(self, -1)
        
        self.lblseed_MBH = wx.StaticText(self, -1, "Seed MBH?")
        self.chkseed_MBH = wx.CheckBox(self, -1)

        self.lblinitial_guess_control_coordinate_system = wx.StaticText(self, -1, "Initial guess control coordinate system")
        control_coordinate_choices = ['Cartesian','Polar']
        self.cmbinitial_guess_control_coordinate_system = wx.ComboBox(self, -1, choices = control_coordinate_choices, style=wx.CB_READONLY)
        
        self.lblinterpolate_initial_guess = wx.StaticText(self, -1, "Interpolate initial guess?")
        self.chkinterpolate_initial_guess = wx.CheckBox(self, -1)
        
        self.lblinitial_guess_num_timesteps = wx.StaticText(self, -1, "Number of timesteps used to create initial guess")
        self.txtinitial_guess_num_timesteps = wx.TextCtrl(self, -1, "initial_guess_num_timesteps")
        
        self.lblinitial_guess_step_size_distribution = wx.StaticText(self, -1, "Initial guess step size distribution")
        initialguessdistributionchoices = ["Uniform","Gaussian","Cauchy"]
        self.cmbinitial_guess_step_size_distribution = wx.ComboBox(self, -1, choices = initialguessdistributionchoices, style=wx.CB_READONLY)

        self.lblinitial_guess_step_size_stdv_or_scale = wx.StaticText(self, -1, "Initial guess scale width/standard deviation")
        self.txtinitial_guess_step_size_stdv_or_scale = wx.TextCtrl(self, -1, "initial_guess_step_size_stdv_or_scale")

        self.lblMBH_zero_control_initial_guess = wx.StaticText(self, -1, "Zero-control initial guess")
        MBH_zero_control_initial_guess_options = ['do not use','zero-control for resets, random perturbations for hops','always use zero-control guess except when seeded']
        self.cmbMBH_zero_control_initial_guess = wx.ComboBox(self, -1, choices = MBH_zero_control_initial_guess_options, style=wx.CB_READONLY)
                
        innerloopgrid.AddMany(  [self.lblInnerLoopSolver, self.cmbInnerLoopSolver,
                                 self.lblNLP_solver_type, self.cmbNLP_solver_type,
                                 self.lblNLP_solver_mode, self.cmbNLP_solver_mode,
                                 self.lblquiet_NLP, self.chkquiet_NLP,
                                 self.lblquiet_MBH, self.chkquiet_MBH,
                                 self.lblMBH_two_step, self.chkMBH_two_step,
                                 self.lblFD_stepsize, self.txtFD_stepsize,
                                 self.lblFD_stepsize_coarse, self.txtFD_stepsize_coarse,
                                 self.lblACE_feasible_point_finder, self.chkACE_feasible_point_finder,
                                self.lblMBH_max_not_improve, self.txtMBH_max_not_improve,
                                self.lblMBH_max_trials, self.txtMBH_max_trials,
                                self.lblMBH_max_run_time, self.txtMBH_max_run_time,
                                self.lblMBH_hop_distribution, self.cmbMBH_hop_distribution,
                                self.lblMBH_max_step_size, self.txtMBH_max_step_size,
                                self.lblMBH_Pareto_alpha, self.txtMBH_Pareto_alpha,
                                self.lblMBH_time_hop_probability, self.txtMBH_time_hop_probability,
                                self.lblsnopt_feasibility_tolerance, self.txtsnopt_feasibility_tolerance,
                                self.lblsnopt_major_iterations, self.txtsnopt_major_iterations,
                                self.lblsnopt_max_run_time, self.txtsnopt_max_run_time,
                                self.lblderivative_type, self.cmbderivative_type,
                                self.lblcheck_derivatives, self.chkcheck_derivatives,
                                self.lblseed_MBH, self.chkseed_MBH,
                                self.lblinitial_guess_control_coordinate_system, self.cmbinitial_guess_control_coordinate_system,
                                self.lblinterpolate_initial_guess, self.chkinterpolate_initial_guess,
                                self.lblinitial_guess_num_timesteps, self.txtinitial_guess_num_timesteps,
                                self.lblinitial_guess_step_size_distribution, self.cmbinitial_guess_step_size_distribution,
                                self.lblinitial_guess_step_size_stdv_or_scale, self.txtinitial_guess_step_size_stdv_or_scale,
                                self.lblMBH_zero_control_initial_guess, self.cmbMBH_zero_control_initial_guess])
                                
        outerloopgrid = wx.GridSizer(12,2,0,0)
        
        self.lblrun_outerloop = wx.StaticText(self, -1, "Outer-Loop Solver")
        outerloop_choices = ["None","Genetic Algorithm"]
        self.cmbrun_outerloop = wx.ComboBox(self, -1, choices=outerloop_choices, style = wx.CB_READONLY)
        
        self.lblouterloop_popsize = wx.StaticText(self, -1, "Population size")
        self.txtouterloop_popsize = wx.TextCtrl(self, -1, "outerloop_popsize")
        
        self.lblouterloop_genmax = wx.StaticText(self, -1, "Maximum number of generations")
        self.txtouterloop_genmax = wx.TextCtrl(self, -1, "outerloop_genmax")
        
        self.lblouterloop_tournamentsize = wx.StaticText(self, -1, "Tournament size")
        self.txtouterloop_tournamentsize = wx.TextCtrl(self, -1, "outerloop_tournamentsize")
        
        self.lblouterloop_CR = wx.StaticText(self, -1, "Crossover ratio")
        self.txtouterloop_CR = wx.TextCtrl(self, -1, "outerloop_CR")
        
        self.lblouterloop_mu = wx.StaticText(self, -1, "Mutation rate")
        self.txtouterloop_mu = wx.TextCtrl(self, -1, "outerloop_mu")
        
        self.lblouterloop_stallmax = wx.StaticText(self, -1, "Maximum stall duration")
        self.txtouterloop_stallmax = wx.TextCtrl(self, -1, "outerloop_stallmax")
                
        self.lblouterloop_tolfit = wx.StaticText(self, -1, "Fitness tolerance")
        self.txtouterloop_tolfit = wx.TextCtrl(self, -1, "outerloop_tolfit")
                
        self.lblouterloop_ntrials = wx.StaticText(self, -1, "Number of trials")
        self.txtouterloop_ntrials = wx.TextCtrl(self, -1, "outerloop_ntrials")
                
        self.lblouterloop_elitecount = wx.StaticText(self, -1, "Number of elite individuals")
        self.txtouterloop_elitecount = wx.TextCtrl(self, -1, "outerloop_elitecount")
                
        self.lblouterloop_useparallel = wx.StaticText(self, -1, "Run outer-loop GA in parallel?")
        self.chkouterloop_useparallel = wx.CheckBox(self, -1)
        
        self.lblouterloop_warmstart = wx.StaticText(self, -1, "Warm-start the outer-loop?")
        self.txtouterloop_warmstart = wx.TextCtrl(self, -1, "outerloop_warmstart")
    
        outerloopgrid.AddMany([self.lblrun_outerloop, self.cmbrun_outerloop,
                               self.lblouterloop_popsize, self.txtouterloop_popsize,
                               self.lblouterloop_genmax, self.txtouterloop_genmax,
                               self.lblouterloop_tournamentsize, self.txtouterloop_tournamentsize,
                               self.lblouterloop_CR, self.txtouterloop_CR,
                               self.lblouterloop_mu, self.txtouterloop_mu,
                               self.lblouterloop_stallmax, self.txtouterloop_stallmax,
                               self.lblouterloop_tolfit, self.txtouterloop_tolfit,
                               self.lblouterloop_ntrials, self.txtouterloop_ntrials,
                               self.lblouterloop_elitecount, self.txtouterloop_elitecount,
                               self.lblouterloop_warmstart, self.txtouterloop_warmstart])

                                
        vboxleft = wx.BoxSizer(wx.VERTICAL)
        vboxright = wx.BoxSizer(wx.VERTICAL)
        lblLeftTitle = wx.StaticText(self, -1, "Inner-Loop Solver Parameters")
        lblRightTitle = wx.StaticText(self, -1, "Outer-Loop Solver Parameters")
        vboxleft.Add(lblLeftTitle)
        vboxleft.Add(innerloopgrid)
        vboxright.Add(lblRightTitle)
        vboxright.Add(outerloopgrid)
        
        font = self.GetFont()
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        lblLeftTitle.SetFont(font)
        lblRightTitle.SetFont(font)

        
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        hbox.Add(vboxleft)
        hbox.AddSpacer(20)
        hbox.Add(vboxright)
        
        
        self.lbltrialX = wx.StaticText(self, -1, "Trial decision vector or initial guess")
        self.txttrialX = wx.TextCtrl(self, -1, style=wx.TE_MULTILINE, size = (700,300))
        self.btntrialX = wx.Button(self, -1, "...")
        trialbox = wx.BoxSizer(wx.HORIZONTAL)
        trialbox.AddMany([self.lbltrialX, self.btntrialX])

        self.mainbox = wx.BoxSizer(wx.VERTICAL)
        self.mainbox.AddMany([hbox, trialbox, self.txttrialX])

        self.SetSizer(self.mainbox)
        self.SetupScrolling()

class PhysicsOptionsPanel(wx.lib.scrolledpanel.ScrolledPanel):
    def __init__(self, parent):
        
        wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent)

        ephemerisgrid = wx.GridSizer(4,2,5,5)
        perturbgrid = wx.GridSizer(4,2,5,5)
        
        self.lblephemeris_source = wx.StaticText(self, -1, "Ephemeris Source")
        ephemeris_source_typestypes = ['Static','SPICE']
        self.cmbephemeris_source = wx.ComboBox(self, -1, choices = ephemeris_source_typestypes, style=wx.CB_READONLY)

        self.lblSPICE_leap_seconds_kernel = wx.StaticText(self, -1, "Leap seconds kernel")
        self.txtSPICE_leap_seconds_kernel = wx.TextCtrl(self, -1, "SPICE_leap_seconds_kernel", size=(200,-1))

        self.lblSPICE_reference_frame_kernel = wx.StaticText(self, -1, "Frame kernel")
        self.txtSPICE_reference_frame_kernel = wx.TextCtrl(self, -1, "SPICE_reference_frame_kernel", size=(200,-1))

        self.lbluniverse_folder = wx.StaticText(self, -1, "Universe folder")
        self.txtuniverse_folder = wx.TextCtrl(self, -1, "universe_folder", size=(400,-1))
        self.btnGetNewUniverseFolder = wx.Button(self, -1, "...")
        self.btnSetDefaultUniverse = wx.Button(self, -1, "Default")
        UniverseButtonSizer = wx.BoxSizer(wx.HORIZONTAL)
        UniverseButtonSizer.AddMany([self.txtuniverse_folder, self.btnGetNewUniverseFolder, self.btnSetDefaultUniverse])

        self.lblperturb_SRP = wx.StaticText(self, -1, "Enable SRP")
        self.chkperturb_SRP = wx.CheckBox(self, -1)

        self.lblperturb_thirdbody = wx.StaticText(self, -1, "Enable third body")
        self.chkperturb_thirdbody = wx.CheckBox(self, -1)

        self.lblspacecraft_area = wx.StaticText(self, -1, "Spacecraft area (in m^2)")
        self.txtspacecraft_area = wx.TextCtrl(self, -1, "spacecraft_area")

        self.lblcoefficient_of_reflectivity = wx.StaticText(self, -1, "Coefficient of reflectivity")
        self.txtcoefficient_of_reflectivity = wx.TextCtrl(self, -1, "coefficient_of_reflectivity")

        ephemerisgrid.AddMany([self.lblephemeris_source, self.cmbephemeris_source,
                              self.lblSPICE_leap_seconds_kernel, self.txtSPICE_leap_seconds_kernel,
                              self.lblSPICE_reference_frame_kernel, self.txtSPICE_reference_frame_kernel,
                              self.lbluniverse_folder, UniverseButtonSizer])
        perturbgrid.AddMany([ self.lblperturb_SRP, self.chkperturb_SRP,
                              self.lblperturb_thirdbody, self.chkperturb_thirdbody,
                              self.lblspacecraft_area, self.txtspacecraft_area,
                              self.lblcoefficient_of_reflectivity, self.txtcoefficient_of_reflectivity])




        lblLeftTitle = wx.StaticText(self, -1, "Ephemeris settings")
        vboxleft = wx.BoxSizer(wx.VERTICAL)
        vboxleft.AddMany([lblLeftTitle, ephemerisgrid])

        lblRightTitle = wx.StaticText(self, -1, "Perturbation settings")
        vboxright = wx.BoxSizer(wx.VERTICAL)
        vboxright.AddMany([lblRightTitle, perturbgrid])

        font = self.GetFont()
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        lblLeftTitle.SetFont(font)
        lblRightTitle.SetFont(font)

        self.mainbox = wx.BoxSizer(wx.HORIZONTAL)
        
        self.mainbox.Add(vboxleft)
        self.mainbox.AddSpacer(20)
        self.mainbox.Add(vboxright)


        spiralgrid = wx.GridSizer(2,2,5,5)
        self.lblspiral_model_type = wx.StaticText(self, -1, "Spiral model type")
        spiral_model_choices = ['Battin','Edelbaum']
        self.cmbspiral_model_type = wx.ComboBox(self, -1, choices = spiral_model_choices, style = wx.CB_READONLY)
        spiralgrid.AddMany([self.lblspiral_model_type, self.cmbspiral_model_type])
        lblBottomTitle = wx.StaticText(self, -1, "Spiral settings")
        lblBottomTitle.SetFont(font)
        vboxspiral = wx.BoxSizer(wx.VERTICAL)
        vboxspiral.AddMany([lblBottomTitle, spiralgrid])

        lambertgrid = wx.GridSizer(2,2,5,5)
        self.lbllambert_type = wx.StaticText(self, -1, "Lambert solver type")
        lambert_choices = ['Arora-Russell','Izzo (not included in open-source)']
        self.cmblambert_type = wx.ComboBox(self, -1, choices = lambert_choices, style = wx.CB_READONLY)
        lambertgrid.AddMany([self.lbllambert_type, self.cmblambert_type])
        lblLambertTitle = wx.StaticText(self, -1, "Lambert settings")
        lblLambertTitle.SetFont(font)
        vboxlambert = wx.BoxSizer(wx.VERTICAL)
        vboxlambert.AddMany([lblLambertTitle, lambertgrid])

        self.mainvbox = wx.BoxSizer(wx.VERTICAL)
        self.mainvbox.Add(self.mainbox)
        self.mainvbox.AddSpacer(20)
        self.mainvbox.AddMany([vboxspiral, vboxlambert])

        self.SetSizer(self.mainvbox)
        self.SetupScrolling()

class OutputOptionsPanel(wx.lib.scrolledpanel.ScrolledPanel):
    def __init__(self, parent):
        
        wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent)

        self.mainbox = wx.FlexGridSizer(20,2,5,5)

        self.lblcreate_GMAT_script = wx.StaticText(self, -1, "Create GMAT scripts")
        self.chkcreate_GMAT_script = wx.CheckBox(self, -1)

        self.lbloutput_units = wx.StaticText(self, -1, "Output units")
        outputchoices = ['km and km/s','LU and LU/day']
        self.cmboutput_units = wx.ComboBox(self, -1, choices=outputchoices, style=wx.CB_READONLY)

        self.mainbox.AddMany([self.lblcreate_GMAT_script, self.chkcreate_GMAT_script,
                           self.lbloutput_units, self.cmboutput_units])

        self.SetSizer(self.mainbox)
        self.SetupScrolling()

        
class OptionsBook(wx.Notebook):
    #class for Options notebook
    def __init__(self, parent):
        wx.Notebook.__init__(self, parent=parent, id=wx.ID_ANY, style=
                             wx.BK_DEFAULT
                             #wx.BK_TOP 
                             #wx.BK_BOTTOM
                             #wx.BK_LEFT
                             #wx.BK_RIGHT
                             )
                             
        font = self.GetFont()
        font.SetPointSize(10)
        self.SetFont(font)


        #create tabs
        self.tabGlobal = GlobalOptionsPanel(self)
        self.AddPage(self.tabGlobal, "Global Mission Options")
        self.tabSpacecraft = SpacecraftOptionsPanel(self)
        self.AddPage(self.tabSpacecraft, "Spacecraft Options")
        self.tabJourney = JourneyOptionsPanel(self)
        self.AddPage(self.tabJourney, "Journey Options")
        self.tabSolver = SolverOptionsPanel(self)
        self.AddPage(self.tabSolver, "Solver Options")
        self.tabPhysics = PhysicsOptionsPanel(self)
        self.AddPage(self.tabPhysics, "Physics Options")
        self.tabOutput = OutputOptionsPanel(self)
        self.AddPage(self.tabOutput, "Output Options")