import wx
import os.path
import JourneyOptions as JO
import PyEMTG_interface as GUI

class MissionOptions(object):
    filename = ""
    success = []
    
    #problem type    
    problem_type = 0
    
    #physical constants
    G = 6.674280e-20
    g0 = 9.806650

    #outer loop solver settings
    run_outerloop = 0 #whether or not to run the outer loop if false
    outerloop_popsize = 200 #population size
    outerloop_genmax = 40 #number of generations
    outerloop_tournamentsize = 4 #tournament size for selection
    outerloop_CR = 0.3 #math::crossover ratio
    outerloop_mu = 0.5 #mutation rate
    outerloop_stallmax = 20 #maximum number of stall generations
    outerloop_tolfit = 1.0e-4 #fitness tolerance
    outerloop_ntrials = 1 #how many times to run the outer loop
    outerloop_elitecount = 1 #how many elite individuals to retain
    outerloop_useparallel = 0 #whether or not to use the parallel outer-loop
    outerloop_warmstart = 0 #if true, read "population.txt" and "solutions.txt"
    outerloop_warm_population = "none"
    outerloop_warm_archive = "none"
    outerloop_reevaluate_full_population = 0

    #outer loop selectable options settings
    outerloop_vary_power = 0
    outerloop_vary_launch_epoch = 0
    outerloop_vary_flight_time_upper_bound = 0
    outerloop_vary_thruster_type = 0
    outerloop_vary_number_of_thrusters = 0
    outerloop_vary_launch_vehicle = 0
    outerloop_vary_departure_C3 = 0
    outerloop_vary_arrival_C3 = 0
    outerloop_power_choices = [10.0]
    outerloop_launch_epoch_choices = [51544.5]
    outerloop_flight_time_upper_bound_choices = [365.25]
    outerloop_thruster_type_choices = [8]
    outerloop_number_of_thrusters_choices = [1]
    outerloop_launch_vehicle_choices = [1]
    outerloop_departure_C3_choices = [25.0]
    outerloop_arrival_C3_choices = [25.0]
    outerloop_restrict_flight_time_lower_bound = 0
    quiet_outerloop = 1#if true, suppress all text outputs except error catches

    #outer-loop point group settings
    outerloop_point_groups_number_to_score = [1]
    outerloop_point_groups_values = [1]
    outerloop_point_groups_members = [[1, 2, 3]]
    
    #outerloop objective settings
    outerloop_objective_function_choices = [2, 6]

    #inner loop solver settings
    NLP_solver_type = 0
    NLP_solver_mode = 1
    quiet_NLP = 0
    quiet_basinhopping = 0
    ACE_feasible_point_finder = 0
    MBH_max_not_improve = 50
    MBH_max_trials = 100000
    MBH_max_run_time = 600
    MBH_max_step_size = 1.0
    MBH_hop_distribution = 1
    MBH_time_hop_probability = 0.05
    MBH_Pareto_alpha = 3.0
    snopt_feasibility_tolerance = 1.0e-6
    snopt_major_iterations = 8000
    snopt_max_run_time = 3600
    derivative_type = 0
    seed_MBH = 0
    interpolate_initial_guess = 0
    initial_guess_num_timesteps = 10
    initial_guess_step_size_distribution = 0 #0: uniform, 1: Gaussian, 2: Cauchy
    initial_guess_step_size_stdv_or_scale = 1.0
    MBH_zero_control_initial_guess = 0
    MBH_two_step = 0 #whether or not to use the 2-step MBH (coarse then fine derivatives)
    FD_stepsize = 1.5e-8 #"fine" finite differencing step size
    FD_stepsize_coarse = 1.0e-4 #"coarse" finite differencing step

    #problem settings set by the user

    #ephemeris data
    ephemeris_source = 1
    SPICE_leap_seconds_kernel = "naif0009.tls"
    SPICE_reference_frame_kernel = "pck00010.tpc"
    universe_folder = "../Universe/"

    #Lambert solver
    LambertSolver = 0 #0: Arora-Russell, 1: Izzo (not included in open-source)

    #low thrust solver parameters
    num_timesteps = 10 #number of timesteps per phase
    control_coordinate_system = 0 #0: cartesian, 1: polar
    initial_guess_control_coordinate_system = 0 #0: cartesian, 1: polar
    step_size_distribution = 0 #0: uniform, 1: Gaussian, 2: Cauchy
    step_size_stdv_or_scale = 1.0
    spiral_model_type = 1#0: Battin, 1: Edelbaum

    #impulsive thrust solver parameters
    maximum_number_of_lambert_revolutions = 0

    #vehicle parameters
    maximum_mass = 1000 #the maximum possible mass of the spacecraft (negative number means use LV max)
    allow_initial_mass_to_vary = 0 #flag on whether or not the solver can choose the initial mass (make the spacecraft wet mass lighter)
    LV_margin = 0.0 #launch vehicle margin
    LV_adapter_mass = 0.0 #launch vehicle adapter mass (kg)
    IspLT = 1 #specific impulse of the engine used for low-thrust maneuvers
    IspLT_minimum = 1 #minimum Isp for VSI systems
    IspChem = 10000 #specific impulse of the engine used for impulsive maneuvers
    IspDS = 1 #specific impulse for the earth departure stage, if applicable
    Thrust = 0.01 #thrust of the spacecraft, in Newtons
    #         -2: custom launch vehicle
    #         -1: burn with departure stage engine
    #         0: fixed initial mass
    #         1: Atlas V (401)	 NLSII
    #         2: Atlas V (411)	NLSII
    #         3: Atlas V (421)	NLSII
    #         4: Atlas V (431)	NLSII
    #         5: Atlas V (501)	NLSII
    #         6: Atlas V (511)	NLSII
    #         7: Atlas V (521)	NLSII
    #         8: Atlas V (531)	NLSII
    #         9: Atlas V (541)	NLSII
    #         10: Atlas V (551)	NLSII
    #         11: Falcon 9 (v1.0)	NLSII
    #         12: Falcon 9 (v1.1)	NLSII
    #         13: Atlas V (551) w/Star 48	NLSI
    #         14: Falcon 9 Heavy
    #         15: Delta IV Heavy	NLSI
    #         16: SLS Block 1
    LV_type = 0

    custom_LV_coefficients = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    custom_LV_C3_bounds = [0.0, 0.0]
    parking_orbit_altitude = 300.0
    parking_orbit_inclination = 28.5
    engine_type = 0
    number_of_engines = 1 #only relevant when modeling an engine
    throttle_logic_mode = 0
    throttle_sharpness = 100.0
    engine_duty_cycle = 1.0 #percentage of time that engine can operate
    power_at_1_AU = 10.0 #in kW
    power_source_type = 0 #0: solar, 1: radioisotope (or other fixed power)
    solar_power_gamma = [0.0, 0.0, 0.0, 0.0, 0.0] #coefficients for solar panels
    power_margin = 0.0 #for propulsion, as a fraction
    engine_input_thrust_coefficients = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    engine_input_mass_flow_rate_coefficients = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    engine_input_power_bounds = [1.0, 5.0]
    user_defined_engine_efficiency = 0.7
    spacecraft_power_coefficients = [0.0, 0.0, 0.0]
    spacecraft_power_model_type = 0
    EP_dry_mass = 0 #in kg
    power_decay_rate = 0.0 #percent per year
        
    #minimum dry mass constraint and related parameters
    minimum_dry_mass = 0 #in kg
    enable_maximum_propellant_mass_constraint = 0
    maximum_propellant_mass = 1000.0
    post_mission_delta_v = 0.0 #in km/s
    post_mission_Isp = 3000.0 #in s
    propellant_margin = 0.0 #percent
        
        
    #perturbation settings
    perturb_SRP = 0
    perturb_thirdbody = 0
    spacecraft_area = 1 #in m^2
    coefficient_of_reflectivity = 1

    #global problem settings
    number_of_journeys = 1
    max_phases_per_journey = 8
    destination_list = [3, 4] #integers
    include_initial_impulse_in_cost = 0.0
    global_timebounded = 1#0: unbounded, 1: bounded total time (note that the global arrival date bound is by definition the same as the last journey"s arrival date bound and is not duplicated
    launch_window_open_date = 60000.0#MJD
    total_flight_time_bounds = [0, 1000.0]#[2] days
    objective_type = 2 #0: minimum deltaV, 1: minimum time, #2: maximum final mass
    DLA_bounds = [-28.5, 28.5] #DLA in degrees
    mission_name = "default"
    mission_type = 2
    initial_V_infinity = [0.0, 0.0, 0.0]
    forced_post_launch_coast = 0.0
    forced_flyby_coast = 0.0
        
    #array of JourneyOptions objects
    Journeys = []
    ActiveJourney = 0
        
    #output format settings
    output_units = 0 #0: km and km/s, 1: LU and LU/day
    create_GMAT_script = 0 #0: no, 1: yes
    generate_initial_guess_file = 0
    mission_type_for_initial_guess_file = 2
    override_working_directory = 0;
    forced_working_directory = "..//EMTG_v8_Results"
    generate_forward_integrated_ephemeris = 0#0 :no, 1: yes

    #debug code
    check_derivatives = 0
    run_inner_loop = 2
    number_of_trial_sequences = 1
    trialX = []

    #************************************************************************************constructor
    def __init__(self, input_file_name):
       self.Journeys = []
       self.trialX = []
       self.parse_options_file(input_file_name)

    #************************************************************************************parse options file
    def parse_options_file(self, input_file_name):
        #Step 1: open the file

        self.filename = input_file_name

        if os.path.isfile(self.filename):
            inputfile = open(input_file_name, "r")
            self.success = 1
        else:
            print "Unable to open", input_file_name, "EMTG Error"
            self.success = 0
            return

        perturb_line_flag = 0
        point_group_members_flag = 0
        sequence_line_flag = 0
        trialX_line_flag = 0
        flyby_choice_line_flag = 0
        destination_choice_line_flag = 0


        #Step 2: scan through the file
        linenumber = 0
        for line in inputfile:
            #strip off the newline character
            line = line.replace("\n","")
            linenumber = linenumber + 1

            if line != "":
                if line[0] != "#":
                    #this is an active line, so it is space delimited
                    linecell = line.split(" ")
                    
                    choice = linecell[0]

                    if perturb_line_flag > 0:
                        perturb_line_flag = perturb_line_flag + 1
                        self.Journeys[j].journey_perturbation_bodies = []
                        for x in linecell:
                            self.Journeys[j].journey_perturbation_bodies.append(int(x))

                    elif point_group_members_flag > 0:
                        point_group_members_flag += 1
                        temp_point_group = []
                        for entry in linecell[1:]:
                            temp_point_group.append(int(entry))
                        self.outerloop_point_groups_members.append(temp_point_group)

                    elif trialX_line_flag > 0:
                        temp_trialX = []
                        for entry in linecell:
                            temp_trialX.append(float(entry))
                        self.trialX.append(temp_trialX)
                        trialX_line_flag = trialX_line_flag + 1
                        

                    elif sequence_line_flag > 0:
                        for j in range(0, self.number_of_journeys):
                            seq = [0] * self.max_phases_per_journey
                            if sequence_line_flag == 1:
                                self.Journeys[j].sequence = []
                                self.Journeys[j].number_of_phases = []
                            for p in range (0, self.max_phases_per_journey):
                                seq[p] = int(linecell[self.max_phases_per_journey*j+p])

                            self.Journeys[j].sequence.append(seq)
                            self.Journeys[j].number_of_phases.append(sum(1 for x in seq if x > 0))
                        sequence_line_flag = sequence_line_flag + 1

                    elif flyby_choice_line_flag > 0:
                        self.Journeys[flyby_choice_line_flag - 1].outerloop_journey_flyby_sequence_choices = []
                        for entry in linecell[1:]:
                            self.Journeys[flyby_choice_line_flag - 1].outerloop_journey_flyby_sequence_choices.append(int(entry))
                        flyby_choice_line_flag = flyby_choice_line_flag + 1

                    elif destination_choice_line_flag > 0:
                        self.Journeys[destination_choice_line_flag - 1].outerloop_journey_destination_choices = []
                        for entry in linecell[1:]:
                            self.Journeys[destination_choice_line_flag - 1].outerloop_journey_destination_choices.append(int(entry))
                        destination_choice_line_flag = destination_choice_line_flag + 1

                    elif choice == "problem_type":
                        self.problem_type = int(linecell[1])

                    #global physical constants
                    elif choice == "G":
                        self.G = float(linecell[1])
                    elif choice == "g0":
                        self.g0 = float(linecell[1])

                    #outer-loop solver settings
                    elif choice == "run_outerloop":
                        self.run_outerloop = int(linecell[1])
                    elif choice == "outerloop_popsize":
                        self.outerloop_popsize = int(linecell[1])
                    elif choice == "outerloop_genmax":
                        self.outerloop_genmax = int(linecell[1])
                    elif choice == "outerloop_tournamentsize":
                        self.outerloop_tournamentsize = int(linecell[1])
                    elif choice == "outerloop_CR":
                        self.outerloop_CR = float(linecell[1])
                    elif choice == "outerloop_mu":
                        self.outerloop_mu = float(linecell[1])
                    elif choice == "outerloop_stallmax":
                        self.outerloop_stallmax = int(linecell[1])
                    elif choice == "outerloop_tolfit":
                        self.outerloop_tolfit = float(linecell[1])
                    elif choice == "outerloop_ntrials":
                        self.outerloop_ntrials = int(linecell[1])
                    elif choice == "outerloop_elitecount":
                        self.outerloop_elitecount = int(linecell[1])
                    elif choice == "outerloop_useparallel":
                        self.outerloop_useparallel = int(linecell[1])
                    elif choice == "outerloop_warmstart":
                        self.outerloop_warmstart = int(linecell[1])
                    elif choice == "outerloop_warm_population":
                        self.outerloop_warm_population = linecell[1]
                    elif choice == "outerloop_warm_archive":
                        self.outerloop_warm_archive = linecell[1]
                    elif choice == "outerloop_reevaluate_full_population":
                        self.outerloop_reevaluate_full_population = int(linecell[1])
                    elif choice == "quiet_outerloop":
                        self.quiet_outerloop = int(linecell[1])

                    #outer loop selectable options settings
                    elif choice == "outerloop_vary_power":
                        self.outerloop_vary_power = int(linecell[1])
                    elif choice == "outerloop_vary_launch_epoch":
                        self.outerloop_vary_launch_epoch = int(linecell[1])
                    elif choice == "outerloop_vary_flight_time_upper_bound":
                        self.outerloop_vary_flight_time_upper_bound = int(linecell[1])
                    elif choice == "outerloop_restrict_flight_time_lower_bound":
                        self.outerloop_restrict_flight_time_lower_bound = int(linecell[1])
                    elif choice == "outerloop_vary_thruster_type":
                        self.outerloop_vary_thruster_type = int(linecell[1])
                    elif choice == "outerloop_vary_number_of_thrusters":
                        self.outerloop_vary_number_of_thrusters = int(linecell[1])
                    elif choice == "outerloop_vary_launch_vehicle":
                        self.outerloop_vary_launch_vehicle = int(linecell[1])
                    elif choice == "outerloop_vary_departure_C3":
                        self.outerloop_vary_departure_C3 = int(linecell[1])
                    elif choice == "outerloop_vary_arrival_C3":
                        self.outerloop_vary_arrival_C3 = int(linecell[1])
                    elif choice == "outerloop_vary_journey_destination":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].outerloop_vary_journey_destination = int(linecell[j+1])
                    elif choice == "outerloop_vary_journey_flyby_sequence":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].outerloop_vary_journey_flyby_sequence = int(linecell[j+1])
                    elif choice == "outerloop_power_choices":
                        self.outerloop_power_choices = []
                        for x in linecell[1:]:
                            self.outerloop_power_choices.append(float(x))
                    elif choice == "outerloop_launch_epoch_choices":
                        self.outerloop_launch_epoch_choices = []
                        for x in linecell[1:]:
                            self.outerloop_launch_epoch_choices.append(float(x))
                    elif choice == "outerloop_flight_time_upper_bound_choices":
                        self.outerloop_flight_time_upper_bound_choices = []
                        for x in linecell[1:]:
                            self.outerloop_flight_time_upper_bound_choices.append(float(x))
                    elif choice == "outerloop_thruster_type_choices":
                        self.outerloop_thruster_type_choices = []
                        for x in linecell[1:]:
                            self.outerloop_thruster_type_choices.append(int(float(x)))
                    elif choice == "outerloop_number_of_thrusters_choices":
                        self.outerloop_number_of_thrusters_choices = []
                        for x in linecell[1:]:
                            self.outerloop_number_of_thrusters_choices.append(int(float(x)))
                    elif choice == "outerloop_launch_vehicle_choices":
                        self.outerloop_launch_vehicle_choices = []
                        for x in linecell[1:]:
                            self.outerloop_launch_vehicle_choices.append(int(float(x)))
                    elif choice == "outerloop_departure_C3_choices":
                        self.outerloop_departure_C3_choices = []
                        for x in linecell[1:]:
                            self.outerloop_departure_C3_choices.append(float(x))
                    elif choice == "outerloop_arrival_C3_choices":
                        self.outerloop_arrival_C3_choices = []
                        for x in linecell[1:]:
                            self.outerloop_arrival_C3_choices.append(float(x))
                    elif choice == "outerloop_journey_flyby_sequence_choices":
                        flyby_choice_line_flag = 1
                    elif choice == "outerloop_journey_destination_choices":
                        destination_choice_line_flag = 1
                    elif choice == "outerloop_journey_maximum_number_of_flybys":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].outerloop_journey_maximum_number_of_flybys = linecell[j+1]

                    #outer-loop point groups settings
                    elif choice == "outerloop_point_groups_values":
                        self.outerloop_point_groups_values = []
                        for x in linecell[1:]:
                            self.outerloop_point_groups_values.append(int(x))

                        #start reading point group members
                        self.outerloop_point_groups_members
                        point_group_members_flag = 1

                    elif choice == "outerloop_point_groups_number_to_score":
                        self.outerloop_point_groups_number_to_score = []
                        for x in linecell[1:]:
                            self.outerloop_point_groups_number_to_score.append(int(x))

                    #outerloop objective settings
                    elif choice == "outerloop_objective_function_choices":
                        self.outerloop_objective_function_choices = []
                        for x in linecell[1:]:
                            self.outerloop_objective_function_choices.append(int(x))

                    
                    #inner loop solver settings
                    elif choice == "NLP_solver_type":
                        self.NLP_solver_type = int(linecell[1])
                    elif choice == "NLP_solver_mode":
                        self.NLP_solver_mode = int(linecell[1])
                    elif choice == "ACE_feasible_point_finder":
                        self.ACE_feasible_point_finder = int(linecell[1])
                    elif choice == "quiet_NLP":
                        self.quiet_NLP = int(linecell[1])
                    elif choice == "quiet_basinhopping":
                        self.quiet_basinhopping = int(linecell[1])
                    elif choice ==  "MBH_max_not_improve":
                        self.MBH_max_not_improve = int(linecell[1])
                    elif choice ==  "MBH_max_trials":
                        self.MBH_max_trials = int(linecell[1])
                    elif choice ==  "MBH_max_run_time":
                        self.MBH_max_run_time = int(linecell[1])
                    elif choice ==  "MBH_max_step_size":
                        self.MBH_max_step_size = float(linecell[1])
                    elif choice == "MBH_hop_distribution":
                        self.MBH_hop_distribution = int(linecell[1])
                    elif choice == "MBH_Pareto_alpha":
                        self.MBH_Pareto_alpha = float(linecell[1])
                    elif choice ==  "MBH_time_hop_probability":
                        self.MBH_time_hop_probability = float(linecell[1])
                    elif choice ==  "snopt_feasibility_tolerance":
                        self.snopt_feasibility_tolerance = float(linecell[1])
                    elif choice ==  "snopt_major_iterations":
                        self.snopt_major_iterations = int(linecell[1])
                    elif choice ==  "snopt_max_run_time":
                        self.snopt_max_run_time = int(linecell[1])
                    elif choice ==  "derivative_type":
                        self.derivative_type = int(linecell[1])
                    elif choice ==  "seed_MBH":
                        self.seed_MBH = int(linecell[1])
                    elif choice ==  "interpolate_initial_guess":
                        self.interpolate_initial_guess = int(linecell[1])
                    elif choice ==  "initial_guess_num_timesteps":
                        self.initial_guess_num_timesteps = int(linecell[1])
                    elif choice ==  "initial_guess_step_size_distribution":
                        self.initial_guess_step_size_distribution = int(linecell[1])
                    elif choice ==  "initial_guess_step_size_stdv_or_scale":
                        self.initial_guess_step_size_stdv_or_scale = float(linecell[1])		
                    elif choice == "MBH_zero_control_initial_guess":
                        self.MBH_zero_control_initial_guess = int(linecell[1])
                    elif choice == "MBH_two_step":
                        self.MBH_two_step = int(linecell[1])
                    elif choice == "FD_stepsize":
                        self.FD_stepsize = float(linecell[1])
                    elif choice == "FD_stepsize_coarse":
                        self.FD_stepsize_coarse = float(linecell[1])

                    #problem settings set by the user
                    elif choice ==  "ephemeris_source":
                        self.ephemeris_source = int(linecell[1])
                    elif choice ==  "SPICE_leap_seconds_kernel":
                        self.SPICE_leap_seconds_kernel = linecell[1]
                    elif choice ==  "SPICE_reference_frame_kernel":
                        self.SPICE_reference_frame_kernel = linecell[1]
                    elif choice ==  "universe_folder":
                        self.universe_folder = linecell[1]

                    elif choice == "LambertSolver":
                        self.LambertSolver = int(linecell[1])

                    #low thrust solver parameters
                    elif choice ==  "num_timesteps":
                        self.num_timesteps = int(linecell[1])
                    elif choice == "control_coordinate_system":
                        self.control_coordinate_system = int(linecell[1])
                    elif choice == "initial_guess_control_coordinate_system":
                        self.initial_guess_control_coordinate_system = int(linecell[1])
                    elif choice ==  "step_size_distribution":
                        self.step_size_distribution = int(linecell[1])
                    elif choice ==  "step_size_stdv_or_scale":
                        self.step_size_stdv_or_scale = float(linecell[1])
                    elif choice == "spiral_model_type":
                        self.spiral_model_type = int(linecell[1])


                    #impulsive thrust solver parameters
                    elif choice == "maximum_number_of_lambert_revolutions":
                        self.maximum_number_of_lambert_revolutions = int(linecell[1])

                    #vehicle parameters
                    elif choice == "maximum_mass":
                        self.maximum_mass = float(linecell[1])
                    elif choice == "enable_maximum_propellant_mass_constraint":
                        self.enable_maximum_propellant_mass_constraint = int(linecell[1])
                    elif choice == "maximum_propellant_mass":
                        self.maximum_propellant_mass = float(linecell[1])
                    elif choice == "LV_margin":
                        self.LV_margin = float(linecell[1])
                    elif choice == "LV_adapter_mass":
                        self.LV_adapter_mass = float(linecell[1])
                    elif choice == "custom_LV_C3_bounds":
                        self.custom_LV_C3_bounds = []
                        for x in linecell[1:]:
                            self.custom_LV_C3_bounds.append(float(x))
                    elif choice == "custom_LV_coefficients":
                        self.custom_LV_coefficients = []
                        for x in linecell[1:]:
                            self.custom_LV_coefficients.append(float(x))
                    elif choice == "IspLT":
                        self.IspLT = float(linecell[1])
                    elif choice == "IspLT_minimum":
                        self.IspLT_minimum = float(linecell[1])
                    elif choice == "IspChem":
                        self.IspChem = float(linecell[1])
                    elif choice == "IspDS":
                        self.IspDS = float(linecell[1])
                    elif choice == "Thrust":
                        self.Thrust = float(linecell[1])
                    elif choice == "LV_type":
                        self.LV_type = float(linecell[1])
                    elif choice == "engine_type":
                        self.engine_type = int(linecell[1])
                    elif choice == "number_of_engines":
                        self.number_of_engines = int(linecell[1])
                    elif choice == "throttle_logic_mode":
                        self.throttle_logic_mode = int(linecell[1])
                    elif choice == "throttle_sharpness":
                        self.throttle_sharpness = float(linecell[1])
                    elif choice == "engine_duty_cycle":
                        self.engine_duty_cycle = float(linecell[1])
                    elif choice == "power_at_1_AU":
                        self.power_at_1_AU = float(linecell[1])
                    elif choice == "power_source_type":
                        self.power_source_type = int(linecell[1])
                    elif choice == "solar_power_gamma":
                        self.solar_power_gamma = []
                        for x in linecell[1:]:
                            self.solar_power_gamma.append(float(x))
                    elif choice == "power_margin":
                        self.power_margin = float(linecell[1])
                    elif choice == "power_decay_rate":
                        self.power_decay_rate = float(linecell[1])
                    elif choice == "spacecraft_power_coefficients":
                        self.spacecraft_power_coefficients = []
                        for x in linecell[1:]:
                            self.spacecraft_power_coefficients.append(float(x))
                    elif choice == "engine_input_thrust_coefficients":
                        self.engine_input_thrust_coefficients = []
                        for x in linecell[1:]:
                            self.engine_input_thrust_coefficients.append(float(x))
                        while len(self.engine_input_thrust_coefficients) < 7:
                            self.engine_input_thrust_coefficients.append(0.0)
                    elif choice == "engine_input_mass_flow_rate_coefficients":
                        self.engine_input_mass_flow_rate_coefficients = []
                        for x in linecell[1:]:
                            self.engine_input_mass_flow_rate_coefficients.append(float(x))
                        while len(self.engine_input_mass_flow_rate_coefficients) < 7:
                            self.engine_input_mass_flow_rate_coefficients.append(0.0)
                    elif choice == "engine_input_power_bounds":
                        self.engine_input_power_bounds = []
                        for x in linecell[1:]:
                            self.engine_input_power_bounds.append(float(x))
                    elif choice == "user_defined_engine_efficiency":
                        self.user_defined_engine_efficiency = float(linecell[1])
                    elif choice == "spacecraft_power_model_type":
                        self.spacecraft_power_model_type = int(linecell[1])
                    elif choice == "EP_dry_mass":
                        self.EP_dry_mass = float(linecell[1])
                    elif choice == "allow_initial_mass_to_vary":
                        self.allow_initial_mass_to_vary = int(linecell[1])
                    elif choice == "parking_orbit_altitude":
                        self.parking_orbit_altitude = float(linecell[1])
                    elif choice == "parking_orbit_inclination":
                        self.parking_orbit_inclination = float(linecell[1])

                    elif choice == "minimum_dry_mass":
                        self.minimum_dry_mass = float(linecell[1])
                    elif choice == "post_mission_delta_v":
                        self.post_mission_delta_v = float(linecell[1])
                    elif choice == "post_mission_Isp":
                        self.post_mission_Isp = float(linecell[1])
                    elif choice == "propellant_margin":
                        self.propellant_margin = float(linecell[1])

                    elif choice == "number_of_journeys":
                        self.number_of_journeys = int(linecell[1])
                        for j in range(0, self.number_of_journeys):
                            temp_JourneyOptions = JO.JourneyOptions(self.mission_type)
                            self.Journeys.append(temp_JourneyOptions)

                    elif choice == "max_phases_per_journey":
                        self.max_phases_per_journey = int(linecell[1])
                    elif choice == "destination_list":
                        for j in range(0,self.number_of_journeys):
                            self.Journeys[j].destination_list = [int(linecell[2*j+1]), int(linecell[2*j+2])]
                    elif choice == "include_initial_impulse_in_cost":
                        self.include_initial_impulse_in_cost = int(float(linecell[1]))
                    elif choice == "global_timebounded":
                        self.global_timebounded = int(linecell[1])
                    elif choice == "launch_window_open_date":
                        self.launch_window_open_date = float(linecell[1])
                    elif choice == "total_flight_time_bounds":
                        self.total_flight_time_bounds = [float(linecell[1]), float(linecell[2])]
                    elif choice == "objective_type":
                        self.objective_type = int(linecell[1])
                    elif choice == "DLA_bounds":
                        self.DLA_bounds = [float(linecell[1]), float(linecell[2])]
                    elif choice == "mission_name":
                        self.mission_name = linecell[1]
                    elif choice == "mission_type":
                        self.mission_type = int(linecell[1])
                    elif choice == "initial_V_infinity":
                        self.initial_V_infinity = [float(linecell[1]), float(linecell[2]), float(linecell[3])]
                    elif choice == "forced_post_launch_coast":
                        self.forced_post_launch_coast = float(linecell[1])
                    elif choice == "forced_flyby_coast":
                        self.forced_flyby_coast = float(linecell[1])

                    #parse all of the journey options and load them
                    #into the journey objects
                    elif choice == "journey_names":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_names = linecell[j+1]

                    elif choice == "journey_starting_mass_increment":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_starting_mass_increment = float(linecell[j+1])
                    
                    elif choice == "journey_variable_mass_increment":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_variable_mass_increment = int(float(linecell[j+1]))
                    
                    elif choice == "journey_timebounded":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_timebounded = int(float(linecell[j+1]))
                    
                    elif choice == "journey_wait_time_bounds":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_wait_time_bounds = [float(linecell[2*j+1]), float(linecell[2*j+2])]
                    
                    elif choice == "journey_flight_time_bounds":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_flight_time_bounds = [float(linecell[2*j+1]), float(linecell[2*j+2])]
                    
                    elif choice == "journey_arrival_date_bounds":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_arrival_date_bounds = [float(linecell[2*j+1]), float(linecell[2*j+2])]
                    
                    elif choice == "journey_initial_impulse_bounds":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_initial_impulse_bounds = [float(linecell[2*j+1]), float(linecell[2*j+2])]
                    
                    elif choice == "journey_departure_type":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_departure_type = int(float(linecell[j+1]))
                    
                    elif choice == "journey_departure_elements_type":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_departure_elements_type = int(float(linecell[j+1]))
                    
                    elif choice == "journey_departure_elements":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_departure_elements = [float(linecell[6*j+1]), float(linecell[6*j+2]), float(linecell[6*j+3]), float(linecell[6*j+4]), float(linecell[6*j+5]), float(linecell[6*j+6])]
                    
                    elif choice == "journey_departure_elements_vary_flag":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_departure_elements_vary_flag = [int(float(linecell[6*j+1])), int(float(linecell[6*j+2])), int(float(linecell[6*j+3])), int(float(linecell[6*j+4])), int(float(linecell[6*j+5])), int(float(linecell[6*j+6]))]
                    
                    elif choice == "journey_departure_elements_bounds":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_departure_elements_bounds = [float(linecell[6*j+1]), float(linecell[6*j+2]), float(linecell[6*j+3]), float(linecell[6*j+4]), float(linecell[6*j+5]), float(linecell[6*j+6]), float(linecell[6*j+7]), float(linecell[6*j+8]), float(linecell[6*j+9]), float(linecell[6*j+10]), float(linecell[6*j+11]), float(linecell[6*j+12])]
                    
                    elif choice == "journey_arrival_type":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_arrival_type = int(float(linecell[j+1]))
                    
                    elif choice == "journey_arrival_elements_type":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_arrival_elements_type = int(float(linecell[j+1]))
                    
                    elif choice == "journey_arrival_elements":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_arrival_elements = [float(linecell[6*j+1]), float(linecell[6*j+2]), float(linecell[6*j+3]), float(linecell[6*j+4]), float(linecell[6*j+5]), float(linecell[6*j+6])]
                    
                    elif choice == "journey_arrival_elements_vary_flag":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_arrival_elements_vary_flag = [int(float(linecell[6*j+1])), int(float(linecell[6*j+2])), int(float(linecell[6*j+3])), int(float(linecell[6*j+4])), int(float(linecell[6*j+5])), int(float(linecell[6*j+6]))]
                    
                    elif choice == "journey_arrival_elements_bounds":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_arrival_elements_bounds = [float(linecell[6*j+1]), float(linecell[6*j+2]), float(linecell[6*j+3]), float(linecell[6*j+4]), float(linecell[6*j+5]), float(linecell[6*j+6]), float(linecell[6*j+7]), float(linecell[6*j+8]), float(linecell[6*j+9]), float(linecell[6*j+10]), float(linecell[6*j+11]), float(linecell[6*j+12])]
                    
                    elif choice == "journey_central_body":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_central_body = linecell[j+1]
                    
                    elif choice == "journey_initial_velocity":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_initial_velocity = [float(linecell[j*3+1]), float(linecell[j*3+2]), float(linecell[j*3+3])]
                        
                    elif choice == "journey_final_velocity":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_final_velocity = [float(linecell[j*3+1]), float(linecell[j*3+2]), float(linecell[j*3+3])]

                    elif choice =="journey_arrival_declination_constraint_flag":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_arrival_declination_constraint_flag = int(float(linecell[j+1]))

                    elif choice == "journey_arrival_declination_bounds":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_arrival_declination_bounds = [float(linecell[j*2+1]), float(linecell[j*2+2])]
                    elif choice == "journey_escape_spiral_starting_radius":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_escape_spiral_starting_radius = float(linecell[j+1])
                    elif choice == "journey_capture_spiral_final_radius":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_capture_spiral_final_radius = float(linecell[j+1])
                    elif choice =="journey_maximum_DSM_magnitude_constraint_flag":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_maximum_DSM_magnitude_constraint_flag = int(float(linecell[j+1]))
                    elif choice == "journey_maximum_DSM_magnitude_constraint":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_maximum_DSM_magnitude_constraint = float(linecell[j+1])
                    
                        
                    #perturbation-related quantities    
                    elif choice == "perturb_SRP":
                        self.perturb_SRP = int(linecell[1])
                    elif choice == "perturb_thirdbody":
                        self.perturb_thirdbody = int(linecell[1])
                    elif choice == "journey_perturbation_bodies":
                        #first get the number of perturbation bodies for each journey
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_number_of_perturbation_bodies = int(linecell[j+1])

                        #next read the bodies line for each journey
                        perturb_line_flag = 1
                    
                    elif choice == "spacecraft_area":
                        self.spacecraft_area = float(linecell[1])
                    elif choice == "coefficient_of_reflectivity":
                        self.coefficient_of_reflectivity = float(linecell[1])
                            
                    #output format settings
                    elif choice == "output_units":
                        self.output_units = int(linecell[1])
                    elif choice == "create_GMAT_script":
                        self.create_GMAT_script = int(linecell[1])
                    elif choice == "generate_initial_guess_file":
                        self.generate_initial_guess_file = int(linecell[1])
                    elif choice == "mission_type_for_initial_guess_file":
                        self.mission_type_for_initial_guess_file = int(linecell[1])
                    elif choice == "override_working_directory":
                        self.override_working_directory = int(linecell[1])
                    elif choice == "forced_working_directory":
                        self.forced_working_directory = linecell[1]
                    elif choice == "generate_forward_integrated_ephemeris":
                        self.generate_forward_integrated_ephemeris = int(linecell[1])
                                
                    #trialX, sequence input, etc
                    elif choice == "check_derivatives":
                        self.check_derivatives = int(linecell[1])
                    elif choice == "run_inner_loop":
                        self.run_inner_loop = int(linecell[1])
                    elif choice == "sequence":
                        self.number_of_trial_sequences = int(linecell[1])
                        sequence_line_flag = 1
                       
                    elif choice == "phase_type":
                        print "phase_type is not currently a valid option"
                    elif choice == "trialX":
                        trialX_line_flag = 1

                    #deprecated options
                    elif choice == "NeuroSpiral_number_of_layers" or choice == "NeuroSpiral_neurons_per_layer" or choice == "journey_capture_spiral_starting_radius" or choice == "lazy_race_tree_allow_duplicates":
                        print choice, " is deprecated."
                                
                    #if option is not recognized
                    else:
                        errorstring = "Option not recognized: " + str(linecell[0]) + " on line " + str(linenumber)
                        print errorstring
                        self.success = 0
                        return
                else:
                    perturb_line_flag = 0
                    point_group_members_flag = 0
                    sequence_line_flag = 0
                    trialX_line_flag = 0
                    flyby_choice_line_flag = 0
                    destination_choice_line_flag = 0
            else:
                perturb_line_flag = 0
                point_group_members_flag = 0
                sequence_line_flag = 0
                trialX_line_flag = 0
                flyby_choice_line_flag = 0
                destination_choice_line_flag = 0
        inputfile.close()

    def write_options_file(self, output_file_name):
        #first make some error-preventing correction
        if (self.IspChem < 1.0):
            self.IspChem = 1.0
            
        #first open the file for writing
        outputfile = open(output_file_name, "w")
        
        outputfile.write("##Options file for EMTG_v8\n")
        outputfile.write("\n")
            
        outputfile.write("##problem type\n")
        outputfile.write("#0: standard EMTG mission\n")
        outputfile.write("problem_type " + str(self.problem_type) + "\n")
        outputfile.write("\n")
            
        outputfile.write("##physical constants\n")
        outputfile.write("#G in km^3/kg/s^2\n")
        outputfile.write("G " + str(self.G) + "\n")
        outputfile.write("#gravity at sea level on Earth in m/s^2\n")
        outputfile.write("g0 " + str(self.g0) + "\n")
        outputfile.write("\n")

        outputfile.write("##outer-loop solver settings\n")
        outputfile.write("#Do you want to run an outer-loop?\n")
        outputfile.write("#0: no\n")
        outputfile.write("#1: Genetic algorithm (number of objective functions determines which GA to run)\n")
        outputfile.write("run_outerloop " + str(self.run_outerloop) + "\n")
        outputfile.write("#outer-loop population size\n")	
        outputfile.write("outerloop_popsize " + str(self.outerloop_popsize) + "\n")
        outputfile.write("#maximum number of outer-loop generations\n")	
        outputfile.write("outerloop_genmax " + str(self.outerloop_genmax) + "\n")
        outputfile.write("#tournament size for selection\n")	
        outputfile.write("outerloop_tournamentsize " + str(self.outerloop_tournamentsize) + "\n")
        outputfile.write("#crossover ratio\n")	
        outputfile.write("outerloop_CR " + str(self.outerloop_CR) + "\n")
        outputfile.write("#mutation rate\n")	
        outputfile.write("outerloop_mu " + str(self.outerloop_mu) + "\n")
        outputfile.write("#maximum number of stall generations\n")	
        outputfile.write("outerloop_stallmax " + str(self.outerloop_stallmax) + "\n")
        outputfile.write("#fitness tolerance for the outer-loop\n")	
        outputfile.write("outerloop_tolfit " + str(self.outerloop_tolfit) + "\n")
        outputfile.write("#how many elite individuals to retain\n")	
        outputfile.write("outerloop_elitecount " + str(self.outerloop_elitecount) + "\n")
        outputfile.write("#how many times to run the outer-loop\n")	
        outputfile.write("outerloop_ntrials " + str(self.outerloop_ntrials) + "\n")
        outputfile.write("#whether or not to use the parallel outer-loop\n")	
        outputfile.write("outerloop_useparallel " + str(self.outerloop_useparallel) + "\n")
        outputfile.write("#whether or not to perform an outer loop warm start\n")
        outputfile.write("outerloop_warmstart " + str(self.outerloop_warmstart) + "\n")
        outputfile.write("#Population file for outerloop warm start (set to none if not warm starting)\n")
        outputfile.write("outerloop_warm_population " + self.outerloop_warm_population + "\n")
        outputfile.write("#Archive file for outerloop warm start (set to none if not warm starting)\n")
        outputfile.write("outerloop_warm_archive " + self.outerloop_warm_archive + "\n")
        outputfile.write("#Re-evaluate the entire outerloop each generation? Otherwise read from the archive.\n")
        outputfile.write("outerloop_reevaluate_full_population " + str(self.outerloop_reevaluate_full_population) + "\n")
        outputfile.write("#Quiet outer-loop?\n")
        outputfile.write("quiet_outerloop " + str(self.quiet_outerloop) + "\n")
        outputfile.write("\n")

        outputfile.write("##inner-loop solver settings\n")
        outputfile.write("#NLP solver type\n")
        outputfile.write("#0: SNOPT\n")
        outputfile.write("#1: WORHP\n")
        outputfile.write("NLP_solver_type " + str(self.NLP_solver_type) + "\n")
        outputfile.write("#NLP solver mode\n")
        outputfile.write("#0: find feasible point only\n")
        outputfile.write("#1: find optimal solution\n")
        outputfile.write("NLP_solver_mode " + str(self.NLP_solver_mode) + "\n")
        outputfile.write("#Quiet NLP solver?\n")
        outputfile.write("quiet_NLP " + str(self.quiet_NLP) + "\n")
        outputfile.write("#Quiet MBH?\n")
        outputfile.write("quiet_basinhopping " + str(self.quiet_basinhopping) + "\n")
        outputfile.write("#Enable ACE feasible point finder?\n")
        outputfile.write("ACE_feasible_point_finder " + str(self.ACE_feasible_point_finder) + "\n")
        outputfile.write("#quantity Max_not_improve for MBH\n")
        outputfile.write("MBH_max_not_improve " + str(self.MBH_max_not_improve) + "\n")
        outputfile.write("#maximum number of trials for MBH\n")
        outputfile.write("MBH_max_trials " + str(self.MBH_max_trials) + "\n")
        outputfile.write("#maximum run time for MBH, in seconds\n")
        outputfile.write("MBH_max_run_time " + str(self.MBH_max_run_time) + "\n")
        outputfile.write("#maximum step size for uniform MBH, or scaling factor for Cauchy MBH\n")
        outputfile.write("MBH_max_step_size " + str(self.MBH_max_step_size) + "\n")
        outputfile.write("#MBH hop probabilty distribution\n")
        outputfile.write("#0: uniform\n")
        outputfile.write("#1: Cauchy\n")
        outputfile.write("#2: Pareto\n")
        outputfile.write("#3: Gaussian\n")
        outputfile.write("MBH_hop_distribution " +str(self.MBH_hop_distribution) + "\n")
        outputfile.write("#Pareto distribution alpha\n")
        outputfile.write("MBH_Pareto_alpha " + str(self.MBH_Pareto_alpha) + "\n");
        outputfile.write("#probability of MBH time hop operation\n")
        outputfile.write("MBH_time_hop_probability " + str(self.MBH_time_hop_probability) + "\n")
        outputfile.write("#feasibility tolerance\n")
        outputfile.write("snopt_feasibility_tolerance " + str(self.snopt_feasibility_tolerance) + "\n")
        outputfile.write("#maximum number of major iterations for SNOPT\n")
        outputfile.write("snopt_major_iterations " + str(self.snopt_major_iterations) + "\n")
        outputfile.write("#Maximum run time, in seconds, for a single call to SNOPT\n")
        outputfile.write("snopt_max_run_time " + str(self.snopt_max_run_time) + "\n")
        outputfile.write("#method of specifying derivatives\n")
        outputfile.write("#0: finite difference\n")
        outputfile.write("#1: analytical flybys and objective function but finite difference the patch points\n")
        outputfile.write("#2: all but time derivatives\n")
        outputfile.write("#3: all but current phase flight time derivatives\n")
        outputfile.write("#4: fully analytical (experimental)\n")
        outputfile.write("derivative_type " + str(self.derivative_type) + "\n")
        outputfile.write("#Will MBH be seeded with an initial point? Otherwise MBH starts from a completely random point.\n")
        outputfile.write("seed_MBH " + str(self.seed_MBH) + "\n")
        outputfile.write("#Will the initial guess be interpolated?\n")
        outputfile.write("#(i.e. are we solving a problem with a different number of time steps than the initial guess?)\n")
        outputfile.write("interpolate_initial_guess " + str(self.interpolate_initial_guess) + "\n")
        outputfile.write("#How many time steps were used to create the initial guess?\n")
        outputfile.write("initial_guess_num_timesteps " + str(self.initial_guess_num_timesteps) + "\n")
        outputfile.write("#Distribution from which the initial guess step sizes were drawn\n")
        outputfile.write("#0: uniform\n")
        outputfile.write("#1: Gaussian\n")
        outputfile.write("#2: Cauchy\n")
        outputfile.write("initial_guess_step_size_distribution " + str(self.initial_guess_step_size_distribution) + "\n")
        outputfile.write("#What scale width (Cauchy) or standard deviation (Gaussian) was used to create the step sizes in the initial guess\n")
        outputfile.write("initial_guess_step_size_stdv_or_scale " + str(self.initial_guess_step_size_stdv_or_scale) + "\n")
        outputfile.write("#Apply zero-control initial guess in MBH?\n")
        outputfile.write("#0: do not use\n")
        outputfile.write("#1: zero-control for resets, random perturbations for hops\n")
        outputfile.write("#2: always use zero-control guess except when seeded\n")
        outputfile.write("MBH_zero_control_initial_guess " + str(self.MBH_zero_control_initial_guess) + "\n")
        outputfile.write("#Enable two-step MBH?\n")
        outputfile.write("MBH_two_step " + str(self.MBH_two_step) + "\n")
        outputfile.write("#'Fine' finite differencing step size\n")
        outputfile.write("FD_stepsize " + str(self.FD_stepsize) + "\n")
        outputfile.write("#'Coarse' finite differencing step size\n")
        outputfile.write("FD_stepsize_coarse " + str(self.FD_stepsize_coarse) + "\n")
        outputfile.write("\n")

        outputfile.write("##low-thrust solver parameters\n")	
        outputfile.write("#number of time steps per phase\n")	
        outputfile.write("num_timesteps " + str(self.num_timesteps) + "\n")
        outputfile.write("#Control coordinate system\n")
        outputfile.write("#0: Cartesian\n")
        outputfile.write("#1: Polar\n")
        outputfile.write("control_coordinate_system " + str(self.control_coordinate_system) + '\n')
        outputfile.write("#Initial guess control coordinate system\n")
        outputfile.write("#0: Cartesian\n")
        outputfile.write("#1: Polar\n")
        outputfile.write("initial_guess_control_coordinate_system " + str(self.initial_guess_control_coordinate_system) + '\n')
        outputfile.write("#Distribution from which to draw the step sizes for each phase\n")
        outputfile.write("#0: uniform\n")
        outputfile.write("#1: Gaussian\n")
        outputfile.write("#2: Cauchy\n")
        outputfile.write("step_size_distribution " + str(self.step_size_distribution) + "\n")
        outputfile.write("#What scale width (Cauchy) or standard deviation (Gaussian) is used to create the step sizes\n")
        outputfile.write("step_size_stdv_or_scale " + str(self.step_size_stdv_or_scale) + "\n")
        outputfile.write("#Spiral model type\n")
        outputfile.write("#0: Battin\n")
        outputfile.write("#1: Edelbaum\n")
        outputfile.write("spiral_model_type " + str(self.spiral_model_type) + "\n")
        outputfile.write("\n")

        outputfile.write("##impulsive-thrust solver parameters\n")
        outputfile.write("#maximum number of revolutions for Lambert's method\n")
        outputfile.write("maximum_number_of_lambert_revolutions " + str(self.maximum_number_of_lambert_revolutions) + "\n")
        outputfile.write("\n")

        outputfile.write("##ephemeris data\n")	
        outputfile.write("#ephemeris source\n")	
        outputfile.write("#0: static\n")	
        outputfile.write("#1: SPICE (default to static if no SPICE file supplied for a body)\n")	
        outputfile.write("ephemeris_source " + str(self.ephemeris_source) + "\n")
        outputfile.write("#Universe folder\n")
        outputfile.write("universe_folder " + str(self.universe_folder) + "\n")
        outputfile.write("#SPICE leap seconds kernel - required for SPICE to work\n")
        outputfile.write("SPICE_leap_seconds_kernel " + str(self.SPICE_leap_seconds_kernel) + "\n")
        outputfile.write("#SPICE_reference_frame_kernel\n")
        outputfile.write("SPICE_reference_frame_kernel " + str(self.SPICE_reference_frame_kernel) + "\n")
        outputfile.write("\n")

        outputfile.write("##lambert solver options\n")
        outputfile.write("#Lambert solver choice\n")
        outputfile.write("#0: Arora-Russell\n")
        outputfile.write("#1: Izzo (not included in open-source package)\n")
        outputfile.write("LambertSolver " + str(self.LambertSolver) + "\n")
        outputfile.write("\n")
            
        outputfile.write("##vehicle parameters\n")	
        outputfile.write("#the maximum possible mass in kg of the spacecraft (negative number means use LV max)\n")	
        outputfile.write("maximum_mass " + str(self.maximum_mass) + "\n")
        outputfile.write("#Launch vehicle type\n")
        outputfile.write("#-2: custom launch vehicle\n")
        outputfile.write("#-1: burn with departure stage engine\n")
        outputfile.write("#0: fixed initial mass\n")
        outputfile.write("#1: Atlas V (401)	 NLSII\n")
        outputfile.write("#2: Atlas V (411)	NLSII\n")
        outputfile.write("#3: Atlas V (421)	NLSII\n")
        outputfile.write("#4: Atlas V (431)	NLSII\n")
        outputfile.write("#5: Atlas V (501)	NLSII\n")
        outputfile.write("#6: Atlas V (511)	NLSII\n")
        outputfile.write("#7: Atlas V (521)	NLSII\n")
        outputfile.write("#8: Atlas V (531)	NLSII\n")
        outputfile.write("#9: Atlas V (541)	NLSII\n")
        outputfile.write("#10: Atlas V (551)	NLSII\n")
        outputfile.write("#11: Falcon 9 (v1.0)	NLSII\n")
        outputfile.write("#12: Falcon 9 (v1.1)	NLSII\n")
        outputfile.write("#13: Atlas V (551) w/Star 48	NLSI\n")
        outputfile.write("#14: Falcon 9 Heavy\n")
        outputfile.write("#15: Delta IV Heavy	NLSI\n")
        outputfile.write("#16: SLS Block 1\n")
        outputfile.write("LV_type " + str(self.LV_type) + "\n")
        outputfile.write("#Launch vehicle margin (0.0 - 1.0)\n")
        outputfile.write("LV_margin " + str(self.LV_margin) + "\n")
        outputfile.write("#Launch vehicle adapter mass (kg)\n")
        outputfile.write("LV_adapter_mass " + str(self.LV_adapter_mass) + "\n")
        outputfile.write("#Custom launch vehicle coefficients (must enter 6 coefficients)\n")
        outputfile.write("#as in a1*C3^5 + a2*C3^4 + a3*C3^3 + a4*C3^2 + a5*C3 + a6\n")
        outputfile.write("custom_LV_coefficients")
        for k in range(0,6):
            outputfile.write(" " + str(self.custom_LV_coefficients[k]))
        outputfile.write("\n")
        outputfile.write("#Custom launch vehicle C3 bounds (two values)\n")
        outputfile.write("custom_LV_C3_bounds " + str(self.custom_LV_C3_bounds[0]) + " " + str(self.custom_LV_C3_bounds[1]) + "\n")
        outputfile.write("#Parking orbit inclination (for use in outputing GMAT scenarios)\n")
        outputfile.write("parking_orbit_inclination " + str(self.parking_orbit_inclination) + "\n")
        outputfile.write("#Parking orbit altitude (for use in outputing GMAT scenarios)\n")
        outputfile.write("parking_orbit_altitude " + str(self.parking_orbit_altitude) + "\n")
        outputfile.write("\n")

        outputfile.write("##parameters that are only relevant for missions that use chemical propulsion\n")	
        outputfile.write("##dummy values should be used if the mission does not use chemical propulsion but are not strictly necessary\n")	
        outputfile.write("#specific impulse in seconds of the engine used for impulsive maneuvers\n")	
        outputfile.write("IspChem " + str(self.IspChem) + "\n")
        outputfile.write("\n")
            
        outputfile.write("##parameters that are only relevant for missions that use a chemical EDS\n")	
        outputfile.write("##dummy values should be used if the mission does not use a chemical EDS but are not strictly necessary\n")	
        outputfile.write("#specific impulse in seconds for the earth departure stage, if applicable\n")	
        outputfile.write("IspDS " + str(self.IspDS) + "\n")
        outputfile.write("\n")
            
        outputfile.write("##parameters that are only relevant for missions that use low-thrust\n")	
        outputfile.write("##dummy values should be used if the mission does not use low-thrust but are not strictly necessary\n")	
        outputfile.write("#specific impulse in seconds of the engine used for low-thrust maneuvers\n")	
        outputfile.write("#for VSI systems, this represents maximum Isp\n")
        outputfile.write("IspLT " + str(self.IspLT) + "\n")
        outputfile.write("#minimum Isp for VSI systems\n")
        outputfile.write("IspLT_minimum " + str(self.IspLT_minimum) + "\n")
        outputfile.write("#thrust of the spacecraft low-thrust motor, in Newtons\n")	
        outputfile.write("Thrust " + str(self.Thrust) + "\n")
        outputfile.write("#low-thrust engine type\n")
        outputfile.write("#0: fixed thrust/Isp\n")
        outputfile.write("#1: constant Isp, efficiency, EMTG computes input power\n")
        outputfile.write("#2: choice of power model, constant efficiency, EMTG chooses Isp\n")
        outputfile.write("#3: choice of power model, constant efficiency and Isp\n")
        outputfile.write("#4: continuously-varying specific impulse\n")
        outputfile.write("#5: custom thrust and mass flow rate polynomial\n")
        outputfile.write("#6: NSTAR\n")
        outputfile.write("#7: XIPS-25\n")
        outputfile.write("#8: BPT-4000 High-Isp\n")
        outputfile.write("#9: BPT-4000 High-Thrust\n")
        outputfile.write("#10: BPT-4000 Ex-High-Isp\n")
        outputfile.write("#11: NEXT high-Isp v9\n")
        outputfile.write("#12: VASIMR (argon, using analytical model)\n")
        outputfile.write("#13: Hall Thruster (Xenon, using analytical model)\n")
        outputfile.write("#14: NEXT high-ISP v10\n")
        outputfile.write("#15: NEXT high-thrust v10\n")
        outputfile.write("#16: BPT-4000 MALTO\n")
        outputfile.write("#17: NEXIS Cardiff 8-15-201\n")
        outputfile.write("#18: H6MS Cardiff 8-15-2013\n")
        outputfile.write("#19: BHT20K Cardiff 8-16-2013\n")
        outputfile.write("#20: Aerojet HiVHAC EM\n")
        outputfile.write("#21: 13 kW STMD Hall high-Isp (not available in open-source)\n")
        outputfile.write("#22: 13 kW STMD Hall high-thrust (not available in open-source)\n")
        outputfile.write("#23: NEXT TT11 High-Thrust\n")
        outputfile.write("#24: NEXT TT11 High-Isp\n")
        outputfile.write("#25: NEXT TT11 Expanded Throttle Table\n")
        outputfile.write("#26: 13 kW STMD Hall high-Isp 9-8-2014 (not available in open-source)\n")
        outputfile.write("#27: 13 kW STMD Hall medium-thrust 9-8-2014 (not available in open-source)\n")
        outputfile.write("#28: 13 kW STMD Hall high-thrust 9-8-2014 (not available in open-source)\n")
        outputfile.write("engine_type " + str(self.engine_type) + "\n")
        outputfile.write("#Custom engine thrust coefficients (T = A + BP + C*P^2 + D*P^3 + E*P^4 + G*P^5 + H*P^6)\n")
        outputfile.write("engine_input_thrust_coefficients")
        for k in range(0,7):
            outputfile.write(" " + str(self.engine_input_thrust_coefficients[k]))
        outputfile.write("\n")
        outputfile.write("#Custom engine mass flow rate coefficients (mdot = A + BP + C*P^2 + D*P^3 + E*P^4 + G*P^5 + H*P^6)\n")
        outputfile.write("engine_input_mass_flow_rate_coefficients")
        for k in range(0,7):
            outputfile.write(" " + str(self.engine_input_mass_flow_rate_coefficients[k]))
        outputfile.write("\n")
        outputfile.write("#Custom engine lower and upper bounds on input power (per engine, in kW)\n")
        outputfile.write("engine_input_power_bounds " + str(self.engine_input_power_bounds[0]) + " " + str(self.engine_input_power_bounds[1]) + "\n")
        outputfile.write("#Custom engine input efficiency\n")
        outputfile.write("user_defined_engine_efficiency " + str(self.user_defined_engine_efficiency) + "\n")
        outputfile.write("#number of low-thrust engines\n")
        outputfile.write("number_of_engines " + str(self.number_of_engines) + "\n")
        outputfile.write("#Throttle logic mode\n")
        outputfile.write("#0: maximum power use\n")
        outputfile.write("#1: maximum thrust\n")
        outputfile.write("#2: maximum Isp\n")
        outputfile.write("#3: maximum efficiency\n")
        outputfile.write("throttle_logic_mode " + str(self.throttle_logic_mode) + "\n")
        outputfile.write("#Throttle sharpness (higher means more precise, lower means smoother)\n")
        outputfile.write("throttle_sharpness " + str(self.throttle_sharpness) + "\n")
        outputfile.write("#engine duty cycle [0,1]\n")
        outputfile.write("engine_duty_cycle " + str(self.engine_duty_cycle) + "\n")
        outputfile.write("#electrical power available at 1 AU (kW)\n")
        outputfile.write("power_at_1_AU " + str(self.power_at_1_AU) + "\n")
        outputfile.write("#power source type, 0: solar, 1: radioisotope (or other fixed power)\n")
        outputfile.write("power_source_type " + str(self.power_source_type) + "\n")
        outputfile.write("#solar power coefficients gamma_1 through gamma_5\n")
        outputfile.write("#if all gamma = [1 0 0 0 0], then solar power is a simple 1/r^2\n")
        outputfile.write("solar_power_gamma")
        for k in range(0,5):
            outputfile.write(" " + str(self.solar_power_gamma[k]))
        outputfile.write("\n")
        outputfile.write("#Power margin (for thrusters, as a fraction)\n")
        outputfile.write("power_margin " + str(self.power_margin) + "\n")
        outputfile.write("#Power system decay rate (per year)\n")
        outputfile.write("power_decay_rate " + str(self.power_decay_rate) + "\n")
        outputfile.write("#spacecraft power coefficients A, B, and C\n")
        outputfile.write("#represent the power requirements of the spacecraft at a distance r from the sun\n")
        outputfile.write("#i.e. heaters, communications, etc\n")
        outputfile.write("spacecraft_power_coefficients " + str(self.spacecraft_power_coefficients[0]) + " " + str(self.spacecraft_power_coefficients[1]) + " " + str(self.spacecraft_power_coefficients[2]) + "\n")
        outputfile.write("#spacecraft power model type\n")
        outputfile.write("#0: P_sc = A + B/r + C/r^2\n")
        outputfile.write("#1: P_sc = A if P > A, A + B(C - P) otherwise\n")
        outputfile.write("spacecraft_power_model_type " + str(self.spacecraft_power_model_type) + "\n")
        outputfile.write("#low-thrust propulsion stage dry mass in kg, will be subtracted before chemical arrival or mid-flight switchover to chemical propulsion\n")	
        outputfile.write("EP_dry_mass " + str(self.EP_dry_mass) + "\n")
        outputfile.write("#Allow initial mass to vary, up to maximum possible mass? (only relevant for MGALT and FBLT)\n")
        outputfile.write("allow_initial_mass_to_vary " + str(self.allow_initial_mass_to_vary) + "\n")
        outputfile.write("#Minimum dry mass\n")
        outputfile.write("minimum_dry_mass " + str(self.minimum_dry_mass) + "\n")
        outputfile.write("#Enable maximum propellant mass constraint?\n")
        outputfile.write("enable_maximum_propellant_mass_constraint " + str(int(self.enable_maximum_propellant_mass_constraint)) + "\n")
        outputfile.write("#Maximum propellant mass (kg)\n")
        outputfile.write("maximum_propellant_mass " + str(self.maximum_propellant_mass) + "\n")
        outputfile.write("#Post-mission delta-v, in km/s (alternatively defined as delta-v margin)\n")
        outputfile.write("post_mission_delta_v " + str(self.post_mission_delta_v) + "\n")
        outputfile.write("#Isp used to compute propellant for post-mission delta-v, in seconds\n")
        outputfile.write("post_mission_Isp " + str(self.post_mission_Isp) + "\n")
        outputfile.write("#Propellant margin, as a fraction of nominal propellant load\n")
        outputfile.write("propellant_margin " + str(self.propellant_margin) + "\n")
        outputfile.write("\n")

        outputfile.write("##Global problem settings\n")
        outputfile.write("#mission name\n")
        outputfile.write("mission_name " + str(self.mission_name) + "\n")
        outputfile.write("#mission type - you can specify MGA, MGA-DSM, MGA-LT, or allow the outer-loop to choose\n")	
        outputfile.write("#0: MGA\n")
        outputfile.write("#1: MGA-DSM\n")
        outputfile.write("#2: MGA-LT\n")
        outputfile.write("#3: FBLT\n")
        outputfile.write("#4: MGA-NDSM\n")
        outputfile.write("#5: DTLT\n")
        outputfile.write("#6: solver chooses (MGA, MGA-DSM)\n")
        outputfile.write("#7: solver chooses (MGA, MGA-LT)\n")
        outputfile.write("#8: solver chooses (MGA-DSM, MGA-LT)\n")
        outputfile.write("#9: solver chooses (MGA, MGA-DSM, MGA-LT)\n")
        outputfile.write("mission_type " + str(self.mission_type) + "\n")
        outputfile.write("#number of journeys (user-defined endpoints)\n")
        outputfile.write("#Each journey has a central body and two boundary points\n")
        outputfile.write("#Each central body has a menu of destinations which is used to choose the boundary points. Every menu is structured:\n")
        outputfile.write("#-1: Boundary at a point in space, either fixed or free\n")
        outputfile.write("#0: Nothing happens. This code is only used to signify ""no flyby"" and should NEVER be coded as a destination.\n")
        outputfile.write("#1: Body 1 (i.e. Mercury, Io, etc)\n")
        outputfile.write("#2: Body 2 (i.e. Venus, Europa, etc)\n")
        outputfile.write("#...\n")
        outputfile.write("#N: Body N\n")
        outputfile.write("number_of_journeys " + str(self.number_of_journeys) + "\n")
        outputfile.write("#maximum number of phases allowed per journey\n")
        outputfile.write("max_phases_per_journey " + str(self.max_phases_per_journey) + "\n")
        outputfile.write("#destination list (number of journeys + 1)\n")
        outputfile.write("destination_list")
        for j in range(0, self.number_of_journeys):
            if j > 0:
                if self.Journeys[j].destination_list[0] != self.Journeys[j-1].destination_list[1]:
                    print "Second entry in destination list for Journey ", j, " does not match first entry in destination list for Journey ", j-1, ". Are you sure that you want to do that?"
            outputfile.write(" " + str(self.Journeys[j].destination_list[0]) + " " + str(self.Journeys[j].destination_list[1]))
        outputfile.write("\n")
        outputfile.write("#the following option is relevant only if optimizing over total deltaV, should the initial impulse be included in the cost?\n")
        outputfile.write("include_initial_impulse_in_cost " + str(self.include_initial_impulse_in_cost) + "\n")
        outputfile.write("#global time bounds\n")
        outputfile.write("#0: unbounded\n")
        outputfile.write("#1: bounded total time (note that the global arrival date bound is by definition the same as the last journey arrival date bound and is not duplicated\n")
        outputfile.write("global_timebounded " + str(self.global_timebounded) + "\n")
        outputfile.write("#MJD of the opening of the launch window\n")
        outputfile.write("launch_window_open_date " + str(self.launch_window_open_date) + "\n")
        outputfile.write("#total flight time bounds, in days\n")
        outputfile.write("total_flight_time_bounds " + str(self.total_flight_time_bounds[0]) + " " + str(self.total_flight_time_bounds[1]) + "\n")
        outputfile.write("#objective function type\n")
        outputfile.write("#0: minimum deltaV\n")
        outputfile.write("#1: minimum time\n")
        outputfile.write("#2: maximum final mass\n")	
        outputfile.write("#3: GTOC 1 asteroid deflection function\n")
        outputfile.write("#4: launch as late as possible in the window\n")
        outputfile.write("#5: launch as early as possible in the window\n")
        outputfile.write("#6: maximize orbit energy\n")
        outputfile.write("#7: minimize launch mass\n")
        outputfile.write("#8: arrive as early as possible\n")
        outputfile.write("#9: arrive as late as possible\n")
        outputfile.write("#10: minimum propellant (not the same as #2)\n")
        outputfile.write("#11: maximum dry/wet ratio\n")
        outputfile.write("#12: maximum arrival kinetic energy\n")
        outputfile.write("#13: minimum BOL power\n")
        outputfile.write("objective_type " + str(self.objective_type) + "\n")
        outputfile.write("#bounds on the DLA, in degrees (typically set to declination of your launch site)\n")	
        outputfile.write("DLA_bounds " + str(self.DLA_bounds[0]) + " " + str(self.DLA_bounds[1]) + "\n")
        outputfile.write("\n")
        outputfile.write("#Initial V-Infinity vector (set to zeros unless starting the mission from periapse of a hyperbolic arrival)\n")	
        outputfile.write("initial_V_infinity " + str(self.initial_V_infinity[0]) + " " + str(self.initial_V_infinity[1]) + " " + str(self.initial_V_infinity[2]) + "\n")
        outputfile.write("#Forced post-launch coast (in days, to be enforced after launch)\n")
        outputfile.write("forced_post_launch_coast " + str(self.forced_post_launch_coast) + "\n")
        outputfile.write("#Forced post flyby/intercept coast (in days, to be enforced before/after each flyby/intercept)\n")
        outputfile.write("forced_flyby_coast " + str(self.forced_flyby_coast) + "\n")
        outputfile.write("\n")

        outputfile.write("##Settings for each journey\n")	
        outputfile.write("##dummy values should be used - they should not be necessary but testing was not exhaustive so please use them\n")
        outputfile.write("#journey names\n")
        outputfile.write("journey_names")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_names))
        outputfile.write("\n")
        outputfile.write("#How much mass to add to the spacecraft at the beginning of the journey (a negative number indicates a mass drop)\n")
        outputfile.write("journey_starting_mass_increment")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_starting_mass_increment))
        outputfile.write("\n")
        outputfile.write("#Is the mass increment variable (i.e. can the optimizer choose how much mass to add)\n")
        outputfile.write("#This option is ignored for journeys with zero or negative mass increment\n")
        outputfile.write("journey_variable_mass_increment")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_variable_mass_increment))
        outputfile.write("\n")
        outputfile.write("#is each journey time bounded (one value per journey)\n")
        outputfile.write("#0: unbounded\n")
        outputfile.write("#1: bounded flight time\n")
        outputfile.write("#2: bounded arrival date\n")
        outputfile.write("#3: bounded aggregate flight time\n")
        outputfile.write("journey_timebounded")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_timebounded))
        outputfile.write("\n")
        outputfile.write("#what are the wait time lower and upper bounds, in days, for each journey (two numbers per journey)\n")	
        outputfile.write("journey_wait_time_bounds")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_wait_time_bounds[0]) + " " + str(self.Journeys[j].journey_wait_time_bounds[1]))
        outputfile.write("\n")
        outputfile.write("#what are the flight time bounds for each journey (two numbers per journey, use dummy values if no flight time bounds)\n")	
        outputfile.write("journey_flight_time_bounds")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_flight_time_bounds[0]) + " " + str(self.Journeys[j].journey_flight_time_bounds[1]))
        outputfile.write("\n")
        outputfile.write("#what are the arrival date bounds for each journey (two numbers per journey, use dummy values if no flight time bounds)\n")	
        outputfile.write("journey_arrival_date_bounds")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_arrival_date_bounds[0]) + " " + str(self.Journeys[j].journey_arrival_date_bounds[1]))
        outputfile.write("\n")
        outputfile.write("#what are the bounds on the initial impulse for each journey in km/s (two numbers per journey)\n")	
        outputfile.write("#you can set a very high upper bound if you are using a launchy vehicle model - the optimizer will find the correct value\n")	
        outputfile.write("journey_initial_impulse_bounds")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_initial_impulse_bounds[0]) + " " + str(self.Journeys[j].journey_initial_impulse_bounds[1]))
        outputfile.write("\n")
        outputfile.write("#journey departure type (one value per journey)\n")	
        outputfile.write("#0: launch or direct insertion\n")
        outputfile.write("#1: depart from parking orbit (you can use this one in place of a launch vehicle model, and the departure burn will be done with the EDS motor)\n")
        outputfile.write("#2: free direct departure, i.e. do not burn to get the departure v_infinity (used for when operations about a small body are not modeled but the departure velocity is known)\n")
        outputfile.write("#3: flyby (only valid for successive journeys)\n")
        outputfile.write("#4: flyby with fixed v-infinity-out (only valid for successive journeys)\n")
        outputfile.write("#5: spiral-out from circular orbit (low-thrust missions only)\n")
        outputfile.write("#6: zero-turn flyby (for small bodies)\n")
        outputfile.write("journey_departure_type")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_departure_type))
        outputfile.write("\n")
        outputfile.write("#journey arrival type (one value per journey)\n")	
        outputfile.write("#0: insertion into parking orbit (use chemical Isp)\n")
        outputfile.write("#1: rendezvous (use chemical Isp)\n")
        outputfile.write("#2: intercept with bounded V_infinity\n")
        outputfile.write("#3: low-thrust rendezvous (does not work if terminal phase is not low-thrust)\n")
        outputfile.write("#4: match final v-infinity vector\n")
        outputfile.write("#5: match final v-infinity vector (low-thrust)\n")
        outputfile.write("#6: escape (E = 0)\n")
        outputfile.write("#7: capture spiral\n")
        outputfile.write("journey_arrival_type")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_arrival_type))
        outputfile.write("\n")
        outputfile.write("#type of orbit elements specified at beginning of journey(0: inertial, 1: COE)\n")
        outputfile.write("journey_departure_elements_type")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_departure_elements_type))
        outputfile.write("\n")
        outputfile.write("#orbit elements at beginning of journey (a(km), e, i, RAAN, AOP, MA) supply angles in degrees\n")
        outputfile.write("journey_departure_elements")
        for j in range (0, self.number_of_journeys):
            for k in range(0, 6):
                outputfile.write(" " + str(self.Journeys[j].journey_departure_elements[k]))
        outputfile.write("\n")
        outputfile.write("#Vary journey departure elements? (one entry per element per journey: 0 means no, 1 means yes)\n")
        outputfile.write("journey_departure_elements_vary_flag")
        for j in range (0, self.number_of_journeys):
            for k in range(0, 6):
                outputfile.write(" " + str(self.Journeys[j].journey_departure_elements_vary_flag[k]))
        outputfile.write("\n")

        outputfile.write("#Lower and upper bounds on journey departure elements (two per element per journey, ignored if vary flag is off for that element)\n")
        outputfile.write("journey_departure_elements_bounds")
        for j in range (0, self.number_of_journeys):
            for k in range(0, 12):
                outputfile.write(" " + str(self.Journeys[j].journey_departure_elements_bounds[k]))
        outputfile.write("\n")    
        outputfile.write("journey_arrival_elements_type")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_arrival_elements_type))
        outputfile.write("\n")
        outputfile.write("#orbit elements at end of journey (a(km), e, i, RAAN, AOP, MA) supply angles in degrees\n")
        outputfile.write("journey_arrival_elements")
        for j in range (0, self.number_of_journeys):
            for k in range(0, 6):
                outputfile.write(" " + str(self.Journeys[j].journey_arrival_elements[k]))
        outputfile.write("\n")
        outputfile.write("#Vary journey arrival elements? (one entry per element per journey: 0 means no, 1 means yes)\n")
        outputfile.write("journey_arrival_elements_vary_flag")
        for j in range (0, self.number_of_journeys):
            for k in range(0, 6):
                outputfile.write(" " + str(self.Journeys[j].journey_arrival_elements_vary_flag[k]))
        outputfile.write("\n")

        outputfile.write("#Lower and upper bounds on journey arrival elements (two per element per journey, ignored if vary flag is off for that element)\n")
        outputfile.write("journey_arrival_elements_bounds")
        for j in range (0, self.number_of_journeys):
            for k in range(0, 12):
                outputfile.write(" " + str(self.Journeys[j].journey_arrival_elements_bounds[k]))
        outputfile.write("\n")
        outputfile.write("#journey central body\n")
        outputfile.write("#use SPICE names, as per http://www-int.stsci.edu/~sontag/spicedocs/req/naif_ids.html\n")
        outputfile.write("journey_central_body")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_central_body))
        outputfile.write("\n")
        outputfile.write("#initial VHP for journeys that begin with flybys, in km/s (three numbers per journey)\n")
        outputfile.write("journey_initial_velocity")
        for j in range(0, self.number_of_journeys):
            for k in range(0,3):
                outputfile.write(" " + str(self.Journeys[j].journey_initial_velocity[k]))
        outputfile.write("\n")
        outputfile.write("#final VHP for journeys that end in intercepts, in km/s (three numbers per journey)\n")
        outputfile.write("journey_final_velocity")
        for j in range(0, self.number_of_journeys):
            for k in range(0,3):
                outputfile.write(" " + str(self.Journeys[j].journey_final_velocity[k]))
        outputfile.write("\n")
        outputfile.write("#Impose arrival declination constraint on each journey?\n")
        outputfile.write("#0: no\n")
        outputfile.write("#1: yes\n")
        outputfile.write("journey_arrival_declination_constraint_flag")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_arrival_declination_constraint_flag))
        outputfile.write("\n")
        outputfile.write("#Arrival declination bounds for each journey\n")
        outputfile.write("#Two numbers per journey, in degrees\n")
        outputfile.write("journey_arrival_declination_bounds")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_arrival_declination_bounds[0]) + " " + str(self.Journeys[j].journey_arrival_declination_bounds[1]))
        outputfile.write("\n")
        outputfile.write("#Starting orbital radius for an escape spiral at the beginning of the journey\n")
        outputfile.write("journey_escape_spiral_starting_radius")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_escape_spiral_starting_radius))
        outputfile.write("\n")
        outputfile.write("#Final orbital radius for a capture spiral at the end of the journey\n")
        outputfile.write("journey_capture_spiral_final_radius")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_capture_spiral_final_radius))
        outputfile.write("\n")
        outputfile.write("#Enable journey maximum DSM magnitude constraint?\n")
        outputfile.write("journey_maximum_DSM_magnitude_constraint_flag")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_maximum_DSM_magnitude_constraint_flag))
        outputfile.write("\n")
        outputfile.write("#Journey maximum DSM magnitude (km/s)\n")
        outputfile.write("journey_maximum_DSM_magnitude_constraint")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_maximum_DSM_magnitude_constraint))
        outputfile.write("\n")
        outputfile.write("\n")
            
        outputfile.write("##Perturbation settings\n")
        outputfile.write("#Enable solar radiation pressure?\n")
        outputfile.write("perturb_SRP " + str(self.perturb_SRP) + "\n")
        outputfile.write("#Enable third-body perturbations?\n")
        outputfile.write("perturb_thirdbody " + str(self.perturb_thirdbody) + "\n")
        outputfile.write("#Journey perturbation bodies. One line per journey. The numbers in the line correspond to\n")
        outputfile.write("#bodies in the journey""s Universe file. If perturbations are off, each line should just have a zero\n")
        outputfile.write("#the numbers in the first line are the number of perturbation bodies for each journey\n")
        outputfile.write("journey_perturbation_bodies")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_number_of_perturbation_bodies))
        outputfile.write("\n")
        for j in range(0, self.number_of_journeys):
            if self.Journeys[j].journey_number_of_perturbation_bodies > 0:
                outputfile.write(str(self.Journeys[j].journey_perturbation_bodies[0]))
                for b in range(1, self.Journeys[j].journey_number_of_perturbation_bodies):
                    outputfile.write(" " + str(self.Journeys[j].journey_perturbation_bodies[b]))
            else:
                outputfile.write("0")
            outputfile.write("\n")
        outputfile.write("#end_journey_perturbation_bodies\n")
        outputfile.write("#Spacecraft area (in m^2)\n")
        outputfile.write("spacecraft_area " + str(self.spacecraft_area) + "\n")
        outputfile.write("#Coefficient of reflectivity\n")
        outputfile.write("#0.0: perfectly translucent\n")
        outputfile.write("#1.0: perfectly absorbing\n")
        outputfile.write("#2.0: perfectly reflecting\n")
        outputfile.write("coefficient_of_reflectivity " + str(self.coefficient_of_reflectivity) + "\n")
        outputfile.write("\n")

        outputfile.write("##Outer-loop selectable options settings\n")
        outputfile.write("#Allow outer-loop to vary power level?\n")
        outputfile.write("outerloop_vary_power " + str(self.outerloop_vary_power) + "\n")
        outputfile.write("#Allow outer-loop to vary launch epoch?\n")
        outputfile.write("outerloop_vary_launch_epoch " + str(self.outerloop_vary_launch_epoch) + "\n")
        outputfile.write("#Allow outer-loop to vary flight time upper bound?\n")
        outputfile.write("outerloop_vary_flight_time_upper_bound " + str(self.outerloop_vary_flight_time_upper_bound) + "\n")
        outputfile.write("#Restrict flight-time lower bound when running outer-loop?\n")
        outputfile.write("outerloop_restrict_flight_time_lower_bound " + str(self.outerloop_restrict_flight_time_lower_bound) + "\n")
        outputfile.write("#Allow outer-loop to vary thruster type?\n")
        outputfile.write("outerloop_vary_thruster_type " + str(self.outerloop_vary_thruster_type) + "\n")
        outputfile.write("#Allow outer-loop to vary number of thrusters?\n")
        outputfile.write("outerloop_vary_number_of_thrusters " + str(self.outerloop_vary_number_of_thrusters) + "\n")
        outputfile.write("#Allow outer-loop to vary launch vehicle?\n")
        outputfile.write("outerloop_vary_launch_vehicle " + str(self.outerloop_vary_launch_vehicle) + "\n")
        outputfile.write("#Allow outer-loop to vary first journey departure C3?\n")
        outputfile.write("outerloop_vary_departure_C3 " + str(self.outerloop_vary_departure_C3) + "\n")
        outputfile.write("#Allow outer-loop to vary last journey arrival C3?\n")
        outputfile.write("outerloop_vary_arrival_C3 " + str(self.outerloop_vary_arrival_C3) + "\n")
        outputfile.write("#Allow outer-loop to vary journey destination? (one value per journey)\n")
        outputfile.write("outerloop_vary_journey_destination")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].outerloop_vary_journey_destination))
        outputfile.write("\n")
        outputfile.write("#Allow outer-loop to vary journey flyby sequence? (one value per journey)\n")
        outputfile.write("outerloop_vary_journey_flyby_sequence")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].outerloop_vary_journey_flyby_sequence))
        outputfile.write("\n")
        outputfile.write("#Outer-loop power at 1 AU choices (in kW)\n")
        outputfile.write("outerloop_power_choices")
        for entry in self.outerloop_power_choices:
            outputfile.write(" " + str(entry))
        outputfile.write("\n")
        outputfile.write("#Outer-loop launch window open epoch choices (in MJD)\n")
        outputfile.write("outerloop_launch_epoch_choices")
        for entry in self.outerloop_launch_epoch_choices:
            outputfile.write(" " + str(entry))
        outputfile.write("\n")
        outputfile.write("#Outer-loop flight time upper bound choices (in days)\n")
        outputfile.write("outerloop_flight_time_upper_bound_choices")
        for entry in self.outerloop_flight_time_upper_bound_choices:
            outputfile.write(" " + str(entry))
        outputfile.write("\n")
        outputfile.write("#Outer-loop thruster type choices (in order of most to least preferable)\n")
        outputfile.write("outerloop_thruster_type_choices")
        for entry in self.outerloop_thruster_type_choices:
            outputfile.write(" " + str(int(entry)))
        outputfile.write("\n")
        outputfile.write("#Outer-loop number of thruster choices\n")
        outputfile.write("outerloop_number_of_thrusters_choices")
        for entry in self.outerloop_number_of_thrusters_choices:
            outputfile.write(" " + str(int(entry)))
        outputfile.write("\n")
        outputfile.write("#Outer-loop launch vehicle choices (in order of most to least preferable)\n")
        outputfile.write("outerloop_launch_vehicle_choices")
        for entry in self.outerloop_launch_vehicle_choices:
            outputfile.write(" " + str(int(entry)))
        outputfile.write("\n")
        outputfile.write("#Outer-loop first journey departure C3 choices\n")
        outputfile.write("outerloop_departure_C3_choices")
        for entry in self.outerloop_departure_C3_choices:
            outputfile.write(" " + str(entry))
        outputfile.write("\n")
        outputfile.write("#Outer-loop last arrival departure C3 choices\n")
        outputfile.write("outerloop_arrival_C3_choices")
        for entry in self.outerloop_arrival_C3_choices:
            outputfile.write(" " + str(entry))
        outputfile.write("\n")
        outputfile.write("#Outer-loop maximum number of flybys (one value for each journey)\n")
        outputfile.write("outerloop_journey_maximum_number_of_flybys")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].outerloop_journey_maximum_number_of_flybys))
        outputfile.write("\n")
        outputfile.write("#Outer-loop journey destination choices (one line for each journey)\n")
        outputfile.write("outerloop_journey_destination_choices\n")
        for j in range(0, self.number_of_journeys):
            for entry in self.Journeys[j].outerloop_journey_destination_choices:
                outputfile.write(" " + str(int(entry)))
            outputfile.write("\n")
        outputfile.write("#Outer-loop flyby sequence choices (one line for each journey)\n")
        outputfile.write("outerloop_journey_flyby_sequence_choices\n")
        for j in range(0, self.number_of_journeys):
            for entry in self.Journeys[j].outerloop_journey_flyby_sequence_choices:
                outputfile.write(" " + str(int(entry)))
            outputfile.write("\n")
        outputfile.write("\n")

        outputfile.write("##Outer-loop objective function settings\n")
        outputfile.write("#Pick as many as you want. The Pareto surface will be generated in these dimensions\n")
        outputfile.write("#0: BOL power at 1 AU (kW)\n")
        outputfile.write("#1: Launch epoch (MJD)\n")
        outputfile.write("#2: Flight time (days)\n")
        outputfile.write("#3: Thruster preference\n")
        outputfile.write("#4: Number of thrusters\n")
        outputfile.write("#5: Launch vehicle preference\n")
        outputfile.write("#6: Delivered mass to final target (kg)\n")
        outputfile.write("#7: Final journey mass increment (for maximizing sample return)\n")
        outputfile.write("#8: First journey departure C3 (km^2/s^2)\n")
        outputfile.write("#9: Final journey arrival C3 (km^2/s^2)\n")
        outputfile.write("#10: Total delta-v (km/s)\n")
        outputfile.write("#11: Inner-loop objective (whatever it was)\n")
        outputfile.write("#12: Point-group value\n")
        outputfile.write("outerloop_objective_function_choices")
        for entry in self.outerloop_objective_function_choices:
            outputfile.write(" " + str(entry))
        outputfile.write("\n")
        outputfile.write("\n")

        outputfile.write("##Outer-loop point group settings\n")
        outputfile.write("#Point group values and members\n")
        outputfile.write("outerloop_point_groups_values")
        for g in range(0, len(self.outerloop_point_groups_values)):
            outputfile.write(" " + str(self.outerloop_point_groups_values[g]))
        outputfile.write("\n")
        for g in range(0, len(self.outerloop_point_groups_values)):
            for m in range(0, len(self.outerloop_point_groups_members[g])):
                outputfile.write(" " + str(self.outerloop_point_groups_members[g][m]))
            outputfile.write("\n")
        outputfile.write("#How many members to score from each point group (additional members add no more points)\n")
        outputfile.write("outerloop_point_groups_number_to_score")
        for g in range(0, len(self.outerloop_point_groups_values)):
            outputfile.write(" " + str(self.outerloop_point_groups_number_to_score[g]))
        outputfile.write("\n")
        outputfile.write("\n")
            
        outputfile.write("##output format settings\n")
        outputfile.write("#output units, 0: km and km/s, 1: LU and LU/day\n")
        outputfile.write("output_units " + str(self.output_units) + "\n")
        outputfile.write("#Output a GMAT script (not compatible with non-body boundary conditions or thruster/power models)\n")
        outputfile.write("create_GMAT_script " + str(self.create_GMAT_script) + "\n")
        outputfile.write("#Generate initial guess file?\n")
        outputfile.write("generate_initial_guess_file " + str(self.generate_initial_guess_file) + "\n")
        outputfile.write("#0: no\n")
        outputfile.write("#1: yes\n")
        outputfile.write("#Mission type for initial guess file (experimental!)\n")
        outputfile.write("#(this is a limited-capability feature and many options will not translate properly)\n")
        outputfile.write("#0: MGA\n")
        outputfile.write("#1: MGADSM\n")
        outputfile.write("#2: MGALT\n")
        outputfile.write("#3: FBLT\n")
        outputfile.write("#4: MGANDSM\n")
        outputfile.write("mission_type_for_initial_guess_file " + str(self.mission_type_for_initial_guess_file) + "\n")
        outputfile.write("#Override working directory?\n")
        outputfile.write("#0: no\n")
        outputfile.write("#1: yes\n")
        outputfile.write("override_working_directory " + str(self.override_working_directory) + "\n")
        outputfile.write("#Custom working directory\n")
        outputfile.write("forced_working_directory " + self.forced_working_directory + "\n")
        outputfile.write("#Generate forward integrated ephemeris (STK compatible)?\n")
        outputfile.write("#0: no\n")
        outputfile.write("#1: yes\n")
        outputfile.write("generate_forward_integrated_ephemeris " + str(self.generate_forward_integrated_ephemeris) + "\n")
        outputfile.write("\n")

        outputfile.write("##debug code\n")	
        outputfile.write("##the purpose of this code is so that you can turn the inner-loop solver on and off, force a sequence of planets and/or phase types\n")	
        outputfile.write("#sequence, must have (max_phases_per_journey) entries for each journey. Use 0 to encode no flyby\n")
        outputfile.write("#integer codes represent planets\n")
        outputfile.write("#this option is NOT used if the outer-loop is turned on\n")
        outputfile.write("#first number is the number of sequences listed, followed by the sequences\n")
        outputfile.write("sequence " + str(self.number_of_trial_sequences) + "\n")
        for entry in range(0, self.number_of_trial_sequences):
            for j in range(0, self.number_of_journeys):
                #first check to make sure that this journey has a sequence
                #line to print
                if j == 0:
                    outputfile.write(str(self.Journeys[j].sequence[entry][0]))
                else:
                    outputfile.write(" " + str(self.Journeys[j].sequence[entry][0]))

                for p in range(1, self.max_phases_per_journey):
                    if p <= len(self.Journeys[j].sequence[entry]):
                        outputfile.write(" " + str(self.Journeys[j].sequence[entry][p]))
                    else:
                        outputfile.write(" 0")
            outputfile.write("\n")
        outputfile.write("#phase type, must have one entry for each phase in the mission\n")
        outputfile.write("#this option allows you to have different phases use different propulsion systems\n")
        outputfile.write("#0: MGA, 1: MGA-DSM, 2: MGA-LT, 3: FBLT, 4: FBLT-S \n")
        outputfile.write("#if mission_type is set to 0, 1, 2 then the following option is ignored\n")
        outputfile.write("#if mission_type > 4 and the outer-loop is ON, then the following option is ignored\n")
        outputfile.write("#the following option is only used if the outer-loop is OFF and mission_type > 4\n")
        outputfile.write("#not specified because either the outer-loop is off or mission_type > 4\n")
        #TODO phase_type is not implemented
        outputfile.write("#Check derivatives against finite differencing?\n")
        outputfile.write("check_derivatives " + str(self.check_derivatives) + "\n")
        outputfile.write("#which inner loop solver to run?\n")	
        outputfile.write("#0: none, evaluate trialX, 1: evaluate a batch of decision vectors, 2: run MBH, 3: run constrained DE, 4: run SNOPT using trialX as initial guess\n")
        outputfile.write("run_inner_loop " + str(self.run_inner_loop) + "\n")
            
        if len(self.trialX) > 0:
            outputfile.write("#trial decision vector\n")
            outputfile.write("#trialX\n")
            outputfile.write("trialX\n")
            for seq in range(0, self.number_of_trial_sequences):
                currentX = self.trialX[seq]
                outputfile.write('%17.20f' % currentX[0])
                for entry in range(1, len(currentX)):
                    outputfile.write(" ")
                    outputfile.write('%17.20f' % currentX[entry])
                outputfile.write("\n")
        outputfile.write("\n")
            
        outputfile.write("#end options file\n")

        outputfile.close()
        
    def update_all_panels(self, optionsnotebook):
        self.update_global_options_panel(optionsnotebook)
        self.update_journey_options_panel(optionsnotebook)
        self.update_spacecraft_options_panel(optionsnotebook)
        self.update_solver_options_panel(optionsnotebook)
        self.update_physics_options_panel(optionsnotebook)
        self.update_output_options_panel(optionsnotebook)
        
    def update_global_options_panel(self, optionsnotebook):

        optionsnotebook.tabGlobal.txtMissionName.SetValue(self.mission_name)
        optionsnotebook.tabGlobal.cmbMissionType.SetSelection(self.mission_type)
        optionsnotebook.tabGlobal.txtmaximum_number_of_lambert_revolutions.SetValue(str(self.maximum_number_of_lambert_revolutions))
        optionsnotebook.tabGlobal.cmbobjective_type.SetSelection(self.objective_type)
        optionsnotebook.tabGlobal.chkinclude_initial_impulse_in_cost.SetValue(self.include_initial_impulse_in_cost)
        optionsnotebook.tabGlobal.txtmax_phases_per_journey.SetValue(str(self.max_phases_per_journey))
        optionsnotebook.tabGlobal.txtlaunch_window_open_date.SetValue(str(self.launch_window_open_date))
        CurrentLaunchDate = wx.DateTimeFromJDN(self.launch_window_open_date + 2400000.5)
        optionsnotebook.tabGlobal.LaunchDateCalendar.SetDate(CurrentLaunchDate.MakeUTC())
        optionsnotebook.tabGlobal.txtnum_timesteps.SetValue(str(self.num_timesteps))
        optionsnotebook.tabGlobal.cmbstep_size_distribution.SetSelection(self.step_size_distribution)
        optionsnotebook.tabGlobal.txtstep_size_stdv_or_scale.SetValue(str(self.step_size_stdv_or_scale))
        optionsnotebook.tabGlobal.cmbcontrol_coordinate_system.SetSelection(self.control_coordinate_system)
        optionsnotebook.tabGlobal.txtDLA_bounds_lower.SetValue(str(self.DLA_bounds[0]))
        optionsnotebook.tabGlobal.txtDLA_bounds_upper.SetValue(str(self.DLA_bounds[1]))
        optionsnotebook.tabGlobal.chkglobal_timebounded.SetValue(self.global_timebounded)
        optionsnotebook.tabGlobal.txttotal_flight_time_bounds_lower.SetValue(str(self.total_flight_time_bounds[0]))
        optionsnotebook.tabGlobal.txttotal_flight_time_bounds_upper.SetValue(str(self.total_flight_time_bounds[1]))
        optionsnotebook.tabGlobal.txtforced_post_launch_coast.SetValue(str(self.forced_post_launch_coast))
        optionsnotebook.tabGlobal.txtforced_flyby_coast.SetValue(str(self.forced_flyby_coast))
        optionsnotebook.tabGlobal.txtinitial_V_infinity_x.SetValue(str(self.initial_V_infinity[0]))
        optionsnotebook.tabGlobal.txtinitial_V_infinity_y.SetValue(str(self.initial_V_infinity[1]))
        optionsnotebook.tabGlobal.txtinitial_V_infinity_z.SetValue(str(self.initial_V_infinity[2]))
        optionsnotebook.tabGlobal.txtminimum_dry_mass.SetValue(str(self.minimum_dry_mass))
        optionsnotebook.tabGlobal.txtpost_mission_delta_v.SetValue(str(self.post_mission_delta_v))

        #for Lambert mission types, show number of Lambert revs
        if self.mission_type == 0 or self.mission_type == 1:
            optionsnotebook.tabGlobal.lblmaximum_number_of_lambert_revolutions.Show(True)
            optionsnotebook.tabGlobal.txtmaximum_number_of_lambert_revolutions.Show(True)
        else:
            optionsnotebook.tabGlobal.lblmaximum_number_of_lambert_revolutions.Show(False)
            optionsnotebook.tabGlobal.txtmaximum_number_of_lambert_revolutions.Show(False)

        #if objective type is delta-v, show include initial impulse in cost
        if self.objective_type == 0:
            optionsnotebook.tabGlobal.lblinclude_initial_impulse_in_cost.Show(True)
            optionsnotebook.tabGlobal.chkinclude_initial_impulse_in_cost.Show(True)
        else:
            optionsnotebook.tabGlobal.lblinclude_initial_impulse_in_cost.Show(False)
            optionsnotebook.tabGlobal.chkinclude_initial_impulse_in_cost.Show(False)

        #if step size distribution is uniform, hide the scale factor box
        if self.step_size_distribution == 0:
            optionsnotebook.tabGlobal.lblstep_size_stdv_or_scale.Show(False)
            optionsnotebook.tabGlobal.txtstep_size_stdv_or_scale.Show(False)
        else:
            optionsnotebook.tabGlobal.lblstep_size_stdv_or_scale.Show(True)
            optionsnotebook.tabGlobal.txtstep_size_stdv_or_scale.Show(True)

        #if global time bounds are off, don't show the bounds fields
        if self.global_timebounded == 1:
            optionsnotebook.tabGlobal.lbltotal_flight_time_bounds.Show(True)
            optionsnotebook.tabGlobal.txttotal_flight_time_bounds_lower.Show(True)
            optionsnotebook.tabGlobal.txttotal_flight_time_bounds_upper.Show(True)
        else:
            optionsnotebook.tabGlobal.lbltotal_flight_time_bounds.Show(False)
            optionsnotebook.tabGlobal.txttotal_flight_time_bounds_lower.Show(False)
            optionsnotebook.tabGlobal.txttotal_flight_time_bounds_upper.Show(False)

        #if the minimum dry mass constraint is not active then make the post-mission delta-v field invisible
        if self.minimum_dry_mass > 0.0:
            optionsnotebook.tabGlobal.lblpost_mission_delta_v.Show(True)
            optionsnotebook.tabGlobal.txtpost_mission_delta_v.Show(True)
        else:
            optionsnotebook.tabGlobal.txtpost_mission_delta_v.Show(False)
            optionsnotebook.tabGlobal.lblpost_mission_delta_v.Show(False)

        #control coordinate system is only shown for low-thrust mission types
        if self.mission_type == 2 or self.mission_type == 3:
            optionsnotebook.tabGlobal.lblcontrol_coordinate_system.Show(True)
            optionsnotebook.tabGlobal.cmbcontrol_coordinate_system.Show(True)
        else:
            optionsnotebook.tabGlobal.lblcontrol_coordinate_system.Show(False)
            optionsnotebook.tabGlobal.cmbcontrol_coordinate_system.Show(False)

        optionsnotebook.tabGlobal.Layout()
        optionsnotebook.tabGlobal.SetupScrolling()

    def update_journey_options_panel(self, optionsnotebook):
        optionsnotebook.tabJourney.Journeylist = []
        for j in range(0, self.number_of_journeys):
            optionsnotebook.tabJourney.Journeylist.append(self.Journeys[j].journey_names)

        optionsnotebook.tabJourney.JourneySelectBox.SetItems(optionsnotebook.tabJourney.Journeylist)
        optionsnotebook.tabJourney.JourneySelectBox.SetSelection(self.ActiveJourney)

        optionsnotebook.tabJourney.txtjourney_names.SetValue(str(self.Journeys[self.ActiveJourney].journey_names))
        optionsnotebook.tabJourney.txtjourney_central_body.SetValue(str(self.Journeys[self.ActiveJourney].journey_central_body))
        optionsnotebook.tabJourney.txtdestination_list.SetValue(str(self.Journeys[self.ActiveJourney].destination_list))
        optionsnotebook.tabJourney.txtjourney_starting_mass_increment.SetValue(str(self.Journeys[self.ActiveJourney].journey_starting_mass_increment))
        optionsnotebook.tabJourney.chkjourney_variable_mass_increment.SetValue(self.Journeys[self.ActiveJourney].journey_variable_mass_increment)
        optionsnotebook.tabJourney.txtjourney_wait_time_bounds_lower.SetValue(str(self.Journeys[self.ActiveJourney].journey_wait_time_bounds[0]))
        optionsnotebook.tabJourney.txtjourney_wait_time_bounds_upper.SetValue(str(self.Journeys[self.ActiveJourney].journey_wait_time_bounds[1]))
        optionsnotebook.tabJourney.cmbjourney_timebounded.SetSelection(self.Journeys[self.ActiveJourney].journey_timebounded)
        optionsnotebook.tabJourney.txtjourney_flight_time_bounds_lower.SetValue(str(self.Journeys[self.ActiveJourney].journey_flight_time_bounds[0]))
        optionsnotebook.tabJourney.txtjourney_flight_time_bounds_upper.SetValue(str(self.Journeys[self.ActiveJourney].journey_flight_time_bounds[1]))
        optionsnotebook.tabJourney.txtjourney_arrival_date_bounds_lower.SetValue(str(self.Journeys[self.ActiveJourney].journey_arrival_date_bounds[0]))
        date = wx.DateTimeFromJDN(self.Journeys[self.ActiveJourney].journey_arrival_date_bounds[0] + 2400000.5)
        optionsnotebook.tabJourney.ArrivalDateLowerCalendar.SetDate(date.MakeUTC())
        optionsnotebook.tabJourney.txtjourney_arrival_date_bounds_upper.SetValue(str(self.Journeys[self.ActiveJourney].journey_arrival_date_bounds[1]))
        date = wx.DateTimeFromJDN(self.Journeys[self.ActiveJourney].journey_arrival_date_bounds[1] + 2400000.5)
        optionsnotebook.tabJourney.ArrivalDateUpperCalendar.SetDate(date.MakeUTC())
        optionsnotebook.tabJourney.txtjourney_initial_impulse_bounds_lower.SetValue(str(self.Journeys[self.ActiveJourney].journey_initial_impulse_bounds[0]))
        optionsnotebook.tabJourney.txtjourney_initial_impulse_bounds_upper.SetValue(str(self.Journeys[self.ActiveJourney].journey_initial_impulse_bounds[1]))
        optionsnotebook.tabJourney.cmbjourney_departure_type.SetSelection(self.Journeys[self.ActiveJourney].journey_departure_type)
        optionsnotebook.tabJourney.txtjourney_escape_spiral_starting_radius.SetValue(str(self.Journeys[self.ActiveJourney].journey_escape_spiral_starting_radius))
        optionsnotebook.tabJourney.txtjourney_initial_velocity0.SetValue(str(self.Journeys[self.ActiveJourney].journey_initial_velocity[0]))
        optionsnotebook.tabJourney.txtjourney_initial_velocity1.SetValue(str(self.Journeys[self.ActiveJourney].journey_initial_velocity[1]))
        optionsnotebook.tabJourney.txtjourney_initial_velocity2.SetValue(str(self.Journeys[self.ActiveJourney].journey_initial_velocity[2]))
        optionsnotebook.tabJourney.cmbjourney_arrival_type.SetSelection(self.Journeys[self.ActiveJourney].journey_arrival_type)
        optionsnotebook.tabJourney.txtjourney_capture_spiral_final_radius.SetValue(str(self.Journeys[self.ActiveJourney].journey_capture_spiral_final_radius))
        optionsnotebook.tabJourney.txtjourney_final_velocity0.SetValue(str(self.Journeys[self.ActiveJourney].journey_final_velocity[0]))
        optionsnotebook.tabJourney.txtjourney_final_velocity1.SetValue(str(self.Journeys[self.ActiveJourney].journey_final_velocity[1]))
        optionsnotebook.tabJourney.txtjourney_final_velocity2.SetValue(str(self.Journeys[self.ActiveJourney].journey_final_velocity[2]))
        optionsnotebook.tabJourney.chkjourney_arrival_declination_constraint_flag.SetValue(self.Journeys[self.ActiveJourney].journey_arrival_declination_constraint_flag)
        if self.Journeys[self.ActiveJourney].journey_arrival_declination_constraint_flag:
            optionsnotebook.tabJourney.txtjourney_arrival_declination_bounds_lower.SetValue(str(self.Journeys[self.ActiveJourney].journey_arrival_declination_bounds[0]))
            optionsnotebook.tabJourney.txtjourney_arrival_declination_bounds_upper.SetValue(str(self.Journeys[self.ActiveJourney].journey_arrival_declination_bounds[1]))
        optionsnotebook.tabJourney.txtsequence.SetValue(str(self.Journeys[self.ActiveJourney].sequence))
        optionsnotebook.tabJourney.txtjourney_perturbation_bodies.SetValue(str(self.Journeys[self.ActiveJourney].journey_perturbation_bodies))
        optionsnotebook.tabJourney.cmbjourney_departure_elements_type.SetSelection(self.Journeys[self.ActiveJourney].journey_departure_elements_type)
        optionsnotebook.tabJourney.chkSMA_departure.SetValue(self.Journeys[self.ActiveJourney].journey_departure_elements_vary_flag[0])
        optionsnotebook.tabJourney.chkECC_departure.SetValue(self.Journeys[self.ActiveJourney].journey_departure_elements_vary_flag[1])
        optionsnotebook.tabJourney.chkINC_departure.SetValue(self.Journeys[self.ActiveJourney].journey_departure_elements_vary_flag[2])
        optionsnotebook.tabJourney.chkRAAN_departure.SetValue(self.Journeys[self.ActiveJourney].journey_departure_elements_vary_flag[3])
        optionsnotebook.tabJourney.chkAOP_departure.SetValue(self.Journeys[self.ActiveJourney].journey_departure_elements_vary_flag[4])
        optionsnotebook.tabJourney.chkMA_departure.SetValue(self.Journeys[self.ActiveJourney].journey_departure_elements_vary_flag[5])
        optionsnotebook.tabJourney.txtSMA_departure.SetValue(str(self.Journeys[self.ActiveJourney].journey_departure_elements[0]))
        optionsnotebook.tabJourney.txtECC_departure.SetValue(str(self.Journeys[self.ActiveJourney].journey_departure_elements[1]))
        optionsnotebook.tabJourney.txtINC_departure.SetValue(str(self.Journeys[self.ActiveJourney].journey_departure_elements[2]))
        optionsnotebook.tabJourney.txtRAAN_departure.SetValue(str(self.Journeys[self.ActiveJourney].journey_departure_elements[3]))
        optionsnotebook.tabJourney.txtAOP_departure.SetValue(str(self.Journeys[self.ActiveJourney].journey_departure_elements[4]))
        optionsnotebook.tabJourney.txtMA_departure.SetValue(str(self.Journeys[self.ActiveJourney].journey_departure_elements[5]))
        optionsnotebook.tabJourney.txtSMA_departure0.SetValue(str(self.Journeys[self.ActiveJourney].journey_departure_elements_bounds[0]))
        optionsnotebook.tabJourney.txtECC_departure0.SetValue(str(self.Journeys[self.ActiveJourney].journey_departure_elements_bounds[2]))
        optionsnotebook.tabJourney.txtINC_departure0.SetValue(str(self.Journeys[self.ActiveJourney].journey_departure_elements_bounds[4]))
        optionsnotebook.tabJourney.txtRAAN_departure0.SetValue(str(self.Journeys[self.ActiveJourney].journey_departure_elements_bounds[6]))
        optionsnotebook.tabJourney.txtAOP_departure0.SetValue(str(self.Journeys[self.ActiveJourney].journey_departure_elements_bounds[8]))
        optionsnotebook.tabJourney.txtMA_departure0.SetValue(str(self.Journeys[self.ActiveJourney].journey_departure_elements_bounds[10]))
        optionsnotebook.tabJourney.txtSMA_departure1.SetValue(str(self.Journeys[self.ActiveJourney].journey_departure_elements_bounds[1]))
        optionsnotebook.tabJourney.txtECC_departure1.SetValue(str(self.Journeys[self.ActiveJourney].journey_departure_elements_bounds[3]))
        optionsnotebook.tabJourney.txtINC_departure1.SetValue(str(self.Journeys[self.ActiveJourney].journey_departure_elements_bounds[5]))
        optionsnotebook.tabJourney.txtRAAN_departure1.SetValue(str(self.Journeys[self.ActiveJourney].journey_departure_elements_bounds[7]))
        optionsnotebook.tabJourney.txtAOP_departure1.SetValue(str(self.Journeys[self.ActiveJourney].journey_departure_elements_bounds[9]))
        optionsnotebook.tabJourney.txtMA_departure1.SetValue(str(self.Journeys[self.ActiveJourney].journey_departure_elements_bounds[11]))
        optionsnotebook.tabJourney.cmbjourney_arrival_elements_type.SetSelection(self.Journeys[self.ActiveJourney].journey_arrival_elements_type)
        optionsnotebook.tabJourney.chkSMA_arrival.SetValue(self.Journeys[self.ActiveJourney].journey_arrival_elements_vary_flag[0])
        optionsnotebook.tabJourney.chkECC_arrival.SetValue(self.Journeys[self.ActiveJourney].journey_arrival_elements_vary_flag[1])
        optionsnotebook.tabJourney.chkINC_arrival.SetValue(self.Journeys[self.ActiveJourney].journey_arrival_elements_vary_flag[2])
        optionsnotebook.tabJourney.chkRAAN_arrival.SetValue(self.Journeys[self.ActiveJourney].journey_arrival_elements_vary_flag[3])
        optionsnotebook.tabJourney.chkAOP_arrival.SetValue(self.Journeys[self.ActiveJourney].journey_arrival_elements_vary_flag[4])
        optionsnotebook.tabJourney.chkMA_arrival.SetValue(self.Journeys[self.ActiveJourney].journey_arrival_elements_vary_flag[5])
        optionsnotebook.tabJourney.txtSMA_arrival.SetValue(str(self.Journeys[self.ActiveJourney].journey_arrival_elements[0]))
        optionsnotebook.tabJourney.txtECC_arrival.SetValue(str(self.Journeys[self.ActiveJourney].journey_arrival_elements[1]))
        optionsnotebook.tabJourney.txtINC_arrival.SetValue(str(self.Journeys[self.ActiveJourney].journey_arrival_elements[2]))
        optionsnotebook.tabJourney.txtRAAN_arrival.SetValue(str(self.Journeys[self.ActiveJourney].journey_arrival_elements[3]))
        optionsnotebook.tabJourney.txtAOP_arrival.SetValue(str(self.Journeys[self.ActiveJourney].journey_arrival_elements[4]))
        optionsnotebook.tabJourney.txtMA_arrival.SetValue(str(self.Journeys[self.ActiveJourney].journey_arrival_elements[5]))
        optionsnotebook.tabJourney.txtSMA_arrival0.SetValue(str(self.Journeys[self.ActiveJourney].journey_arrival_elements_bounds[0]))
        optionsnotebook.tabJourney.txtECC_arrival0.SetValue(str(self.Journeys[self.ActiveJourney].journey_arrival_elements_bounds[2]))
        optionsnotebook.tabJourney.txtINC_arrival0.SetValue(str(self.Journeys[self.ActiveJourney].journey_arrival_elements_bounds[4]))
        optionsnotebook.tabJourney.txtRAAN_arrival0.SetValue(str(self.Journeys[self.ActiveJourney].journey_arrival_elements_bounds[6]))
        optionsnotebook.tabJourney.txtAOP_arrival0.SetValue(str(self.Journeys[self.ActiveJourney].journey_arrival_elements_bounds[8]))
        optionsnotebook.tabJourney.txtMA_arrival0.SetValue(str(self.Journeys[self.ActiveJourney].journey_arrival_elements_bounds[10]))
        optionsnotebook.tabJourney.txtSMA_arrival1.SetValue(str(self.Journeys[self.ActiveJourney].journey_arrival_elements_bounds[1]))
        optionsnotebook.tabJourney.txtECC_arrival1.SetValue(str(self.Journeys[self.ActiveJourney].journey_arrival_elements_bounds[3]))
        optionsnotebook.tabJourney.txtINC_arrival1.SetValue(str(self.Journeys[self.ActiveJourney].journey_arrival_elements_bounds[5]))
        optionsnotebook.tabJourney.txtRAAN_arrival1.SetValue(str(self.Journeys[self.ActiveJourney].journey_arrival_elements_bounds[7]))
        optionsnotebook.tabJourney.txtAOP_arrival1.SetValue(str(self.Journeys[self.ActiveJourney].journey_arrival_elements_bounds[9]))
        optionsnotebook.tabJourney.txtMA_arrival1.SetValue(str(self.Journeys[self.ActiveJourney].journey_arrival_elements_bounds[11]))
        optionsnotebook.tabJourney.chkjourney_maximum_DSM_magnitude_flag.SetValue(self.Journeys[self.ActiveJourney].journey_maximum_DSM_magnitude_constraint_flag)
        optionsnotebook.tabJourney.txtjourney_maximum_DSM_magnitude.SetValue(str(self.Journeys[self.ActiveJourney].journey_maximum_DSM_magnitude_constraint))

        #if there is only one journey in the list then disable delete, up, and down
        if self.number_of_journeys == 1:
            optionsnotebook.tabJourney.btnDeleteJourney.Disable()
            optionsnotebook.tabJourney.btnMoveJourneyUp.Disable()
            optionsnotebook.tabJourney.btnMoveJourneyDown.Disable()
        else:
            optionsnotebook.tabJourney.btnDeleteJourney.Enable()

            #if the first journey in the list is active then you cannot move up
            if self.ActiveJourney == 0:
                optionsnotebook.tabJourney.btnMoveJourneyUp.Disable()
            else:
                optionsnotebook.tabJourney.btnMoveJourneyUp.Enable()

            #if the last journey in the list is active then you cannot move down
            if self.ActiveJourney == self.number_of_journeys - 1:
                optionsnotebook.tabJourney.btnMoveJourneyDown.Disable()
            else:
                optionsnotebook.tabJourney.btnMoveJourneyDown.Enable()

        #hide or show flight time and arrival date bounds
        if self.Journeys[self.ActiveJourney].journey_timebounded == 0:
            optionsnotebook.tabJourney.lbljourney_flight_time_bounds.Show(False)
            optionsnotebook.tabJourney.txtjourney_flight_time_bounds_lower.Show(False)
            optionsnotebook.tabJourney.txtjourney_flight_time_bounds_upper.Show(False)
            optionsnotebook.tabJourney.lbljourney_arrival_date_bounds.Show(False)
            optionsnotebook.tabJourney.txtjourney_arrival_date_bounds_lower.Show(False)
            optionsnotebook.tabJourney.txtjourney_arrival_date_bounds_upper.Show(False)
            optionsnotebook.tabJourney.ArrivalDateLowerCalendar.Show(False)
            optionsnotebook.tabJourney.ArrivalDateUpperCalendar.Show(False)
        elif self.Journeys[self.ActiveJourney].journey_timebounded == 1 or self.Journeys[self.ActiveJourney].journey_timebounded == 3:
            optionsnotebook.tabJourney.lbljourney_flight_time_bounds.Show(True)
            optionsnotebook.tabJourney.txtjourney_flight_time_bounds_lower.Show(True)
            optionsnotebook.tabJourney.txtjourney_flight_time_bounds_upper.Show(True)
            optionsnotebook.tabJourney.lbljourney_arrival_date_bounds.Show(False)
            optionsnotebook.tabJourney.txtjourney_arrival_date_bounds_lower.Show(False)
            optionsnotebook.tabJourney.txtjourney_arrival_date_bounds_upper.Show(False)
            optionsnotebook.tabJourney.ArrivalDateLowerCalendar.Show(False)
            optionsnotebook.tabJourney.ArrivalDateUpperCalendar.Show(False)
        elif self.Journeys[self.ActiveJourney].journey_timebounded == 2:
            optionsnotebook.tabJourney.lbljourney_flight_time_bounds.Show(False)
            optionsnotebook.tabJourney.txtjourney_flight_time_bounds_lower.Show(False)
            optionsnotebook.tabJourney.txtjourney_flight_time_bounds_upper.Show(False)
            optionsnotebook.tabJourney.lbljourney_arrival_date_bounds.Show(True)
            optionsnotebook.tabJourney.txtjourney_arrival_date_bounds_lower.Show(True)
            optionsnotebook.tabJourney.txtjourney_arrival_date_bounds_upper.Show(True)
            optionsnotebook.tabJourney.ArrivalDateLowerCalendar.Show(True)
            optionsnotebook.tabJourney.ArrivalDateUpperCalendar.Show(True)

        #hide or show DSM magnitude constraint
        if (self.mission_type == 1): #only relevant for MGA-DSM
            optionsnotebook.tabJourney.lbljourney_maximum_DSM_magnitude_flag.Show(True)
            optionsnotebook.tabJourney.chkjourney_maximum_DSM_magnitude_flag.Show(True)
            if self.Journeys[self.ActiveJourney].journey_maximum_DSM_magnitude_constraint_flag:
                optionsnotebook.tabJourney.lbljourney_maximum_DSM_magnitude.Show(True)
                optionsnotebook.tabJourney.txtjourney_maximum_DSM_magnitude.Show(True)
            else:
                optionsnotebook.tabJourney.lbljourney_maximum_DSM_magnitude.Show(False)
                optionsnotebook.tabJourney.txtjourney_maximum_DSM_magnitude.Show(False)
        else:
            optionsnotebook.tabJourney.lbljourney_maximum_DSM_magnitude_flag.Show(False)
            optionsnotebook.tabJourney.chkjourney_maximum_DSM_magnitude_flag.Show(False)
            optionsnotebook.tabJourney.lbljourney_maximum_DSM_magnitude.Show(False)
            optionsnotebook.tabJourney.txtjourney_maximum_DSM_magnitude.Show(False)

        #enable or disable the orbit elements selection boxes as appropriate
        if self.Journeys[self.ActiveJourney].destination_list[0] == -1:
            #enable departure orbit elements box
            optionsnotebook.tabJourney.boxjourney_departure_elements.Show(True)
            optionsnotebook.tabJourney.lbljourney_departure_elements_type.Show(True)
            optionsnotebook.tabJourney.cmbjourney_departure_elements_type.Show(True)
            optionsnotebook.tabJourney.lblvarydepartureelements.Show(True)
            optionsnotebook.tabJourney.lbldepartureelementsvalue.Show(True)
            optionsnotebook.tabJourney.lbldepartureelementslower.Show(True)
            optionsnotebook.tabJourney.lbldepartureelementsupper.Show(True)
            optionsnotebook.tabJourney.lblSMA_departure.Show(True)
            optionsnotebook.tabJourney.lblECC_departure.Show(True)
            optionsnotebook.tabJourney.lblINC_departure.Show(True)
            optionsnotebook.tabJourney.lblRAAN_departure.Show(True)
            optionsnotebook.tabJourney.lblAOP_departure.Show(True)
            optionsnotebook.tabJourney.lblMA_departure.Show(True)
            optionsnotebook.tabJourney.chkSMA_departure.Show(True)
            optionsnotebook.tabJourney.chkECC_departure.Show(True)
            optionsnotebook.tabJourney.chkINC_departure.Show(True)
            optionsnotebook.tabJourney.chkRAAN_departure.Show(True)
            optionsnotebook.tabJourney.chkAOP_departure.Show(True)
            optionsnotebook.tabJourney.chkMA_departure.Show(True)
            optionsnotebook.tabJourney.txtSMA_departure.Show(True)
            optionsnotebook.tabJourney.txtECC_departure.Show(True)
            optionsnotebook.tabJourney.txtINC_departure.Show(True)
            optionsnotebook.tabJourney.txtRAAN_departure.Show(True)
            optionsnotebook.tabJourney.txtAOP_departure.Show(True)
            optionsnotebook.tabJourney.txtMA_departure.Show(True)
            optionsnotebook.tabJourney.txtSMA_departure0.Show(True)
            optionsnotebook.tabJourney.txtECC_departure0.Show(True)
            optionsnotebook.tabJourney.txtINC_departure0.Show(True)
            optionsnotebook.tabJourney.txtRAAN_departure0.Show(True)
            optionsnotebook.tabJourney.txtAOP_departure0.Show(True)
            optionsnotebook.tabJourney.txtMA_departure0.Show(True)
            optionsnotebook.tabJourney.txtSMA_departure1.Show(True)
            optionsnotebook.tabJourney.txtECC_departure1.Show(True)
            optionsnotebook.tabJourney.txtINC_departure1.Show(True)
            optionsnotebook.tabJourney.txtRAAN_departure1.Show(True)
            optionsnotebook.tabJourney.txtAOP_departure1.Show(True)
            optionsnotebook.tabJourney.txtMA_departure1.Show(True)

            if self.Journeys[self.ActiveJourney].journey_departure_elements_type == 0:
                #inertial elements
                optionsnotebook.tabJourney.lblSMA_departure.SetLabel("x (km)")
                optionsnotebook.tabJourney.lblECC_departure.SetLabel("y (km)")
                optionsnotebook.tabJourney.lblINC_departure.SetLabel("z (km/s)")
                optionsnotebook.tabJourney.lblRAAN_departure.SetLabel("vx (km/s)")
                optionsnotebook.tabJourney.lblAOP_departure.SetLabel("vy (km/s)")
                optionsnotebook.tabJourney.lblMA_departure.SetLabel("vz (km/s)")
            else:
                #classical orbit elements
                optionsnotebook.tabJourney.lblSMA_departure.SetLabel("SMA (km)")
                optionsnotebook.tabJourney.lblECC_departure.SetLabel("ECC")
                optionsnotebook.tabJourney.lblINC_departure.SetLabel("INC (deg)")
                optionsnotebook.tabJourney.lblRAAN_departure.SetLabel("RAAN (deg)")
                optionsnotebook.tabJourney.lblAOP_departure.SetLabel("AOP (deg)")
                optionsnotebook.tabJourney.lblMA_departure.SetLabel("MA (deg)")

            if self.Journeys[self.ActiveJourney].journey_departure_elements_vary_flag[0] == 0:
                optionsnotebook.tabJourney.txtSMA_departure.Enable()
                optionsnotebook.tabJourney.txtSMA_departure0.Disable()
                optionsnotebook.tabJourney.txtSMA_departure1.Disable()
            else:
                optionsnotebook.tabJourney.txtSMA_departure.Disable()
                optionsnotebook.tabJourney.txtSMA_departure0.Enable()
                optionsnotebook.tabJourney.txtSMA_departure1.Enable()

            if self.Journeys[self.ActiveJourney].journey_departure_elements_vary_flag[1] == 0:
                optionsnotebook.tabJourney.txtECC_departure.Enable()
                optionsnotebook.tabJourney.txtECC_departure0.Disable()
                optionsnotebook.tabJourney.txtECC_departure1.Disable()
            else:
                optionsnotebook.tabJourney.txtECC_departure.Disable()
                optionsnotebook.tabJourney.txtECC_departure0.Enable()
                optionsnotebook.tabJourney.txtECC_departure1.Enable()

            if self.Journeys[self.ActiveJourney].journey_departure_elements_vary_flag[2] == 0:
                optionsnotebook.tabJourney.txtINC_departure.Enable()
                optionsnotebook.tabJourney.txtINC_departure0.Disable()
                optionsnotebook.tabJourney.txtINC_departure1.Disable()
            else:
                optionsnotebook.tabJourney.txtINC_departure.Disable()
                optionsnotebook.tabJourney.txtINC_departure0.Enable()
                optionsnotebook.tabJourney.txtINC_departure1.Enable()

            if self.Journeys[self.ActiveJourney].journey_departure_elements_vary_flag[3] == 0:
                optionsnotebook.tabJourney.txtRAAN_departure.Enable()
                optionsnotebook.tabJourney.txtRAAN_departure0.Disable()
                optionsnotebook.tabJourney.txtRAAN_departure1.Disable()
            else:
                optionsnotebook.tabJourney.txtRAAN_departure.Disable()
                optionsnotebook.tabJourney.txtRAAN_departure0.Enable()
                optionsnotebook.tabJourney.txtRAAN_departure1.Enable()
            
            if self.Journeys[self.ActiveJourney].journey_departure_elements_vary_flag[4] == 0:
                optionsnotebook.tabJourney.txtAOP_departure.Enable()
                optionsnotebook.tabJourney.txtAOP_departure0.Disable()
                optionsnotebook.tabJourney.txtAOP_departure1.Disable()
            else:
                optionsnotebook.tabJourney.txtAOP_departure.Disable()
                optionsnotebook.tabJourney.txtAOP_departure0.Enable()
                optionsnotebook.tabJourney.txtAOP_departure1.Enable()

            if self.Journeys[self.ActiveJourney].journey_departure_elements_vary_flag[5] == 0:
                optionsnotebook.tabJourney.txtMA_departure.Enable()
                optionsnotebook.tabJourney.txtMA_departure0.Disable()
                optionsnotebook.tabJourney.txtMA_departure1.Disable()
            else:
                optionsnotebook.tabJourney.txtMA_departure.Disable()
                optionsnotebook.tabJourney.txtMA_departure0.Enable()
                optionsnotebook.tabJourney.txtMA_departure1.Enable()
        else:
            optionsnotebook.tabJourney.boxjourney_departure_elements.Show(False)
            optionsnotebook.tabJourney.lbljourney_departure_elements_type.Show(False)
            optionsnotebook.tabJourney.cmbjourney_departure_elements_type.Show(False)
            optionsnotebook.tabJourney.lblvarydepartureelements.Show(False)
            optionsnotebook.tabJourney.lbldepartureelementsvalue.Show(False)
            optionsnotebook.tabJourney.lbldepartureelementslower.Show(False)
            optionsnotebook.tabJourney.lbldepartureelementsupper.Show(False)
            optionsnotebook.tabJourney.lblSMA_departure.Show(False)
            optionsnotebook.tabJourney.lblECC_departure.Show(False)
            optionsnotebook.tabJourney.lblINC_departure.Show(False)
            optionsnotebook.tabJourney.lblRAAN_departure.Show(False)
            optionsnotebook.tabJourney.lblAOP_departure.Show(False)
            optionsnotebook.tabJourney.lblMA_departure.Show(False)
            optionsnotebook.tabJourney.chkSMA_departure.Show(False)
            optionsnotebook.tabJourney.chkECC_departure.Show(False)
            optionsnotebook.tabJourney.chkINC_departure.Show(False)
            optionsnotebook.tabJourney.chkRAAN_departure.Show(False)
            optionsnotebook.tabJourney.chkAOP_departure.Show(False)
            optionsnotebook.tabJourney.chkMA_departure.Show(False)
            optionsnotebook.tabJourney.txtSMA_departure.Show(False)
            optionsnotebook.tabJourney.txtECC_departure.Show(False)
            optionsnotebook.tabJourney.txtINC_departure.Show(False)
            optionsnotebook.tabJourney.txtRAAN_departure.Show(False)
            optionsnotebook.tabJourney.txtAOP_departure.Show(False)
            optionsnotebook.tabJourney.txtMA_departure.Show(False)
            optionsnotebook.tabJourney.txtSMA_departure0.Show(False)
            optionsnotebook.tabJourney.txtECC_departure0.Show(False)
            optionsnotebook.tabJourney.txtINC_departure0.Show(False)
            optionsnotebook.tabJourney.txtRAAN_departure0.Show(False)
            optionsnotebook.tabJourney.txtAOP_departure0.Show(False)
            optionsnotebook.tabJourney.txtMA_departure0.Show(False)
            optionsnotebook.tabJourney.txtSMA_departure1.Show(False)
            optionsnotebook.tabJourney.txtECC_departure1.Show(False)
            optionsnotebook.tabJourney.txtINC_departure1.Show(False)
            optionsnotebook.tabJourney.txtRAAN_departure1.Show(False)
            optionsnotebook.tabJourney.txtAOP_departure1.Show(False)
            optionsnotebook.tabJourney.txtMA_departure1.Show(False)



        if self.Journeys[self.ActiveJourney].destination_list[1] == -1 or self.Journeys[self.ActiveJourney].journey_arrival_type == 0:
            #enable arrival orbit elements box
            optionsnotebook.tabJourney.boxjourney_arrival_elements.Show(True)
            optionsnotebook.tabJourney.lbljourney_arrival_elements_type.Show(True)
            optionsnotebook.tabJourney.cmbjourney_arrival_elements_type.Show(True)
            optionsnotebook.tabJourney.lblvaryarrivalelements.Show(True)
            optionsnotebook.tabJourney.lblarrivalelementsvalue.Show(True)
            optionsnotebook.tabJourney.lblarrivalelementslower.Show(True)
            optionsnotebook.tabJourney.lblarrivalelementsupper.Show(True)
            optionsnotebook.tabJourney.lblSMA_arrival.Show(True)
            optionsnotebook.tabJourney.lblECC_arrival.Show(True)
            optionsnotebook.tabJourney.lblINC_arrival.Show(True)
            optionsnotebook.tabJourney.lblRAAN_arrival.Show(True)
            optionsnotebook.tabJourney.lblAOP_arrival.Show(True)
            optionsnotebook.tabJourney.lblMA_arrival.Show(True)
            optionsnotebook.tabJourney.chkSMA_arrival.Show(True)
            optionsnotebook.tabJourney.chkECC_arrival.Show(True)
            optionsnotebook.tabJourney.chkINC_arrival.Show(True)
            optionsnotebook.tabJourney.chkRAAN_arrival.Show(True)
            optionsnotebook.tabJourney.chkAOP_arrival.Show(True)
            optionsnotebook.tabJourney.chkMA_arrival.Show(True)
            optionsnotebook.tabJourney.txtSMA_arrival.Show(True)
            optionsnotebook.tabJourney.txtECC_arrival.Show(True)
            optionsnotebook.tabJourney.txtINC_arrival.Show(True)
            optionsnotebook.tabJourney.txtRAAN_arrival.Show(True)
            optionsnotebook.tabJourney.txtAOP_arrival.Show(True)
            optionsnotebook.tabJourney.txtMA_arrival.Show(True)
            optionsnotebook.tabJourney.txtSMA_arrival0.Show(True)
            optionsnotebook.tabJourney.txtECC_arrival0.Show(True)
            optionsnotebook.tabJourney.txtINC_arrival0.Show(True)
            optionsnotebook.tabJourney.txtRAAN_arrival0.Show(True)
            optionsnotebook.tabJourney.txtAOP_arrival0.Show(True)
            optionsnotebook.tabJourney.txtMA_arrival0.Show(True)
            optionsnotebook.tabJourney.txtSMA_arrival1.Show(True)
            optionsnotebook.tabJourney.txtECC_arrival1.Show(True)
            optionsnotebook.tabJourney.txtINC_arrival1.Show(True)
            optionsnotebook.tabJourney.txtRAAN_arrival1.Show(True)
            optionsnotebook.tabJourney.txtAOP_arrival1.Show(True)
            optionsnotebook.tabJourney.txtMA_arrival1.Show(True)

            if self.Journeys[self.ActiveJourney].journey_arrival_elements_type == 0:
                #inertial elements
                optionsnotebook.tabJourney.lblSMA_arrival.SetLabel("x (km)")
                optionsnotebook.tabJourney.lblECC_arrival.SetLabel("y (km)")
                optionsnotebook.tabJourney.lblINC_arrival.SetLabel("z (km/s)")
                optionsnotebook.tabJourney.lblRAAN_arrival.SetLabel("vx (km/s)")
                optionsnotebook.tabJourney.lblAOP_arrival.SetLabel("vy (km/s)")
                optionsnotebook.tabJourney.lblMA_arrival.SetLabel("vz (km/s)")
            else:
                #classical orbit elements
                optionsnotebook.tabJourney.lblSMA_arrival.SetLabel("SMA (km)")
                optionsnotebook.tabJourney.lblECC_arrival.SetLabel("ECC")
                optionsnotebook.tabJourney.lblINC_arrival.SetLabel("INC (deg)")
                optionsnotebook.tabJourney.lblRAAN_arrival.SetLabel("RAAN (deg)")
                optionsnotebook.tabJourney.lblAOP_arrival.SetLabel("AOP (deg)")
                optionsnotebook.tabJourney.lblMA_arrival.SetLabel("MA (deg)")

            if self.Journeys[self.ActiveJourney].journey_arrival_elements_vary_flag[0] == 0:
                optionsnotebook.tabJourney.txtSMA_arrival.Enable()
                optionsnotebook.tabJourney.txtSMA_arrival0.Disable()
                optionsnotebook.tabJourney.txtSMA_arrival1.Disable()
            else:
                optionsnotebook.tabJourney.txtSMA_arrival.Disable()
                optionsnotebook.tabJourney.txtSMA_arrival0.Enable()
                optionsnotebook.tabJourney.txtSMA_arrival1.Enable()

            if self.Journeys[self.ActiveJourney].journey_arrival_elements_vary_flag[1] == 0:
                optionsnotebook.tabJourney.txtECC_arrival.Enable()
                optionsnotebook.tabJourney.txtECC_arrival0.Disable()
                optionsnotebook.tabJourney.txtECC_arrival1.Disable()
            else:
                optionsnotebook.tabJourney.txtECC_arrival.Disable()
                optionsnotebook.tabJourney.txtECC_arrival0.Enable()
                optionsnotebook.tabJourney.txtECC_arrival1.Enable()

            if self.Journeys[self.ActiveJourney].journey_arrival_elements_vary_flag[2] == 0:
                optionsnotebook.tabJourney.txtINC_arrival.Enable()
                optionsnotebook.tabJourney.txtINC_arrival0.Disable()
                optionsnotebook.tabJourney.txtINC_arrival1.Disable()
            else:
                optionsnotebook.tabJourney.txtINC_arrival.Disable()
                optionsnotebook.tabJourney.txtINC_arrival0.Enable()
                optionsnotebook.tabJourney.txtINC_arrival1.Enable()

            if self.Journeys[self.ActiveJourney].journey_arrival_elements_vary_flag[3] == 0:
                optionsnotebook.tabJourney.txtRAAN_arrival.Enable()
                optionsnotebook.tabJourney.txtRAAN_arrival0.Disable()
                optionsnotebook.tabJourney.txtRAAN_arrival1.Disable()
            else:
                optionsnotebook.tabJourney.txtRAAN_arrival.Disable()
                optionsnotebook.tabJourney.txtRAAN_arrival0.Enable()
                optionsnotebook.tabJourney.txtRAAN_arrival1.Enable()
            
            if self.Journeys[self.ActiveJourney].journey_arrival_elements_vary_flag[4] == 0:
                optionsnotebook.tabJourney.txtAOP_arrival.Enable()
                optionsnotebook.tabJourney.txtAOP_arrival0.Disable()
                optionsnotebook.tabJourney.txtAOP_arrival1.Disable()
            else:
                optionsnotebook.tabJourney.txtAOP_arrival.Disable()
                optionsnotebook.tabJourney.txtAOP_arrival0.Enable()
                optionsnotebook.tabJourney.txtAOP_arrival1.Enable()

            if self.Journeys[self.ActiveJourney].journey_arrival_elements_vary_flag[5] == 0:
                optionsnotebook.tabJourney.txtMA_arrival.Enable()
                optionsnotebook.tabJourney.txtMA_arrival0.Disable()
                optionsnotebook.tabJourney.txtMA_arrival1.Disable()
            else:
                optionsnotebook.tabJourney.txtMA_arrival.Disable()
                optionsnotebook.tabJourney.txtMA_arrival0.Enable()
                optionsnotebook.tabJourney.txtMA_arrival1.Enable()
        else:
            optionsnotebook.tabJourney.boxjourney_arrival_elements.Show(False)
            optionsnotebook.tabJourney.lbljourney_arrival_elements_type.Show(False)
            optionsnotebook.tabJourney.cmbjourney_arrival_elements_type.Show(False)
            optionsnotebook.tabJourney.lblvaryarrivalelements.Show(False)
            optionsnotebook.tabJourney.lblarrivalelementsvalue.Show(False)
            optionsnotebook.tabJourney.lblarrivalelementslower.Show(False)
            optionsnotebook.tabJourney.lblarrivalelementsupper.Show(False)
            optionsnotebook.tabJourney.lblSMA_arrival.Show(False)
            optionsnotebook.tabJourney.lblECC_arrival.Show(False)
            optionsnotebook.tabJourney.lblINC_arrival.Show(False)
            optionsnotebook.tabJourney.lblRAAN_arrival.Show(False)
            optionsnotebook.tabJourney.lblAOP_arrival.Show(False)
            optionsnotebook.tabJourney.lblMA_arrival.Show(False)
            optionsnotebook.tabJourney.chkSMA_arrival.Show(False)
            optionsnotebook.tabJourney.chkECC_arrival.Show(False)
            optionsnotebook.tabJourney.chkINC_arrival.Show(False)
            optionsnotebook.tabJourney.chkRAAN_arrival.Show(False)
            optionsnotebook.tabJourney.chkAOP_arrival.Show(False)
            optionsnotebook.tabJourney.chkMA_arrival.Show(False)
            optionsnotebook.tabJourney.txtSMA_arrival.Show(False)
            optionsnotebook.tabJourney.txtECC_arrival.Show(False)
            optionsnotebook.tabJourney.txtINC_arrival.Show(False)
            optionsnotebook.tabJourney.txtRAAN_arrival.Show(False)
            optionsnotebook.tabJourney.txtAOP_arrival.Show(False)
            optionsnotebook.tabJourney.txtMA_arrival.Show(False)
            optionsnotebook.tabJourney.txtSMA_arrival0.Show(False)
            optionsnotebook.tabJourney.txtECC_arrival0.Show(False)
            optionsnotebook.tabJourney.txtINC_arrival0.Show(False)
            optionsnotebook.tabJourney.txtRAAN_arrival0.Show(False)
            optionsnotebook.tabJourney.txtAOP_arrival0.Show(False)
            optionsnotebook.tabJourney.txtMA_arrival0.Show(False)
            optionsnotebook.tabJourney.txtSMA_arrival1.Show(False)
            optionsnotebook.tabJourney.txtECC_arrival1.Show(False)
            optionsnotebook.tabJourney.txtINC_arrival1.Show(False)
            optionsnotebook.tabJourney.txtRAAN_arrival1.Show(False)
            optionsnotebook.tabJourney.txtAOP_arrival1.Show(False)
            optionsnotebook.tabJourney.txtMA_arrival1.Show(False)

        #only show the journey v-infinity boxes if doing a bounded v-infinity-in flyby
        if self.Journeys[self.ActiveJourney].journey_departure_type == 4:
            optionsnotebook.tabJourney.lbljourney_initial_velocity.Show(True)
            optionsnotebook.tabJourney.txtjourney_initial_velocity0.Show(True)
            optionsnotebook.tabJourney.txtjourney_initial_velocity1.Show(True)
            optionsnotebook.tabJourney.txtjourney_initial_velocity2.Show(True)
        else:
            optionsnotebook.tabJourney.lbljourney_initial_velocity.Show(False)
            optionsnotebook.tabJourney.txtjourney_initial_velocity0.Show(False)
            optionsnotebook.tabJourney.txtjourney_initial_velocity1.Show(False)
            optionsnotebook.tabJourney.txtjourney_initial_velocity2.Show(False)

        if self.Journeys[self.ActiveJourney].journey_departure_type == 3 or self.Journeys[self.ActiveJourney].journey_departure_type == 4 or self.Journeys[self.ActiveJourney].journey_departure_type == 6:
            optionsnotebook.tabJourney.lbljourney_wait_time_bounds.Show(False)
            optionsnotebook.tabJourney.txtjourney_wait_time_bounds_lower.Show(False)
            optionsnotebook.tabJourney.txtjourney_wait_time_bounds_upper.Show(False)
            optionsnotebook.tabJourney.lbljourney_initial_impulse_bounds.Show(False)
            optionsnotebook.tabJourney.txtjourney_initial_impulse_bounds_lower.Show(False)
            optionsnotebook.tabJourney.txtjourney_initial_impulse_bounds_upper.Show(False)
        else:
            optionsnotebook.tabJourney.lbljourney_wait_time_bounds.Show(True)
            optionsnotebook.tabJourney.txtjourney_wait_time_bounds_lower.Show(True)
            optionsnotebook.tabJourney.txtjourney_wait_time_bounds_upper.Show(True)
            optionsnotebook.tabJourney.lbljourney_initial_impulse_bounds.Show(True)
            optionsnotebook.tabJourney.txtjourney_initial_impulse_bounds_lower.Show(True)
            optionsnotebook.tabJourney.txtjourney_initial_impulse_bounds_upper.Show(True)

        if self.Journeys[self.ActiveJourney].journey_arrival_type == 4 or self.Journeys[self.ActiveJourney].journey_arrival_type == 5:
            optionsnotebook.tabJourney.lbljourney_final_velocity.Show(True)
            optionsnotebook.tabJourney.lbljourney_final_velocity.SetLabel("Journey final velocity vector")
            optionsnotebook.tabJourney.txtjourney_final_velocity0.Show(True)
            optionsnotebook.tabJourney.txtjourney_final_velocity1.Show(True)
            optionsnotebook.tabJourney.txtjourney_final_velocity2.Show(True)
        elif self.Journeys[self.ActiveJourney].journey_arrival_type == 2:
            optionsnotebook.tabJourney.lbljourney_final_velocity.Show(True)
            optionsnotebook.tabJourney.lbljourney_final_velocity.SetLabel("Journey final velocity bounds")
            optionsnotebook.tabJourney.txtjourney_final_velocity0.Show(True)
            optionsnotebook.tabJourney.txtjourney_final_velocity1.Show(True)
            optionsnotebook.tabJourney.txtjourney_final_velocity2.Show(False)
        else:
            optionsnotebook.tabJourney.lbljourney_final_velocity.Show(False)
            optionsnotebook.tabJourney.txtjourney_final_velocity0.Show(False)
            optionsnotebook.tabJourney.txtjourney_final_velocity1.Show(False)
            optionsnotebook.tabJourney.txtjourney_final_velocity2.Show(False)

        #only show the declination constraint flag if this is a bounded v-infinity intercept or orbit insertion
        if self.Journeys[self.ActiveJourney].journey_arrival_type == 0 or self.Journeys[self.ActiveJourney].journey_arrival_type == 2:
            optionsnotebook.tabJourney.lbljourney_arrival_declination_constraint_flag.Show(True)
            optionsnotebook.tabJourney.chkjourney_arrival_declination_constraint_flag.Show(True)
            if self.Journeys[self.ActiveJourney].journey_arrival_declination_constraint_flag:
                optionsnotebook.tabJourney.lbljourney_arrival_declination_bounds.Show(True)
                optionsnotebook.tabJourney.txtjourney_arrival_declination_bounds_lower.Show(True)
                optionsnotebook.tabJourney.txtjourney_arrival_declination_bounds_upper.Show(True)
            else:
                optionsnotebook.tabJourney.lbljourney_arrival_declination_bounds.Show(False)
                optionsnotebook.tabJourney.txtjourney_arrival_declination_bounds_lower.Show(False)
                optionsnotebook.tabJourney.txtjourney_arrival_declination_bounds_upper.Show(False)
        else:
            optionsnotebook.tabJourney.lbljourney_arrival_declination_constraint_flag.Show(False)
            optionsnotebook.tabJourney.chkjourney_arrival_declination_constraint_flag.Show(False)
            optionsnotebook.tabJourney.lbljourney_arrival_declination_bounds.Show(False)
            optionsnotebook.tabJourney.txtjourney_arrival_declination_bounds_lower.Show(False)
            optionsnotebook.tabJourney.txtjourney_arrival_declination_bounds_upper.Show(False)

        #options for an escape spiral
        if self.Journeys[self.ActiveJourney].journey_departure_type == 5:
            optionsnotebook.tabJourney.lbljourney_initial_impulse_bounds.Show(False)
            optionsnotebook.tabJourney.txtjourney_initial_impulse_bounds_lower.Show(False)
            optionsnotebook.tabJourney.txtjourney_initial_impulse_bounds_upper.Show(False)
            optionsnotebook.tabJourney.lbljourney_escape_spiral_starting_radius.Show(True)
            optionsnotebook.tabJourney.txtjourney_escape_spiral_starting_radius.Show(True)
        #free direct departure
        elif self.Journeys[self.ActiveJourney].journey_departure_type == 2:
            optionsnotebook.tabJourney.lbljourney_initial_impulse_bounds.Show(False)
            optionsnotebook.tabJourney.txtjourney_initial_impulse_bounds_lower.Show(False)
            optionsnotebook.tabJourney.txtjourney_initial_impulse_bounds_upper.Show(False)
            optionsnotebook.tabJourney.lbljourney_escape_spiral_starting_radius.Show(False)
            optionsnotebook.tabJourney.txtjourney_escape_spiral_starting_radius.Show(False)
        else:
            #hide the initial v-infinity options
            optionsnotebook.tabJourney.lbljourney_initial_impulse_bounds.Show(True)
            optionsnotebook.tabJourney.txtjourney_initial_impulse_bounds_lower.Show(True)
            optionsnotebook.tabJourney.txtjourney_initial_impulse_bounds_upper.Show(True)
            optionsnotebook.tabJourney.lbljourney_escape_spiral_starting_radius.Show(False)
            optionsnotebook.tabJourney.txtjourney_escape_spiral_starting_radius.Show(False)

        #options for a capture spiral
        if self.Journeys[self.ActiveJourney].journey_arrival_type == 7:
            optionsnotebook.tabJourney.lbljourney_capture_spiral_final_radius.Show(True)
            optionsnotebook.tabJourney.txtjourney_capture_spiral_final_radius.Show(True)
        else:
            optionsnotebook.tabJourney.lbljourney_capture_spiral_final_radius.Show(False)
            optionsnotebook.tabJourney.txtjourney_capture_spiral_final_radius.Show(False)



        if self.perturb_thirdbody == 1:
            optionsnotebook.tabJourney.lbljourney_perturbation_bodies.Show(True)
            optionsnotebook.tabJourney.txtjourney_perturbation_bodies.Show(True)
            optionsnotebook.tabJourney.btnjourney_perturbation_bodies.Show(True)
        else:
            optionsnotebook.tabJourney.lbljourney_perturbation_bodies.Show(False)
            optionsnotebook.tabJourney.txtjourney_perturbation_bodies.Show(False)
            optionsnotebook.tabJourney.btnjourney_perturbation_bodies.Show(False)

        optionsnotebook.tabJourney.Layout()
        optionsnotebook.tabJourney.SetupScrolling()

    def update_spacecraft_options_panel(self, optionsnotebook):
        optionsnotebook.tabSpacecraft.txtmaximum_mass.SetValue(str(self.maximum_mass))
        optionsnotebook.tabSpacecraft.chkallow_initial_mass_to_vary.SetValue(self.allow_initial_mass_to_vary)
        optionsnotebook.tabSpacecraft.txtEP_dry_mass.SetValue(str(self.EP_dry_mass))
        optionsnotebook.tabSpacecraft.cmbLV_type.SetSelection(self.LV_type + 2)
        optionsnotebook.tabSpacecraft.txtLV_adapter_mass.SetValue(str(self.LV_adapter_mass))
        optionsnotebook.tabSpacecraft.txtIspDS.SetValue(str(self.IspDS))
        optionsnotebook.tabSpacecraft.txtIspChem.SetValue(str(self.IspChem))
        optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients0.SetValue(str(self.custom_LV_coefficients[0]))
        optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients1.SetValue(str(self.custom_LV_coefficients[1]))
        optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients2.SetValue(str(self.custom_LV_coefficients[2]))
        optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients3.SetValue(str(self.custom_LV_coefficients[3]))
        optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients4.SetValue(str(self.custom_LV_coefficients[4]))
        optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients5.SetValue(str(self.custom_LV_coefficients[5]))
        optionsnotebook.tabSpacecraft.txtcustom_LV_C3_bounds_lower.SetValue(str(self.custom_LV_C3_bounds[0]))
        optionsnotebook.tabSpacecraft.txtcustom_LV_C3_bounds_upper.SetValue(str(self.custom_LV_C3_bounds[1]))
        optionsnotebook.tabSpacecraft.txtparking_orbit_altitude.SetValue(str(self.parking_orbit_altitude))
        optionsnotebook.tabSpacecraft.txtparking_orbit_inclination.SetValue(str(self.parking_orbit_inclination))
        optionsnotebook.tabSpacecraft.txtpost_mission_Isp.SetValue(str(self.post_mission_Isp))
        optionsnotebook.tabSpacecraft.txtpropellant_margin.SetValue(str(self.propellant_margin))
        optionsnotebook.tabSpacecraft.txtpower_margin.SetValue(str(self.power_margin))
        optionsnotebook.tabSpacecraft.txtLV_margin.SetValue(str(self.LV_margin))
        optionsnotebook.tabSpacecraft.chkenable_propellant_mass_constraint.SetValue(self.enable_maximum_propellant_mass_constraint)
        optionsnotebook.tabSpacecraft.txtmaximum_propellant_mass.SetValue(str(self.maximum_propellant_mass))
        optionsnotebook.tabSpacecraft.cmbengine_type.SetSelection(self.engine_type)
        optionsnotebook.tabSpacecraft.txtnumber_of_engines.SetValue(str(self.number_of_engines))
        optionsnotebook.tabSpacecraft.cmbthrottle_logic_mode.SetSelection(self.throttle_logic_mode)
        optionsnotebook.tabSpacecraft.txtthrottle_sharpness.SetValue(str(self.throttle_sharpness))
        optionsnotebook.tabSpacecraft.txtengine_duty_cycle.SetValue(str(self.engine_duty_cycle))
        optionsnotebook.tabSpacecraft.txtThrust.SetValue(str(self.Thrust))
        optionsnotebook.tabSpacecraft.txtIspLT.SetValue(str(self.IspLT))
        optionsnotebook.tabSpacecraft.txtIspLT_minimum.SetValue(str(self.IspLT_minimum))
        optionsnotebook.tabSpacecraft.txtuser_defined_engine_efficiency.SetValue(str(self.user_defined_engine_efficiency))
        optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients0.SetValue(str(self.engine_input_thrust_coefficients[0]))
        optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients1.SetValue(str(self.engine_input_thrust_coefficients[1]))
        optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients2.SetValue(str(self.engine_input_thrust_coefficients[2]))
        optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients3.SetValue(str(self.engine_input_thrust_coefficients[3]))
        optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients4.SetValue(str(self.engine_input_thrust_coefficients[4]))
        optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients5.SetValue(str(self.engine_input_thrust_coefficients[5]))
        optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients6.SetValue(str(self.engine_input_thrust_coefficients[6]))
        optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients0.SetValue(str(self.engine_input_mass_flow_rate_coefficients[0]))
        optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients1.SetValue(str(self.engine_input_mass_flow_rate_coefficients[1]))
        optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients2.SetValue(str(self.engine_input_mass_flow_rate_coefficients[2]))
        optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients3.SetValue(str(self.engine_input_mass_flow_rate_coefficients[3]))
        optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients4.SetValue(str(self.engine_input_mass_flow_rate_coefficients[4]))
        optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients5.SetValue(str(self.engine_input_mass_flow_rate_coefficients[5]))
        optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients6.SetValue(str(self.engine_input_mass_flow_rate_coefficients[6]))
        optionsnotebook.tabSpacecraft.txtengine_input_power_bounds_lower.SetValue(str(self.engine_input_power_bounds[0]))
        optionsnotebook.tabSpacecraft.txtengine_input_power_bounds_upper.SetValue(str(self.engine_input_power_bounds[1]))
        optionsnotebook.tabSpacecraft.txtpower_at_1_AU.SetValue(str(self.power_at_1_AU))
        optionsnotebook.tabSpacecraft.cmbpower_source_type.SetSelection(self.power_source_type)
        optionsnotebook.tabSpacecraft.txtsolar_power_gamma0.SetValue(str(self.solar_power_gamma[0]))
        optionsnotebook.tabSpacecraft.txtsolar_power_gamma1.SetValue(str(self.solar_power_gamma[1]))
        optionsnotebook.tabSpacecraft.txtsolar_power_gamma2.SetValue(str(self.solar_power_gamma[2]))
        optionsnotebook.tabSpacecraft.txtsolar_power_gamma3.SetValue(str(self.solar_power_gamma[3]))
        optionsnotebook.tabSpacecraft.txtsolar_power_gamma4.SetValue(str(self.solar_power_gamma[4]))
        optionsnotebook.tabSpacecraft.cmbspacecraft_power_model_type.SetSelection(self.spacecraft_power_model_type)
        optionsnotebook.tabSpacecraft.txtspacecraft_power_coefficients0.SetValue(str(self.spacecraft_power_coefficients[0]))
        optionsnotebook.tabSpacecraft.txtspacecraft_power_coefficients1.SetValue(str(self.spacecraft_power_coefficients[1]))
        optionsnotebook.tabSpacecraft.txtspacecraft_power_coefficients2.SetValue(str(self.spacecraft_power_coefficients[2]))
        optionsnotebook.tabSpacecraft.txtpower_decay_rate.SetValue(str(self.power_decay_rate))

        #if the minimum dry mass constraint is not active then make the post-mission delta-v field invisible
        if self.minimum_dry_mass > 0.0 or self.enable_maximum_propellant_mass_constraint > 0:
            optionsnotebook.tabSpacecraft.lblpost_mission_Isp.Show(True)
            optionsnotebook.tabSpacecraft.txtpost_mission_Isp.Show(True)
        else:
            optionsnotebook.tabSpacecraft.lblpost_mission_Isp.Show(False)
            optionsnotebook.tabSpacecraft.txtpost_mission_Isp.Show(False)

        if self.enable_maximum_propellant_mass_constraint:
            optionsnotebook.tabSpacecraft.lblmaximum_propellant_mass.Show(True)
            optionsnotebook.tabSpacecraft.txtmaximum_propellant_mass.Show(True)
        else:
            optionsnotebook.tabSpacecraft.lblmaximum_propellant_mass.Show(False)
            optionsnotebook.tabSpacecraft.txtmaximum_propellant_mass.Show(False)

        #switch various spacecraft fields visible and invisible depending on the mission type
        #launch vehicle types
        if self.LV_type == -2:
            #custom launch vehicle
            optionsnotebook.tabSpacecraft.lblIspDS.Show(False)
            optionsnotebook.tabSpacecraft.lblcustom_LV_coefficients.Show(True)
            optionsnotebook.tabSpacecraft.lblcustom_LV_C3_bounds.Show(True)
            optionsnotebook.tabSpacecraft.lblparking_orbit_altitude.Show(False)
            optionsnotebook.tabSpacecraft.lblparking_orbit_inclination.Show(False)
            optionsnotebook.tabSpacecraft.txtIspDS.Show(False)
            optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients0.Show(True)
            optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients1.Show(True)
            optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients2.Show(True)
            optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients3.Show(True)
            optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients4.Show(True)
            optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients5.Show(True)
            optionsnotebook.tabSpacecraft.txtcustom_LV_C3_bounds_lower.Show(True)
            optionsnotebook.tabSpacecraft.txtcustom_LV_C3_bounds_upper.Show(True)
            optionsnotebook.tabSpacecraft.txtparking_orbit_altitude.Show(False)
            optionsnotebook.tabSpacecraft.txtparking_orbit_inclination.Show(False)
            optionsnotebook.tabSpacecraft.lblLV_margin.Show(True)
            optionsnotebook.tabSpacecraft.txtLV_margin.Show(True)
            optionsnotebook.tabSpacecraft.lblLV_adapter_mass.Show(True)
            optionsnotebook.tabSpacecraft.txtLV_adapter_mass.Show(True)
        elif self.LV_type == -1:
            #burn with departure stage
            optionsnotebook.tabSpacecraft.lblIspDS.Show(True)
            optionsnotebook.tabSpacecraft.lblcustom_LV_coefficients.Show(False)
            optionsnotebook.tabSpacecraft.lblcustom_LV_C3_bounds.Show(False)
            optionsnotebook.tabSpacecraft.lblparking_orbit_altitude.Show(False)
            optionsnotebook.tabSpacecraft.lblparking_orbit_inclination.Show(False)
            optionsnotebook.tabSpacecraft.txtIspDS.Show(True)
            optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients0.Show(False)
            optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients1.Show(False)
            optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients2.Show(False)
            optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients3.Show(False)
            optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients4.Show(False)
            optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients5.Show(False)
            optionsnotebook.tabSpacecraft.txtcustom_LV_C3_bounds_lower.Show(False)
            optionsnotebook.tabSpacecraft.txtcustom_LV_C3_bounds_upper.Show(False)
            optionsnotebook.tabSpacecraft.txtparking_orbit_altitude.Show(False)
            optionsnotebook.tabSpacecraft.txtparking_orbit_inclination.Show(False)
            optionsnotebook.tabSpacecraft.lblLV_margin.Show(False)
            optionsnotebook.tabSpacecraft.txtLV_margin.Show(False)
            optionsnotebook.tabSpacecraft.lblLV_adapter_mass.Show(False)
            optionsnotebook.tabSpacecraft.txtLV_adapter_mass.Show(False)
        else:
            #fixed initial mass or hard-coded launch vehicle
            optionsnotebook.tabSpacecraft.lblIspDS.Show(False)
            optionsnotebook.tabSpacecraft.lblcustom_LV_coefficients.Show(False)
            optionsnotebook.tabSpacecraft.lblcustom_LV_C3_bounds.Show(False)
            optionsnotebook.tabSpacecraft.lblparking_orbit_altitude.Show(False)
            optionsnotebook.tabSpacecraft.lblparking_orbit_inclination.Show(False)
            optionsnotebook.tabSpacecraft.txtIspDS.Show(False)
            optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients0.Show(False)
            optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients1.Show(False)
            optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients2.Show(False)
            optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients3.Show(False)
            optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients4.Show(False)
            optionsnotebook.tabSpacecraft.txtcustom_LV_coefficients5.Show(False)
            optionsnotebook.tabSpacecraft.txtcustom_LV_C3_bounds_lower.Show(False)
            optionsnotebook.tabSpacecraft.txtcustom_LV_C3_bounds_upper.Show(False)
            optionsnotebook.tabSpacecraft.txtparking_orbit_altitude.Show(False)
            optionsnotebook.tabSpacecraft.txtparking_orbit_inclination.Show(False)

            #LV margin and adapter mass are not applicable to fixed-initial mass missions
            if self.LV_type > 0:
                optionsnotebook.tabSpacecraft.lblLV_margin.Show(True)
                optionsnotebook.tabSpacecraft.txtLV_margin.Show(True)
                optionsnotebook.tabSpacecraft.lblLV_adapter_mass.Show(True)
                optionsnotebook.tabSpacecraft.txtLV_adapter_mass.Show(True)
            else:
                optionsnotebook.tabSpacecraft.lblLV_margin.Show(False)
                optionsnotebook.tabSpacecraft.txtLV_margin.Show(False)
                optionsnotebook.tabSpacecraft.lblLV_adapter_mass.Show(False)
                optionsnotebook.tabSpacecraft.txtLV_adapter_mass.Show(False)

        #impulsive vs low-thrust missions
        if self.mission_type == 0 or self.mission_type == 1 or self.mission_type == 4:
            #impulsive mission
            optionsnotebook.tabSpacecraft.powergridtitle.Show(False)
            optionsnotebook.tabSpacecraft.lblIspChem.Show(True)
            optionsnotebook.tabSpacecraft.lblpower_margin.Show(False)
            optionsnotebook.tabSpacecraft.lblengine_type.Show(False)
            optionsnotebook.tabSpacecraft.lblnumber_of_engines.Show(False)
            optionsnotebook.tabSpacecraft.lblengine_duty_cycle.Show(False)
            optionsnotebook.tabSpacecraft.lblThrust.Show(False)
            optionsnotebook.tabSpacecraft.lblIspLT.Show(False)
            optionsnotebook.tabSpacecraft.lblIspLT_minimum.Show(False)
            optionsnotebook.tabSpacecraft.lbluser_defined_engine_efficiency.Show(False)
            optionsnotebook.tabSpacecraft.lblengine_input_thrust_coefficients.Show(False)
            optionsnotebook.tabSpacecraft.lblengine_input_mass_flow_rate_coefficients.Show(False)
            optionsnotebook.tabSpacecraft.lblengine_input_power_bounds.Show(False)
            optionsnotebook.tabSpacecraft.lblpower_at_1_AU.Show(False)
            optionsnotebook.tabSpacecraft.lblpower_source_type.Show(False)
            optionsnotebook.tabSpacecraft.lblsolar_power_gamma.Show(False)
            optionsnotebook.tabSpacecraft.lblspacecraft_power_model_type.Show(False)
            optionsnotebook.tabSpacecraft.lblspacecraft_power_coefficients.Show(False)
            optionsnotebook.tabSpacecraft.lblpower_decay_rate.Show(False)
            optionsnotebook.tabSpacecraft.txtIspChem.Show(True)
            optionsnotebook.tabSpacecraft.txtpower_margin.Show(False)
            optionsnotebook.tabSpacecraft.cmbengine_type.Show(False)
            optionsnotebook.tabSpacecraft.txtnumber_of_engines.Show(False)
            optionsnotebook.tabSpacecraft.txtengine_duty_cycle.Show(False)
            optionsnotebook.tabSpacecraft.txtThrust.Show(False)
            optionsnotebook.tabSpacecraft.txtIspLT.Show(False)
            optionsnotebook.tabSpacecraft.txtIspLT_minimum.Show(False)
            optionsnotebook.tabSpacecraft.txtuser_defined_engine_efficiency.Show(False)
            optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients0.Show(False)
            optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients1.Show(False)
            optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients2.Show(False)
            optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients3.Show(False)
            optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients4.Show(False)
            optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients5.Show(False)
            optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients6.Show(False)
            optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients0.Show(False)
            optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients1.Show(False)
            optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients2.Show(False)
            optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients3.Show(False)
            optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients4.Show(False)
            optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients5.Show(False)
            optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients6.Show(False)
            optionsnotebook.tabSpacecraft.txtengine_input_power_bounds_lower.Show(False)
            optionsnotebook.tabSpacecraft.txtengine_input_power_bounds_upper.Show(False)
            optionsnotebook.tabSpacecraft.txtpower_at_1_AU.Show(False)
            optionsnotebook.tabSpacecraft.cmbpower_source_type.Show(False)
            optionsnotebook.tabSpacecraft.txtsolar_power_gamma0.Show(False)
            optionsnotebook.tabSpacecraft.txtsolar_power_gamma1.Show(False)
            optionsnotebook.tabSpacecraft.txtsolar_power_gamma2.Show(False)
            optionsnotebook.tabSpacecraft.txtsolar_power_gamma3.Show(False)
            optionsnotebook.tabSpacecraft.txtsolar_power_gamma4.Show(False)
            optionsnotebook.tabSpacecraft.cmbspacecraft_power_model_type.Show(False)
            optionsnotebook.tabSpacecraft.txtspacecraft_power_coefficients0.Show(False)
            optionsnotebook.tabSpacecraft.txtspacecraft_power_coefficients1.Show(False)
            optionsnotebook.tabSpacecraft.txtspacecraft_power_coefficients2.Show(False)
            optionsnotebook.tabSpacecraft.txtpower_decay_rate.Show(False)
            optionsnotebook.tabSpacecraft.lblthrottle_logic_mode.Show(False)
            optionsnotebook.tabSpacecraft.cmbthrottle_logic_mode.Show(False)
            optionsnotebook.tabSpacecraft.lblthrottle_sharpness.Show(False)
            optionsnotebook.tabSpacecraft.txtthrottle_sharpness.Show(False)
        else:
            #low-thrust missions
            optionsnotebook.tabSpacecraft.lblengine_type.Show(True)
            optionsnotebook.tabSpacecraft.cmbengine_type.Show(True)

            if self.engine_type == 0:
                #fixed thrust/Isp, no power information required
                optionsnotebook.tabSpacecraft.powergridtitle.Show(False)
                optionsnotebook.tabSpacecraft.lblIspChem.Show(False)
                optionsnotebook.tabSpacecraft.lblpower_margin.Show(False)
                optionsnotebook.tabSpacecraft.lblnumber_of_engines.Show(False)
                optionsnotebook.tabSpacecraft.lblengine_duty_cycle.Show(True)
                optionsnotebook.tabSpacecraft.lblThrust.Show(True)
                optionsnotebook.tabSpacecraft.lblIspLT.Show(True)
                optionsnotebook.tabSpacecraft.lblIspLT_minimum.Show(False)
                optionsnotebook.tabSpacecraft.lbluser_defined_engine_efficiency.Show(False)
                optionsnotebook.tabSpacecraft.lblengine_input_thrust_coefficients.Show(False)
                optionsnotebook.tabSpacecraft.lblengine_input_mass_flow_rate_coefficients.Show(False)
                optionsnotebook.tabSpacecraft.lblengine_input_power_bounds.Show(False)
                optionsnotebook.tabSpacecraft.lblpower_at_1_AU.Show(False)
                optionsnotebook.tabSpacecraft.lblpower_source_type.Show(False)
                optionsnotebook.tabSpacecraft.lblsolar_power_gamma.Show(False)
                optionsnotebook.tabSpacecraft.lblspacecraft_power_model_type.Show(False)
                optionsnotebook.tabSpacecraft.lblspacecraft_power_coefficients.Show(False)
                optionsnotebook.tabSpacecraft.lblpower_decay_rate.Show(False)
                optionsnotebook.tabSpacecraft.txtIspChem.Show(False)
                optionsnotebook.tabSpacecraft.txtpower_margin.Show(False)
                optionsnotebook.tabSpacecraft.txtnumber_of_engines.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_duty_cycle.Show(True)
                optionsnotebook.tabSpacecraft.txtThrust.Show(True)
                optionsnotebook.tabSpacecraft.txtIspLT.Show(True)
                optionsnotebook.tabSpacecraft.txtIspLT_minimum.Show(False)
                optionsnotebook.tabSpacecraft.txtuser_defined_engine_efficiency.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients0.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients1.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients2.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients3.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients4.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients5.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients6.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients0.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients1.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients2.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients3.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients4.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients5.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients6.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_power_bounds_lower.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_power_bounds_upper.Show(False)
                optionsnotebook.tabSpacecraft.txtpower_at_1_AU.Show(False)
                optionsnotebook.tabSpacecraft.cmbpower_source_type.Show(False)
                optionsnotebook.tabSpacecraft.txtsolar_power_gamma0.Show(False)
                optionsnotebook.tabSpacecraft.txtsolar_power_gamma1.Show(False)
                optionsnotebook.tabSpacecraft.txtsolar_power_gamma2.Show(False)
                optionsnotebook.tabSpacecraft.txtsolar_power_gamma3.Show(False)
                optionsnotebook.tabSpacecraft.txtsolar_power_gamma4.Show(False)
                optionsnotebook.tabSpacecraft.cmbspacecraft_power_model_type.Show(False)
                optionsnotebook.tabSpacecraft.txtspacecraft_power_coefficients0.Show(False)
                optionsnotebook.tabSpacecraft.txtspacecraft_power_coefficients1.Show(False)
                optionsnotebook.tabSpacecraft.txtspacecraft_power_coefficients2.Show(False)
                optionsnotebook.tabSpacecraft.txtpower_decay_rate.Show(False)
                optionsnotebook.tabSpacecraft.lblthrottle_logic_mode.Show(False)
                optionsnotebook.tabSpacecraft.cmbthrottle_logic_mode.Show(False)
            elif self.engine_type == 1:
                #constant Isp, efficiency, EMTG computes input power
                #do not need anything except Isp, efficiency, and duty cycle
                optionsnotebook.tabSpacecraft.powergridtitle.Show(False)
                optionsnotebook.tabSpacecraft.lblIspChem.Show(False)
                optionsnotebook.tabSpacecraft.lblpower_margin.Show(False)
                optionsnotebook.tabSpacecraft.lblnumber_of_engines.Show(False)
                optionsnotebook.tabSpacecraft.lblengine_duty_cycle.Show(True)
                optionsnotebook.tabSpacecraft.lblThrust.Show(False)
                optionsnotebook.tabSpacecraft.lblIspLT.Show(True)
                optionsnotebook.tabSpacecraft.lblIspLT_minimum.Show(False)
                optionsnotebook.tabSpacecraft.lbluser_defined_engine_efficiency.Show(True)
                optionsnotebook.tabSpacecraft.lblengine_input_thrust_coefficients.Show(False)
                optionsnotebook.tabSpacecraft.lblengine_input_mass_flow_rate_coefficients.Show(False)
                optionsnotebook.tabSpacecraft.lblengine_input_power_bounds.Show(False)
                optionsnotebook.tabSpacecraft.lblpower_at_1_AU.Show(False)
                optionsnotebook.tabSpacecraft.lblpower_source_type.Show(False)
                optionsnotebook.tabSpacecraft.lblsolar_power_gamma.Show(False)
                optionsnotebook.tabSpacecraft.lblspacecraft_power_model_type.Show(False)
                optionsnotebook.tabSpacecraft.lblspacecraft_power_coefficients.Show(False)
                optionsnotebook.tabSpacecraft.lblpower_decay_rate.Show(False)
                optionsnotebook.tabSpacecraft.txtIspChem.Show(False)
                optionsnotebook.tabSpacecraft.txtpower_margin.Show(False)
                optionsnotebook.tabSpacecraft.txtnumber_of_engines.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_duty_cycle.Show(True)
                optionsnotebook.tabSpacecraft.txtThrust.Show(False)
                optionsnotebook.tabSpacecraft.txtIspLT.Show(True)
                optionsnotebook.tabSpacecraft.txtIspLT_minimum.Show(False)
                optionsnotebook.tabSpacecraft.txtuser_defined_engine_efficiency.Show(True)
                optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients0.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients1.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients2.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients3.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients4.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients5.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients6.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients0.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients1.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients2.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients3.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients4.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients5.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients6.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_power_bounds_lower.Show(False)
                optionsnotebook.tabSpacecraft.txtengine_input_power_bounds_upper.Show(False)
                optionsnotebook.tabSpacecraft.txtpower_at_1_AU.Show(False)
                optionsnotebook.tabSpacecraft.cmbpower_source_type.Show(False)
                optionsnotebook.tabSpacecraft.txtsolar_power_gamma0.Show(False)
                optionsnotebook.tabSpacecraft.txtsolar_power_gamma1.Show(False)
                optionsnotebook.tabSpacecraft.txtsolar_power_gamma2.Show(False)
                optionsnotebook.tabSpacecraft.txtsolar_power_gamma3.Show(False)
                optionsnotebook.tabSpacecraft.txtsolar_power_gamma4.Show(False)
                optionsnotebook.tabSpacecraft.cmbspacecraft_power_model_type.Show(False)
                optionsnotebook.tabSpacecraft.txtspacecraft_power_coefficients0.Show(False)
                optionsnotebook.tabSpacecraft.txtspacecraft_power_coefficients1.Show(False)
                optionsnotebook.tabSpacecraft.txtspacecraft_power_coefficients2.Show(False)
                optionsnotebook.tabSpacecraft.txtpower_decay_rate.Show(False)
                optionsnotebook.tabSpacecraft.lblthrottle_logic_mode.Show(False)
                optionsnotebook.tabSpacecraft.cmbthrottle_logic_mode.Show(False)
            else:
                #engine types 2 and greater all require power information but NOT chemical Isp information
                optionsnotebook.tabSpacecraft.powergridtitle.Show(True)
                optionsnotebook.tabSpacecraft.lblIspChem.Show(False)
                optionsnotebook.tabSpacecraft.lblpower_margin.Show(True)
                optionsnotebook.tabSpacecraft.lblpower_at_1_AU.Show(True)
                optionsnotebook.tabSpacecraft.lblpower_source_type.Show(True)
                optionsnotebook.tabSpacecraft.lblspacecraft_power_model_type.Show(True)
                optionsnotebook.tabSpacecraft.lblspacecraft_power_coefficients.Show(True)
                optionsnotebook.tabSpacecraft.lblpower_decay_rate.Show(True)
                optionsnotebook.tabSpacecraft.txtIspChem.Show(False)
                optionsnotebook.tabSpacecraft.txtpower_margin.Show(True)
                optionsnotebook.tabSpacecraft.txtpower_at_1_AU.Show(True)
                optionsnotebook.tabSpacecraft.cmbpower_source_type.Show(True)
                if self.power_source_type == 0:
                    optionsnotebook.tabSpacecraft.lblpower_at_1_AU.SetLabel("Power at BOL, 1 AU (kW)")
                    optionsnotebook.tabSpacecraft.lblsolar_power_gamma.Show(True)
                    optionsnotebook.tabSpacecraft.txtsolar_power_gamma0.Show(True)
                    optionsnotebook.tabSpacecraft.txtsolar_power_gamma1.Show(True)
                    optionsnotebook.tabSpacecraft.txtsolar_power_gamma2.Show(True)
                    optionsnotebook.tabSpacecraft.txtsolar_power_gamma3.Show(True)
                    optionsnotebook.tabSpacecraft.txtsolar_power_gamma4.Show(True)
                else:
                    optionsnotebook.tabSpacecraft.lblpower_at_1_AU.SetLabel("Power at BOL (kW)")
                    optionsnotebook.tabSpacecraft.lblsolar_power_gamma.Show(False)
                    optionsnotebook.tabSpacecraft.txtsolar_power_gamma0.Show(False)
                    optionsnotebook.tabSpacecraft.txtsolar_power_gamma1.Show(False)
                    optionsnotebook.tabSpacecraft.txtsolar_power_gamma2.Show(False)
                    optionsnotebook.tabSpacecraft.txtsolar_power_gamma3.Show(False)
                    optionsnotebook.tabSpacecraft.txtsolar_power_gamma4.Show(False)

                optionsnotebook.tabSpacecraft.cmbspacecraft_power_model_type.Show(True)
                optionsnotebook.tabSpacecraft.txtspacecraft_power_coefficients0.Show(True)
                optionsnotebook.tabSpacecraft.txtspacecraft_power_coefficients1.Show(True)
                optionsnotebook.tabSpacecraft.txtspacecraft_power_coefficients2.Show(True)
                optionsnotebook.tabSpacecraft.txtpower_decay_rate.Show(True)

                if self.engine_type == 2:
                    #choice of power model, constant efficiency, EMTG chooses Isp
                    #all other options off
                    optionsnotebook.tabSpacecraft.lblnumber_of_engines.Show(False)
                    optionsnotebook.tabSpacecraft.lblengine_duty_cycle.Show(True)
                    optionsnotebook.tabSpacecraft.lblThrust.Show(False)
                    optionsnotebook.tabSpacecraft.lblIspLT.Show(True)
                    optionsnotebook.tabSpacecraft.lblIspLT_minimum.Show(True)
                    optionsnotebook.tabSpacecraft.lbluser_defined_engine_efficiency.Show(True)
                    optionsnotebook.tabSpacecraft.lblengine_input_thrust_coefficients.Show(False)
                    optionsnotebook.tabSpacecraft.lblengine_input_mass_flow_rate_coefficients.Show(False)
                    optionsnotebook.tabSpacecraft.lblengine_input_power_bounds.Show(True)
                    optionsnotebook.tabSpacecraft.txtnumber_of_engines.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_duty_cycle.Show(True)
                    optionsnotebook.tabSpacecraft.txtThrust.Show(False)
                    optionsnotebook.tabSpacecraft.txtIspLT.Show(True)
                    optionsnotebook.tabSpacecraft.txtIspLT_minimum.Show(True)
                    optionsnotebook.tabSpacecraft.txtuser_defined_engine_efficiency.Show(True)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients0.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients1.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients2.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients3.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients4.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients5.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients6.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients0.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients1.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients2.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients3.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients4.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients5.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients6.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_power_bounds_lower.Show(True)
                    optionsnotebook.tabSpacecraft.txtengine_input_power_bounds_upper.Show(True)
                    optionsnotebook.tabSpacecraft.lblthrottle_logic_mode.Show(False)
                    optionsnotebook.tabSpacecraft.cmbthrottle_logic_mode.Show(False)  
                    optionsnotebook.tabSpacecraft.lblthrottle_sharpness.Show(False)
                    optionsnotebook.tabSpacecraft.txtthrottle_sharpness.Show(False)

                elif self.engine_type == 3:
                    #choice of power model, constant efficiency and Isp
                    #all other options off
                    optionsnotebook.tabSpacecraft.lblnumber_of_engines.Show(False)
                    optionsnotebook.tabSpacecraft.lblengine_duty_cycle.Show(True)
                    optionsnotebook.tabSpacecraft.lblThrust.Show(False)
                    optionsnotebook.tabSpacecraft.lblIspLT.Show(True)
                    optionsnotebook.tabSpacecraft.lblIspLT_minimum.Show(False)
                    optionsnotebook.tabSpacecraft.lbluser_defined_engine_efficiency.Show(True)
                    optionsnotebook.tabSpacecraft.lblengine_input_thrust_coefficients.Show(False)
                    optionsnotebook.tabSpacecraft.lblengine_input_mass_flow_rate_coefficients.Show(False)
                    optionsnotebook.tabSpacecraft.lblengine_input_power_bounds.Show(True)
                    optionsnotebook.tabSpacecraft.txtnumber_of_engines.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_duty_cycle.Show(True)
                    optionsnotebook.tabSpacecraft.txtThrust.Show(False)
                    optionsnotebook.tabSpacecraft.txtIspLT.Show(True)
                    optionsnotebook.tabSpacecraft.txtIspLT_minimum.Show(False)
                    optionsnotebook.tabSpacecraft.txtuser_defined_engine_efficiency.Show(True)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients0.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients1.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients2.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients3.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients4.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients5.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients6.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients0.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients1.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients2.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients3.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients4.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients5.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients6.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_power_bounds_lower.Show(True)
                    optionsnotebook.tabSpacecraft.txtengine_input_power_bounds_upper.Show(True)
                    optionsnotebook.tabSpacecraft.lblthrottle_logic_mode.Show(False)
                    optionsnotebook.tabSpacecraft.cmbthrottle_logic_mode.Show(False)
                    optionsnotebook.tabSpacecraft.lblthrottle_sharpness.Show(False)
                    optionsnotebook.tabSpacecraft.txtthrottle_sharpness.Show(False)
                
                elif self.engine_type == 4:
                    #continuously-varying specific impulse
                    #requires efficiency, Isp min and max, duty cycle
                    optionsnotebook.tabSpacecraft.lblnumber_of_engines.Show(False)
                    optionsnotebook.tabSpacecraft.lblengine_duty_cycle.Show(True)
                    optionsnotebook.tabSpacecraft.lblThrust.Show(False)
                    optionsnotebook.tabSpacecraft.lblIspLT.Show(True)
                    optionsnotebook.tabSpacecraft.lblIspLT_minimum.Show(True)
                    optionsnotebook.tabSpacecraft.lbluser_defined_engine_efficiency.Show(True)
                    optionsnotebook.tabSpacecraft.lblengine_input_thrust_coefficients.Show(False)
                    optionsnotebook.tabSpacecraft.lblengine_input_mass_flow_rate_coefficients.Show(False)
                    optionsnotebook.tabSpacecraft.lblengine_input_power_bounds.Show(False)
                    optionsnotebook.tabSpacecraft.txtnumber_of_engines.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_duty_cycle.Show(True)
                    optionsnotebook.tabSpacecraft.txtThrust.Show(False)
                    optionsnotebook.tabSpacecraft.txtIspLT.Show(True)
                    optionsnotebook.tabSpacecraft.txtIspLT_minimum.Show(True)
                    optionsnotebook.tabSpacecraft.txtuser_defined_engine_efficiency.Show(True)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients0.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients1.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients2.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients3.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients4.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients5.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients6.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients0.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients1.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients2.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients3.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients4.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients5.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients6.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_power_bounds_lower.Show(True)
                    optionsnotebook.tabSpacecraft.txtengine_input_power_bounds_upper.Show(True)
                    optionsnotebook.tabSpacecraft.lblthrottle_logic_mode.Show(False)
                    optionsnotebook.tabSpacecraft.cmbthrottle_logic_mode.Show(False)
                    optionsnotebook.tabSpacecraft.lblthrottle_sharpness.Show(False)
                    optionsnotebook.tabSpacecraft.txtthrottle_sharpness.Show(False)

                elif self.engine_type == 5:
                    #custom thrust and mass flow rate polynomial
                    optionsnotebook.tabSpacecraft.lblnumber_of_engines.Show(True)
                    optionsnotebook.tabSpacecraft.lblengine_duty_cycle.Show(True)
                    optionsnotebook.tabSpacecraft.lblThrust.Show(False)
                    optionsnotebook.tabSpacecraft.lblIspLT.Show(False)
                    optionsnotebook.tabSpacecraft.lblIspLT_minimum.Show(False)
                    optionsnotebook.tabSpacecraft.lbluser_defined_engine_efficiency.Show(False)
                    optionsnotebook.tabSpacecraft.lblengine_input_thrust_coefficients.Show(True)
                    optionsnotebook.tabSpacecraft.lblengine_input_mass_flow_rate_coefficients.Show(True)
                    optionsnotebook.tabSpacecraft.lblengine_input_power_bounds.Show(True)
                    optionsnotebook.tabSpacecraft.txtnumber_of_engines.Show(True)
                    optionsnotebook.tabSpacecraft.txtengine_duty_cycle.Show(True)
                    optionsnotebook.tabSpacecraft.txtThrust.Show(False)
                    optionsnotebook.tabSpacecraft.txtIspLT.Show(False)
                    optionsnotebook.tabSpacecraft.txtIspLT_minimum.Show(False)
                    optionsnotebook.tabSpacecraft.txtuser_defined_engine_efficiency.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients0.Show(True)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients1.Show(True)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients2.Show(True)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients3.Show(True)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients4.Show(True)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients5.Show(True)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients6.Show(True)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients0.Show(True)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients1.Show(True)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients2.Show(True)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients3.Show(True)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients4.Show(True)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients5.Show(True)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients6.Show(True)
                    optionsnotebook.tabSpacecraft.txtengine_input_power_bounds_lower.Show(True)
                    optionsnotebook.tabSpacecraft.txtengine_input_power_bounds_upper.Show(True)
                    optionsnotebook.tabSpacecraft.lblthrottle_logic_mode.Show(True)
                    optionsnotebook.tabSpacecraft.cmbthrottle_logic_mode.Show(True)
                    optionsnotebook.tabSpacecraft.lblthrottle_sharpness.Show(True)
                    optionsnotebook.tabSpacecraft.txtthrottle_sharpness.Show(True)
                else:
                    #hard-coded thrust model
                    optionsnotebook.tabSpacecraft.lblnumber_of_engines.Show(True)
                    optionsnotebook.tabSpacecraft.lblengine_duty_cycle.Show(True)
                    optionsnotebook.tabSpacecraft.lblThrust.Show(False)
                    optionsnotebook.tabSpacecraft.lblIspLT.Show(False)
                    optionsnotebook.tabSpacecraft.lblIspLT_minimum.Show(False)
                    optionsnotebook.tabSpacecraft.lbluser_defined_engine_efficiency.Show(False)
                    optionsnotebook.tabSpacecraft.lblengine_input_thrust_coefficients.Show(False)
                    optionsnotebook.tabSpacecraft.lblengine_input_mass_flow_rate_coefficients.Show(False)
                    optionsnotebook.tabSpacecraft.lblengine_input_power_bounds.Show(False)
                    optionsnotebook.tabSpacecraft.txtnumber_of_engines.Show(True)
                    optionsnotebook.tabSpacecraft.txtengine_duty_cycle.Show(True)
                    optionsnotebook.tabSpacecraft.txtThrust.Show(False)
                    optionsnotebook.tabSpacecraft.txtIspLT.Show(False)
                    optionsnotebook.tabSpacecraft.txtIspLT_minimum.Show(False)
                    optionsnotebook.tabSpacecraft.txtuser_defined_engine_efficiency.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients0.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients1.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients2.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients3.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients4.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients5.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_thrust_coefficients6.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients0.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients1.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients2.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients3.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients4.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients5.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_mass_flow_rate_coefficients6.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_power_bounds_lower.Show(False)
                    optionsnotebook.tabSpacecraft.txtengine_input_power_bounds_upper.Show(False)
                    optionsnotebook.tabSpacecraft.lblthrottle_logic_mode.Show(True)
                    optionsnotebook.tabSpacecraft.cmbthrottle_logic_mode.Show(True)
                    optionsnotebook.tabSpacecraft.lblthrottle_sharpness.Show(True)
                    optionsnotebook.tabSpacecraft.txtthrottle_sharpness.Show(True)

        #re-size the panel
        optionsnotebook.tabSpacecraft.Layout()
        optionsnotebook.tabSpacecraft.SetupScrolling()




    def update_solver_options_panel(self, optionsnotebook):

        if self.snopt_max_run_time > self.MBH_max_run_time:
            self.snopt_max_run_time = self.MBH_max_run_time - 1

        #inner-loop solver options
        optionsnotebook.tabSolver.cmbInnerLoopSolver.SetSelection(self.run_inner_loop)
        optionsnotebook.tabSolver.cmbNLP_solver_type.SetSelection(self.NLP_solver_type)
        optionsnotebook.tabSolver.cmbNLP_solver_mode.SetSelection(self.NLP_solver_mode)
        optionsnotebook.tabSolver.chkquiet_NLP.SetValue(self.quiet_NLP)
        optionsnotebook.tabSolver.chkquiet_MBH.SetValue(self.quiet_basinhopping)
        optionsnotebook.tabSolver.chkMBH_two_step.SetValue(self.MBH_two_step)
        optionsnotebook.tabSolver.txtFD_stepsize.SetValue(str(self.FD_stepsize))
        optionsnotebook.tabSolver.txtFD_stepsize_coarse.SetValue(str(self.FD_stepsize_coarse))
        optionsnotebook.tabSolver.chkACE_feasible_point_finder.SetValue(self.ACE_feasible_point_finder)
        optionsnotebook.tabSolver.txtMBH_max_not_improve.SetValue(str(self.MBH_max_not_improve))
        optionsnotebook.tabSolver.txtMBH_max_trials.SetValue(str(self.MBH_max_trials))
        optionsnotebook.tabSolver.txtMBH_max_run_time.SetValue(str(self.MBH_max_run_time))
        optionsnotebook.tabSolver.txtMBH_max_step_size.SetValue(str(self.MBH_max_step_size))
        optionsnotebook.tabSolver.cmbMBH_hop_distribution.SetSelection(self.MBH_hop_distribution)
        optionsnotebook.tabSolver.txtMBH_Pareto_alpha.SetValue(str(self.MBH_Pareto_alpha))
        optionsnotebook.tabSolver.txtMBH_time_hop_probability.SetValue(str(self.MBH_time_hop_probability))
        optionsnotebook.tabSolver.txtsnopt_feasibility_tolerance.SetValue(str(self.snopt_feasibility_tolerance))
        optionsnotebook.tabSolver.txtsnopt_major_iterations.SetValue(str(self.snopt_major_iterations))
        optionsnotebook.tabSolver.txtsnopt_max_run_time.SetValue(str(self.snopt_max_run_time))
        optionsnotebook.tabSolver.cmbderivative_type.SetSelection(self.derivative_type)
        optionsnotebook.tabSolver.chkcheck_derivatives.SetValue(self.check_derivatives)
        optionsnotebook.tabSolver.chkseed_MBH.SetValue(self.seed_MBH)
        optionsnotebook.tabSolver.cmbinitial_guess_control_coordinate_system.SetSelection(self.initial_guess_control_coordinate_system)
        optionsnotebook.tabSolver.chkinterpolate_initial_guess.SetValue(self.interpolate_initial_guess)
        optionsnotebook.tabSolver.txtinitial_guess_num_timesteps.SetValue(str(self.initial_guess_num_timesteps))
        optionsnotebook.tabSolver.cmbinitial_guess_step_size_distribution.SetSelection(self.initial_guess_step_size_distribution)
        optionsnotebook.tabSolver.txtinitial_guess_step_size_stdv_or_scale.SetValue(str(self.initial_guess_step_size_stdv_or_scale))
        optionsnotebook.tabSolver.cmbMBH_zero_control_initial_guess.SetSelection(self.MBH_zero_control_initial_guess)
        optionsnotebook.tabSolver.txttrialX.SetValue(str(self.trialX))

        if self.run_inner_loop == 2:
            optionsnotebook.tabSolver.lblMBH_two_step.Show(True)
            optionsnotebook.tabSolver.chkMBH_two_step.Show(True)

            if self.MBH_two_step:
                optionsnotebook.tabSolver.lblFD_stepsize_coarse.Show(True)
                optionsnotebook.tabSolver.txtFD_stepsize_coarse.Show(True)
            else:
                optionsnotebook.tabSolver.lblFD_stepsize_coarse.Show(False)
                optionsnotebook.tabSolver.txtFD_stepsize_coarse.Show(False)
        else:
            optionsnotebook.tabSolver.lblMBH_two_step.Show(False)
            optionsnotebook.tabSolver.chkMBH_two_step.Show(False)
            optionsnotebook.tabSolver.lblFD_stepsize_coarse.Show(False)
            optionsnotebook.tabSolver.txtFD_stepsize_coarse.Show(False)

        if self.derivative_type < 3 and (self.run_inner_loop == 2 or self.run_inner_loop == 4):
            optionsnotebook.tabSolver.lblFD_stepsize.Show(True)
            optionsnotebook.tabSolver.txtFD_stepsize.Show(True)
        else:
            optionsnotebook.tabSolver.lblFD_stepsize.Show(False)
            optionsnotebook.tabSolver.txtFD_stepsize.Show(False)
        
        if self.run_inner_loop == 0: #trialX
            optionsnotebook.tabSolver.lblMBH_max_not_improve.Show(False)
            optionsnotebook.tabSolver.lblMBH_max_trials.Show(False)
            optionsnotebook.tabSolver.lblMBH_max_run_time.Show(False)
            optionsnotebook.tabSolver.lblMBH_max_step_size.Show(False)
            optionsnotebook.tabSolver.lblMBH_time_hop_probability.Show(False)
            optionsnotebook.tabSolver.lblsnopt_feasibility_tolerance.Show(False)
            optionsnotebook.tabSolver.lblsnopt_major_iterations.Show(False)
            optionsnotebook.tabSolver.lblsnopt_max_run_time.Show(False)
            optionsnotebook.tabSolver.lblderivative_type.Show(False)
            optionsnotebook.tabSolver.lblcheck_derivatives.Show(True)
            optionsnotebook.tabSolver.lblseed_MBH.Show(False)
            optionsnotebook.tabSolver.lblMBH_hop_distribution.Show(False)
            optionsnotebook.tabSolver.lblNLP_solver_type.Show(False)
            optionsnotebook.tabSolver.lblNLP_solver_mode.Show(False)
            optionsnotebook.tabSolver.lblquiet_NLP.Show(False)
            optionsnotebook.tabSolver.lblACE_feasible_point_finder.Show(False)
            
            optionsnotebook.tabSolver.txtMBH_max_not_improve.Show(False)
            optionsnotebook.tabSolver.txtMBH_max_trials.Show(False)
            optionsnotebook.tabSolver.txtMBH_max_run_time.Show(False)
            optionsnotebook.tabSolver.txtMBH_max_step_size.Show(False)
            optionsnotebook.tabSolver.txtMBH_time_hop_probability.Show(False)
            optionsnotebook.tabSolver.txtsnopt_feasibility_tolerance.Show(False)
            optionsnotebook.tabSolver.txtsnopt_major_iterations.Show(False)
            optionsnotebook.tabSolver.txtsnopt_max_run_time.Show(False)
            optionsnotebook.tabSolver.cmbderivative_type.Show(False)
            optionsnotebook.tabSolver.chkcheck_derivatives.Show(True)
            optionsnotebook.tabSolver.chkseed_MBH.Show(False)
            optionsnotebook.tabSolver.cmbMBH_hop_distribution.Show(False)
            optionsnotebook.tabSolver.cmbNLP_solver_type.Show(False)
            optionsnotebook.tabSolver.cmbNLP_solver_mode.Show(False)
            optionsnotebook.tabSolver.chkquiet_NLP.Show(False)
            optionsnotebook.tabSolver.chkACE_feasible_point_finder.Show(False)

            optionsnotebook.tabSolver.lblMBH_Pareto_alpha.Show(False)
            optionsnotebook.tabSolver.txtMBH_Pareto_alpha.Show(False)

            optionsnotebook.tabSolver.lblMBH_zero_control_initial_guess.Show(False)
            optionsnotebook.tabSolver.cmbMBH_zero_control_initial_guess.Show(False)

            optionsnotebook.tabSolver.lblquiet_MBH.Show(False)
            optionsnotebook.tabSolver.chkquiet_MBH.Show(False)

        elif self.run_inner_loop == 1: #batch trialX
            optionsnotebook.tabSolver.lblMBH_max_not_improve.Show(False)
            optionsnotebook.tabSolver.lblMBH_max_trials.Show(False)
            optionsnotebook.tabSolver.lblMBH_max_run_time.Show(False)
            optionsnotebook.tabSolver.lblMBH_max_step_size.Show(False)
            optionsnotebook.tabSolver.lblMBH_time_hop_probability.Show(False)
            optionsnotebook.tabSolver.lblsnopt_feasibility_tolerance.Show(False)
            optionsnotebook.tabSolver.lblsnopt_major_iterations.Show(False)
            optionsnotebook.tabSolver.lblsnopt_max_run_time.Show(False)
            optionsnotebook.tabSolver.lblderivative_type.Show(False)
            optionsnotebook.tabSolver.lblcheck_derivatives.Show(False)
            optionsnotebook.tabSolver.lblseed_MBH.Show(False)
            optionsnotebook.tabSolver.lblMBH_hop_distribution.Show(False)
            optionsnotebook.tabSolver.lblNLP_solver_type.Show(False)
            optionsnotebook.tabSolver.lblNLP_solver_mode.Show(False)
            optionsnotebook.tabSolver.lblquiet_NLP.show(False)
            optionsnotebook.tabSolver.lblACE_feasible_point_finder.Show(False)
            
            optionsnotebook.tabSolver.txtMBH_max_not_improve.Show(False)
            optionsnotebook.tabSolver.txtMBH_max_trials.Show(False)
            optionsnotebook.tabSolver.txtMBH_max_run_time.Show(False)
            optionsnotebook.tabSolver.txtMBH_max_step_size.Show(False)
            optionsnotebook.tabSolver.txtMBH_time_hop_probability.Show(False)
            optionsnotebook.tabSolver.txtsnopt_feasibility_tolerance.Show(False)
            optionsnotebook.tabSolver.txtsnopt_major_iterations.Show(False)
            optionsnotebook.tabSolver.txtsnopt_max_run_time.Show(False)
            optionsnotebook.tabSolver.cmbderivative_type.Show(False)
            optionsnotebook.tabSolver.chkcheck_derivatives.Show(False)
            optionsnotebook.tabSolver.chkseed_MBH.Show(False)
            optionsnotebook.tabSolver.cmbMBH_hop_distribution.Show(False)
            optionsnotebook.tabSolver.cmbNLP_solver_type.Show(False)
            optionsnotebook.tabSolver.cmbNLP_solver_mode.Show(False)
            optionsnotebook.tabSolver.chkquiet_NLP.Show(False)
            optionsnotebook.tabSolver.chkACE_feasible_point_finder.Show(False)

            optionsnotebook.tabSolver.lblMBH_Pareto_alpha.Show(False)
            optionsnotebook.tabSolver.txtMBH_Pareto_alpha.Show(False)

            optionsnotebook.tabSolver.lblMBH_zero_control_initial_guess.Show(False)
            optionsnotebook.tabSolver.cmbMBH_zero_control_initial_guess.Show(False)

            optionsnotebook.tabSolver.lblquiet_MBH.Show(False)
            optionsnotebook.tabSolver.chkquiet_MBH.Show(False)
            
        elif self.run_inner_loop == 2: #MBH
            optionsnotebook.tabSolver.lblMBH_max_not_improve.Show(True)
            optionsnotebook.tabSolver.lblMBH_max_trials.Show(True)
            optionsnotebook.tabSolver.lblMBH_max_run_time.Show(True)
            optionsnotebook.tabSolver.lblMBH_max_step_size.Show(True)
            optionsnotebook.tabSolver.lblMBH_time_hop_probability.Show(True)
            optionsnotebook.tabSolver.lblsnopt_feasibility_tolerance.Show(True)
            optionsnotebook.tabSolver.lblsnopt_major_iterations.Show(True)
            optionsnotebook.tabSolver.lblsnopt_max_run_time.Show(True)
            optionsnotebook.tabSolver.lblderivative_type.Show(True)
            optionsnotebook.tabSolver.lblcheck_derivatives.Show(True)
            optionsnotebook.tabSolver.lblseed_MBH.Show(True)
            optionsnotebook.tabSolver.lblMBH_hop_distribution.Show(True)
            optionsnotebook.tabSolver.lblNLP_solver_type.Show(True)
            optionsnotebook.tabSolver.lblNLP_solver_mode.Show(True)
            optionsnotebook.tabSolver.lblquiet_NLP.Show(True)
            optionsnotebook.tabSolver.lblACE_feasible_point_finder.Show(True)
            
            optionsnotebook.tabSolver.txtMBH_max_not_improve.Show(True)
            optionsnotebook.tabSolver.txtMBH_max_trials.Show(True)
            optionsnotebook.tabSolver.txtMBH_max_run_time.Show(True)
            optionsnotebook.tabSolver.txtMBH_max_step_size.Show(True)
            optionsnotebook.tabSolver.txtMBH_time_hop_probability.Show(True)
            optionsnotebook.tabSolver.txtsnopt_feasibility_tolerance.Show(True)
            optionsnotebook.tabSolver.txtsnopt_major_iterations.Show(True)
            optionsnotebook.tabSolver.txtsnopt_max_run_time.Show(True)
            optionsnotebook.tabSolver.cmbderivative_type.Show(True)
            optionsnotebook.tabSolver.chkcheck_derivatives.Show(True)
            optionsnotebook.tabSolver.chkseed_MBH.Show(True)
            optionsnotebook.tabSolver.cmbMBH_hop_distribution.Show(True)
            optionsnotebook.tabSolver.cmbNLP_solver_type.Show(True)
            optionsnotebook.tabSolver.cmbNLP_solver_mode.Show(True)
            optionsnotebook.tabSolver.chkquiet_NLP.Show(True)
            optionsnotebook.tabSolver.chkACE_feasible_point_finder.Show(True)

            optionsnotebook.tabSolver.lblMBH_zero_control_initial_guess.Show(True)
            optionsnotebook.tabSolver.cmbMBH_zero_control_initial_guess.Show(True)

            optionsnotebook.tabSolver.lblquiet_MBH.Show(True)
            optionsnotebook.tabSolver.chkquiet_MBH.Show(True)

            #change the available parameters and labels based on which distribution is selected
            if self.MBH_hop_distribution == 0: #uniform
                optionsnotebook.tabSolver.lblMBH_max_step_size.SetLabel("MBH uniform hop ball size")
                optionsnotebook.tabSolver.lblMBH_Pareto_alpha.Show(False)
                optionsnotebook.tabSolver.txtMBH_Pareto_alpha.Show(False)
            elif self.MBH_hop_distribution == 1: #Cauchy
                optionsnotebook.tabSolver.lblMBH_max_step_size.SetLabel("MBH hop scale factor")
                optionsnotebook.tabSolver.lblMBH_Pareto_alpha.Show(False)
                optionsnotebook.tabSolver.txtMBH_Pareto_alpha.Show(False)
            elif self.MBH_hop_distribution == 2: #Pareto
                optionsnotebook.tabSolver.lblMBH_max_step_size.SetLabel("MBH hop scale factor")
                optionsnotebook.tabSolver.lblMBH_Pareto_alpha.Show(True)
                optionsnotebook.tabSolver.txtMBH_Pareto_alpha.Show(True)
            elif self.MBH_hop_distribution == 3: #Gaussian
                optionsnotebook.tabSolver.lblMBH_max_step_size.SetLabel("MBH hop standard deviation")
                optionsnotebook.tabSolver.lblMBH_Pareto_alpha.Show(False)
                optionsnotebook.tabSolver.txtMBH_Pareto_alpha.Show(False)
            
        elif self.run_inner_loop == 3: #ACDE
            optionsnotebook.tabSolver.lblMBH_max_not_improve.Show(True)
            optionsnotebook.tabSolver.lblMBH_max_trials.Show(True)
            optionsnotebook.tabSolver.lblMBH_max_run_time.Show(True)
            optionsnotebook.tabSolver.lblMBH_max_step_size.Show(False)
            optionsnotebook.tabSolver.lblMBH_time_hop_probability.Show(False)
            optionsnotebook.tabSolver.lblsnopt_feasibility_tolerance.Show(True)
            optionsnotebook.tabSolver.lblsnopt_major_iterations.Show(False)
            optionsnotebook.tabSolver.lblsnopt_max_run_time.Show(False)
            optionsnotebook.tabSolver.lblderivative_type.Show(False)
            optionsnotebook.tabSolver.lblcheck_derivatives.Show(False)
            optionsnotebook.tabSolver.lblseed_MBH.Show(False)
            optionsnotebook.tabSolver.lblMBH_hop_distribution.Show(False)
            optionsnotebook.tabSolver.lblNLP_solver_type.Show(False)
            optionsnotebook.tabSolver.lblNLP_solver_mode.Show(False)
            optionsnotebook.tabSolver.lblquiet_NLP.Show(False)
            optionsnotebook.tabSolver.lblACE_feasible_point_finder.Show(False)
            
            optionsnotebook.tabSolver.txtMBH_max_not_improve.Show(True)
            optionsnotebook.tabSolver.txtMBH_max_trials.Show(True)
            optionsnotebook.tabSolver.txtMBH_max_run_time.Show(True)
            optionsnotebook.tabSolver.txtMBH_max_step_size.Show(False)
            optionsnotebook.tabSolver.txtMBH_time_hop_probability.Show(False)
            optionsnotebook.tabSolver.txtsnopt_feasibility_tolerance.Show(True)
            optionsnotebook.tabSolver.txtsnopt_major_iterations.Show(False)
            optionsnotebook.tabSolver.txtsnopt_max_run_time.Show(False)
            optionsnotebook.tabSolver.cmbderivative_type.Show(False)
            optionsnotebook.tabSolver.chkcheck_derivatives.Show(False)
            optionsnotebook.tabSolver.chkseed_MBH.Show(False)
            optionsnotebook.tabSolver.cmbMBH_hop_distribution.Show(False)
            optionsnotebook.tabSolver.cmbNLP_solver_type.Show(False)
            optionsnotebook.tabSolver.cmbNLP_solver_mode.Show(False)
            optionsnotebook.tabSolver.chkquiet_NLP.Show(False)
            optionsnotebook.tabSolver.chkACE_feasible_point_finder.Show(False)

            optionsnotebook.tabSolver.lblMBH_Pareto_alpha.Show(False)
            optionsnotebook.tabSolver.txtMBH_Pareto_alpha.Show(False)

            optionsnotebook.tabSolver.lblMBH_zero_control_initial_guess.Show(False)
            optionsnotebook.tabSolver.cmbMBH_zero_control_initial_guess.Show(False)

            optionsnotebook.tabSolver.lblquiet_MBH.Show(False)
            optionsnotebook.tabSolver.chkquiet_MBH.Show(False)
            
        elif self.run_inner_loop == 4: #SNOPT
            optionsnotebook.tabSolver.lblMBH_max_not_improve.Show(False)
            optionsnotebook.tabSolver.lblMBH_max_trials.Show(False)
            optionsnotebook.tabSolver.lblMBH_max_run_time.Show(False)
            optionsnotebook.tabSolver.lblMBH_max_step_size.Show(False)
            optionsnotebook.tabSolver.lblMBH_time_hop_probability.Show(False)
            optionsnotebook.tabSolver.lblsnopt_feasibility_tolerance.Show(True)
            optionsnotebook.tabSolver.lblsnopt_major_iterations.Show(True)
            optionsnotebook.tabSolver.lblsnopt_max_run_time.Show(True)
            optionsnotebook.tabSolver.lblderivative_type.Show(True)
            optionsnotebook.tabSolver.lblcheck_derivatives.Show(True)
            optionsnotebook.tabSolver.lblseed_MBH.Show(False)
            optionsnotebook.tabSolver.lblMBH_hop_distribution.Show(False)
            optionsnotebook.tabSolver.lblNLP_solver_type.Show(True)
            optionsnotebook.tabSolver.lblNLP_solver_mode.Show(True)
            optionsnotebook.tabSolver.lblquiet_NLP.Show(True)
            optionsnotebook.tabSolver.lblACE_feasible_point_finder.Show(False)
            
            optionsnotebook.tabSolver.txtMBH_max_not_improve.Show(False)
            optionsnotebook.tabSolver.txtMBH_max_trials.Show(False)
            optionsnotebook.tabSolver.txtMBH_max_run_time.Show(False)
            optionsnotebook.tabSolver.txtMBH_max_step_size.Show(False)
            optionsnotebook.tabSolver.txtMBH_time_hop_probability.Show(False)
            optionsnotebook.tabSolver.txtsnopt_feasibility_tolerance.Show(True)
            optionsnotebook.tabSolver.txtsnopt_major_iterations.Show(True)
            optionsnotebook.tabSolver.txtsnopt_max_run_time.Show(True)
            optionsnotebook.tabSolver.cmbderivative_type.Show(True)
            optionsnotebook.tabSolver.chkcheck_derivatives.Show(True)
            optionsnotebook.tabSolver.chkseed_MBH.Show(False)
            optionsnotebook.tabSolver.cmbMBH_hop_distribution.Show(False)
            optionsnotebook.tabSolver.cmbNLP_solver_type.Show(True)
            optionsnotebook.tabSolver.cmbNLP_solver_mode.Show(True)
            optionsnotebook.tabSolver.chkquiet_NLP.Show(True)
            optionsnotebook.tabSolver.chkACE_feasible_point_finder.Show(False)

            optionsnotebook.tabSolver.lblMBH_Pareto_alpha.Show(False)
            optionsnotebook.tabSolver.txtMBH_Pareto_alpha.Show(False)

            optionsnotebook.tabSolver.lblMBH_zero_control_initial_guess.Show(False)
            optionsnotebook.tabSolver.cmbMBH_zero_control_initial_guess.Show(False)

            optionsnotebook.tabSolver.lblquiet_MBH.Show(False)
            optionsnotebook.tabSolver.chkquiet_MBH.Show(False)
            
        if (self.run_inner_loop == 2 and self.seed_MBH == 1) or (self.run_inner_loop == 4):
            optionsnotebook.tabSolver.lblinterpolate_initial_guess.Show(True)
            optionsnotebook.tabSolver.chkinterpolate_initial_guess.Show(True)
            optionsnotebook.tabSolver.lbltrialX.Show(True)
            optionsnotebook.tabSolver.txttrialX.Show(True)
            optionsnotebook.tabSolver.btntrialX.Show(True)
        else:
            optionsnotebook.tabSolver.lblinterpolate_initial_guess.Show(False)
            optionsnotebook.tabSolver.chkinterpolate_initial_guess.Show(False)
            
            if self.run_inner_loop == 0 or self.run_inner_loop == 1:
                optionsnotebook.tabSolver.lbltrialX.Show(True)
                optionsnotebook.tabSolver.txttrialX.Show(True)
                optionsnotebook.tabSolver.btntrialX.Show(True)
            else:
                optionsnotebook.tabSolver.lbltrialX.Show(False)
                optionsnotebook.tabSolver.txttrialX.Show(False)
                optionsnotebook.tabSolver.btntrialX.Show(False)
            
        if self.interpolate_initial_guess:
            optionsnotebook.tabSolver.lblinitial_guess_num_timesteps.Show(True)
            optionsnotebook.tabSolver.lblinitial_guess_step_size_distribution.Show(True)
            optionsnotebook.tabSolver.txtinitial_guess_num_timesteps.Show(True)
            optionsnotebook.tabSolver.cmbinitial_guess_step_size_distribution.Show(True)

            #if step size distribution is uniform, hide the scale factor box
            if self.initial_guess_step_size_distribution == 0:
                optionsnotebook.tabSolver.lblinitial_guess_step_size_stdv_or_scale.Show(False)
                optionsnotebook.tabSolver.txtinitial_guess_step_size_stdv_or_scale.Show(False)
            else:
                optionsnotebook.tabSolver.lblinitial_guess_step_size_stdv_or_scale.Show(True)
                optionsnotebook.tabSolver.txtinitial_guess_step_size_stdv_or_scale.Show(True)
            
        else:
            optionsnotebook.tabSolver.lblinitial_guess_num_timesteps.Show(False)
            optionsnotebook.tabSolver.lblinitial_guess_step_size_distribution.Show(False)
            optionsnotebook.tabSolver.txtinitial_guess_num_timesteps.Show(False)
            optionsnotebook.tabSolver.cmbinitial_guess_step_size_distribution.Show(False)
            
            #if step size distribution is uniform, hide the scale factor box
            if self.initial_guess_step_size_distribution == 0:
                optionsnotebook.tabSolver.lblinitial_guess_step_size_stdv_or_scale.Show(False)
                optionsnotebook.tabSolver.txtinitial_guess_step_size_stdv_or_scale.Show(False)
            else:
                optionsnotebook.tabSolver.lblinitial_guess_step_size_stdv_or_scale.Show(True)
                optionsnotebook.tabSolver.txtinitial_guess_step_size_stdv_or_scale.Show(True)
        

        #control coordinate system is only shown for low-thrust mission types with MBH or SNOPT
        if (self.mission_type == 2 or self.mission_type == 3) and ( (self.run_inner_loop == 2 and self.seed_MBH == 1) or self.run_inner_loop == 4):
            optionsnotebook.tabSolver.lblinitial_guess_control_coordinate_system.Show(True)
            optionsnotebook.tabSolver.cmbinitial_guess_control_coordinate_system.Show(True)
        else:
            optionsnotebook.tabSolver.lblinitial_guess_control_coordinate_system.Show(False)
            optionsnotebook.tabSolver.cmbinitial_guess_control_coordinate_system.Show(False)


        #outer-loop solver options

                                                                                
        optionsnotebook.tabSolver.cmbrun_outerloop.SetSelection(self.run_outerloop)
        optionsnotebook.tabSolver.txtouterloop_popsize.SetValue(str(self.outerloop_popsize))
        optionsnotebook.tabSolver.txtouterloop_genmax.SetValue(str(self.outerloop_genmax))
        optionsnotebook.tabSolver.txtouterloop_tournamentsize.SetValue(str(self.outerloop_tournamentsize))
        optionsnotebook.tabSolver.txtouterloop_CR.SetValue(str(self.outerloop_CR))
        optionsnotebook.tabSolver.txtouterloop_mu.SetValue(str(self.outerloop_mu))
        optionsnotebook.tabSolver.txtouterloop_stallmax.SetValue(str(self.outerloop_stallmax))
        optionsnotebook.tabSolver.txtouterloop_tolfit.SetValue(str(self.outerloop_tolfit))
        optionsnotebook.tabSolver.txtouterloop_ntrials.SetValue(str(self.outerloop_ntrials))
        optionsnotebook.tabSolver.txtouterloop_elitecount.SetValue(str(self.outerloop_elitecount))
        optionsnotebook.tabSolver.txtouterloop_warmstart.SetValue(str(self.outerloop_warmstart))
        if self.run_outerloop == 1:
            optionsnotebook.tabSolver.txtouterloop_popsize.Show(True)
            optionsnotebook.tabSolver.txtouterloop_genmax.Show(True)
            optionsnotebook.tabSolver.txtouterloop_tournamentsize.Show(True)
            optionsnotebook.tabSolver.txtouterloop_CR.Show(True)
            optionsnotebook.tabSolver.txtouterloop_mu.Show(True)
            optionsnotebook.tabSolver.txtouterloop_stallmax.Show(True)
            optionsnotebook.tabSolver.txtouterloop_tolfit.Show(True)
            optionsnotebook.tabSolver.txtouterloop_ntrials.Show(True)
            optionsnotebook.tabSolver.txtouterloop_elitecount.Show(True)
            optionsnotebook.tabSolver.txtouterloop_warmstart.Show(True)
            optionsnotebook.tabSolver.lblouterloop_popsize.Show(True)
            optionsnotebook.tabSolver.lblouterloop_genmax.Show(True)
            optionsnotebook.tabSolver.lblouterloop_tournamentsize.Show(True)
            optionsnotebook.tabSolver.lblouterloop_CR.Show(True)
            optionsnotebook.tabSolver.lblouterloop_mu.Show(True)
            optionsnotebook.tabSolver.lblouterloop_stallmax.Show(True)
            optionsnotebook.tabSolver.lblouterloop_tolfit.Show(True)
            optionsnotebook.tabSolver.lblouterloop_ntrials.Show(True)
            optionsnotebook.tabSolver.lblouterloop_elitecount.Show(True)
            optionsnotebook.tabSolver.lblouterloop_warmstart.Show(True)
        elif self.run_outerloop == 2:
            optionsnotebook.tabSolver.txtouterloop_popsize.Show(False)
            optionsnotebook.tabSolver.txtouterloop_genmax.Show(False)
            optionsnotebook.tabSolver.txtouterloop_tournamentsize.Show(False)
            optionsnotebook.tabSolver.txtouterloop_CR.Show(False)
            optionsnotebook.tabSolver.txtouterloop_mu.Show(False)
            optionsnotebook.tabSolver.txtouterloop_stallmax.Show(False)
            optionsnotebook.tabSolver.txtouterloop_tolfit.Show(False)
            optionsnotebook.tabSolver.txtouterloop_ntrials.Show(False)
            optionsnotebook.tabSolver.txtouterloop_elitecount.Show(False)
            optionsnotebook.tabSolver.txtouterloop_warmstart.Show(False)
            optionsnotebook.tabSolver.lblouterloop_popsize.Show(False)
            optionsnotebook.tabSolver.lblouterloop_genmax.Show(False)
            optionsnotebook.tabSolver.lblouterloop_tournamentsize.Show(False)
            optionsnotebook.tabSolver.lblouterloop_CR.Show(False)
            optionsnotebook.tabSolver.lblouterloop_mu.Show(False)
            optionsnotebook.tabSolver.lblouterloop_stallmax.Show(False)
            optionsnotebook.tabSolver.lblouterloop_tolfit.Show(False)
            optionsnotebook.tabSolver.lblouterloop_ntrials.Show(False)
            optionsnotebook.tabSolver.lblouterloop_elitecount.Show(False)
            optionsnotebook.tabSolver.lblouterloop_warmstart.Show(False)
        else:
            optionsnotebook.tabSolver.txtouterloop_popsize.Show(False)
            optionsnotebook.tabSolver.txtouterloop_genmax.Show(False)
            optionsnotebook.tabSolver.txtouterloop_tournamentsize.Show(False)
            optionsnotebook.tabSolver.txtouterloop_CR.Show(False)
            optionsnotebook.tabSolver.txtouterloop_mu.Show(False)
            optionsnotebook.tabSolver.txtouterloop_stallmax.Show(False)
            optionsnotebook.tabSolver.txtouterloop_tolfit.Show(False)
            optionsnotebook.tabSolver.txtouterloop_ntrials.Show(False)
            optionsnotebook.tabSolver.txtouterloop_elitecount.Show(False)
            optionsnotebook.tabSolver.txtouterloop_warmstart.Show(False)
            optionsnotebook.tabSolver.lblouterloop_popsize.Show(False)
            optionsnotebook.tabSolver.lblouterloop_genmax.Show(False)
            optionsnotebook.tabSolver.lblouterloop_tournamentsize.Show(False)
            optionsnotebook.tabSolver.lblouterloop_CR.Show(False)
            optionsnotebook.tabSolver.lblouterloop_mu.Show(False)
            optionsnotebook.tabSolver.lblouterloop_stallmax.Show(False)
            optionsnotebook.tabSolver.lblouterloop_tolfit.Show(False)
            optionsnotebook.tabSolver.lblouterloop_ntrials.Show(False)
            optionsnotebook.tabSolver.lblouterloop_elitecount.Show(False)
            optionsnotebook.tabSolver.lblouterloop_warmstart.Show(False)

        #re-size the panel
        optionsnotebook.tabSolver.Layout()
        optionsnotebook.tabSolver.SetupScrolling()


    def update_physics_options_panel(self, optionsnotebook):

        optionsnotebook.tabPhysics.cmbephemeris_source.SetSelection(self.ephemeris_source)
        optionsnotebook.tabPhysics.txtSPICE_leap_seconds_kernel.SetValue(str(self.SPICE_leap_seconds_kernel))
        optionsnotebook.tabPhysics.txtSPICE_reference_frame_kernel.SetValue(str(self.SPICE_reference_frame_kernel))
        optionsnotebook.tabPhysics.txtuniverse_folder.SetValue(self.universe_folder)
        optionsnotebook.tabPhysics.chkperturb_SRP.SetValue(self.perturb_SRP)
        optionsnotebook.tabPhysics.chkperturb_thirdbody.SetValue(self.perturb_thirdbody)
        optionsnotebook.tabPhysics.txtspacecraft_area.SetValue(str(self.spacecraft_area))
        optionsnotebook.tabPhysics.txtcoefficient_of_reflectivity.SetValue(str(self.coefficient_of_reflectivity))
        optionsnotebook.tabPhysics.cmbspiral_model_type.SetSelection(self.spiral_model_type)
        optionsnotebook.tabPhysics.cmblambert_type.SetSelection(self.LambertSolver)

        #if SRP is disabled, make the options associated with it invisible
        if self.perturb_SRP == 1:
            optionsnotebook.tabPhysics.lblspacecraft_area.Show(True)
            optionsnotebook.tabPhysics.lblcoefficient_of_reflectivity.Show(True)
            optionsnotebook.tabPhysics.txtspacecraft_area.Show(True)
            optionsnotebook.tabPhysics.txtcoefficient_of_reflectivity.Show(True)
        else:
            optionsnotebook.tabPhysics.lblspacecraft_area.Show(False)
            optionsnotebook.tabPhysics.lblcoefficient_of_reflectivity.Show(False)
            optionsnotebook.tabPhysics.txtspacecraft_area.Show(False)
            optionsnotebook.tabPhysics.txtcoefficient_of_reflectivity.Show(False)

        #re-size the panel
        optionsnotebook.tabPhysics.Layout()
        optionsnotebook.tabPhysics.SetupScrolling()


    def update_output_options_panel(self, optionsnotebook):
        optionsnotebook.tabOutput.chkcreate_GMAT_script.SetValue(self.create_GMAT_script)
        optionsnotebook.tabOutput.cmboutput_units.SetSelection(self.output_units)
        optionsnotebook.tabOutput.chkgenerate_initial_guess_file.SetValue(self.generate_initial_guess_file)
        optionsnotebook.tabOutput.cmbmission_type_for_initial_guess_file.SetSelection(self.mission_type_for_initial_guess_file)
        optionsnotebook.tabOutput.chkoverride_working_directory.SetValue(self.override_working_directory)
        optionsnotebook.tabOutput.txtforced_working_directory.SetValue(self.forced_working_directory)
        optionsnotebook.tabOutput.chkgenerate_forward_integrated_ephemeris.SetValue(self.generate_forward_integrated_ephemeris)

        if self.generate_initial_guess_file:
            optionsnotebook.tabOutput.lblmission_type_for_initial_guess_file.Show(True)
            optionsnotebook.tabOutput.cmbmission_type_for_initial_guess_file.Show(True)
        else:
            optionsnotebook.tabOutput.lblmission_type_for_initial_guess_file.Show(False)
            optionsnotebook.tabOutput.cmbmission_type_for_initial_guess_file.Show(False)

        if self.override_working_directory:
            optionsnotebook.tabOutput.lblforced_working_directory.Show(True)
            optionsnotebook.tabOutput.txtforced_working_directory.Show(True)
            optionsnotebook.tabOutput.btnforced_working_directory.Show(True)
        else:
            optionsnotebook.tabOutput.lblforced_working_directory.Show(False)
            optionsnotebook.tabOutput.txtforced_working_directory.Show(False)
            optionsnotebook.tabOutput.btnforced_working_directory.Show(False)