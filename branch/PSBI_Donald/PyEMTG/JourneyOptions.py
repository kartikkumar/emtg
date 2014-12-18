#Python class file for EMTG JourneyOptions
class JourneyOptions(object):
    journey_names = 'New_Journey'
    journey_timebounded = 0#0: unbounded, 1: bounded flight time, 2: bounded arrival date
    journey_wait_time_bounds = [0.0]*2#days
    journey_flight_time_bounds = [0.0]*2
    journey_arrival_date_bounds = [0.0]*2
    journey_initial_impulse_bounds = [0.0]*2 #in km/s
    journey_arrival_type = 1 #0: orbit insertion (use chemical Isp), 1: rendezvous (use chemical Isp), 2: flyby with bounded VHP, 3: low-thrust rendezvous (does not work if terminal phase is not low-thrust), 4: match v-infinity, 5: match v-infinity low-thrust, 6: E=0, 7: capture spiral
    journey_departure_type = 0 #0: launch or direct insertion, 1: depart from parking orbit (you can use this one in place of a launch vehicle model, and the departure burn will be done with the EDS motor), 2: 'free' direct departure, i.e. do not burn to get the departure v_infinity, 3: Start from Sphere of Influence (use SOI angles chosen by previous journey's endpoint, i.e. after a spiral-out or fully modeled departure from parking orbit), 3: flyby (only valid for successive journeys), 4: flyby with fixed VHP, 5: escape spiral
    journey_arrival_elements_type = 1 #0: cartesian, 1: COE
    journey_arrival_elements = [0.0]*6 #a(km), e, i, RAAN, AOP, TA
    journey_arrival_elements_bounds = [0.0]*12
    journey_arrival_elements_vary_flag = [0]*6
    journey_departure_elements_type = 1 #0: cartesian, 1: COE
    journey_departure_elements = [0.0]*6 #a(km), e, i, RAAN, AOP, TA
    journey_departure_elements_bounds = [0.0]*12
    journey_departure_elements_vary_flag = [0]*6
    journey_central_body = 'Sun' #spice names
    journey_initial_velocity = [0.0]*3 #in km/s
    journey_final_velocity = [0.0]*3 #in km/s
    journey_starting_mass_increment = 0 #in kg
    journey_variable_mass_increment = 0 #whether or not the optimizer can choose the mass increment (ignored for non-positive mass increments)
    journey_arrival_declination_constraint_flag = 0
    journey_arrival_declination_bounds = [0.0]*2#in degrees
    journey_number_of_perturbation_bodies = 1
    journey_perturbation_bodies = [0]
    journey_escape_spiral_starting_radius = 6678#in km
    journey_capture_spiral_final_radius = 6678#in km
    journey_maximum_DSM_magnitude_constraint_flag = 0
    journey_maximum_DSM_magnitude_constraint = 2.0 #in km/s
    journey_distance_constraint_number_of_bodies = 0
    journey_distance_constraint_bodies = []
    journey_distance_constraint_bounds = []
        
    #sequence information
    number_of_phases = []
    sequence = []
    destination_list = [1, 1]
    phase_type = []

    #outer loop selectable options settings
    outerloop_vary_journey_destination = 0
    outerloop_vary_journey_flyby_sequence = 0
    outerloop_journey_destination_choices = [1]
    outerloop_journey_flyby_sequence_choices = [1]
    outerloop_journey_maximum_number_of_flybys = 8

    #************************************************************************************constructor
    def __init__(self, mission_type):
        self.sequence = []
        self.journey_number_of_perturbation_bodies = 1
        self.journey_perturbation_bodies = [0]

        if mission_type == 2 | mission_type == 3 | mission_type == 4:
            self.journey_arrival_type = 3



