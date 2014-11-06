import Body
import os
import wx
import copy

class Universe(object):
    #fields
    central_body_name = 'Sun'
    central_body_SPICE_ID = 10
    central_body_radius = 4.379e+6
    mu = 132712440017.99
    LU = 1.49597870691e+8
    reference_angles = [-90.0, 0.0, 66.560709000000003, 0.0, 84.176, 14.1844000]
    r_SOI = 747989353500
    minimum_safe_distance = 2.0e+7

    #list of bodies
    bodies = []
    flyby_menu = []
    flyby_indices = []
    perturbation_menu = []
    perturbation_indices = []
    destination_menu = []
    destination_indices = []
    number_of_bodies = 0

    #success flag
    success = 1

    def __init__(self, input_file_name):
        self.parse_universe_file(input_file_name)
        self.generate_flyby_menu()
        self.generate_perturbation_menu()
        self.generate_destination_menu()

    def parse_universe_file(self, input_file_name):
        #Step 1: open the file
        if os.path.isfile(input_file_name):
            inputfile = open(input_file_name, "r")
            self.success = 1
        else:
            dlg = wx.MessageDialog(self, "Unable to open" + input_file_name, "EMTG Error", wx.OK)
            dlg.ShowModal()
            dlg.Destroy()
            self.success = 0
            return

        body_list_line = 0

        #Step 2: parse the file
        linenumber = 0
        for line in inputfile:
            #strip off the newline character
            line = line.replace("\n","")
            linenumber = linenumber + 1

            if line != "":
                if line[0] != "#":
                    #this is an active line, so it is space delimited
                    #first we need to check if this line is part of the body list
                    if body_list_line > 0:
                        if line == "end_body_list":
                           body_list_line = 0
                        else:
                            self.bodies.append(Body.Body(line))
                            self.number_of_bodies += 1

                    else:
                        linecell = line.split(" ")
                    
                        choice = linecell[0]

                        if choice == "central_body_name":
                            self.central_body_name = linecell[1]

                        elif choice == "central_body_SPICE_ID":
                            self.central_body_SPICE_ID = eval(linecell[1])

                        elif choice == "central_body_radius":
                            self.central_body_radius = eval(linecell[1])

                        elif choice == "mu":
                            self.mu = eval(linecell[1])

                        elif choice == "LU":
                            self.LU = eval(linecell[1])

                        elif choice == "reference_angles":
                            self.reference_angles = [linecell[1], linecell[2], linecell[3], linecell[4], linecell[5], linecell[6]]

                        elif choice == "r_SOI":
                            self.r_SOI = eval(linecell[1])

                        elif choice == "minimum_safe_distance":
                            self.minimum_safe_distance = eval(linecell[1])

                        elif choice == "begin_body_list":
                            self.bodies = []
                            body_list_line = 1

                        #if option is not recognized
                        else:
                            print "Option not recognized: ",str(linecell[0]), " on line " ,str(linenumber)
                            self.success = 0
                            inputfile.close()
                            return
        inputfile.close()
    
    def generate_flyby_menu(self):
        self.flyby_menu = []
        self.flyby_indices = []

        for b in range(0, self.number_of_bodies):
            if self.bodies[b].minimum_flyby_altitude > 0.0:
                self.flyby_indices.append(b)
                self.flyby_menu.append(self.bodies[b].name)

    def generate_perturbation_menu(self):
        self.perturbation_menu = []
        self.flyby_indices = []

        for b in range(0, self.number_of_bodies):
            if self.bodies[b].mass > 1.0e+10:
                self.perturbation_indices.append(b)
                self.perturbation_menu.append(self.bodies[b].name)

    def generate_destination_menu(self):
        self.destination_menu = []
        self.destination_indices = []
        
        #add the free point destination option to the beginning of the list
        self.destination_menu.append('Free point')
        self.destination_indices.append(-1)

        for b in range(0, self.number_of_bodies):
            self.destination_indices.append(b)
            self.destination_menu.append(self.bodies[b].name)

    def write_universe_file(self, output_file_name):
        #write out the universe file
        outputfile = open(output_file_name, "w")

        outputfile.write("#universe file\n")
        outputfile.write("\n")

        outputfile.write("#Central body name\n")
        outputfile.write("central_body_name " + str(self.central_body_name) + "\n")
        outputfile.write("#Central body SPICE ID\n")
        outputfile.write("central_body_SPICE_ID " + str(self.central_body_SPICE_ID) + "\n")
        outputfile.write("#central body radius\n")
        outputfile.write("central_body_radius " + str(self.central_body_radius) + "\n")
        outputfile.write("#gravitational constant of central body, in km^3/s^2\n")
        outputfile.write("mu " + str(self.mu) + "\n")
        outputfile.write("#characteristic length unit, in km\n")
        outputfile.write("LU " + str(self.LU) + "\n")
        outputfile.write("#angles defining the local reference frame relative to ICRF, given in degrees\n")
        outputfile.write("#alpha0, alphadot, delta0, deltadot, W, Wdot\n")
        outputfile.write("reference_angles " + str(self.reference_angles[0]) + " " + str(self.reference_angles[1]) + " " + str(self.reference_angles[2]) + " " + str(self.reference_angles[3]) + " " + str(self.reference_angles[4]) + " " + str(self.reference_angles[5]) + "\n")
        outputfile.write("#radius of the central body's sphere of influence\n")
        outputfile.write("r_SOI " + str(self.r_SOI) + "\n")
        outputfile.write("#minimum safe distance from the central body\n")
        outputfile.write("minimum_safe_distance " + str(self.minimum_safe_distance) + "\n")
        outputfile.write("\n")

        outputfile.write("#menu of bodies\n")
        outputfile.write("#name shortname number SPICE_ID minimum_flyby_altitude mass radius ephemeris_epoch alpha0 alphadot delta0 deltadot W Wdot SMA ECC INC RAAN AOP MA\n")
        outputfile.write("#SMA, radius, and minimum_flyby_altitude in km, angles in degrees, mass in kg, \n")
        outputfile.write("#if minimum_flyby_altitude <= 0, then this object is not placed on the flyby menu\n")
        outputfile.write("#orbit elements are for MJD 51544 / January 1, 2000 00:00 CT, (""Coordinate Time"" as per JPL HORIZONS documentation) in ICRF frame\n")
        outputfile.write("begin_body_list\n")
        for body in self.bodies:
            outputfile.write(body.body_line() + "\n")
        outputfile.write("end_body_list")

        outputfile.close()