class Body(object):
    #class to contain all data for a body
    name = 'ChinchillaLand'
    shortname = 'Fuzzy'
    number = 0
    SPICE_ID = 10000000000
    minimum_flyby_altitude = -1
    mass = 1000.0
    radius = 1000.0
    ephemeris_epoch = 51544
    alpha0 = 0.0
    alphadot = 0.0
    delta0 = 0.0
    deltadot = 0.0
    W = 0.0
    Wdot = 0.0
    SMA = 1.5e+9
    ECC = 0.0
    INC = 1.0e-6
    RAAN = 1.0e-6
    AOP = 1.0e-6
    MA = 0.0

    def __init__(self, linestring = []):
        if linestring != []:
            self.parse_input_line(linestring)

    def parse_input_line(self, linestring):
        linecell = linestring.split(' ')
        self.name = linecell[0]
        self.shortname = linecell[1]
        self.number = eval(linecell[2])
        self.SPICE_ID = eval(linecell[3])
        self.minimum_flyby_altitude = float(linecell[4])
        self.mass = float(linecell[5])
        self.radius = float(linecell[6])
        self.ephemeris_epoch = float(linecell[7])
        self.alpha0 = float(linecell[8])
        self.alphadot = float(linecell[9])
        self.delta0 = float(linecell[10])
        self.deltadot = float(linecell[11])
        self.W = float(linecell[12])
        self.Wdot = float(linecell[13])
        self.SMA = float(linecell[14])
        self.ECC = float(linecell[15])
        self.INC = float(linecell[16])
        self.RAAN = float(linecell[17])
        self.AOP = float(linecell[18])
        self.MA = float(linecell[19])

    def body_line(self):

        body_line = ""
        body_line += str(self.name) + " "
        body_line += str(self.shortname) + " "
        body_line += str(self.number) + " "
        body_line += str(self.SPICE_ID) + " "
        body_line += str(self.minimum_flyby_altitude) + " "
        body_line += str(self.mass) + " "
        body_line += str(self.radius) + " "
        body_line += str(self.ephemeris_epoch) + " "
        body_line += str(self.alpha0) + " "
        body_line += str(self.alphadot) + " "
        body_line += str(self.delta0) + " "
        body_line += str(self.deltadot) + " "
        body_line += str(self.W) + " "
        body_line += str(self.Wdot) + " "
        body_line += str(self.SMA) + " "
        body_line += str(self.ECC) + " "
        body_line += str(self.INC) + " "
        body_line += str(self.RAAN) + " "
        body_line += str(self.AOP) + " "
        body_line += str(self.MA)

        return body_line