#class file for LRTS sequence archives

class LRTS_sequence(object):
    def __init__(self):
        self.clear()

    def __init__(self, line):
        self.clear()
        self.readline(line)

    def clear(self):
        self.CaseName = []
        self.InitialMass = []
        self.FinalMass = []
        self.StartingEpoch = []
        self.FinalEpoch = []
        self.EpochDifference = []
        self.TotalScore = []
        self.BodyList = []
        self.BodyArrivalEpochs = []

    def readline(self, line):
        print line
        linecell = line.split(',')
        self.CaseName = linecell[0]
        self.InitialMass = float(linecell[1])
        self.FinalMass = float(linecell[2])
        self.StartingEpoch = float(linecell[3])
        self.FinalEpoch = float(linecell[4])
        self.EpochDifference = self.FinalEpoch - self.StartingEpoch
        self.TotalScore = int(linecell[6])

        for b in range(0, self.TotalScore):
            self.BodyList.append(int(linecell[6 + 2*b + 1]))
            self.BodyArrivalEpochs.append(int(linecell[6 + 2*b + 2]))

class LRTS_archive(object):
    def __init__(self):
        self.clear()

    def __init__(self, filename):
        self.clear()
        self.readfile(filename)

    def clear(self):
        self.number_of_sequences = []
        self.sequences = []
        self.filename = []

    def readfile(self, filename):
        linenumber = 0

        inputfile = open(filename, 'r')

        self.filename = filename

        for line in inputfile:
            linenumber += 1
            if linenumber > 2:
                self.sequences.append(LRTS_sequence(line))

        self.number_of_sequences = len(self.sequences)