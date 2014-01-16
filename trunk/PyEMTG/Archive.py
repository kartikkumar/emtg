import numpy as np
import pylab
from matplotlib import ticker
import matplotlib.pyplot as plt
import os
import datetime

import astropy.time

#EMTG archive item class
class ArchiveItem(object):
    def __init__(self, namestring):
        self.clear()

        if namestring != []:
            self.setname(namestring)

    def clear(self):
        self.values = []
        self.name = []

    def setname(self, namestring):
        self.name = namestring

    def append(self, value):
        self.values.append(value)

    def pop(self, popvalue):
        self.values.pop(popvalue)

    def remove(self, removevalue):
        self.values.remove(removevalue)

#EMTG archive class
class Archive(object):
    def __init__(self, input_file_name):
        self.clear()
        if input_file_name != []:
            self.parse_archive(input_file_name)

    def clear(self):
        self.ArchiveItems = []

    def parse_archive(self, input_file_name):
        #Step 1: open the file
        self.filename = input_file_name

        if os.path.isfile(self.filename):
            inputfile = open(input_file_name, "r")
            self.success = 1
        else:
            print "Unable to open", input_file_name, "EMTG Error"
            self.success = 0
            return

        #Step 2: scan through the file
        linenumber = 0
        for line in inputfile:
            #strip off the newline character
            line = line.replace("\n","")
            linenumber = linenumber + 1

            #the first line of the archive defines everything else
            if linenumber == 1:
               linecell = line.split(',')

               for entry in linecell:
                   self.ArchiveItems.append(ArchiveItem(entry))

            else:
                linecell = line.split(',')
                entryindex = 0
                for entry in linecell:
                    self.ArchiveItems[entryindex].append(float(entry))
                    entryindex += 1

        #Close the file
        inputfile.close()

    def write_archive(self, output_file_name):
        outputfile = open(output_file_name, "w")

        outputfile.write(self.ArchiveItems[0].name)
        for entry in self.ArchiveItems[1:]:
            outputfile.write(',' + entry.name)
        outputfile.write('\n')

        for linenumber in range(0,len(self.ArchiveItems[0].values)):
            outputfile.write(str(self.ArchiveItems[0].values[linenumber]))
            for item in self.ArchiveItems[1:]:
                outputfile.write(',' + str(item.values[linenumber]))
            outputfile.write('\n')

        outputfile.close()

    def plot_objective_vs_arrival_date(self):
        #first we need to create an array of arrival dates

        ArrivalDates = []
        for linenumber in range(0,len(self.ArchiveItems[0].values)):
            currentArrivalDate = 0.0
            for item in self.ArchiveItems:
                
                if "flight time" in item.name or "launch epoch" in item.name or "stay time" in item.name or "wait time" in item.name:
                    currentArrivalDate += item.values[linenumber]
            ArrivalDates.append(currentArrivalDate)

        date_vector = []
        for date in ArrivalDates:
            datestring = astropy.time.Time(date, format='mjd', scale='tdb', out_subfmt='date').iso
            date_vector.append(datetime.datetime.strptime(datestring,'%Y-%m-%d').date())
        
        self.DataFigure = plt.figure()
        self.DataAxes = self.DataFigure.add_axes([0.1, 0.1, 0.8, 0.8])
        self.DataAxes.scatter(date_vector, -np.array(self.ArchiveItems[-1].values))
        self.DataAxes.set_xlabel('Arrival Epoch')

        def format_date(x, pos=None):
            return pylab.num2date(x).strftime('%m-%d-%Y')

        self.DataAxes.xaxis.set_major_formatter(ticker.FuncFormatter(format_date))
        self.DataFigure.autofmt_xdate()
        plt.show()