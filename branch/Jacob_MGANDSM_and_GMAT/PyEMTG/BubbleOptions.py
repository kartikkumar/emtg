#class file for BubbleOptions class
#this class contains settings for bubble search
#bubble search is used to find flybys of opportunity from an existing trajectory

class BubbleOptions(object):
    def __init__(self):
        self.initialize()

    def __init__(self, smallbodyfile):
        self.initialize()
        self.smallbodyfile = smallbodyfile

    def initialize(self):
        self.smallbodyfile = './MainBeltBright.SmallBody'
        self.LU = 149597870.691
        self.mu = 1.32712440018e+11
        self.RelativePositionFilterMagnitude = 0.1 * self.LU # km
        self.RelativeVelocityFilterMagnitude = 2.0 # km/s
        self.MaximumMagnitude = 14.0

    def update_mission_panel(self, missionpanel):
        missionpanel.txtBubbleSearchFile.SetValue(self.smallbodyfile)
        missionpanel.txtLU.SetValue(str(self.LU))
        missionpanel.txtmu.SetValue(str(self.mu))
        missionpanel.txtRelativePositionFilterMagnitude.SetValue(str(self.RelativePositionFilterMagnitude))
        missionpanel.txtRelativeVelocityFilterMagnitude.SetValue(str(self.RelativeVelocityFilterMagnitude))
        missionpanel.txtMaximumMagnitude.SetValue(str(self.MaximumMagnitude))