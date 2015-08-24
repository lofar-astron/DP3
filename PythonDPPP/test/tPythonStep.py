# tPythonStep.py: Test python DPPP class
# Copyright (C) 2015
# ASTRON (Netherlands Institute for Radio Astronomy)
# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
#
# This file is part of the LOFAR software suite.
# The LOFAR software suite is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The LOFAR software suite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
#
# $Id: __init__.py 23074 2012-12-03 07:51:29Z diepen $


from lofar.pythondppp import DPStep
from lofar.parameterset import parameterset

class tPythonStep(DPStep):
    def __init__(self, parsetDict):
        # The constructor gets the subset of the NDPPP parset containing
        # all keys-value pairs for this step.
        # Note: the superclass constructor MUST be called.
        DPStep.__init__(self, parsetDict)
        parset = parameterset(parsetDict)
        self.itsIncr = parset.getDouble('incr', 1)

    def updateInfo(self, dpinfo):
        # This function must be implemented.
        self.itsInfo = dpinfo
        # Make the arrays that will get the input buffer data from
        # the getData, etc. calls in the process function.
        self.itsData = self.makeArrayDataIn()
        self.itsFlags = self.makeArrayFlagsIn()
        self.itsWeights = self.makeArrayWeightsIn()
        self.itsUVW = self.makeArrayUVWIn()
        # Return the dict with info fields that change in this step.
        return {};

    def process(self, time, exposure):
        # This function must be implemented.
        # First get the data arrays needed by this step.
        self.getData (self.itsData);
        self.getFlags (self.itsFlags);
        self.getWeights (self.itsWeights);
        self.getUVW (self.itsUVW);
        # Process the data.
        print "process tPythonStep", time-4.47203e9, exposure, self.itsData.sum(), self.itsFlags.sum(), self.itsWeights.sum(), self.itsUVW.sum()
        # Execute the next step in the DPPP pipeline. TIME,UVW are changed.
        return self.processNext ({'TIME': time+self.itsIncr, 'UVW': self.itsUVW+self.itsIncr})

    def finish(self):
        # Finish the step as needed.
        # This function does not need to be implemented.
        # Note: finish of the next step is called by the C++ layer.
        print "finish tPythonStep"

    def showCounts(self):
        # Show the counts of this test.
        # This function does not need to be implemented.
        return "   **showcounttest**"

    def addToMS(self, msname):
        # Add some info the the output MeasurementSet.
        # This function does not need to be implemented.
        print "addToMS tPythonStep", msname
