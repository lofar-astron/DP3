# __init__.py: Top level .py file for python DPPP step
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

from lofar.pythondppp._pythondppp import _DPStepBase
import numpy as np


class DPStep(_DPStepBase):
    """
    The superclass for a python DPPP step.
    An DPPP step implemented in Python must derive from this class.
    An example can be seen in CEP/DP3/PythonDPPP/test/tPythonStep.py.
    """

    def __init__(self, parset):
        """ This constructor must be called by the subclass """
        _DPStepBase.__init__(self)
        self.itsParset = parset

    def _updateInfo(self, dpinfo):
        """ Private function (called by C++ layer) to update the info.

        The dpinfo argument is a dict containing all fields in the
        C++ DPInfo object. The names of the keys in the dict are the
        DPInfo member names without their 'its' prefix.

        This function extracts the information telling the buffer sizes.
        Thereafter it calls updateInfo to let a subclass extract and
        optionally set the info.
        Finally it extracts the information again to get the output sizes.

        """
        # Get the input sizes.
        self.itsNCorrIn = dpinfo['NCorr']
        self.itsNChanIn = dpinfo['NChan']
        self.itsNBlIn   = len(dpinfo['Ant1'])
        # Default output size is the input size.
        self.itsNCorrOut = self.itsNCorrIn
        self.itsNChanOut = self.itsNChanIn
        self.itsNBlOut   = self.itsNBlIn
        # Let the subclass extract and set info.
        infoOut = self.updateInfo(dpinfo)
        # Extract the output sizes (if defined).
        if 'NCorr' in infoOut:
            self.itsNCorrOut = infoOut['NCorr']
        if 'NChan' in infoOut:
            self.itsNChanOut = infoOut['NChan']
        if 'Ant1' in infoOut:
            self.itsNBlOut   = len(infoOut['Ant1'])
        return infoOut

    # The following functions can be overwritten in a subclass.
    def updateInfo(self, dpinfo):
        """ Extract and optionally set the DPInfo fields.

        This function must be implemented in a subclass.

        The dpinfo argument is a dict containing all fields in the
        C++ DPInfo object. The names of the keys in the dict are the
        DPInfo member names without their 'its' prefix.

        It must return a (possibly empty) dict containing values of the
        same keys. Only the keys that have changed must be part of the
        dict; other keys can be part of it.

        """
        raise ValueError("A class derived from DPStep must implement updateInfo")

    def needVisData(self):
        """ Does the subclass need the visibility data?

        This function only needs to be implemented in a subclass
        if it does not need the visibility data.

        """
        return True

    def needWrite(self):
        """ Does the subclass change data to be written by DPPP?

        This function needs to be implemented in a subclass
        if it changes data (visibility data, flags, weights, and/or UVW)
        that needs to be written (or possibly updated) in the MS.

        """
        return False

    def process(self, time, exposure):
        """ Process the data of the given time slot.

        This function must be implemented in a subclass.
        it will need the functions getData, etc. to get the visibility data,
        flags, etc..
        It should call processNext to process the next DPPP step on the
        output of this step.

        """        
        raise ValueError("A class derived from DPStep must implement process")

    def finish(self):
        """ Finish the processing.

        This function must be implemented in a subclass if data were
        buffered and the last buffers have not been (fully) processed yet.
        In that case finish also has to call processNext.

        """        
        pass
    
    def show(self):
        """ Show the parset parameters.
 
        It should return a string (with newlines) that will be
        printed by the C++ layer.
        An empty string is not printed.

        This function does not need to be implemented.
        The default implementation shows all parset keys.
        """
        s = ''
        for k,v in self.itsParset.iteritems():
            if k not in ['type', 'python.class', 'python.module']:
                s += '  %-15s %s\n' % (k+':', v)
        return s

    def showCounts(self):
        """ Show possible counts (e.g., nr of flags set).

        It should return a string (with newlines) that will be
        printed by the C++ layer.
        An empty string is not printed.

        This function does not need to be implemented.
        """
        return ''

    def showTimings(self, elapsedDuration):
        """ Show possible timings.

        It should return a string (with newlines) that will be
        printed by the C++ layer.
        An empty string is not printed.

        This function does not need to be implemented.
        """
        return ''

    def addToMS(self, msname):
        """ Add information to the output MeasurementSet.

        This function will only be needed in very special cases
        where dedicated info needs to be added to the MS.
        
        """
        pass


    # The following functions are to be called from Python.
    def processNext(self, arraysDict):
        """ Let the next step in the DPPP pipeline execute its process.

        If data (time, exposure, visibility, flags, weights, and/or UVW)
        have been changed by the process function, they should be passed
        using the arraysDict argument.
        This is a dict that can contain 6 fields: TIME, EXPOSURE, DATA,
        FLAGS, WEIGHTS, and UVW.

        If process does not buffer data, the dict only needs to contain
        the changed data items. But if buffered and pocessNext is called
        later, all data items have to be part of the dict.
        """
        return self._processNext(arraysDict)

    def getData(self, nparray):
        """ Get the visibility data into the given numpy array.

        The array must be a contiguous array of the complex64 data type and
        the correct shape.
        The array can be created with the makeDataArrayIn function.

        """
        if not nparray.flags.c_contiguous  or  nparray.size == 0:
            raise ValueError("getData argument 'nparray' has to be a contiguous numpy a\
rray")
        return self._getData (nparray)

    def getFlags(self, nparray):
        """ Get the flags into the given numpy array.

        The array must be a contiguous array of the boolean data type and
        the correct shape.
        The array can be created with the makeFlagsArrayIn function.

        """
        if not nparray.flags.c_contiguous  or  nparray.size == 0:
            raise ValueError("getFlags argument 'nparray' has to be a contiguous numpy a\
rray")
        return self._getFlags (nparray)

    def getWeights(self, nparray):
        """ Get the weights into the given numpy array.

        The array must be a contiguous array of the float32 data type and
        the correct shape.
        The array can be created with the makeWeightsArrayIn function.

        """
        if not nparray.flags.c_contiguous  or  nparray.size == 0:
            raise ValueError("getWeights argument 'nparray' has to be a contiguous numpy a\
rray")
        return self._getWeights (nparray)

    def getUVW(self, nparray):
        """ Get the UVW coordinates into the given numpy array.

        The array must be a contiguous array of the float64 data type and
        the correct shape.
        The array can be created with the makeUVWArrayIn function.

        """
        if not nparray.flags.c_contiguous  or  nparray.size == 0:
            raise ValueError("getUVW argument 'nparray' has to be a contiguous numpy a\
rray")
        return self._getUVW (nparray)

    def getModelData(self, nparray):
        """ Get the model data into the given numpy array.

        The array must be a contiguous array of the complex64 data type and
        the correct shape.
        The array can be created with the makeDataArrayIn function.

        """
        if not nparray.flags.c_contiguous  or  nparray.size == 0:
            raise ValueError("getModelData argument 'nparray' has to be a contiguous numpy a\
rray")
        return self._getModelData (nparray)

    def makeArrayDataIn(self):
        """ Make a numpy array for the visibility input data. """
        return np.empty([self.itsNBlIn, self.itsNChanIn, self.itsNCorrIn], dtype='complex64')

    def makeArrayFlagsIn(self):
        """ Make a numpy array for the input flags. """
        return np.empty([self.itsNBlIn, self.itsNChanIn, self.itsNCorrIn], dtype='bool')

    def makeArrayWeightsIn(self):
        """ Make a numpy array for the input weights. """
        return np.empty([self.itsNBlIn, self.itsNChanIn, self.itsNCorrIn], dtype='float32')

    def makeArrayUVWIn(self):
        """ Make a numpy array for the input UVW coordinates. """
        return np.empty([self.itsNBlIn, 3], dtype='float64')

