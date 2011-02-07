# __init__.py: Top level .py file for DPPP flagging results plotting
#
# Copyright (C) 2007
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
# $Id$


# Plot the NDPPP count results by frequency or baseline.
# Frequencies are concatenated; baselines are averaged.

""" Module to plot NDPPP count results. """

import pyrap.tables as pt
import numpy
import pylab

def plotflags (tabnames):
    """Plot NDPPP Count results

    A flagging or count step in NDPPP can save the flagging percentages per
    frequency or station. They are saved in a table (per subband) with the
    extension ''.flagfreq'' or ''.flagstat''.
    The flag percentages of a subband can be plotted by giving the name of
    the table containing the results.

    It is also possible to plot the results of multiple subbands by giving
    a list of table names. Frequency results will be sorted in order of
    frequency, while station results are averaged over the subbands.

    """
    t = pt.table(tabnames)
    if 'Frequency' in t.colnames():
        t1 = t.sort ('Frequency')
        pylab.plot (t1.getcol('Frequency'), t1.getcol('Percentage'))
    elif 'Station' in t.colnames():
        percs = []
        names = []
        for t1 in t.iter ('Station'):
            percs.append (t1.getcol('Percentage').mean())
            names.append (t1.getcell('Name', 0))
        pylab.plot (numpy.array(percs), '+')
    else:
        raise RuntimeError('Table appears not to be a NDPPP Count result; it does not contain a Frequency or Station column')

