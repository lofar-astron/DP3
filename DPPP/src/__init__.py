# Plot the NDPPP count results by frequency or baseline.
# Frequencies are concatenated; baselines are averaged.

import pyrap.tables as pt
import numpy
import pylab

def plotflags (tabnames):
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
        raise 'Table appears not to be a NDPPP Count result; it does not contain a Frequency or Station column'

