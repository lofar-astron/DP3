//# ProgressMeter.h: Visual indication of a tasks progress.
//# Copyright (C) 1997,2000
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id$

#ifndef LOFAR_COMMON_PROGRESSMETER_H
#define LOFAR_COMMON_PROGRESSMETER_H

//# Includes
#include <Common/lofar_string.h>

namespace LOFAR {

// @ingroup NDPPP

// This class shows the progress of a task.
// It is copied from casacore.
//
// It shows the progress on stdout using a line with dots and percentages.
//
// It is possible to attach the progressmeter to function in, say, a GUI
// that can show the progress in a more visual way.
//
// The progress meter will usually be removed from the screen once the maximum
// value is set, so you should not reuse the ProgressMeter after that has
// happened. It is harmless, but it will not result in any visual feedback for
// the user.
//
// While the "min" is usually less than "max", if in fact it is greater than
// "max" the progress meter will count down correctly.
//
// For example:
// <srcblock>
// void calculate(uint n) {
//   int skip = n / 200;
//   ProgressMeter meter(0, n, "Title", "Subtitle", "", "", true, skip);
//   for (uint i=0; i<n; i++) {
//       ... calculate ...
//       meter.update(i);
//   }
// }
// </srcblock>

class ProgressMeter
{
public:
    // Makes a null progress meter, i.e. no updates to the screen are
    // generated.
    ProgressMeter();

    // Create a progress meter with the given min and max values and labels.
    // if <tt>estimateTime=true</tt>, an estimate of the
    // time remaining will be made for the user. This estimate assumes that
    // the remaining portion will compute at the same rate as the portion
    // completed so far, so the time should not be estimated for processes
    // which are non-linear.
    //
    // Any labels which are set to the empty string will have sensible defaults
    // supplied. In particular, <tt>minlabel</tt> and <tt>maxlabel</tt>
    // will be set to the display the minimum and maximum values.
    //
    // Normally the progress bar will be updated with every call to
    // <tt>update()</tt>. If however you will be sending many events
    // then you might want to update the GUI every <tt>updateEvery</tt>'th
    // event for efficiency. Generally there's no point updating more than
    // a couple of hundred times since the eye can't distinguish differences
    // in the progress bar position at that level. If updateEvery is <=0, it
    // is set to 1 for you.
    ProgressMeter(double min, double max, 
		  const string& title, const string& subtitle,
		  const string& minlabel, const string& maxlabel,
		  bool estimateTime=true, int updateEvery=1);

    // The destruction of the meter will cause an update to be sent with the
    // maximum value. This will usually cause the GUI window to be removed
    // from the screen. Thus the progress meter should generally live as long
    // as the calculation it is tracking.
    ~ProgressMeter();

    void update(double value, bool force=false);

    // Get the min and max values of the progress meter.
    // <group>
    double min() const
      { return min_p; }
    double max() const
      { return max_p; }
    // </group>

    friend class ObjectController;
private:
    int id_p;
    double min_p, max_p;
    int update_every_p, update_count_p;
    
    // These are set by ObjectController for executables that have the tasking
    // system in them, otherwise they are null and this class just does no-ops.
    static int (*creation_function_p)(double, double, 
                                      const string&, const string&,
                                      const string&, const string&,
                                      bool);
    static void (*update_function_p)(int, double);

    // Undefined and inaccessible
    ProgressMeter(const ProgressMeter&);
    ProgressMeter& operator=(const ProgressMeter&);
};


} //# NAMESPACE LOFAR END

#endif
