//# Box.cc:
//#
//# Copyright (C) 2008
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
//# $Id: Box.cc 16977 2010-12-20 08:40:36Z diepen $

#include "Box.h"

#include <iostream>

namespace DP3 {
namespace BBS {

  Box::Box (double x1, double x2, double y1, double y2, bool asStartEnd)
  {
    if (!asStartEnd) {
      x1 -= x2 * 0.5;
      x2 += x1;
      y1 -= y2 * 0.5;
      y2 += y1;
    }
    itsStart = Point(x1,y1);
    itsEnd   = Point(x2,y2);
  }

  Box unite (const Box& lhs, const Box& rhs)
  {
    Point start(std::min(lhs.lowerX(), rhs.lowerX()),
                std::min(lhs.lowerY(), rhs.lowerY()));
    Point end(std::max(lhs.upperX(), rhs.upperX()),
              std::max(lhs.upperY(), rhs.upperY()));
    return Box(start, end);
  }


  Box intersect (const Box& lhs, const Box& rhs)
  {
    Point start(std::max(lhs.lowerX(), rhs.lowerX()),
                std::max(lhs.lowerY(), rhs.lowerY()));
    Point end(std::min(lhs.upperX(), rhs.upperX()),
              std::min(lhs.upperY(), rhs.upperY()));
    if(start.first < end.first
       && !casacore::near(start.first, end.first)
       && start.second < end.second
       && !casacore::near(start.second, end.second))
      {
        return Box(start, end);
      }
    return Box();
  }

  Box::Box (const std::vector<double>& values)
  {
    double stx = -1e30;
    double sty = -1e30;
    double enx =  1e30;
    double eny =  1e30;
    int sz = values.size();
    if (sz > 4) sz=4;
    switch (sz) {
    case 4:
      eny = values[3];
    case 3:
      enx = values[2];
    case 2:
      sty = values[1];
    case 1:
      stx = values[0];
      break;
    default:
      break;
    }
    assert (stx <= enx);
    assert (sty <= eny);
    itsStart = Point(stx, sty);
    itsEnd   = Point(enx, eny);
  }

  void Box::print() const
  {
    std::cout << itsStart.first << "\t" << itsStart.second << "\t" 
         << itsEnd.first << "\t" << itsEnd.second << '\n';
  }

} //# namespace BBS
} //# namespace LOFAR
