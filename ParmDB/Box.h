//# Box.h: Class representing a 2-dim box
//#
//# Copyright (C) 2007
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
//# $Id: Box.h 16977 2010-12-20 08:40:36Z diepen $

// @file
// @brief Class representing a 2-dim box
// @author Ger van Diepen (diepen AT astron nl)

#ifndef LOFAR_PARMDB_BOX_H
#define LOFAR_PARMDB_BOX_H

#include <casacore/casa/BasicMath/Math.h>

#include <cassert>
#include <utility>
#include <vector>

namespace DP3 {
namespace BBS {

  // Point: A point in a 2-D space.
  typedef std::pair<double, double>    Point;


  // @ingroup ParmDB
  // @{

  // @brief Class representing a 2-dim box
  class Box;
  Box unite(const Box& lhs, const Box& rhs);
  Box intersect(const Box& lhs, const Box& rhs);


  // A Box is a rectangular region defined by two Point objects defining
  // the bottom-left and top-rigth corner of the box.
  // The bounding box of a Grid is a Box object.
  class Box
  {
  public:
    // Default constructor creates an empty box.
    Box()
      : itsStart(0, 0),
        itsEnd  (0, 0)
    {}

    // Create a box from the bottom-left and top-right corner.
    // both coordinates in start must be <= end.
    Box (const Point& start, const Point& end)
      : itsStart(start),
        itsEnd  (end)
    {
      assert(start.first <= end.first && start.second <= end.second);
    }            

    // Create a box from a vector which must be ordered as stx,sty,endx,endy.
    // Trailing values may be omitted and default to -1e30 or 1e30.
    Box (const std::vector<double>&);

    // Create from start/end or center/width.
    Box (double x1, double x2, double y1, double y2, bool asStartEnd=false);

    // Test if boxes are exactly the same.
    // <group>
    bool operator== (const Box& that) const
      { return itsStart==that.itsStart  &&  itsEnd==that.itsEnd; }
    bool operator!= (const Box& that) const
      { return ! operator== (that); }
    // </group>

    // Get start and end values.
    // <group>
    const Point& lower() const
      { return itsStart; }
    const Point& upper() const
      { return itsEnd; }
    double lowerX() const
      { return itsStart.first; }
    double lowerY() const
      { return itsStart.second; }
    double upperX() const
      { return itsEnd.first; }
    double upperY() const
      { return itsEnd.second; }
    // </group>

    // Get widths.
    // <group>
    double widthX() const
      { return itsEnd.first - itsStart.first; }
    double widthY() const
      { return itsEnd.second - itsStart.second; }
    // </group>

    // Box A only intersects box B if there is at least one point within or
    // on the border of A that falls within B (excluding its border).
    bool intersects (const Box& other) const
    {
      return (other.itsStart.first < itsEnd.first
              && !casacore::near(other.itsStart.first, itsEnd.first)
              && other.itsEnd.first > itsStart.first
              && !casacore::near(other.itsEnd.first, itsStart.first)
              && other.itsStart.second < itsEnd.second
              && !casacore::near(other.itsStart.second, itsEnd.second)
              && other.itsEnd.second > itsStart.second
              && !casacore::near(other.itsEnd.second, itsStart.second));
    }

    // A box A contains a box B if all points within or on the border of B
    // fall within or on the border of A.
    bool contains (const Box& other) const
    {
      return ((other.itsStart.first >= itsStart.first
               || casacore::near(other.itsStart.first, itsStart.first))
              && (other.itsEnd.first <= itsEnd.first
                  || casacore::near(other.itsEnd.first, itsEnd.first))
              && (other.itsStart.second >= itsStart.second
                  || casacore::near(other.itsStart.second, itsStart.second))
              && (other.itsEnd.second <= itsEnd.second
                  || casacore::near(other.itsEnd.second, itsEnd.second)));
    }

    // Check if the box is empty.
    bool empty() const
    {
      return (casacore::near(itsStart.first, itsEnd.first)
              || casacore::near(itsStart.second, itsEnd.second));
    }

    // Return the intersection of this and that box. An empty box is
    // returned if the boxes are disjoint.
    // Note that the operator has a low precedence, so it is advised to
    // enclose the expression in parentheses.
    Box operator& (const Box& that) const
      { return intersect(*this, that); }

    // Return the union of this and that box.
    // The union also contains the points between disjoint boxes.
    // Note that the operator has a low precedence, so it is advised to
    // enclose the expression in parentheses.
    Box operator|(const Box& that) const
      { return unite(*this, that); }

    // Define an ordering functions to be able to sort boxes.
    // The ordering is on startY,startX.
    // <group>
    bool operator< (const Box& that) const
      { return itsStart.second < that.itsStart.second
          ||  (itsStart.second == that.itsStart.second  &&
               itsStart.first < that.itsStart.first); }
    bool operator> (const Box& that) const
      { return itsStart.second > that.itsStart.second
          ||  (itsStart.second == that.itsStart.second  &&
               itsStart.first > that.itsStart.first); }
    // </group>

    // Output the start and end point coordinates of the Box
    void print() const;

  private:
    Point itsStart;
    Point itsEnd;
  };

  // @}

} //# namespace BBS
} //# namespace LOFAR

#endif
