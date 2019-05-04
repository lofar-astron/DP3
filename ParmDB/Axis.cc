//# Axis.cc: Classes representing a regular or irregular axis.
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
//# $Id: Axis.cc 20771 2012-04-19 12:04:48Z diepen $

#include "Axis.h"

#include "../Blob/BlobIStream.h"
#include "../Blob/BlobOStream.h"
#include "../Blob/BlobSTL.h"

#include "../Common/StreamUtil.h" 

#include <casacore/casa/BasicMath/Math.h>

namespace DP3 {
namespace BBS {

  using DP3::operator<<;

  // Initialize static.
  unsigned int Axis::theirId = 0;


  // Register with the BlobStreamableFactory. Use an anonymous namespace. This
  // ensures that the 'dummy' variables get their own private storage area and
  // are only visible in this compilation unit.
  namespace
  {
    bool dummy1 = BlobStreamableFactory::instance().registerClass<RegularAxis>("RegularAxis");
    bool dummy2 = BlobStreamableFactory::instance().registerClass<OrderedAxis>("OrderedAxis");
  }


  Axis::Axis()
  {
    itsId = theirId++;
  }

  void Axis::setup (double start, double width, unsigned int count)
  {
    itsIsRegular = true;
    assert(width > 0  &&  count > 0);
    itsCenter.resize (count);
    itsWidth.resize  (count);
    itsUpper.resize  (count);
    itsLower.resize  (count);
    for (unsigned int i=0; i<count; ++i) {
      itsWidth[i]  = width;
      itsCenter[i] = start + width*0.5;
      itsLower[i]  =  start;
      start += width;
      itsUpper[i]  = start;
    }
    // Note: small rounding errors will be made and the last Upper value
    // will not be exactly the same as the given v2.
    // An alternative would be to smear out the rest as shown below,
    // but that would mean the width of each cell can be slightly different.
//     double end = count*width;
//     unsigned int   nr  = count;
//     for (unsigned int i=0; i<count; ++i) {
//       itsLower[i]  = start;
//       itsWidth[i]  = (end-start)/nr;
//       start += itsWidth[i];
//       --nr;
//       itsUpper[i]  = start;
//       itsCenter[i] = (itsUpper[i] + itsLower[i]) * 0.5;
//     }
  }

  void Axis::setup (const std::vector<double>& v1, const std::vector<double>& v2,
                    bool asStartEnd)
  {
    itsIsRegular = false;
    assert(v1.size() == v2.size());
    assert(v1.size() > 0);
    unsigned int nr = v1.size();
    if (!asStartEnd) {
      itsCenter = v1;
      itsWidth  = v2;
      itsLower.resize (nr);
      itsUpper.resize (nr);
      for (unsigned int i=0; i<nr; ++i) {
        itsLower[i] = itsCenter[i] - itsWidth[i] * 0.5;
        itsUpper[i] = itsLower[i] + itsWidth[i];
      }
    } else {
      itsLower = v1;
      itsUpper = v2;
      itsCenter.resize (nr);
      itsWidth.resize  (nr);
      for (unsigned int i=0; i<nr; ++i) {
        itsCenter[i] = (v1[i] + v2[i]) * 0.5;
        itsWidth[i]  = v2[i] - v1[i];
      }
    }
    for (unsigned int i=0; i<nr; ++i) {
      assert(itsWidth[i] > 0);
      if (i > 0) {
        assert(itsUpper[i-1] <= itsLower[i]  ||
               casacore::near(itsUpper[i-1],itsLower[i]));
      }
    }
  }

  bool Axis::operator== (const Axis& that) const
  {
    if (isRegular()  &&  that.isRegular()) {
      return start()==that.start() && end()==that.end() && size()==that.size();
    }
    return itsCenter==that.itsCenter && itsWidth==that.itsWidth;
  }

  bool Axis::checkIntervals (const Axis& that) const
  {
    size_t index;
    Axis::ShPtr ax1 (subset(that.start(), that.end(), index));
    Axis::ShPtr ax2 (that.subset(start(), end(), index));
    unsigned int nr = ax1->size();
    if (ax2->size() != nr) {
      return false;
    }
    for (unsigned int i=0; i<nr; ++i) {
      double low1 = ax1->lower(i);
      double low2 = ax2->lower(i);
      if (!casacore::near(low1, low2)) return false;
      double upp1 = ax1->upper(i);
      double upp2 = ax2->upper(i);
      if (!casacore::near(upp1, upp2)) return false;
    }
    return true;
  }

  pair<size_t,bool> Axis::find (double x, bool biasRight, size_t start) const
  {
    // The locate function searches in a linear way because usually
    // domains are looked up in an ordered way. The start argument can
    // contain the index of the previous interval, so most of the time the
    // current or next interval is the desired one.
    // It is important to note that due to rounding errors the start of
    // an interval could be slightly before the end of the previous interval.
    // The near function needs to be used to check for these cases.
    size_t nr = itsCenter.size();
    // Start searching at the given start position (if possible).
    if (start >= nr) {
      // At the end, so restart at the beginning.
      start = 0;
    } else if (casacore::near(x, itsLower[start])) {
      // At the left edge, so use it if a bias to the right.
      if (biasRight) {
        return pair<size_t,bool> (start, true);
      }
      // Otherwise start at beginning (note that the start of this interval
      // does not need to be the end of the previous one).
      start = 0;
    } else if (x < itsLower[start]) {
      // Before the interval, so start at the beginning.
      start = 0;
    }
    // Search until found.
    bool fnd = true;
    while (true) {
      // Not found if at the end.
      if (start >= nr) {
        fnd = false;
        break;
      }
      double s = itsLower[start];
      double e = itsUpper[start];
      if (casacore::near(x,s)) {
        // On the left edge; take it if biased to the right.
        // Otherwise the value is before the interval (thus nothing found).
        if (!biasRight) {
          fnd = false;
        }
        break;
       } else if (casacore::near(x,e)) {
        // On the right edge; take it if biased to the left.
        // Otherwise continue searching.
        if (!biasRight) {
          break;
        }
      } else if (x < s) {
        // The value is before the interval (thus nothing found).
        fnd = false;
        break;
      } else if (x < e) {
        // Inside the interval, so take it.
        break;
      }
      // Past the current interval, so continue searching if possible.
      ++start;
    }
    return pair<size_t,bool> (start,fnd);
  }

  void Axis::throwNotFound (double x) const
  {
    throw std::runtime_error("Axis::locate: cell " + std::to_string(x) + " not found");
  }

  Axis::ShPtr Axis::combine (const Axis& that,
                             int& s1, int& e1, int& s2, int& e2) const
  {
    pair<double,double> range1 = range();
    pair<double,double> range2 = that.range();
    if (range1.second <= range2.first  ||
        casacore::near(range1.second, range2.first)) {
      // this is fully left of that.
      Axis::ShPtr newAxis (add (that));
      s1 = 0;
      e1 = size();
      e2 = newAxis->size();
      s2 = e2 - that.size();
      return newAxis;
    }
    if (range2.second <= range1.first  ||
        casacore::near(range2.second, range1.first)) {
      // that is fully left of this.
      Axis::ShPtr newAxis (that.add (*this));
      e1 = newAxis->size();
      s1 = e1 - size();
      s2 = 0;
      e2 = that.size();
      return newAxis;
    }
    // Full or partial overlap.
    int nr1 = size();
    int nr2 = that.size();
    double sc1 = center(0);
    double ec1 = center(nr1-1);
    double sc2 = that.center(0);
    double ec2 = that.center(nr2-1);
    // Find out where the range starts are.
    if (range1.first < range2.first) {
      s1 = 0;
      s2 = find(sc2).first;
    } else {
      s1 = that.find(sc1).first;
      s2 = 0;
    }
    // Find out how many intervals are part of the overlap.
    int nr;
    if (range1.second > range2.second) {
      nr = find(ec2).first - s2;
    } else {
      nr = that.find(ec1).first - s1;
    }
    ++nr;
    // Check if the overlapping parts match.
    for (int i=0; i<nr; ++i) {
      double low1 = lower(s2+i);
      double upp1 = upper(s2+i);
      double low2 = that.lower(s1+i);
      double upp2 = that.upper(s1+i);
      if(!((low1==low2 || casacore::near(low1,low2))  &&
                 (upp1==upp2 || casacore::near(upp1,upp2))) )
				throw std::runtime_error(
                 "Axis::combine: interval [" + std::to_string(low1) + ',' + std::to_string(upp1)
                 + "] mismatches [" + std::to_string(low2) + ',' + std::to_string(upp2) + ']');
    }
    // If this fully covers that, return this axis.
    if (s1 == 0  &&  nr == nr2) {
      e1 = nr1;
      e2 = s2+nr;
      return Axis::ShPtr (this->clone());
    }
    // If that fully covers this, return that axis and set which parts
    // of that are not contained in this.
    if (s2 == 0  &&  nr == nr1) {
      e1 = s1+nr;
      e2 = nr2;
      return Axis::ShPtr (that.clone());
    }
    // Partial overlap, so make a new axis.
    std::vector<double> low, upp;
    if (s1 == 0) {
      // this starts before that.
      e1 = nr1;
      e2 = s2+nr2;
      low.reserve (e2);
      upp.reserve (e2);
      for (int i=0; i<s2; ++i) {
        low.push_back (lower(i));
        upp.push_back (upper(i));
      }
      for (int i=0; i<nr2; ++i) {
        low.push_back (that.lower(i));
        upp.push_back (that.upper(i));
      }
    } else {
      // that starts before this.
      e1 = s1+nr1;
      e2 = nr2;
      low.reserve (e1);
      upp.reserve (e1);
      for (int i=0; i<s1; ++i) {
        low.push_back (that.lower(i));
        upp.push_back (that.upper(i));
      }
      for (int i=0; i<nr1; ++i) {
        low.push_back (lower(i));
        upp.push_back (upper(i));
      }
    }
    return makeAxis (low, upp);
  }

  Axis::ShPtr Axis::add (const Axis& that) const
  {
    // That axis will be appended to this one.
    pair<double,double> range1 = range();
    pair<double,double> range2 = that.range();
    int nr1 = size();
    int nr2 = that.size();
    std::vector<double> low, upp;
    low.reserve (nr1+nr2+1);
    upp.reserve (nr1+nr2+1);
    // Copy the first axis bounds.
    for (int i=0; i<nr1; ++i) {
      low.push_back (lower(i));
      upp.push_back (upper(i));
    }
    // See if there is a hole between the axes.
    // If so, create an extra interval for the hole.
    if (! casacore::near (range1.second, range2.first)) {
      // Check this is before that.
      assert(range1.second < range2.first);
      low.push_back (upp[nr1-1]);
      upp.push_back (that.lower(0));
    }
    // Copy the second axis bounds.
    for (int i=0; i<nr2; ++i) {
      low.push_back (that.lower(i));
      upp.push_back (that.upper(i));
    }
    return makeAxis (low, upp);
  }

  Axis::ShPtr Axis::makeAxis (const std::vector<double>& low,
                              const std::vector<double>& upp)
  {
    assert(low.size() == upp.size()  &&  low.size() > 0);
    // Check if the width is constant, thus if the result is a regular axis.
    double width = upp[0] - low[0];
    for (unsigned int i=1; i<low.size(); ++i) {
      if (!casacore::near (width, upp[i]-low[i])) {
        return Axis::ShPtr (new OrderedAxis (low, upp, true));
      }
    }
    return Axis::ShPtr (new RegularAxis (low[0], width, low.size()));
  }

  Axis::ShPtr Axis::subset (double start, double end, size_t& index) const
  {
    pair<double,double> rng = range();
    int sinx = 0;
    int einx = size() - 1;
    if (start > rng.first) {
      sinx = find(start, true).first;
    }
    if (end < rng.second) {
      pair<size_t,bool> res = find(end, false);
      einx = res.first;
      if (einx == 0  &&  !res.second) {
        // Entire interval before the range, so make sinx>einx to return
        // a default axis.
        sinx = 1;
      }
    }

    index = sinx;
    return subset (size_t(sinx), size_t(einx));
  }

  Axis::ShPtr Axis::subset (double start, double end) const
  {
    size_t index;
    return subset (start, end, index);
  }


  RegularAxis::RegularAxis()
    : itsStart (-1e30),
      itsWidth ( 2e30),
      itsCount (1)
  {
    setup (itsStart, itsWidth, itsCount);
  }     
    
  RegularAxis::RegularAxis (double start, double width, unsigned int count,
                            bool asStartEnd)
    : itsStart (start),
      itsWidth (width),
      itsCount (count)
  {
    if (asStartEnd) {
      itsWidth = (itsWidth - itsStart) / itsCount;
    }
    setup (itsStart, itsWidth, itsCount);
  }

  RegularAxis::~RegularAxis()
  {}

  RegularAxis* RegularAxis::clone() const
  {
    return new RegularAxis(*this);
  }

//   size_t RegularAxis::locate (double x, bool biasRight, size_t) const
//   {
//     // Find the cell that contains x.
//     // A value not in the axis domain gets the first or last cell.
//     double inxd = (x - itsBegin) / itsWidth;
//     int inx = int(inxd);
//     int last = int(itsCount) - 1;
//     if (inx < 0) return 0;
//     if (inx > last) return last;
//     // If near the border of the cell, use left or right cell depending
//     // on the biasRight argument.
//     if (biasRight) {
//       if (inx < last) {
//         if (casacore::near (double(inx+1), inxd)) ++inx;
//       }
//     } else {
//       if (inx > 0) {
//         if (casacore::near(double(inx), inxd)) --inx;
//       }
//     }
//     return inx;
//   }

  Axis::ShPtr RegularAxis::doSubset (size_t start, size_t end) const
  {
    if (end >= size()) {
      end = size() - 1;
    }
    if (start > end) {
      return Axis::ShPtr (new RegularAxis());
    }
    return Axis::ShPtr (new RegularAxis (itsStart + start*itsWidth,
                                         itsWidth, end - start + 1));
  }

  Axis::ShPtr RegularAxis::compress (size_t factor) const
  {
    // Is the resulting axis still regular?
    if (itsCount % factor == 0) {
      return Axis::ShPtr(new RegularAxis(itsStart, itsWidth * factor,
                                         itsCount / factor));
    }
    std::vector<double> centers(itsCount / factor + 1);
    for (size_t i = 0; i < itsCount / factor; ++i) {
      centers[i] = itsStart + (i + 0.5) * factor * itsWidth;
    }
    centers.back() = lower(itsCount - (itsCount % factor)) +
      0.5 * (itsCount % factor) * itsWidth;
    std::vector<double> widths(itsCount / factor + 1, factor * itsWidth);
    widths.back() = (itsCount % factor) * itsWidth;
    return Axis::ShPtr(new OrderedAxis(centers, widths));
  }


  void RegularAxis::write (BlobOStream& bos) const
  {
    bos << itsStart << itsWidth << itsCount;
  }

  void RegularAxis::read (BlobIStream& bis)
  {
    bis >> itsStart >> itsWidth >> itsCount;
    setup (itsStart, itsWidth, itsCount);
  }

  const std::string& RegularAxis::classType() const
  {
    static std::string type("RegularAxis");
    return type;
  }


  OrderedAxis::OrderedAxis()
  {
    setup (-1e30, 2e30, 1);
  }

  OrderedAxis::OrderedAxis (const std::vector<double>& starts,
                            const std::vector<double>& ends,
                            bool asStartEnd)
  {
    setup (starts, ends, asStartEnd);
  }

  OrderedAxis::~OrderedAxis()
  {}
    
  OrderedAxis* OrderedAxis::clone() const
  {
    return new OrderedAxis(*this);
  }

  Axis::ShPtr OrderedAxis::doSubset(size_t start, size_t end) const
  {
    if (end >= size()) {
      end = size() - 1;
    }
    if (start > end) {
      return Axis::ShPtr(new RegularAxis());
    }
    std::vector<double> centers(itsCenter.begin()+start, itsCenter.begin()+end+1);
    std::vector<double> widths (itsWidth.begin()+start,  itsWidth.begin()+end+1);
    return Axis::ShPtr(new OrderedAxis(centers, widths));
  }

  Axis::ShPtr OrderedAxis::compress (size_t factor) const
  {
    assert(factor > 0);
    size_t count = static_cast<size_t>(std::ceil(static_cast<double>(size())
        / factor));
    std::vector<double> centers(count), widths(count);
    for (size_t i=0; i<count; ++i) {
      size_t start = i * factor;
      size_t end   = std::min(start + factor, size()) - 1;
      centers[i] = 0.5 * (upper(end) + lower(start));
      widths[i] = upper(end) - lower(start);
    }
    return Axis::ShPtr(new OrderedAxis(centers, widths));
  }
    
  const std::string& OrderedAxis::classType() const
  {
    static std::string type("OrderedAxis");
    return type;
  }

  void OrderedAxis::write (BlobOStream& bos) const
  {
    bos << itsCenter << itsWidth;
  }

  void OrderedAxis::read (BlobIStream& bis)
  {
    bis >> itsCenter >> itsWidth;
    setup (itsCenter, itsWidth, false);
  }


} // namespace BBS
} // namespace LOFAR
