//# Axis.h: Classes representing a regular or irregular axis.
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
//# $Id: Axis.h 20771 2012-04-19 12:04:48Z diepen $

// @file
// @brief Classes representing a regular or irregular axis.
// @author Ger van Diepen (diepen AT astron nl)

#ifndef LOFAR_PARMDB_AXIS_H
#define LOFAR_PARMDB_AXIS_H

#include "../Blob/BlobStreamable.h"

#include <memory>

namespace DP3 {
namespace BBS {

  using std::pair;
  using std::make_pair;

  // @ingroup ParmDB
  // @{

  // @brief Abstract base class for a cell centered axis.
  class Axis: public BlobStreamable
  {
  public:
    // Define a shared_ptr for this class.
    typedef std::shared_ptr<Axis> ShPtr;

    // The constructor sets the unique id.
    Axis();

    virtual ~Axis()
    {}

    // Clone the object.
    virtual Axis* clone() const = 0;

    // Check if two axes are equal. They are if they the same type and values.
    // <group>
    bool operator== (const Axis& that) const;
    bool operator!= (const Axis& that) const
      { return ! operator== (that); }
    // </group>

    // Get the unique axis id.
    uint getId() const
      { return itsId; }

    // Get the center, etc. of the i-th cell.
    // <group>
    double center(size_t n) const
      { return itsCenter[n]; }
    double width(size_t n) const
      { return itsWidth[n]; }
    double upper(size_t n) const
      { return itsUpper[n]; }
    double lower(size_t n) const
      { return itsLower[n]; }
    // </group>

    // Get all centers, etc.
    // <group>
    const std::vector<double>& centers() const
      { return itsCenter; }
    const std::vector<double>& widths() const
      { return itsWidth; }
    const std::vector<double>& uppers() const
      { return itsUpper; }
    const std::vector<double>& lowers() const
      { return itsLower; }
    // </group>

    // Is the axis regular?
    bool isRegular() const
      { return itsIsRegular; }

    // Get nr of cells.
    size_t size() const
      { return itsCenter.size(); }

    // Get the start and end value.
    // <group>
    double start() const
      { return itsLower[0]; }
    double end() const
      { return itsUpper[itsUpper.size()-1]; }
    // </group>

    // Get the total range of the axis.
    pair<double,double> range() const
      { return pair<double,double> (start(), end()); }

    // Get the cellnr of the cell containing value x.
    // If x is right on the edge, biasRight tells if the left or right cell
    // is taken.
    // As a search hint one can tell where to start the search (usually the
    // result of the previous locate).
    // If the cell is not found, it returns the cellnr of the next higher cell
    // and sets the bool to false.
    pair<size_t,bool> find (double x, bool biasRight = true,
                            size_t start=0) const;

    // Get the cellnr as above, but throw an exception if not found.
    size_t locate (double x, bool biasRight = true, size_t start=0) const
      { pair<size_t,bool> res = find (x, biasRight, start);
        if (!res.second) throwNotFound (x);
        return res.first;
      }

    // Check if the corresponding intervals in this and that axis are the same.
    bool checkIntervals (const Axis& that) const;

    // Make a subset of the axis for the given start/end value.
    // It fills the index of the starting point of the subset on the axis.
    // <group>
    Axis::ShPtr subset (double start, double end, size_t& index) const;
    Axis::ShPtr subset (double start, double end) const;
    // </group>

    // Make a subset of the axis for the given start/end index.
    Axis::ShPtr subset (size_t start, size_t end) const
      { return doSubset (start, end); }

    // Compress the axis.
    virtual Axis::ShPtr compress(size_t factor) const = 0;

    // Return the union of this and that axis.
    // If checks if matching intervals are the same.
    // It fills s1,e1 with the first and last index of this axis in the new
    // one. Similarly s2,e2 are filled for that axis.
    // Note the e1 and e2 are one past the end.
    // The returned object is a RegularAxis if all intervals in the result
    // are consecutive and have the same width, otherwise the result
    // is an OrderedAxis.
    Axis::ShPtr combine (const Axis& that,
                         int& s1, int& e1, int& s2, int& e2) const;

    // Return the type of \c *this as a string.
    virtual const std::string& classType() const = 0;

    // Write the contents of \c *this into the blob output stream \a bos.
    virtual void write (BlobOStream& bos) const = 0;

    // Read the contents from the blob input stream \a bis into \c *this.
    virtual void read (BlobIStream& bis) = 0;

    // Make an Axis object from the intervals defined by the low/upp values.
    // If all intervals have the same width, a RegularAxis object is made.
    // Otherwise an OrderedAxis object.
    // The intervals must be consecutive.
    // <br>It checks if both vectors have equal length.
    static Axis::ShPtr makeAxis (const std::vector<double>& low,
                                 const std::vector<double>& high);

  private:
    // Add this and that axis, where this axis must be before that axis.
    Axis::ShPtr add (const Axis& that) const;

    virtual Axis::ShPtr doSubset (size_t start, size_t end) const = 0;

    // Throw an exception for a non-found cell.
    void throwNotFound (double x) const;

protected:
    // Set up the object for a regular axis.
    void setup (double start, double width, uint count);
    // Set up the object for an irregular axis.
    void setup (const std::vector<double>& v1, const std::vector<double>& v2,
                bool asStartEnd);

    //# Unique seqnr of an Axis object. Used in class AxisMapping.
    static uint theirId;
    uint        itsId;
    bool        itsIsRegular;
    std::vector<double> itsCenter;
    std::vector<double> itsWidth;
    std::vector<double> itsLower;
    std::vector<double> itsUpper;
  };


  // @brief Regularly strided cell centered axis.
  class RegularAxis: public Axis
  {
  public:
    // Default constructor creates one cell from -1e30 till 1e30.
    RegularAxis();

    // Construct giving the beginning of the axis and the width of each cell.
    RegularAxis(double begin, double cellWidth, uint count,
                bool asStartEnd=false);

    virtual ~RegularAxis();

    // Clone the object.
    virtual RegularAxis* clone() const;

    virtual Axis::ShPtr doSubset (size_t start, size_t end) const;
    virtual Axis::ShPtr compress(size_t factor) const;

    // Write the contents of \c *this into the blob output stream \a bos.
    virtual void write (BlobOStream& bos) const;

    // Read the contents from the blob input stream \a bis into \c *this.
    virtual void read (BlobIStream& bis);

    // Return the type of \c *this as a string.
    virtual const std::string& classType() const;
    
  private:
    double itsStart;
    double itsWidth;
    uint32_t itsCount;
  };


  // @brief Ordered irregularly strided cell centered axis.
  // The cells are ordered and disjoint, but gaps may be present.
  // \todo Implementation needs carefull inspection if gaps are to be allowed.
  class OrderedAxis: public Axis
  {
  public:
    // Default constructor creates one cell from -1e30 till 1e30.
    OrderedAxis();

    // Specify the intervals defined by v1/v2 as center/width or start/end.
    // The vectors must have equal sizes. The intervals must be in ascending
    // order and they have to be disjoint. However, they do not need to be
    // consecutive. 
    OrderedAxis (const std::vector<double>& v1, const std::vector<double>& v2,
                 bool asStartEnd=false);
    
    virtual ~OrderedAxis();

    // Clone the object.
    virtual OrderedAxis* clone() const;

    virtual Axis::ShPtr doSubset (size_t start, size_t end) const;
    virtual Axis::ShPtr compress (size_t factor) const;

    // Write the contents of \c *this into the blob output stream \a bos.
    virtual void write (BlobOStream& bos) const;

    // Read the contents from the blob input stream \a bis into \c *this.
    virtual void read (BlobIStream& bis);

    // Return the type of \c *this as a string.
    virtual const std::string& classType() const;
  };


  // @}

} //# namespace BBS
} //# namespace LOFAR

#endif
