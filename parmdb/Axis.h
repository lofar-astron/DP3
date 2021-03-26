// Axis.h: Classes representing a regular or irregular axis.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Classes representing a regular or irregular axis.
/// @author Ger van Diepen (diepen AT astron nl)

#ifndef LOFAR_PARMDB_AXIS_H
#define LOFAR_PARMDB_AXIS_H

#include "../blob/BlobStreamable.h"

#include <memory>
#include <utility>
#include <vector>

namespace dp3 {
namespace parmdb {

/// @ingroup ParmDB
/// @brief Classes representing a regular or irregular axis.

/// @{

/// @brief Abstract base class for a cell centered axis.
class Axis : public blob::BlobStreamable {
 public:
  /// Define a shared_ptr for this class.
  typedef std::shared_ptr<Axis> ShPtr;

  /// The constructor sets the unique id.
  Axis();

  virtual ~Axis() {}

  /// Clone the object.
  virtual Axis::ShPtr clone() const = 0;

  /// Check if two axes are equal. They are if they the same type and values.
  ///@{
  bool operator==(const Axis& that) const;
  bool operator!=(const Axis& that) const { return !operator==(that); }
  ///@}

  /// Get the unique axis id.
  unsigned int getId() const { return itsId; }

  /// Get the center, etc. of the i-th cell.
  ///@{
  double center(size_t n) const { return itsCenter[n]; }
  double width(size_t n) const { return itsWidth[n]; }
  double upper(size_t n) const { return itsUpper[n]; }
  double lower(size_t n) const { return itsLower[n]; }
  ///@}

  /// Get all centers, etc.
  ///@{
  const std::vector<double>& centers() const { return itsCenter; }
  const std::vector<double>& widths() const { return itsWidth; }
  const std::vector<double>& uppers() const { return itsUpper; }
  const std::vector<double>& lowers() const { return itsLower; }
  ///@}

  /// Is the axis regular?
  bool isRegular() const { return itsIsRegular; }

  /// Get nr of cells.
  size_t size() const { return itsCenter.size(); }

  /// Get the start and end value.
  ///@{
  double start() const { return itsLower[0]; }
  double end() const { return itsUpper[itsUpper.size() - 1]; }
  ///@}

  /// Get the total range of the axis.
  std::pair<double, double> range() const {
    return std::make_pair(start(), end());
  }

  /// Get the cellnr of the cell containing value x.
  /// If x is right on the edge, biasRight tells if the left or right cell
  /// is taken.
  /// As a search hint one can tell where to start the search (usually the
  /// result of the previous locate).
  /// If the cell is not found, it returns the cellnr of the next higher cell
  /// and sets the bool to false.
  std::pair<size_t, bool> find(double x, bool biasRight = true,
                               size_t start = 0) const;

  /// Get the cellnr as above, but throw an exception if not found.
  size_t locate(double x, bool biasRight = true, size_t start = 0) const {
    std::pair<size_t, bool> res = find(x, biasRight, start);
    if (!res.second) throwNotFound(x);
    return res.first;
  }

  /// Check if the corresponding intervals in this and that axis are the same.
  bool checkIntervals(const Axis& that) const;

  /// Make a subset of the axis for the given start/end value.
  /// It fills the index of the starting point of the subset on the axis.
  ///@{
  Axis::ShPtr subset(double start, double end, size_t& index) const;
  Axis::ShPtr subset(double start, double end) const;
  ///@}

  /// Make a subset of the axis for the given start/end index.
  Axis::ShPtr subset(size_t start, size_t end) const {
    return doSubset(start, end);
  }

  /// Compress the axis.
  virtual Axis::ShPtr compress(size_t factor) const = 0;

  /// Return the union of this and that axis.
  /// If checks if matching intervals are the same.
  /// It fills s1,e1 with the first and last index of this axis in the new
  /// one. Similarly s2,e2 are filled for that axis.
  /// Note the e1 and e2 are one past the end.
  /// The returned object is a RegularAxis if all intervals in the result
  /// are consecutive and have the same width, otherwise the result
  /// is an OrderedAxis.
  Axis::ShPtr combine(const Axis& that, int& s1, int& e1, int& s2,
                      int& e2) const;

  /// Return the type of \c *this as a string.
  virtual const std::string& classType() const = 0;

  /// Write the contents of \c *this into the blob output stream \a bos.
  virtual void write(blob::BlobOStream& bos) const = 0;

  /// Read the contents from the blob input stream \a bis into \c *this.
  virtual void read(blob::BlobIStream& bis) = 0;

  /// Make an Axis object from the intervals defined by the low/upp values.
  /// If all intervals have the same width, a RegularAxis object is made.
  /// Otherwise an OrderedAxis object.
  /// The intervals must be consecutive.
  /// <br>It checks if both vectors have equal length.
  static Axis::ShPtr makeAxis(const std::vector<double>& low,
                              const std::vector<double>& high);

 private:
  /// Add this and that axis, where this axis must be before that axis.
  Axis::ShPtr add(const Axis& that) const;

  virtual Axis::ShPtr doSubset(size_t start, size_t end) const = 0;

  /// Throw an exception for a non-found cell.
  void throwNotFound(double x) const;

 protected:
  /// Set up the object for a regular axis.
  void setup(double start, double width, unsigned int count);
  /// Set up the object for an irregular axis.
  void setup(const std::vector<double>& v1, const std::vector<double>& v2,
             bool asStartEnd);

  /// Unique seqnr of an Axis object. Used in class AxisMapping.
  static unsigned int theirId;

  unsigned int itsId;
  bool itsIsRegular;
  std::vector<double> itsCenter;
  std::vector<double> itsWidth;
  std::vector<double> itsLower;
  std::vector<double> itsUpper;
};

/// @brief Regularly strided cell centered axis.
class RegularAxis : public Axis {
 public:
  /// Default constructor creates one cell from -1e30 till 1e30.
  RegularAxis();

  /// Construct giving the beginning of the axis and the width of each cell.
  RegularAxis(double begin, double cellWidth, unsigned int count,
              bool asStartEnd = false);

  virtual ~RegularAxis();

  /// Clone the object.
  Axis::ShPtr clone() const override;

  Axis::ShPtr doSubset(size_t start, size_t end) const override;
  Axis::ShPtr compress(size_t factor) const override;

  /// Write the contents of \c *this into the blob output stream \a bos.
  void write(blob::BlobOStream& bos) const override;

  /// Read the contents from the blob input stream \a bis into \c *this.
  void read(blob::BlobIStream& bis) override;

  /// Return the type of \c *this as a string.
  const std::string& classType() const override;

 private:
  double itsStart;
  double itsWidth;
  uint32_t itsCount;
};

/// @brief Ordered irregularly strided cell centered axis.
/// The cells are ordered and disjoint, but gaps may be present.
/// \todo Implementation needs carefull inspection if gaps are to be allowed.
class OrderedAxis : public Axis {
 public:
  /// Default constructor creates one cell from -1e30 till 1e30.
  OrderedAxis();

  /// Specify the intervals defined by v1/v2 as center/width or start/end.
  /// The vectors must have equal sizes. The intervals must be in ascending
  /// order and they have to be disjoint. However, they do not need to be
  /// consecutive.
  OrderedAxis(const std::vector<double>& v1, const std::vector<double>& v2,
              bool asStartEnd = false);

  virtual ~OrderedAxis();

  /// Clone the object.
  Axis::ShPtr clone() const override;

  Axis::ShPtr doSubset(size_t start, size_t end) const override;
  Axis::ShPtr compress(size_t factor) const override;

  /// Write the contents of \c *this into the blob output stream \a bos.
  void write(blob::BlobOStream& bos) const override;

  /// Read the contents from the blob input stream \a bis into \c *this.
  void read(blob::BlobIStream& bis) override;

  /// Return the type of \c *this as a string.
  const std::string& classType() const override;
};

/// @}

}  // namespace parmdb
}  // namespace dp3

#endif
