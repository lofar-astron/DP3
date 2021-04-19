// ParameterSet.h: Implements a map of Key-Value pairs.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LOFAR_COMMON_PARAMETERSET_H
#define LOFAR_COMMON_PARAMETERSET_H

// Never #include <config.h> or #include <lofar_config.h> in a header file!
#include "ParameterSetImpl.h"
#include "KVpair.h"

namespace dp3 {
namespace common {

/// \ingroup Common
/// \brief Implements a map of Key-Value pairs.

/// @{

/// The ParameterSet class is a key-value implementation of the type
/// map<string, string>.
/// This means that values are stored as a string which allows easy merging and
/// splitting of ParameterSets because no conversions have to be done.
/// A couple of getXxx routines are provided to convert the strings to the
/// desired type.
///
class ParameterSet {
 public:
  typedef ParameterSetImpl::iterator iterator;
  typedef ParameterSetImpl::const_iterator const_iterator;

  /// \name Construction and Destruction
  /// A ParameterSet can be constructed as empty collection, can be
  /// read from a file or copied from another collection.
  /// @{

  /// Create an empty collection. The optional argument \a mode
  /// determines how keys should be compared.
  explicit ParameterSet(KeyCompare::Mode mode = KeyCompare::NORMAL);

  /// Create an empty collection.
  /// Tell if keys have to be compared case-insenstitive.
  explicit ParameterSet(bool caseInsensitive);

  /// Destroy the contents.
  ~ParameterSet();

  /// Construct a ParameterSet from the contents of \a theFilename. The
  /// optional argument \a mode determines how keys should be compared.
  /// @{
  explicit ParameterSet(const std::string& theFilename,
                        KeyCompare::Mode = KeyCompare::NORMAL);
  /// This one is needed to avoid problems with the bool constructor above.
  explicit ParameterSet(const char* theFilename,
                        KeyCompare::Mode = KeyCompare::NORMAL);
  /// @}

  /// Construct a ParameterSet from the contents of \a theFilename.
  /// Tell if keys have to be compared case-insenstitive.
  explicit ParameterSet(const std::string& theFilename, bool caseInsensitive);

  /// Copying is allowed.
  ParameterSet(const ParameterSet& that);

  /// Copying is allowed.
  ParameterSet& operator=(const ParameterSet& that);
  //@}

  /// Is the set empty?
  bool empty() const;

  /// Get the number of parameters.
  int size() const;

  /// Iteration.
  //@{
  iterator begin();
  iterator end();
  const_iterator begin() const;
  const_iterator end() const;
  //@}

  /// Get the ParameterValue.
  /// @{
  const ParameterValue& get(const std::string& aKey) const;
  const ParameterValue& operator[](const std::string& aKey) const;
  /// @}

  /// Key comparison mode.
  KeyCompare::Mode keyCompareMode() const;

  /// Clear the set.
  void clear();

  /// \name Merging or appending collections
  /// An existing collection can be extended/merged with another collection.
  /// @{

  /// Adds the Key-Values pair in the given file to the current
  /// collection. Each key will be prefixed with the optional argument \a
  /// thePrefix.
  /// <br> Note that no dot is added to the prefix.
  void adoptFile(const std::string& theFilename,
                 const std::string& thePrefix = "");

  /// Adds the Key-Values pair in the given buffer to the current
  /// collection. Each key will be prefixed with the optional argument \a
  /// thePrefix.
  /// <br> Note that no dot is added to the prefix.
  void adoptBuffer(const std::string& theBuffer,
                   const std::string& thePrefix = "");

  /// Adds the Key-Values pair in the given collection to the current
  /// collection. Each key will be prefixed with the optional argument \a
  /// thePrefix.
  /// <br> Note that no dot is added to the prefix.
  /// <br>If theCollection is this collection, nothing will be done.
  void adoptCollection(const ParameterSet& theCollection,
                       const std::string& thePrefix = "");

  /// Adds the Key-Values pairs in the argument list.
  /// It ignores arguments not having the Key=Value syntax.
  void adoptArgv(int nr, char const* const argv[]);
  /// @}

  /// \name Saving the collection
  /// The map of key-value pair can be saved in a file or a string.
  /// @{

  /// Writes the Key-Values pair from the current ParCollection to the file.
  void writeFile(const std::string& theFilename, bool append = false) const;

  /// Writes the Key-Values pair from the current ParCollection to the
  /// string buffer.
  void writeBuffer(std::string& theBuffer) const;

  /// Writes the Key-Values pair from the current ParCollection to the
  /// output stream.
  void writeStream(std::ostream& os) const;
  //@}

  /// \name Handle subsets
  /// A subset from the current collection can be made based on the prefix
  /// of the keys in the collection.
  /// @{

  /// Creates a subset from the current ParameterSet containing all the
  /// parameters that start with the given baseKey. The baseKey is cut off
  /// from the Keynames in the created subset, the optional prefix is put
  /// before all keys in the subset.
  /// <br>It is important to note that no dot is added to the prefix, so
  /// it has to be given explicitly. So when giving a prefix like "pre",
  /// a key "pre.key" gets ".key" and "prefix.key" get "fix.key".
  ParameterSet makeSubset(const std::string& baseKey,
                          const std::string& prefix = "") const;

  /// Subtract a subset from the current ParameterSet. Every parameter
  /// whose key starts with the given name will be removed from the
  /// ParameterSet.
  /// <br> Similarly to makeSubset, no dot is added to the prefix.
  void subtractSubset(const std::string& fullPrefix);
  /// @}

  /// \name Handling single key-value pairs
  /// Single key-value pairs can ofcourse be added, replaced or removed from
  /// a collection.
  /// @{

  /// Add the given pair to the collection. When the \p aKey already exist
  /// in the collection an exception is thrown.
  void add(const std::string& aKey, const std::string& aValue);
  void add(const KVpair& aPair);

  /// Replaces the given pair in the collection. If \p aKey does not exist in
  /// the collection the pair is just added to the collection.
  void replace(const std::string& aKey, const std::string& aValue);
  void replace(const KVpair& aPair);

  /// Removes the pair with the given key. Removing a non-existing key is ok.
  void remove(const std::string& aKey);
  /// @}

  /// \name Searching and retrieving
  /// The following functions support searching the collection for existance
  /// of given keys an the retrieval of the corresponding value. In the getXxx
  /// retrieve functions the stored string-value is converted to the wanted
  /// type.
  /// @{

  /// Find a key.
  iterator find(const std::string& searchKey);
  const_iterator find(const std::string& searchKey) const;

  /// Checks if the given Key is defined in the ParameterSet.
  bool isDefined(const std::string& searchKey) const;

  /// Searches for a module whose name end in the given modulename.
  /// e.g: a.b.c.d.param=xxxx --> locateModule(d)-->a.b.c.
  std::string locateModule(const std::string& shortName) const;

  /// Searches the module name or module hierarchy and returns its fullposition.
  /// e.g: a.b.c.d.param=xxxx --> fullModuleName(d)-->a.b.c.d
  /// e.g: a.b.c.d.param=xxxx --> fullModuleName(b.c)-->a.b.c
  std::string fullModuleName(const std::string& shortName) const;

  /// Return the value of the key as a vector of values.
  /// This can only be done if the value is enclosed in square brackets.
  std::vector<ParameterValue> getVector(const std::string& aKey) const;

  /// Return the value of the key as a parameter record.
  /// This can only be done if the value is enclosed in curly braces.
  ParameterRecord getRecord(const std::string& aKey) const;

  /// Return scalar value.
  /// @{
  bool getBool(const std::string& aKey) const;
  bool getBool(const std::string& aKey, bool aValue) const;
  int getInt(const std::string& aKey) const;
  int getInt(const std::string& aKey, int aValue) const;
  unsigned int getUint(const std::string& aKey) const;
  unsigned int getUint(const std::string& aKey, unsigned int aValue) const;
  int16_t getInt16(const std::string& aKey) const;
  int16_t getInt16(const std::string& aKey, int16_t aValue) const;
  uint16_t getUint16(const std::string& aKey) const;
  uint16_t getUint16(const std::string& aKey, uint16_t aValue) const;
  int32_t getInt32(const std::string& aKey) const;
  int32_t getInt32(const std::string& aKey, int32_t aValue) const;
  uint32_t getUint32(const std::string& aKey) const;
  uint32_t getUint32(const std::string& aKey, uint32_t aValue) const;
  int64_t getInt64(const std::string& aKey) const;
  int64_t getInt64(const std::string& aKey, int64_t aValue) const;
  uint64_t getUint64(const std::string& aKey) const;
  uint64_t getUint64(const std::string& aKey, uint64_t aValue) const;
  float getFloat(const std::string& aKey) const;
  float getFloat(const std::string& aKey, float aValue) const;
  double getDouble(const std::string& aKey) const;
  double getDouble(const std::string& aKey, double aValue) const;
  std::string getString(const std::string& aKey) const;
  std::string getString(const std::string& aKey,
                        const std::string& aValue) const;
  /// Returns the value as a time value (seconds since 1970).
  /// @{
  time_t getTime(const std::string& aKey) const;
  time_t getTime(const std::string& aKey, const time_t& aValue) const;
  /// @}
  /// @}

  /// Return vector of values.
  /// @{
  std::vector<bool> getBoolVector(const std::string& aKey,
                                  bool expandable = false) const;
  std::vector<bool> getBoolVector(const std::string& aKey,
                                  const std::vector<bool>& aValue,
                                  bool expandable = false) const;
  std::vector<int> getIntVector(const std::string& aKey,
                                bool expandable = false) const;
  std::vector<int> getIntVector(const std::string& aKey,
                                const std::vector<int>& aValue,
                                bool expandable = false) const;
  std::vector<unsigned int> getUintVector(const std::string& aKey,
                                          bool expandable = false) const;
  std::vector<unsigned int> getUintVector(
      const std::string& aKey, const std::vector<unsigned int>& aValue,
      bool expandable = false) const;
  std::vector<int16_t> getInt16Vector(const std::string& aKey,
                                      bool expandable = false) const;
  std::vector<int16_t> getInt16Vector(const std::string& aKey,
                                      const std::vector<int16_t>& aValue,
                                      bool expandable = false) const;
  std::vector<uint16_t> getUint16Vector(const std::string& aKey,
                                        bool expandable = false) const;
  std::vector<uint16_t> getUint16Vector(const std::string& aKey,
                                        const std::vector<uint16_t>& aValue,
                                        bool expandable = false) const;
  std::vector<int32_t> getInt32Vector(const std::string& aKey,
                                      bool expandable = false) const;
  std::vector<int32_t> getInt32Vector(const std::string& aKey,
                                      const std::vector<int32_t>& aValue,
                                      bool expandable = false) const;
  std::vector<uint32_t> getUint32Vector(const std::string& aKey,
                                        bool expandable = false) const;
  std::vector<uint32_t> getUint32Vector(const std::string& aKey,
                                        const std::vector<uint32_t>& aValue,
                                        bool expandable = false) const;
  std::vector<int64_t> getInt64Vector(const std::string& aKey,
                                      bool expandable = false) const;
  std::vector<int64_t> getInt64Vector(const std::string& aKey,
                                      const std::vector<int64_t>& aValue,
                                      bool expandable = false) const;
  std::vector<uint64_t> getUint64Vector(const std::string& aKey,
                                        bool expandable = false) const;
  std::vector<uint64_t> getUint64Vector(const std::string& aKey,
                                        const std::vector<uint64_t>& aValue,
                                        bool expandable = false) const;
  std::vector<float> getFloatVector(const std::string& aKey,
                                    bool expandable = false) const;
  std::vector<float> getFloatVector(const std::string& aKey,
                                    const std::vector<float>& aValue,
                                    bool expandable = false) const;
  std::vector<double> getDoubleVector(const std::string& aKey,
                                      bool expandable = false) const;
  std::vector<double> getDoubleVector(const std::string& aKey,
                                      const std::vector<double>& aValue,
                                      bool expandable = false) const;
  std::vector<std::string> getStringVector(const std::string& aKey,
                                           bool expandable = false) const;
  std::vector<std::string> getStringVector(
      const std::string& aKey, const std::vector<std::string>& aValue,
      bool expandable = false) const;
  std::vector<time_t> getTimeVector(const std::string& aKey,
                                    bool expandable = false) const;
  std::vector<time_t> getTimeVector(const std::string& aKey,
                                    const std::vector<time_t>& aValue,
                                    bool expandable = false) const;
  /// @}

  /// @}

  /// Get all unused parameter names, thus the names of parameters
  /// that have not been asked for.
  std::vector<std::string> unusedKeys() const;

  /// \name Printing
  /// Mostly for debug purposes the collection can be printed.
  /// @{

  /// Allow printing the whole parameter collection.
  friend std::ostream& operator<<(std::ostream& os, const ParameterSet& thePS);
  /// @}

 private:
  /// Construct from an existing impl object.
  explicit ParameterSet(const std::shared_ptr<ParameterSetImpl>& set) {
    itsSet = set;
  }

  std::shared_ptr<ParameterSetImpl> itsSet;
};

//
// ---------- inline functions ----------
//
inline bool ParameterSet::empty() const { return itsSet->empty(); }
inline int ParameterSet::size() const { return itsSet->size(); }
inline ParameterSet::iterator ParameterSet::begin() { return itsSet->begin(); }
inline ParameterSet::iterator ParameterSet::end() { return itsSet->end(); }
inline ParameterSet::const_iterator ParameterSet::begin() const {
  return itsSet->begin();
}
inline ParameterSet::const_iterator ParameterSet::end() const {
  return itsSet->end();
}

inline const ParameterValue& ParameterSet::get(const std::string& aKey) const {
  return itsSet->get(aKey);
}
inline const ParameterValue& ParameterSet::operator[](
    const std::string& aKey) const {
  return itsSet->get(aKey);
}

inline KeyCompare::Mode ParameterSet::keyCompareMode() const {
  return itsSet->keyCompareMode();
}

inline void ParameterSet::clear() { itsSet->clear(); }

inline void ParameterSet::adoptFile(const std::string& theFilename,
                                    const std::string& thePrefix) {
  itsSet->adoptFile(theFilename, thePrefix);
}

inline void ParameterSet::adoptBuffer(const std::string& theBuffer,
                                      const std::string& thePrefix) {
  itsSet->adoptBuffer(theBuffer, thePrefix);
}

inline void ParameterSet::adoptCollection(const ParameterSet& theCollection,
                                          const std::string& thePrefix) {
  itsSet->adoptCollection(*theCollection.itsSet, thePrefix);
}

inline void ParameterSet::adoptArgv(int nr, const char* const argv[]) {
  itsSet->adoptArgv(nr, argv);
}

inline void ParameterSet::writeFile(const std::string& theFilename,
                                    bool append) const {
  itsSet->writeFile(theFilename, append);
}

inline void ParameterSet::writeBuffer(std::string& theBuffer) const {
  itsSet->writeBuffer(theBuffer);
}

inline void ParameterSet::writeStream(std::ostream& os) const {
  itsSet->writeStream(os);
}

inline ParameterSet ParameterSet::makeSubset(const std::string& baseKey,
                                             const std::string& prefix) const {
  return ParameterSet(itsSet->makeSubset(baseKey, prefix));
}

inline void ParameterSet::subtractSubset(const std::string& fullPrefix) {
  itsSet->subtractSubset(fullPrefix);
}

inline void ParameterSet::add(const std::string& aKey,
                              const std::string& aValue) {
  itsSet->add(aKey, ParameterValue(aValue, false));
}
inline void ParameterSet::add(const KVpair& aPair) {
  add(aPair.first, aPair.second);
}

inline void ParameterSet::replace(const std::string& aKey,
                                  const std::string& aValue) {
  itsSet->replace(aKey, ParameterValue(aValue, false));
}
inline void ParameterSet::replace(const KVpair& aPair) {
  replace(aPair.first, aPair.second);
}

inline void ParameterSet::remove(const std::string& aKey) {
  itsSet->remove(aKey);
}

inline ParameterSet::iterator ParameterSet::find(const std::string& searchKey) {
  return itsSet->find(searchKey);
}

inline ParameterSet::const_iterator ParameterSet::find(
    const std::string& searchKey) const {
  return itsSet->find(searchKey);
}

inline bool ParameterSet::isDefined(const std::string& searchKey) const {
  return itsSet->isDefined(searchKey);
}

inline std::string ParameterSet::locateModule(
    const std::string& shortName) const {
  return (itsSet->locateModule(shortName));
}

inline std::string ParameterSet::fullModuleName(
    const std::string& shortName) const {
  return (itsSet->fullModuleName(shortName));
}

inline std::vector<ParameterValue> ParameterSet::getVector(
    const std::string& aKey) const {
  return get(aKey).getVector();
}

// getBool(key)
inline bool ParameterSet::getBool(const std::string& aKey) const {
  return itsSet->getBool(aKey);
}

// getBool(key, value)
inline bool ParameterSet::getBool(const std::string& aKey, bool aValue) const {
  return itsSet->getBool(aKey, aValue);
}

// getInt(key)
inline int ParameterSet::getInt(const std::string& aKey) const {
  return itsSet->getInt(aKey);
}

// getInt(key, value)
inline int ParameterSet::getInt(const std::string& aKey, int aValue) const {
  return itsSet->getInt(aKey, aValue);
}

// getUint(key)
inline unsigned int ParameterSet::getUint(const std::string& aKey) const {
  return itsSet->getUint(aKey);
}

// getUint(key, value)
inline unsigned int ParameterSet::getUint(const std::string& aKey,
                                          unsigned int aValue) const {
  return itsSet->getUint(aKey, aValue);
}

// getInt16(key)
inline int16_t ParameterSet::getInt16(const std::string& aKey) const {
  return itsSet->getInt16(aKey);
}

// getInt16(key, value)
inline int16_t ParameterSet::getInt16(const std::string& aKey,
                                      int16_t aValue) const {
  return itsSet->getInt16(aKey, aValue);
}

// getUint16(key)
inline uint16_t ParameterSet::getUint16(const std::string& aKey) const {
  return itsSet->getUint16(aKey);
}

// getUint16(key, value)
inline uint16_t ParameterSet::getUint16(const std::string& aKey,
                                        uint16_t aValue) const {
  return itsSet->getUint16(aKey, aValue);
}

// getInt32(key)
inline int32_t ParameterSet::getInt32(const std::string& aKey) const {
  return itsSet->getInt32(aKey);
}

// getInt32(key, value)
inline int32_t ParameterSet::getInt32(const std::string& aKey,
                                      int32_t aValue) const {
  return itsSet->getInt32(aKey, aValue);
}

// getUint32(key)
inline uint32_t ParameterSet::getUint32(const std::string& aKey) const {
  return itsSet->getUint32(aKey);
}

// getUint32(key, value)
inline uint32_t ParameterSet::getUint32(const std::string& aKey,
                                        uint32_t aValue) const {
  return itsSet->getUint32(aKey, aValue);
}

// getInt64(key)
inline int64_t ParameterSet::getInt64(const std::string& aKey) const {
  return itsSet->getInt64(aKey);
}

// getInt64(key, value)
inline int64_t ParameterSet::getInt64(const std::string& aKey,
                                      int64_t aValue) const {
  return itsSet->getInt64(aKey, aValue);
}

// getUint64(key)
inline uint64_t ParameterSet::getUint64(const std::string& aKey) const {
  return itsSet->getUint64(aKey);
}

// getUint64(key, value)
inline uint64_t ParameterSet::getUint64(const std::string& aKey,
                                        uint64_t aValue) const {
  return itsSet->getUint64(aKey, aValue);
}

// getFloat(key)
inline float ParameterSet::getFloat(const std::string& aKey) const {
  return itsSet->getFloat(aKey);
}

// getFloat(key, value)
inline float ParameterSet::getFloat(const std::string& aKey,
                                    float aValue) const {
  return itsSet->getFloat(aKey, aValue);
}

// getDouble(key)
inline double ParameterSet::getDouble(const std::string& aKey) const {
  return itsSet->getDouble(aKey);
}

// getDouble(key, value)
inline double ParameterSet::getDouble(const std::string& aKey,
                                      double aValue) const {
  return itsSet->getDouble(aKey, aValue);
}

// getString(key)
inline std::string ParameterSet::getString(const std::string& aKey) const {
  return itsSet->getString(aKey);
}

// getString(key, value)
inline std::string ParameterSet::getString(const std::string& aKey,
                                           const std::string& aValue) const {
  return itsSet->getString(aKey, aValue);
}

// getTime(key)
inline time_t ParameterSet::getTime(const std::string& aKey) const {
  return itsSet->getTime(aKey);
}

// getTime(key, value)
inline time_t ParameterSet::getTime(const std::string& aKey,
                                    const time_t& aValue) const {
  return itsSet->getTime(aKey, aValue);
}

// getBoolVector(key)
inline std::vector<bool> ParameterSet::getBoolVector(const std::string& aKey,
                                                     bool expandable) const {
  return itsSet->getBoolVector(aKey, expandable);
}

// getBoolVector(key, value)
inline std::vector<bool> ParameterSet::getBoolVector(
    const std::string& aKey, const std::vector<bool>& aValue,
    bool expandable) const {
  return itsSet->getBoolVector(aKey, aValue, expandable);
}

// getIntVector(key)
inline std::vector<int> ParameterSet::getIntVector(const std::string& aKey,
                                                   bool expandable) const {
  return itsSet->getIntVector(aKey, expandable);
}

// getIntVector(key, value)
inline std::vector<int> ParameterSet::getIntVector(
    const std::string& aKey, const std::vector<int>& aValue,
    bool expandable) const {
  return itsSet->getIntVector(aKey, aValue, expandable);
}

// getUintVector(key)
inline std::vector<unsigned int> ParameterSet::getUintVector(
    const std::string& aKey, bool expandable) const {
  return itsSet->getUintVector(aKey, expandable);
}

// getUintVector(key, value)
inline std::vector<unsigned int> ParameterSet::getUintVector(
    const std::string& aKey, const std::vector<unsigned int>& aValue,
    bool expandable) const {
  return itsSet->getUintVector(aKey, aValue, expandable);
}

// getInt16Vector(key)
inline std::vector<int16_t> ParameterSet::getInt16Vector(
    const std::string& aKey, bool expandable) const {
  return itsSet->getInt16Vector(aKey, expandable);
}

// getInt16Vector(key, value)
inline std::vector<int16_t> ParameterSet::getInt16Vector(
    const std::string& aKey, const std::vector<int16_t>& aValue,
    bool expandable) const {
  return itsSet->getInt16Vector(aKey, aValue, expandable);
}

// getUint16Vector(key)
inline std::vector<uint16_t> ParameterSet::getUint16Vector(
    const std::string& aKey, bool expandable) const {
  return itsSet->getUint16Vector(aKey, expandable);
}

// getUint16Vector(key, value)
inline std::vector<uint16_t> ParameterSet::getUint16Vector(
    const std::string& aKey, const std::vector<uint16_t>& aValue,
    bool expandable) const {
  return itsSet->getUint16Vector(aKey, aValue, expandable);
}

// getInt32Vector(key)
inline std::vector<int32_t> ParameterSet::getInt32Vector(
    const std::string& aKey, bool expandable) const {
  return itsSet->getInt32Vector(aKey, expandable);
}

// getInt32Vector(key, value)
inline std::vector<int32_t> ParameterSet::getInt32Vector(
    const std::string& aKey, const std::vector<int32_t>& aValue,
    bool expandable) const {
  return itsSet->getInt32Vector(aKey, aValue, expandable);
}

// getUint32Vector(key)
inline std::vector<uint32_t> ParameterSet::getUint32Vector(
    const std::string& aKey, bool expandable) const {
  return itsSet->getUint32Vector(aKey, expandable);
}

// getUint32Vector(key, value)
inline std::vector<uint32_t> ParameterSet::getUint32Vector(
    const std::string& aKey, const std::vector<uint32_t>& aValue,
    bool expandable) const {
  return itsSet->getUint32Vector(aKey, aValue, expandable);
}

// getInt64Vector(key)
inline std::vector<int64_t> ParameterSet::getInt64Vector(
    const std::string& aKey, bool expandable) const {
  return itsSet->getInt64Vector(aKey, expandable);
}

// getInt64Vector(key, value)
inline std::vector<int64_t> ParameterSet::getInt64Vector(
    const std::string& aKey, const std::vector<int64_t>& aValue,
    bool expandable) const {
  return itsSet->getInt64Vector(aKey, aValue, expandable);
}

// getUint64Vector(key)
inline std::vector<uint64_t> ParameterSet::getUint64Vector(
    const std::string& aKey, bool expandable) const {
  return itsSet->getUint64Vector(aKey, expandable);
}

// getUint64Vector(key, value)
inline std::vector<uint64_t> ParameterSet::getUint64Vector(
    const std::string& aKey, const std::vector<uint64_t>& aValue,
    bool expandable) const {
  return itsSet->getUint64Vector(aKey, aValue, expandable);
}

// getFloatVector(key)
inline std::vector<float> ParameterSet::getFloatVector(const std::string& aKey,
                                                       bool expandable) const {
  return itsSet->getFloatVector(aKey, expandable);
}

// getFloatVector(key, value)
inline std::vector<float> ParameterSet::getFloatVector(
    const std::string& aKey, const std::vector<float>& aValue,
    bool expandable) const {
  return itsSet->getFloatVector(aKey, aValue, expandable);
}

// getDoubleVector(key)
inline std::vector<double> ParameterSet::getDoubleVector(
    const std::string& aKey, bool expandable) const {
  return itsSet->getDoubleVector(aKey, expandable);
}
// getDoubleVector(key, value)
inline std::vector<double> ParameterSet::getDoubleVector(
    const std::string& aKey, const std::vector<double>& aValue,
    bool expandable) const {
  return itsSet->getDoubleVector(aKey, aValue, expandable);
}

// getStringVector(key)
inline std::vector<std::string> ParameterSet::getStringVector(
    const std::string& aKey, bool expandable) const {
  return itsSet->getStringVector(aKey, expandable);
}

// getStringVector(key, value)
inline std::vector<std::string> ParameterSet::getStringVector(
    const std::string& aKey, const std::vector<std::string>& aValue,
    bool expandable) const {
  return itsSet->getStringVector(aKey, aValue, expandable);
}

// getTimeVector(key)
inline std::vector<time_t> ParameterSet::getTimeVector(const std::string& aKey,
                                                       bool expandable) const {
  return itsSet->getTimeVector(aKey, expandable);
}

// getTimeVector(key, value)
inline std::vector<time_t> ParameterSet::getTimeVector(
    const std::string& aKey, const std::vector<time_t>& aValue,
    bool expandable) const {
  return itsSet->getTimeVector(aKey, aValue, expandable);
}

inline std::vector<std::string> ParameterSet::unusedKeys() const {
  return itsSet->unusedKeys();
}

}  // namespace common
}  // namespace dp3

#endif
