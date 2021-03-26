// ParameterSetImpl.h: Implements a map of Key-Value pairs.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LOFAR_COMMON_PARAMETERSETIMPL_H
#define LOFAR_COMMON_PARAMETERSETIMPL_H

// Never #include <config.h> or #include <lofar_config.h> in a header file!
#include "ParameterValue.h"
#include "StringTools.h"

#include <map>
#include <memory>
#include <mutex>
#include <set>

namespace dp3 {
namespace common {

/// \ingroup Common
/// \brief Implements a map of Key-Value pairs.

/// @{

const char PC_QUAL_STABLE[] = {"stable"};
const char PC_QUAL_TEST[] = {"test"};
const char PC_QUAL_DEVELOP[] = {"development"};
const char PC_KEY_VERSIONNR[] = {"versionnr"};
const char PC_KEY_QUAL[] = {"qualification"};

typedef common::stringtools::Compare KeyCompare;

/// A key/value map is defined as a map of strings. The third template
/// parameter, \c KeyCompare, is the string comparison functor that will be
/// used to compare keys.
typedef std::map<std::string, ParameterValue, KeyCompare> KVMap;

/// @brief Implements a map of Key-Value pairs.
/// Description of class.
/// The ParameterSetImpl class is a key-value implementation of the type
/// map<string, string, KeyCompare>.
/// This means that values are stored as a string which allows easy merging and
/// splitting of ParameterSetImpls because no conversions have to be done.
/// A couple of getXxx routines are provided to convert the strings to the
/// desired type.
///
/// It keeps track of the keys being asked for. Unused keys might mean
/// that their names were misspelled.
/// The class is fully thread-safe.
class ParameterSetImpl : public KVMap {
 public:
  typedef KVMap::iterator iterator;
  typedef KVMap::const_iterator const_iterator;

  /// \name Construction and Destruction
  /// A ParameterSetImpl can be constructed as empty collection, can be
  /// read from a file or copied from another collection.
  /// @{

  /// Create an empty collection. The argument \a mode determines how
  /// keys should be compared.
  explicit ParameterSetImpl(KeyCompare::Mode mode);

  /// Destroy the contents.
  ~ParameterSetImpl();

  /// Construct a ParameterSet from the contents of \a theFilename. The
  /// argument \a mode determines how keys should be compared.
  explicit ParameterSetImpl(const std::string& theFilename,
                            KeyCompare::Mode mode);
  /// @}

  /// Get the ParameterValue.
  /// @{
  const ParameterValue& get(const std::string& aKey) const;

  /// Return the key comparison mode.
  KeyCompare::Mode keyCompareMode() const { return itsMode; }

  /// \name Merging or appending collections
  /// An existing collection can be extended/merged with another collection.
  /// @{

  /// Adds the Key-Values pair in the given file to the current
  /// collection. Each key will be prefixed with the optional argument
  /// \p thePrefix.
  void adoptFile(const std::string& theFilename,
                 const std::string& thePrefix = "");

  /// Adds the Key-Values pair in the given buffer to the current
  /// collection. Each key will be prefixed with the optional argument
  /// \p thePrefix.
  void adoptBuffer(const std::string& theBuffer,
                   const std::string& thePrefix = "");

  /// Adds the Key-Values pair in the given collection to the current
  /// collection. Each key will be prefixed with the optional argument
  /// \p thePrefix.
  void adoptCollection(const ParameterSetImpl& theCollection,
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

  /// Write the Key-Value pairs from the current ParCollection to the
  /// output stream.
  void writeStream(std::ostream& os) const;
  ///@}

  /// \name Handle subsets
  /// A subset from the current collection can be made based on the prefix
  /// of the keys in the collection.
  /// @{

  /// Creates a subset from the current ParameterSetImpl containing all
  /// parameters starting with the given baseKey. The baseKey is cut off
  /// from the Keynames in the created subset, the optional prefix is put
  /// before all keys in the subset.
  std::shared_ptr<ParameterSetImpl> makeSubset(
      const std::string& baseKey, const std::string& prefix = "") const;

  /// Subtract a subset from the current ParameterSet. Every parameter
  /// whose key starts with the given name will be removed from the
  /// ParameterSet.
  void subtractSubset(const std::string& fullPrefix);
  /// @}

  /// \name Handling single key-value pairs
  /// Single key-value pairs can ofcourse be added, replaced or removed from
  /// a collection.
  /// @{

  /// Add the given pair to the collection. When the \p aKey already exist
  /// in the collection an exception is thrown.
  void add(const std::string& aKey, const ParameterValue& aValue);

  /// Replaces the given pair in the collection. If \p aKey does not exist in
  /// the collection the pair is just added to the collection.
  void replace(const std::string& aKey, const ParameterValue& aValue);

  /// Removes the pair with the given key. Removing a non-existing key is ok.
  void remove(const std::string& aKey);
  /// @}

  /// \name Searching and retrieving
  /// The following functions support searching the collection for existance
  /// of given keys an the retrieval of the corresponding value. In the getXxx
  /// retrieve functions the stored string-value is converted to the wanted
  /// type.
  /// @{

  /// Checks if the given Key is defined in the ParameterSetImpl.
  bool isDefined(const std::string& searchKey) const {
    return (find(searchKey) != end());
  };

  /// Searches for a module whose name end in the given modulename.
  /// e.g: a.b.c.d.param=xxx ; locateModule('d') --> 'a.b.c.'
  std::string locateModule(const std::string& shortName) const;

  /// Searches the module name or module hierarchy and returns its fullposition.
  /// e.g: a.b.c.d.param=xxxx --> fullModuleName(d)-->a.b.c.d
  /// e.g: a.b.c.d.param=xxxx --> fullModuleName(b.c)-->a.b.c
  std::string fullModuleName(const std::string& shortName) const;

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
                                  bool expandable) const;
  std::vector<bool> getBoolVector(const std::string& aKey,
                                  const std::vector<bool>& aValue,
                                  bool expandable) const;
  std::vector<int> getIntVector(const std::string& aKey, bool expandable) const;
  std::vector<int> getIntVector(const std::string& aKey,
                                const std::vector<int>& aValue,
                                bool expandable) const;
  std::vector<unsigned int> getUintVector(const std::string& aKey,
                                          bool expandable) const;
  std::vector<unsigned int> getUintVector(
      const std::string& aKey, const std::vector<unsigned int>& aValue,
      bool expandable) const;
  std::vector<int16_t> getInt16Vector(const std::string& aKey,
                                      bool expandable) const;
  std::vector<int16_t> getInt16Vector(const std::string& aKey,
                                      const std::vector<int16_t>& aValue,
                                      bool expandable) const;
  std::vector<uint16_t> getUint16Vector(const std::string& aKey,
                                        bool expandable) const;
  std::vector<uint16_t> getUint16Vector(const std::string& aKey,
                                        const std::vector<uint16_t>& aValue,
                                        bool expandable) const;
  std::vector<int32_t> getInt32Vector(const std::string& aKey,
                                      bool expandable) const;
  std::vector<int32_t> getInt32Vector(const std::string& aKey,
                                      const std::vector<int32_t>& aValue,
                                      bool expandable) const;
  std::vector<uint32_t> getUint32Vector(const std::string& aKey,
                                        bool expandable) const;
  std::vector<uint32_t> getUint32Vector(const std::string& aKey,
                                        const std::vector<uint32_t>& aValue,
                                        bool expandable) const;
  std::vector<int64_t> getInt64Vector(const std::string& aKey,
                                      bool expandable) const;
  std::vector<int64_t> getInt64Vector(const std::string& aKey,
                                      const std::vector<int64_t>& aValue,
                                      bool expandable) const;
  std::vector<uint64_t> getUint64Vector(const std::string& aKey,
                                        bool expandable) const;
  std::vector<uint64_t> getUint64Vector(const std::string& aKey,
                                        const std::vector<uint64_t>& aValue,
                                        bool expandable) const;
  std::vector<float> getFloatVector(const std::string& aKey,
                                    bool expandable) const;
  std::vector<float> getFloatVector(const std::string& aKey,
                                    const std::vector<float>& aValue,
                                    bool expandable) const;
  std::vector<double> getDoubleVector(const std::string& aKey,
                                      bool expandable) const;
  std::vector<double> getDoubleVector(const std::string& aKey,
                                      const std::vector<double>& aValue,
                                      bool expandable) const;
  std::vector<std::string> getStringVector(const std::string& aKey,
                                           bool expandable) const;
  std::vector<std::string> getStringVector(
      const std::string& aKey, const std::vector<std::string>& aValue,
      bool expandable) const;
  std::vector<time_t> getTimeVector(const std::string& aKey,
                                    bool expandable) const;
  std::vector<time_t> getTimeVector(const std::string& aKey,
                                    const std::vector<time_t>& aValue,
                                    bool expandable) const;
  /// @}

  /// @}

  /// \name Printing
  /// Mostly for debug purposes the collection can be printed.
  /// @{

  // Allow printing the whole parameter collection.
  friend std::ostream& operator<<(std::ostream& os,
                                  const ParameterSetImpl& thePS);
  /// @}

  /// Get all unused parameter names.
  std::vector<std::string> unusedKeys() const;

 private:
  /// Copying is not needed, thus not allowed.
  /// @{
  ParameterSetImpl(const ParameterSetImpl& that);
  ParameterSetImpl& operator=(const ParameterSetImpl& that);
  /// @}

  /// \name Implementation of the 'adopt' methods
  /// The 'adopt' methods are implemented in the readStream method. The 'read'
  /// methods do some preprocessing so the 'adopt' method can use the
  /// \c readStream method.
  /// @{
  void readFile(const std::string& theFile, const std::string& prefix,
                const bool merge);
  void readBuffer(const std::string& theFile, const std::string& prefix,
                  const bool merge);
  void readStream(std::istream& inputStream, const std::string& prefix,
                  const bool merge);
  /// @}

  /// Find the key \p aKey. If \p doThrow \c == \c true (the default) an
  /// exception is thrown when \p aKey is not found. Otherwise, end() is
  /// returned.
  const_iterator findKV(const std::string& aKey, bool doThrow = true) const;

  /// Merge in a key/value. A warning is logged if already existing
  /// and merge=false.
  void addMerge(const std::string& key, const std::string& value, bool merge);

  /// For internal use, add a key without locking.
  void addUnlocked(const std::string& aKey, const ParameterValue& aValue);

  /// For internal use, replace a key without locking.
  void replaceUnlocked(const std::string& aKey, const ParameterValue& aValue);

  /// Key comparison mode.
  const KeyCompare::Mode itsMode;
  /// The set of keys that have been asked.
  mutable std::set<std::string> itsAskedParms;
  /// Mutex to make access to parset thread-safe.
  mutable std::mutex itsMutex;
};

// -------------------- Global functions --------------------
/// Checks if the given string is a valid versionnumber (x.y.z)
bool isValidVersionNr(const std::string& versionNr);

/// Checks if the given string is a valid versionnumber reference. This may be
/// of the form \c x.y.z or the words \c stable, \c test or \c development
/// (defined as \c PC_QUAL_STABLE, \c PC_QUAL_TEST and \c PC_QUAL_DEVELOP).
bool isValidVersionNrRef(const std::string& versionNr);

// Returns the value of the given string or 0 if it is not a valid seqnr
// uint32	seqNr(const string& aString);

// When a hierarchical keyname is passed to \c fullKeyName the methods returns
/// the last part of the keyname. For example:
/// \code
/// moduleName("base.sub.leaf")
/// \endcode
/// returns \c "leaf". When a keyname without dots is passed the whole key
/// is returned.<br>
/// \c keyName is a kind of \c dirname function for keys.
std::string keyName(const std::string& fullKeyName);

// When a hierarchical keyname is passed to \c moduleName the methods returns
/// all but the last part of the keyname. For example:
/// \code
/// moduleName("base.sub.leaf")
/// \endcode
/// returns \c "base.sub". When a keyname without dots is passed and empty
/// string is returned.<br> \c moduleName is a kind of \c basename function for
/// keys.
std::string moduleName(const std::string& fullKeyName);

/// Returns the raw keypart of a parameterline that contains a key-value pair.
/// The returned string is \e not trimmed for whitespace.
std::string keyPart(const std::string& parameterLine);

/// Returns the raw value-part of a parameterline that contains a key-value
/// pair. This means that the string is \e not trimmed for whitespace and that
/// comments at the end of the line are also returned.<br> It simply returns
/// everything behind the first \c = sign.
std::string valuePart(const std::string& parameterLine);

/// Returns the value of the index if the string contains an index otherwise
/// 0 is returned. The \c indexMarker argument must be used to pass the two
/// chars that are used to delimeter the index. The index must be a literal
/// value not an expression. For example: \code
///  indexValue("label{25}andmore", "{}");
/// \endcode
/// returns the value 25. When more indexdelimiters are found in the string the
/// last pair is used.
int32_t indexValue(const std::string& label, const char indexMarker[2]);

/// @} addgroup

inline const ParameterValue& ParameterSetImpl::get(
    const std::string& aKey) const {
  return findKV(aKey)->second;
}

inline bool ParameterSetImpl::getBool(const std::string& aKey) const {
  return findKV(aKey)->second.getBool();
}

inline bool ParameterSetImpl::getBool(const std::string& aKey,
                                      bool aValue) const {
  const_iterator it = findKV(aKey, false);
  if (it == end()) return aValue;
  return it->second.getBool();
}

inline int ParameterSetImpl::getInt(const std::string& aKey) const {
  return findKV(aKey)->second.getInt();
}

inline int ParameterSetImpl::getInt(const std::string& aKey, int aValue) const {
  const_iterator it = findKV(aKey, false);
  if (it == end()) return aValue;
  return it->second.getInt();
}

inline unsigned int ParameterSetImpl::getUint(const std::string& aKey) const {
  return findKV(aKey)->second.getUint();
}

inline unsigned int ParameterSetImpl::getUint(const std::string& aKey,
                                              unsigned int aValue) const {
  const_iterator it = findKV(aKey, false);
  if (it == end()) return aValue;
  return it->second.getUint();
}

inline int16_t ParameterSetImpl::getInt16(const std::string& aKey) const {
  return findKV(aKey)->second.getInt16();
}

inline int16_t ParameterSetImpl::getInt16(const std::string& aKey,
                                          int16_t aValue) const {
  const_iterator it = findKV(aKey, false);
  if (it == end()) return aValue;
  return it->second.getInt16();
}

inline uint16_t ParameterSetImpl::getUint16(const std::string& aKey) const {
  return findKV(aKey)->second.getUint16();
}

inline uint16_t ParameterSetImpl::getUint16(const std::string& aKey,
                                            uint16_t aValue) const {
  const_iterator it = findKV(aKey, false);
  if (it == end()) return aValue;
  return it->second.getUint16();
}

inline int32_t ParameterSetImpl::getInt32(const std::string& aKey) const {
  return findKV(aKey)->second.getInt32();
}

inline int32_t ParameterSetImpl::getInt32(const std::string& aKey,
                                          int32_t aValue) const {
  const_iterator it = findKV(aKey, false);
  if (it == end()) return aValue;
  return it->second.getInt32();
}

inline uint32_t ParameterSetImpl::getUint32(const std::string& aKey) const {
  return findKV(aKey)->second.getUint32();
}

inline uint32_t ParameterSetImpl::getUint32(const std::string& aKey,
                                            uint32_t aValue) const {
  const_iterator it = findKV(aKey, false);
  if (it == end()) return aValue;
  return it->second.getUint32();
}

inline int64_t ParameterSetImpl::getInt64(const std::string& aKey) const {
  return findKV(aKey)->second.getInt64();
}

inline int64_t ParameterSetImpl::getInt64(const std::string& aKey,
                                          int64_t aValue) const {
  const_iterator it = findKV(aKey, false);
  if (it == end()) return aValue;
  return it->second.getInt64();
}

inline uint64_t ParameterSetImpl::getUint64(const std::string& aKey) const {
  return findKV(aKey)->second.getUint64();
}

inline uint64_t ParameterSetImpl::getUint64(const std::string& aKey,
                                            uint64_t aValue) const {
  const_iterator it = findKV(aKey, false);
  if (it == end()) return aValue;
  return it->second.getUint64();
}

inline float ParameterSetImpl::getFloat(const std::string& aKey) const {
  return findKV(aKey)->second.getFloat();
}

inline float ParameterSetImpl::getFloat(const std::string& aKey,
                                        float aValue) const {
  const_iterator it = findKV(aKey, false);
  if (it == end()) return aValue;
  return it->second.getFloat();
}

inline double ParameterSetImpl::getDouble(const std::string& aKey) const {
  return findKV(aKey)->second.getDouble();
}

inline double ParameterSetImpl::getDouble(const std::string& aKey,
                                          double aValue) const {
  const_iterator it = findKV(aKey, false);
  if (it == end()) return aValue;
  return it->second.getDouble();
}

inline std::string ParameterSetImpl::getString(const std::string& aKey) const {
  return findKV(aKey)->second.getString();
}

inline std::string ParameterSetImpl::getString(
    const std::string& aKey, const std::string& aValue) const {
  const_iterator it = findKV(aKey, false);
  if (it == end()) return aValue;
  return it->second.getString();
}

inline time_t ParameterSetImpl::getTime(const std::string& aKey) const {
  return findKV(aKey)->second.getTime();
}

inline time_t ParameterSetImpl::getTime(const std::string& aKey,
                                        const time_t& aValue) const {
  const_iterator it = findKV(aKey, false);
  if (it == end()) return aValue;
  return it->second.getTime();
}

}  // namespace common
}  // namespace dp3

#endif
