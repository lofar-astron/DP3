//# ParameterSet.cc: Implements a map of Key-Value pairs.
//#
//# Copyright (C) 2002-2003
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
//# $Id: ParameterSetImpl.cc 23459 2013-01-08 08:49:52Z diepen $

//# Always #include <lofar_config.h> first!
#include "ParameterSetImpl.h"

#include "StreamUtil.h"
#include "StringUtil.h"

#include <boost/algorithm/string/case_conv.hpp>

#include <fstream>
#include <sstream>
#include <iostream>
#include <memory>
#include <string>
#include <stdexcept>
#include <vector>

typedef std::runtime_error APSException;

namespace DP3 {

//-------------------------- creation and destroy ---------------------------
//
// Default constructor
//
ParameterSetImpl::ParameterSetImpl(KeyCompare::Mode	mode)
	: KVMap(mode),
	  itsMode(mode)
{}

//
// Construction by reading a parameter file.
//
ParameterSetImpl::ParameterSetImpl(const std::string&	theFilename,
				   KeyCompare::Mode	mode)
	: KVMap(mode), 
	  itsMode(mode)
{
	readFile(theFilename, "", false);
}

//
//	Destructor
//
ParameterSetImpl::~ParameterSetImpl()
{}

//
// operator<<
//
std::ostream&	operator<< (std::ostream& os, const ParameterSetImpl &thePS)
{
	thePS.writeStream(os);
	return os;
}

//-------------------------- merge and split --------------------------------
//
// makeSubset(baseKey [, prefix])
//
// Creates a Subset from the current ParameterSetImpl containing all the 
// parameters that start with the given baseKey. 
// The baseKey is cut off from the Keynames in the created subset, the 
// optional prefix is put before the keynames.
//
std::shared_ptr<ParameterSetImpl>
ParameterSetImpl::makeSubset(const std::string& baseKey, 
                             const std::string& prefix) const
{
  // Thread-safety.
  std::lock_guard<std::mutex> locker(itsMutex);
  // Convert \a baseKey to lowercase, if we need to do case insensitve compare.
  std::string   base = (itsMode == KeyCompare::NOCASE) ? boost::to_lower_copy(baseKey) : baseKey;
  std::shared_ptr<ParameterSetImpl> subSet (new ParameterSetImpl(itsMode));
  iterator pos  = subSet->begin();

  // Start scanning at the point where \a base might first occur.
  for (const_iterator it = lower_bound(base); it != end(); ++it) {

    bool match = (itsMode == KeyCompare::NOCASE) ?
      boost::to_lower_copy(it->first).compare(0, base.size(), base) == 0 :
      it->first.compare(0, base.size(), base) == 0;

    // We can stop scanning once \a match becomes false, since keys are sorted.
    if (!match) break;

    // cut off baseString and copy to subset
    pos = subSet->insert(pos, make_pair(prefix + it->first.substr(base.size()),
                                        it->second));
    itsAskedParms.insert (it->first);
  }
  
  return (subSet);
}

//
// subtractSubset(fullPrefix)
//
// Removes all the keys the start with the given prefix from the parset.
//
void ParameterSetImpl::subtractSubset(const std::string& fullPrefix) 
{
        // Thread-safety.
        std::lock_guard<std::mutex> locker(itsMutex);

	// Convert \a baseKey to lowercase, if we need to do case insensitve compare.
	std::string	prefix = (itsMode == KeyCompare::NOCASE) ?
	boost::to_lower_copy(fullPrefix) : fullPrefix;
	int		length = prefix.size();

	// Loop over parset and delete the matching keys.
	iterator	endIter  = end();
	iterator	iter     = lower_bound(prefix);
	while (iter != endIter) {
		bool match = (itsMode == KeyCompare::NOCASE) ?
					boost::to_lower_copy(iter->first).compare(0, length, prefix) == 0 :
					iter->first.compare(0, length, prefix) == 0;

		// We can stop scanning once a match becomes false, since keys are sorted.
		if (!match) {
			return;
		}

		// erase the matching element
		iterator	old_iter = iter;
		iter++;
		erase (old_iter);
	}
}

//
// adoptFile
//
// Adds the parameters from the file to the current ParameterSetImpl.
//
void ParameterSetImpl::adoptFile(const std::string&	theFilename,
				 const std::string&	thePrefix)
{
	readFile(theFilename, thePrefix, true);
}

//
// adoptBuffer
//
// Adds the parameters from the string to the current ParameterSetImpl.
//
void ParameterSetImpl::adoptBuffer(const std::string&	theBuffer,
				   const std::string&	thePrefix)
{
	readBuffer(theBuffer, thePrefix, true);
}

//
// adoptCollection
//
// Adds the parameters from the ParColl. to the current ParameterSetImpl.
//
void ParameterSetImpl::adoptCollection(const ParameterSetImpl& theCollection,
				       const std::string&	thePrefix)
{
  // Thread-safety.
  std::lock_guard<std::mutex> locker(itsMutex);
  // Cannot adopt itself.
  if (&theCollection != this) {
    for (const_iterator iter = theCollection.begin();
         iter != theCollection.end(); ++iter) {
      replaceUnlocked(thePrefix+iter->first, iter->second);
    }
  } else if (! thePrefix.empty()) {
    // However, adopt itself if a prefix is given.
    // Iterate on a copy, otherwise an endless loop occurs.
    KVMap tmp(theCollection);
    for (const_iterator iter = tmp.begin();
         iter != tmp.end(); ++iter) {
      replaceUnlocked(thePrefix+iter->first, iter->second);
    }
  }
}

void ParameterSetImpl::adoptArgv (int nr, char const * const argv[])
{
  // Thread-safety.
  std::lock_guard<std::mutex> locker(itsMutex);
  for (int i=0; i<nr; ++i) {
    std::string arg(argv[i]);
    // Only add arguments containing an =-sign.
    std::string::size_type eqs = arg.find('=');
    if (eqs != std::string::npos) {
      replaceUnlocked(arg.substr(0, eqs), ParameterValue(arg.substr(eqs+1)));
    }
  }
}

//
// readFile
// (private)
//
// Disentangles the file and adds the Key-Values pair to the current ParameterSetImpl.
//
void ParameterSetImpl::readFile(const	std::string&	theFilename, 
				const	std::string&	prefix,
				const	bool	merge)
{
	std::ifstream		paramFile;

	// Try to pen the file
	paramFile.open(theFilename.c_str(), std::ifstream::in);
	if (!paramFile) {
		throw APSException( 
		       formatString("Unable to open file %s", theFilename.c_str()));
	}

	if (paramFile.eof()) {
		throw APSException( 
		       formatString("file %s is empty", theFilename.c_str()));
	}

	readStream(paramFile, prefix, merge);

	paramFile.close();
}

//
// readBuffer
// (private)
//
// Disentangles the file and adds the Key-Values pair to the current ParameterSetImpl.
//
void ParameterSetImpl::readBuffer(const	std::string&	theBuffer, 
				  const	std::string&	prefix,
				  const	bool	merge)
{
	std::istringstream		iss(theBuffer, std::istringstream::in);
	readStream(iss, prefix, merge);
}

void ParameterSetImpl::readStream (std::istream& inputStream, 
                                   const std::string& prefix,
                                   bool	merge)
{
  // Thread-safety.
  std::lock_guard<std::mutex> locker(itsMutex);
  // Define key and value.
  std::string value;
  std::string key;
  // Read first line.
  std::string line;
  getline (inputStream, line);
  while (inputStream) {
    // Skip leading and trailing whitespace.
    uint st = lskipws (line, 0, line.size());
    if (line[st] != '#') {                         // skip if only comment
      uint end = rskipws (line, st, line.size());
      if (st < end) {                              // skip empty line
        uint nonbl = st;               // Position of last non-blank character
        uint stval = st;               // Start of value
        while (st<end) {
          if (line[st] == '"'  ||  line[st] == '\'') {
            st = skipQuoted (line, st);
            nonbl = st-1;                       // last char of quoted string
          } else if (line[st] == '#') {
            end = rskipws(line, stval, st);     // A comment ends the line
          } else {
            if (line[st] == '=') {
              if (! key.empty()) {
                // Add previous key/value to map.
                addMerge (key, value, merge);
                value.erase();
              }
              // Use ParameterValue to get the key without possible quotes.
              ParameterValue pvkey(line.substr(stval, nonbl-stval+1), false);
              key = pvkey.getString();
              if (key.empty())
								throw APSException("Empty key given in line " + line);
              key = prefix + key;
              stval = st+1;
            } else if (line[st] != ' '  &&  line[st] != '\t') {
              nonbl = st;                        // Position of last non-blank
            }
            st++;
          }
        }
        // Skip possible whitespace before the value.
        // A trailing backslash is processed for backward compatibility.
        if (line[end-1] == '\\') {
          end = rskipws(line, stval, end-1);
        }
        stval = lskipws(line, stval, end);
        if (stval < end) {
          // Append the line's non-empty value to the value.
          if (value.empty()) {
            value = line.substr(stval, end-stval);
          } else {
            // Add a blank if the continuation line is not quoted.
            if (line[stval] != '"'  &&  line[stval] != '\'') {
              value += ' ';
            }
            value += line.substr(stval, end-stval);
          }
        }
      }
    }
    // Get next line.
    getline (inputStream, line);
  }
  // Add last key.
  if (! key.empty()) {
    addMerge (key, value, merge);
  }
}


void ParameterSetImpl::addMerge (const std::string& key,
                                 const std::string& value,
                                 bool merge)
{
  // remove any existing value and insert this value
  if ((erase(key) > 0)  &&  !merge) {
    std::cout << "Key " + key + " is defined twice; ignoring first value";
  }
  addUnlocked (key, ParameterValue(value));
}

//------------------------- single pair functions ----------------------------
//
// findKV(key) [private]
//
ParameterSetImpl::const_iterator
ParameterSetImpl::findKV(const std::string& aKey, bool doThrow) const
{
        // Thread-safety.
        std::lock_guard<std::mutex> locker(itsMutex);

	const_iterator	iter = find(aKey);

	if (iter == end()) {
          if (doThrow) {
		throw APSException(formatString("Key %s unknown", aKey.c_str()));
          }
        } else {
          itsAskedParms.insert (aKey);          \
	}

	return (iter);
}

//
// add (key, value)
//
void ParameterSetImpl::add(const std::string& aKey, const ParameterValue& aValue)
{
  std::lock_guard<std::mutex> locker(itsMutex);
  addUnlocked (aKey, aValue);
}

void ParameterSetImpl::addUnlocked(const std::string& aKey,
                                   const ParameterValue& aValue)
{
  if (!insert(make_pair(aKey, aValue)).second) {
    throw APSException("add: Key " + aKey + " double defined?"); 
  }
}

//
// replace (key, value)
//
void ParameterSetImpl::replace(const std::string& aKey, const ParameterValue& aValue)
{
  std::lock_guard<std::mutex> locker(itsMutex);
  replaceUnlocked (aKey, aValue);
}

void ParameterSetImpl::replaceUnlocked(const std::string& aKey,
                                       const ParameterValue& aValue)
{
  (*this)[aKey] = aValue;
}

//
// remove (key)
//
void ParameterSetImpl::remove(const std::string& aKey)
{
  std::lock_guard<std::mutex> locker(itsMutex);
  // remove any existed value
  erase(aKey);
}

//
//-------------------------- retrieve functions -----------------------------
// getName
//
//string	ParameterSetImpl::getName() const
//{
//	string fullKeyName = begin()->first;
//	char*	firstPoint = strchr(fullKeyName.c_str(), '.');
//
//	return(fullKeyName.substr(0, firstPoint - fullKeyName.c_str()));
//}

//
// getVersionNr
//
//string	ParameterSetImpl::getVersionNr() const
//{
//	const_iterator	iter = find(getName()+"."+PC_KEY_VERSIONNR);
//
//	if (iter != end()) {
//		return (iter->second);
//	}
//	
//	return("");
//}

#define PARAMETERSETIMPL_GETVECTOR(TPC,TPL) \
std::vector<TPL> ParameterSetImpl::get##TPC##Vector(const std::string& aKey, \
                                               bool expandable) const \
{ \
  ParameterValue value (findKV(aKey)->second); \
  if (expandable) value = value.expand(); \
  return value.get##TPC##Vector(); \
} \
 \
std::vector<TPL> ParameterSetImpl::get##TPC##Vector(const std::string& aKey, \
                                               const std::vector<TPL>& aValue, \
                                               bool expandable) const   \
{ \
  const_iterator it = findKV(aKey,false); \
  if (it == end()) return aValue; \
  ParameterValue value (it->second); \
  if (expandable) value = value.expand(); \
  return value.get##TPC##Vector(); \
}

PARAMETERSETIMPL_GETVECTOR (Bool, bool)
PARAMETERSETIMPL_GETVECTOR (Int, int)
PARAMETERSETIMPL_GETVECTOR (Uint, uint)
PARAMETERSETIMPL_GETVECTOR (Int16, int16_t)
PARAMETERSETIMPL_GETVECTOR (Uint16, uint16_t)
PARAMETERSETIMPL_GETVECTOR (Int32, int32_t)
PARAMETERSETIMPL_GETVECTOR (Uint32, uint32_t)
PARAMETERSETIMPL_GETVECTOR (Int64, int64_t)
PARAMETERSETIMPL_GETVECTOR (Uint64, uint64_t)
PARAMETERSETIMPL_GETVECTOR (Float, float)
PARAMETERSETIMPL_GETVECTOR (Double, double)
PARAMETERSETIMPL_GETVECTOR (String, std::string)
PARAMETERSETIMPL_GETVECTOR (Time, time_t)


//---------------------------- save functions -------------------------------
//
// writeFile
//
// Writes the Key-Values pair from the current ParameterSetImpl to the given file
// thereby overwritting any file contents.
//
void ParameterSetImpl::writeFile(const std::string&	theFilename,
								    bool			append) const
{
	std::ofstream		paramFile;

	// Try to open the file
	paramFile.open(theFilename.c_str(), 
				   std::ofstream::out | (append ? std::ofstream::app : std::ofstream::trunc));
	if (!paramFile) {
		throw APSException(formatString("Unable to open file %s", theFilename.c_str()));
	}

	// Write all the pairs to the file
	writeStream(paramFile);

	// Close the file
	paramFile.close();
}

//
// writeBuffer
//
// Writes the Key-Values pair from the current ParameterSetImpl to the given 
// string.
//
void ParameterSetImpl::writeBuffer(std::string&	aBuffer) const
{
	std::ostringstream oss;
	writeStream(oss);
	aBuffer = oss.str();
}

//
// writeStream
//
// Writes the Key-Value pairs from the current ParameterSetImpl to the given
// output stream.
//
void ParameterSetImpl::writeStream(std::ostream&	os) const
{
        std::lock_guard<std::mutex> locker(itsMutex);
	// Write all the pairs to the file
	const_iterator		curPair = begin();
	while (curPair != end()) {
		// Key can always be written.
		os << curPair->first << "=";
                os << curPair->second.get() << '\n';
		curPair++;
	}
}


//
// isValidVersionNr(versionNr)
//
// Check format of versionnumber.
//
bool isValidVersionNr   (const std::string& versionNr)
{
	int		release, update, patch;
	char	toomuch[20];

	return (sscanf(versionNr.c_str(), "%d.%d.%d%10s", 
					&release, &update, &patch, toomuch) == 3);
}

//
// isValidVersionNrRef(versionNr)
//
// Check format of versionnumber, used as a reference.
//
bool isValidVersionNrRef(const std::string& versionNr)
{
	return (isValidVersionNr(versionNr) || (versionNr == PC_QUAL_STABLE) || 
		    (versionNr == PC_QUAL_TEST) || (versionNr == PC_QUAL_DEVELOP));

}

#if 0
//
// seqNr(aString)
//
// Check is given string is a valid sequencenumber
//
uint32	seqNr(const std::string& aString)
{
	int32	theNumber;

	sscanf(aString.c_str(), "%d", &theNumber);

	if (theNumber <= 0) {
		return (0);
	}

	return (theNumber);
}
#endif

//
// keyName(fullKeyName)
//
// Returns the real name of the key (without the module hierachy)
//
std::string keyName(const std::string& fullKeyName)
{
	std::string::size_type lastPoint = fullKeyName.rfind('.');
	if (lastPoint == std::string::npos) 
		return fullKeyName;
	else 
		return fullKeyName.substr(lastPoint+1);
}


//
// moduleName(fullKeyName)
//
// Returns the module hierarchy of key
//
std::string moduleName(const std::string& fullKeyName)
{
	std::string::size_type lastPoint = fullKeyName.rfind('.');
	if (lastPoint == std::string::npos) 
		return "";
	else 
		return fullKeyName.substr(0, lastPoint);
}

//
// keyPart(parameterline)
//
// Returns the key part of a parameter line.
//
std::string	keyPart	  (const std::string& parameterLine)
{
	std::string::size_type firstEqual = parameterLine.find('=');
	if (firstEqual == std::string::npos)
		return parameterLine;
	else
		return parameterLine.substr(0, firstEqual);
}

// Returns the value of a parameterline
std::string	valuePart   (const std::string& parameterLine)
{
	std::string::size_type firstEqual = parameterLine.find('=');
	if (firstEqual == std::string::npos)
		return parameterLine;
	else
		return parameterLine.substr(firstEqual+1);
}

// Returns the value of the index if the string contains an index otherwise
// 0 is returned. The second string contains the opening and closing chars
// that are used to indicate the index. The index must be a literal value
// not an expression.
int32_t 	indexValue (const std::string&	label, const char	indexMarker[2])
{
	std::string::size_type	start = label.find_last_of(indexMarker[0]);
	if (start == std::string::npos) {
		return (0);
	}

	std::string::size_type	end = label.find(indexMarker[1], start);
	if (end == std::string::npos) {
		return(0);
	}

	return (strtol(label.data()+start+1, 0 ,0));

}

//
// locateModule(shortKey)
//
// Searches for a key ending in the given 'shortkey' and returns it full name.
// e.g: a.b.c.d.param=xxxx --> locateModule(d)-->a.b.c.
std::string	ParameterSetImpl::locateModule(const std::string&	shortKey) const
{
	const_iterator		iter = begin();
	const_iterator		eom  = end();
	while ((iter != eom)) {
		if (keyName(moduleName(iter->first)) == shortKey) {
			std::string prefix = moduleName(moduleName((iter->first)));
			if (prefix.length() > 0) {
				prefix += ".";
			}
			return (prefix);
		}
		iter++;
	}
	return ("");
}

//
// fullModuleName(shortKey)
//
// Searches for a key ending in the given 'shortkey' and returns it full name.
// e.g: a.b.c.d.param=xxxx --> fullModuleName(d)      --> a.b.c.d
// e.g: a.b.c.d.param=xxxx --> fullModuleName(b.c)    --> a.b.c
// e.g: a.b.c.d.param=xxxx --> fullModuleName(d.param)-->
std::string	ParameterSetImpl::fullModuleName(const std::string&	shortKey) const
{
	const_iterator		iter = begin();
	const_iterator		eom  = end();
	while ((iter != eom)) {
		std::string::size_type	start = moduleName(iter->first).rfind(shortKey);
		if (start != std::string::npos) {
			std::string::size_type	keyLen = shortKey.length();
			std::string::size_type	last = start+keyLen;
			std::string::size_type	iterLen = (iter->first).length();
			if ((last == iterLen || (last < iterLen && (iter->first)[last] == '.')) &&
				(start == 0 || (start > 0 && (iter->first)[start-1] == '.'))) {
				std::string prefix = iter->first.substr(0, start);
				return (prefix+shortKey);
			}
		}
		iter++;
	}

	return ("");
}

std::vector<std::string> ParameterSetImpl::unusedKeys() const
{
  std::lock_guard<std::mutex> locker(itsMutex);
  std::vector<std::string> vec;
  for (const_iterator iter = begin(); iter != end(); ++iter) {
    if (itsAskedParms.find (iter->first) == itsAskedParms.end()) {
      vec.push_back (iter->first);
    }
  }
  return vec;
}

} // namespace LOFAR
