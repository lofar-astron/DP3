//# makesourcedb.cc: Fill a SourceDB from an ASCII file
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
//# $Id: makesourcedb.cc 37340 2017-05-11 12:39:06Z dijkema $

// This program writes patch and source information into a SourceDB.
// The input is read from an ASCII file that can be formatted in various ways.
// By default the ra/dec/flux of a patch is determined from the sources it
// contains by adding the fluxes and taking the flux-weighted average of the
// ra/dec in xyz-coordinates.
//
// The program can be run as:
//   makesourcedb in=inname out=outname format="fmt" append=true/false (or 1/0)
//                average=true/false
// in      gives the input file name
// out     gives the sourcedb name which will be created or appended
// format  defines the format of the input file
// append  defines if the sourcedb is created or appended (default is appended)
// average defines if the average patch ra/dec are calculated (default is true)
// center  defines the field center (ra and dec) of a search cone or box
// radius  defines the radius if searching using a cone
// width   defines the widths in ra and dec if searching using a box
//
//// Also minflux (skip if < minflux), maxflux (always if >= maxflux),
//// beammodel (only WSRT for time being); if beammodel, convert to app flux. 
//
// The format string can be given as an argument. If its value starts with a <
// it means that it is read from the file following it. No file means from the
// text file (headers) itself. A format string in a file is like a normal
// format string, but preceeded with 'format='. For backward compatibility
// it can also be preceeded by '#(' and followed by ')=format'.
// Whitespace is allowed between these characters.
//
// A format string looks like:
//    name type ra dec cat
// It defines the fields in the input file and which separator is used between
// the fields (in this example all fields are separated by whitespace).
// Note that both in the format string and input file whitespace is ignored
// unless used as separator.
// The default format string is:
//     "Name,Type,Ra,Dec,I,Q,U,V,SpectralIndex,MajorAxis,MinorAxis,Orientation"
// thus all fields are separated by commas.
// If the format string contains only one character, the default format is used
// with the given character as separator.
// A field name can consists of alphanumeric characters, underscores and colons.
// However, a colon can not be used as the first character.
// In this way a colon can be used as separator as long as it is surrounded by
// whitespace (otherwise it would be part of the field name).
//
// It is possible to define default values in the format string by giving a
// value to a field. The default value will be used if the field in the
// input file is empty (see below). The value must be enclosed in single
// or double quotes. For example:
//     "Name,Type,Ra,Dec,I,Q,U,V,SpectralIndex='1',ReferenceFrequency='1e9'"
// If no default value is given for empty field values, an appropriate default
// is used (which is usually 0).
//
// In a similar way it is possible to define fixed values by predeeding the
// value with the string 'fixed'. For example:
//     "Name,Type,Ra,Dec,I,Q,U,V,Major,Minor,Phi, Category=fixed'2'"
// It means that Category is not a field in the input file, but the given
// value is used for all patches and/or sources in the input file.
// So this example will make all the patches/sources Cat2.
//
// It is possible to ignore a column in the input file by giving an empty name
// or the name 'dummy' in the format string. For example:
//     "Name,Type,Ra,Dec,,dummy,I,Q,U,V
// will ignore the two columns in the input file between Dec and I.

// Each line in the input file represents a patch or source.
// Lines starting with a # (possibly preceeded by whitespace) are ignored.
// The order of the fields in the input file must (of course) match the
// order of the fields in the format string. Fields with fixed values can be
// put anywhere in the format string.
// Each value in the input can be enclosed in single or double quotes. Quotes
// must be used if the value contains the separator character.
// A field in the input can be empty. If whitespace is used as separator, an
// empty field must be represented by "" (or '').
// An input line can contain less values than there are fields in the format
// string. Missing values are empty.
//
// An input line can define a patch or a source.
// A patch is defined by having an empty source name, otherwise it is a source.
// Thus an empty source name and empty patch name is invalid.
//
// Currently 3 source types are supported (Point, Gaussian, and Shapelet).
//
// Ra can be specified in various forms.
// If it is a single value, it must be in the Casacore MVAngle format which is:
// - a value followed by a unit (e.g. 1.2rad). Default unit is deg.
// - a string like hh:mm:ss.s (or dd.mm.ss.s). H and m (or d and m) can
//   be used instead of : (or .).
// However, some catalogues use whitespace to separate hh, mm, and ss.
// They can be seen as individual fields (with whitespace as separator).
// Therefore some extra format fields exist which are Rah, Ram, and Ras.
// They define the hh, mm, and ss parts. Instead of Rah one can use Rad which
// defines it as degrees instead of hours.
// <br>The same is true for Dec (which extra fields Dech, Decd, Decm, and Decs).
// Please note that in a value like '-10 23 45.6' the mm and ss parts are
// also treated negative, thus it is handled as -10:23:45.6
//
// It is possible to specify a reference type string for the source position.
// If such a column is not given, it defaults to J2000. The reference type
// given must be a valid casacore measures type (which is case-insensitive).

// See the various test/tmakesourcedb files for an example.

#include "SourceDB.h"

#include "../Common/StringUtil.h"
#include "../Common/StreamUtil.h"

#include <string>                //# for getline
#include <iostream>
#include <fstream>
#include <sstream>

#include <casacore/casa/OS/Path.h>
#include <casacore/casa/Quanta/MVAngle.h>
#include <casacore/casa/Inputs/Input.h>
#include <casacore/casa/BasicSL/Constants.h>

#include <boost/algorithm/string/case_conv.hpp>

using namespace std;
using namespace casacore;
using namespace DP3;
using namespace BBS;

// Define the sequence nrs of the various fields.
enum FieldNr {
  // First the standard fields.
  NameNr, TypeNr, RefTypeNr, RaNr, DecNr,
  INr, QNr, UNr, VNr, SpInxNr, LogSINr, RefFreqNr,
  MajorNr, MinorNr, OrientNr, RotMeasNr, PolFracNr, PolAngNr, RefWavelNr,
  IShapeletNr, QShapeletNr, UShapeletNr, VShapeletNr,
  NrKnownFields,
  // Now other fields
  CatNr=NrKnownFields, PatchNr,
  RahNr, RadNr, RamNr, RasNr, DechNr, DecdNr, DecmNr, DecsNr,
  // Nr of fields.
  NrFields
};

vector<string> fillKnown()
{
  // The order in which the names are pushed must match exactly with the
  // enum above.
  vector<string> names;
  names.reserve (NrFields);
  names.push_back ("Name");
  names.push_back ("Type");
  names.push_back ("RefType");
  names.push_back ("Ra");
  names.push_back ("Dec");
  names.push_back ("I");
  names.push_back ("Q");
  names.push_back ("U");
  names.push_back ("V");
  names.push_back ("SpectralIndex");
  names.push_back ("LogarithmicSI");
  names.push_back ("ReferenceFrequency");
  names.push_back ("MajorAxis");
  names.push_back ("MinorAxis");
  names.push_back ("Orientation");
  names.push_back ("RotationMeasure");
  names.push_back ("PolarizedFraction");
  names.push_back ("PolarizationAngle");
  names.push_back ("ReferenceWavelength");
  names.push_back ("IShapelet");
  names.push_back ("QShapelet");
  names.push_back ("UShapelet");
  names.push_back ("VShapelet");
  names.push_back ("Category");
  names.push_back ("Patch");
  names.push_back ("rah");
  names.push_back ("rad");
  names.push_back ("ram");
  names.push_back ("ras");
  names.push_back ("dech");
  names.push_back ("decd");
  names.push_back ("decm");
  names.push_back ("decs");
  assert (names.size() == NrFields);
  return names;
}
// Define all field names.
vector<string> theFieldNames = fillKnown();

enum FieldType {
  KNOWNFIELD = 1,
  FIXEDVALUE = 2,
  SKIPFIELD  = 4,
  DEFAULTVALUE=8
};

struct SdbFormat
{
  // Define fieldnr of the known fields.
  // -1 means not given.
  vector<int> fieldNrs;
  // Define the separators.
  vector<char> sep;
  // Define the field names.
  vector<string> names;
  // Define the field type.
  vector<int> types;
  // Define the fixed values.
  vector<string> values;
};

// Read a line and remove a possible carriage-return at the end.
void getInLine (istream& infile, string& line)
{
  getline (infile, line);
  int sz = line.size();
  if (sz > 0  &&  line[sz-1] == '\r') {
    line = line.substr(0,sz-1);
  }
}

void checkRaDec (const SdbFormat& sdbf, int nr,
                 int hnr, int dnr, int mnr, int snr)
{
  int  v = sdbf.fieldNrs[nr];
  int hv = sdbf.fieldNrs[hnr];
  int dv = sdbf.fieldNrs[dnr];
  int mv = sdbf.fieldNrs[mnr];
  int sv = sdbf.fieldNrs[snr];
  if (hv>=0 && dv>=0)
		throw std::runtime_error(
             "rah and rad cannot be used both (same for dec)");
  if (v >= 0) {
    if(!(hv<0 && dv<0 && mv<0 && sv<0))
			throw std::runtime_error(
              "rah/d/m/s cannot be used if ra is given (same for dec)");
  }
  if(!(v>=0 || hv>=0 || dv>=0 || mv>=0 || sv>=0))
		throw std::runtime_error(
            "No Ra or Dec info given");
}

uint ltrim (const string& value, uint st, uint end)
{
  for (; st<end; ++st) {
    if (value[st] != ' '  &&  value[st] != '\t') {
      break;
    }
  }
  return st;
}
  
uint rtrim (const string& value, uint st, uint end)
{
  for (; end>st; --end) {
    if (value[end-1] != ' '  &&  value[end-1] != '\t') {
      break;
    }
  }
  return end;
}

// Get the next value by looking for the separator.
// The separator is ignored in parts enclosed in quotes or square brackets.
// Square brackets can be nested (they indicate arrays).
uint nextValue (const string& str, char sep, uint st, uint end)
{
  uint posbracket = 0;
  int nbracket = 0;
  while (st < end) {
    if (str[st] == '\''  ||  str[st] == '"') {
      // Ignore a quoted part.
      string::size_type pos = str.find (str[st], st+1);
      if (pos == string::npos)
				throw std::runtime_error("Unbalanced quoted string at position "
                 + to_string(st) + " in " + str);
      st = pos;
    } else if (str[st] == '[') {
      nbracket++;
      posbracket = st;
    } else if (str[st] == ']') {
      if (nbracket<=0)
				throw std::runtime_error("Unbalanced square brackets at position "
             + to_string(st) + " in " + str);
      nbracket--;
    } else if (nbracket == 0) {
      if (str[st] == sep) {
        return st;
      }
    }
    ++st;
  }
  if (nbracket!=0)
				throw std::runtime_error("Unbalanced square brackets at position "
             + to_string(posbracket) + " in " + str);
  return end;
}

SdbFormat getFormat (const string& format)
{
  // Fill the map with known names.
  map<string,int> nameMap;
  for (uint i=0; i<theFieldNames.size(); ++i) {
    nameMap[boost::to_lower_copy(theFieldNames[i])] = i;
  }
  // Skip possible left and right whitespace.
  uint end = format.size();
  uint st  = ltrim (format, 0, end);
  end      = rtrim (format, st, end);
  // Initialize the format.
  SdbFormat sdbf;
  sdbf.fieldNrs.resize (NrFields);
  for (uint i=0; i<NrFields; ++i) {
    sdbf.fieldNrs[i] = -1;
  }
  // Use default if the string is empty.
  if (st >= end-1) {
    char sep = ',';
    if (st == end-1) sep = format[st];
    for (uint i=0; i<NrKnownFields; ++i) {
      sdbf.fieldNrs[i] = i;
      sdbf.sep.push_back (sep);
      sdbf.names.push_back (theFieldNames[i]);
      sdbf.types.push_back (KNOWNFIELD);
      sdbf.values.push_back ("");
    }
    return sdbf;
  }
  // Parse the format string.
  uint nr = 0;
  uint i = st;
  while (i < end) {
    if ((format[i] >= 'a'  &&  format[i] <= 'z')  ||
        (format[i] >= 'A'  &&  format[i] <= 'Z')  ||
        (format[i] >= '0'  &&  format[i] <= '9')  ||
        format[i] == '_'  ||  (format[i] == ':'  &&  i > st)) {
      ++i;    // part of name
    } else {
      // End of name
      string name  = format.substr(st, i-st);
      string lname = boost::to_lower_copy(name);
      int fieldType = 0;
      if (lname.empty()  ||  lname == "dummy") {
        fieldType = SKIPFIELD;
      } else {
        /// Remove ASSERT branch once we're sure SpectralIndexDegree is not used anymore.
        if (lname == "spectralindexdegree")
					throw std::runtime_error(
                   "Use SpectralIndex=[v1,v2,...] instead of SpectralIndexDegree "
                   "and SpectralIndex:i");
        map<string,int>::const_iterator namepos = nameMap.find(lname);
        // Fill in fieldnr of a known field.
        if (namepos != nameMap.end()) {
          fieldType = KNOWNFIELD;
          sdbf.fieldNrs[namepos->second] = nr;
        }
      }
      // See if a default or fixed value is given.
      string fixedValue;
      i = ltrim(format, i, end);
      if (i < end  &&  format[i] == '=') {
        i = ltrim(format, i+1, end);
        // See if it is a fixed value.
        bool isDefault = True;
        if (i+5 < end   &&  boost::to_lower_copy(format.substr(i,5)) == "fixed") {
          isDefault = False;
          i += 5;
        }
        if (i>=end)
					throw std::runtime_error("No value given after " + name + '=');
        if (format[i] != '"' && format[i] != '\'')
                   throw std::runtime_error("value after " + name + "= must start with a quote");
        // Find the ending quote.
        string::size_type pos = format.find (format[i], i+1);
        if (pos==string::npos)
                   throw std::runtime_error("No closing value quote given after " + name + '=');
        fixedValue = format.substr(i+1, pos-i-1);
        i = ltrim(format, pos+1, end);
        if (isDefault) {
          fieldType |= DEFAULTVALUE;
        } else {
          fieldType |= FIXEDVALUE;
        }
      }
      // Now look for a separator.
      char sep = ' ';
      if (!((format[i] >= 'a'  &&  format[i] <= 'z')  ||
            (format[i] >= 'A'  &&  format[i] <= 'Z')  ||
            (format[i] >= '0'  &&  format[i] <= '9')  ||
            format[i] == '_')) {
        sep = format[i];
        if (sep=='"' || sep=='\'')
					throw std::runtime_error("A quote is found as separator; "
                   "probably a quote around a value in the format string is "
                   "missing");
        i = ltrim(format, i+1, end);
      }
      sdbf.sep.push_back (sep);
      sdbf.names.push_back (name);
      sdbf.types.push_back (fieldType);
      sdbf.values.push_back (fixedValue);
      nr++;
      st = i;
    }
  }
  if (st < end) {
    // The last item was just a name which has to be processed.
    string name = format.substr(st, end-st);
    string lname = boost::to_lower_copy(name);
    int fieldType = 0;
    if (lname.empty()  ||  lname == "dummy") {
      fieldType = SKIPFIELD;
    } else {
      map<string,int>::const_iterator namepos = nameMap.find(lname);
      // Fill in fieldnr of a known field.
      if (namepos != nameMap.end()) {
        sdbf.fieldNrs[namepos->second] = nr;
        fieldType = KNOWNFIELD;
      }
    }
    sdbf.sep.push_back (' ');
    sdbf.names.push_back (name);
    sdbf.types.push_back (fieldType);
    sdbf.values.push_back ("");
  }
  // Make sure Ra and Dec are given correctly.
  checkRaDec (sdbf, RaNr, RahNr, RadNr, RamNr, RasNr);
  checkRaDec (sdbf, DecNr, DechNr, DecdNr, DecmNr, DecsNr);
  return sdbf;
}

SourceInfo::Type string2type (const string& str)
{
  string s = boost::to_lower_copy(str);
  if (s == "point"  ||  s.empty()) {
    return SourceInfo::POINT;
  } else if (s == "gaussian") {
    return SourceInfo::GAUSSIAN;
  } else if (s == "disk") {
    return SourceInfo::DISK;
  } else if (s == "shapelet") {
    return SourceInfo::SHAPELET;
  }
  throw std::runtime_error(str + " is an invalid source type");
}

string unquote (const string& value)
{
  string res(value);
  if (res.size() > 1) {
    int last = res.size() - 1;
    if (last >= 1  &&  ((res[0] == '"'  && res[last] == '"')  ||
                        (res[0] == '\'' && res[last] == '\''))) {
      res = res.substr(1,last-1);
    }
  }
  return res;
}

string getValue (const vector<string>& values, int nr,
                 const string& defVal=string())
{
  if (nr < 0) {
    return defVal;
  }
  return unquote (values[nr]);
}

bool string2bool (const vector<string>& values, int nr, bool defVal)
{
  string value = getValue (values, nr);
  if (value.empty()) {
    return defVal;
  }
  return strToBool (value);
}

int string2int (const vector<string>& values, int nr, int defVal)
{
  string value = getValue (values, nr);
  if (value.empty()) {
    return defVal;
  }
  return strToInt (value);
}

double string2real (const string& value, double defVal)
{
  if (value.empty()) {
    return defVal;
  }
  return strToDouble (value);
}

double string2real (const vector<string>& values, int nr, double defVal)
{
  return string2real (getValue (values, nr), defVal);
}

// Convert values in a possibly bracketed string to a vector of strings
// taking quoted or bracketed values into account.
// A comma isused as separator.
vector<string> string2vector (const string& value, const vector<string>& defVal)
{
  vector<string> result;
  // Test if anything is given.
  if (value.empty()) {
    result = defVal;
  } else {
    uint end = value.size();
    // If no brackets given, it is a single value.
    if (value.size() < 2  ||  value[0] != '['  ||  value[end-1] != ']') {
      result.push_back (value);
    } else {
      // Skip opening and closing bracket and possible whitespace.
      uint st = ltrim(value, 1, end-1);
      end = rtrim(value, st, end-1);
      while (st < end) {
        uint pos = nextValue (value, ',' , st, end);
        result.push_back (value.substr(st, rtrim(value, st, pos) - st));
        st = ltrim (value, pos+1, end);
      }
    }
  }
  return result;
}

vector<string> string2vector (const vector<string>& values, int nr,
                              const vector<string>& defVal)
{
  return string2vector (getValue (values, nr), defVal);
}

vector<double> vector2real (const vector<string>& values, double defVal)
{
  vector<double> result;
  result.reserve (values.size());
  for (uint i=0; i<values.size(); ++i) {
    result.push_back (string2real (values[i], defVal));
  }
  return result;
}

double string2pos (const vector<string>& values, int pnr, int hnr, int dnr,
                   int mnr, int snr, bool canUseColon)
{
  double deg = 0;
  bool fnd = false;
  if (pnr >= 0) {
    string value = getValue(values, pnr);
    if (! value.empty()) {
      if (!canUseColon) {
        if (value.find(':') != string::npos)
					throw std::runtime_error(
                   "Colons cannot be used in declination value " + value);
      }
      Quantity q;
      if (!MVAngle::read (q, values[pnr]))
				throw std::runtime_error("Error in reading position " + values[pnr]);
      deg = q.getValue ("deg");
      fnd = true;
    }
  } else {
    if (hnr >= 0) {
      string value = getValue(values, hnr);
      if (! value.empty()) {
        deg = string2real (values[hnr], 0);
        fnd = true;
      }
    } else if (dnr >= 0) {
      string value = getValue(values, dnr);
      if (! value.empty()) {
        deg = string2real (values[dnr], 0);
        fnd = true;
      }
    }
    double ms = 0;
    if (mnr >= 0) {
      string value = getValue(values, mnr);
      if (! value.empty()) {
        ms = string2real (values[mnr], 0);
        fnd = true;
      }
    }
    if (snr >= 0) {
      string value = getValue(values, snr);
      if (! value.empty()) {
        ms += string2real (values[snr], 0) / 60;
        fnd = true;
      }
    }
    if (deg < 0) {
      deg -= ms/60;
    } else {
      deg += ms/60;
    }
    if (hnr >= 0) {
      deg *= 15;
    }
  }
  if (fnd) {
    Quantity q (deg, "deg");
    return q.getValue ("rad");
  }
  return 1e-9;
}

void checkRefType (const string& refType)
{
  string type = boost::to_upper_copy(refType);
  if (type != "J2000"   &&  type != "B1950"    &&
      type != "SUN"     &&  type != "MOON"     &&  type != "VENUS"   &&
      type != "MARS"    &&  type != "JUPITER"  &&  type != "SATURN"  &&
      type != "URANUS"  &&  type != "NEPTUNE"  &&  type != "MERCURY") {
    throw std::runtime_error("Reference type " + refType + " is incorrect");
  }
}

struct SearchInfo
{
  double ra, dec, sinDec, cosDec, cosRadius, raStart, raEnd, decStart, decEnd;
  bool   search;        // false no search info given
  bool   asCone;        // true is search in a cone, otherwise a box
};

// Get the search cone or box values.
SearchInfo getSearchInfo (const string& center, const string& radius,
                          const string& width)
{
  SearchInfo searchInfo;
  if (center.empty()) {
    searchInfo.search = false;
  } else {
    searchInfo.search = true;
    vector<string> pos = StringUtil::split (center, ',');
    if(pos.size() != 2) throw std::runtime_error("center not specified as ra,dec");
    searchInfo.ra  = string2pos (pos, 0, -1, -1, -1, -1, true);
    searchInfo.dec = string2pos (pos, 1, -1, -1, -1, -1, false);
    searchInfo.sinDec = sin(searchInfo.dec);
    searchInfo.cosDec = cos(searchInfo.dec);
    if(radius.empty() == width.empty())
               throw std::runtime_error("radius OR width must be given if center is given (not both)");
    if (radius.empty()) {
      double raw, decw;
      searchInfo.asCone = false;
      pos = StringUtil::split(width, ',');
      if (pos.size() != 1 && pos.size() != 2)
        throw std::runtime_error("width should be specified as 1 or 2 values");
      raw = string2pos (pos, 0, -1, -1, -1, -1, true);
      if (pos.size() > 1) {
        decw = string2pos (pos, 1, -1, -1, -1, -1, false);
      } else {
        decw = raw;
      }
      searchInfo.raStart  = searchInfo.ra  - raw/2;
      searchInfo.raEnd    = searchInfo.ra  + raw/2;
      searchInfo.decStart = searchInfo.dec - decw/2;
      searchInfo.decEnd   = searchInfo.dec + decw/2;
    } else {
      searchInfo.asCone = true;
      pos[0] = radius;
      searchInfo.cosRadius = cos (string2pos (pos, 0, -1, -1, -1, -1, false));
    }
  }
  return searchInfo;
}

bool matchSearchInfo (double ra, double dec, const SearchInfo& si)
{
  if (!si.search) {
    return true;
  }
  bool match = false;
  if (si.asCone) {
    match = (si.cosRadius <=
             si.sinDec * sin(dec) + si.cosDec * cos(dec) * cos(si.ra - ra));
  } else {
    // Ra can be around 0 or 360 degrees, so make sure all cases are handled.
    ra -= C::_2pi;
    for (int i=0; i<4; ++i) {
      if (ra >= si.raStart  &&  ra <= si.raEnd) {
        match = true;
        break;
      }
      ra += C::_2pi;
    }
    if (match) {
      match = (dec >= si.decStart  &&  dec <= si.decEnd);
    }
  }
  return match;
}

void addValue (ParmMap& fieldValues, const string& name, double value)
{
  fieldValues.define (name, ParmValueSet(ParmValue(value)));
}

void add (ParmMap& fieldValues, FieldNr field, double value)
{
  fieldValues.define (theFieldNames[field], ParmValueSet(ParmValue(value)));
}

void addSpInx (ParmMap& fieldValues, const vector<double>& spinx,
               double refFreq)
{
  if (spinx.size() > 0) {
    if (refFreq<=0)
			throw std::runtime_error("SpectralIndex given, but no ReferenceFrequency");
    /// Remove the following lines if not needed anymore for BBS.
    addValue (fieldValues, "SpectralIndexDegree", int(spinx.size()) - 1);
    for (uint i=0; i<spinx.size(); ++i) {
      ostringstream ostr;
      ostr << "SpectralIndex:" << i;
      addValue (fieldValues, ostr.str(), spinx[i]);
    }
  }
}

void readShapelet (const string& fileName, Array<double>& coeff,
                   double& scale)
{
  ifstream file(fileName.c_str());
  if (!file)
		throw std::runtime_error("Shapelet file " + fileName + " could not be opened");
  string line;
  getInLine (file, line);    // ra dec
  getInLine (file, line);    // order scale
  vector<string> parts = StringUtil::split (line, ' ');
  if (parts.size() != 2)
		throw std::runtime_error("Expected 2 values in shapelet line " + line);
  int order = string2int (parts, 0, 0);
  scale = string2real (parts, 1, 0.);
  if (order <= 0)
		throw std::runtime_error("Invalid order in shapelet line " + line);
  coeff.resize (IPosition(2, order, order));
  double* coeffData = coeff.data();
  for (uint i=0; i<coeff.size(); ++i) {
    getInLine (file, line);    // index coeff
    vector<string> parts = StringUtil::split (line, ' ');
    if (parts.size() != 2)
			throw std::runtime_error("Expected 2 values in shapelet line " + line);
    if (string2int(parts, 0, -1) != int(i))
			throw std::runtime_error("Expected shapelet line with index " + to_string(i));
    *coeffData++ = string2real (parts, 1, 0.);
  }
}

void fillShapelet (SourceInfo& srcInfo,
                   const string& shpI,
                   const string& shpQ,
                   const string& shpU,
                   const string& shpV)
{
  double scaleI = 0;
  double scaleQ = 0;
  double scaleU = 0;
  double scaleV = 0;
  Array<double> coeffI, coeffQ, coeffU, coeffV;
  readShapelet (shpI, coeffI, scaleI);
  if (shpQ.empty()) {
    coeffQ = coeffI;
    scaleQ = scaleI;
  } else {
    readShapelet (shpQ, coeffQ, scaleQ);
  }
  if (shpU.empty()) {
    coeffU = coeffI;
    scaleU = scaleI;
  } else {
    readShapelet (shpU, coeffU, scaleU);
  }
  if (shpV.empty()) {
    coeffV = coeffI;
    scaleV = scaleI;
  } else {
    readShapelet (shpV, coeffV, scaleV);
  }
  srcInfo.setShapeletCoeff (coeffI, coeffQ, coeffU, coeffV);
  srcInfo.setShapeletScale (scaleI, scaleQ, scaleU, scaleV);
}

// Calculate the polarization angle and polarized fraction given Q and U
// for a given reference wavelength.
// A spectral index can be used to calculate Stokes I.
void calcRMParam (double& polfrac, double& polang,
                  double fluxi0, double fluxq, double fluxu,
                  const vector<double>& spinx, double rm,
                  double refFreq, double rmRefWavel)
{
  // polfrac = sqrt(q^2 + u^2) / i
  // where i = i(0) * spinx
  // Compute spectral index for the RM reference wavelength as:
  // (v / v0) ^ (c0 + c1 * log10(v / v0) + c2 * log10(v / v0)^2 + ...)
  // Where v is the RM frequency and v0 is the spinx reference frequency.
  double si = 1;
  if (spinx.size() > 0) {
    if (rmRefWavel <= 0)
			throw std::runtime_error("No RM reference wavelength given");
    double rmFreq = C::c / rmRefWavel;
    double vv0 = rmFreq / refFreq;
    double factor = 1;
    double sum = 0;
    for (uint i=0; i<spinx.size(); ++i) {
      sum += factor * spinx[i];
      factor *= log10(vv0);
    }
    si = std::pow(vv0, sum);
  }
  double fluxi = fluxi0 * si;
  polfrac = sqrt(fluxq*fluxq + fluxu*fluxu) / fluxi;
  // Calculate polang(0) from Q and U given for the reference lambda.
  // polang = atan2(u,q) / 2
  // polang(lambda) = polang(0) + lambda^2 * rm
  // Scale between 0 and pi.
  double pa = 0.5 * atan2(fluxu, fluxq) - rmRefWavel*rmRefWavel*rm;
  polang = fmod(pa, C::pi);
  if (polang < 0) {
    polang += C::pi;
  }
}

void process (const string& line, SourceDB& pdb, const SdbFormat& sdbf,
              const string& prefix, const string& suffix,
              bool check, int& nrpatch, int& nrsource,
              int& nrpatchfnd, int& nrsourcefnd,
              map<string,PatchSumInfo>& patchSumInfo,
              const SearchInfo& searchInfo)
{
  //  cout << line << endl;
  // Hold the values.
  ParmMap fieldValues;
  vector<string> values;
  // Process the line.
  uint end = line.size();
  uint st  = ltrim(line, 0, end);
  for (uint i=0; i<sdbf.names.size(); ++i) {
    string value;
    if ((sdbf.types[i] & FIXEDVALUE) == FIXEDVALUE) {
      value = sdbf.values[i];
    } else if (st < end) {
      uint pos = nextValue (line, sdbf.sep[i], st, end);
      value = line.substr(st, rtrim(line, st, pos) - st);
      st = ltrim (line, pos+1, end);
    }
    if (value.empty()  &&  (sdbf.types[i] & DEFAULTVALUE) == DEFAULTVALUE) {
      value = sdbf.values[i];
    }
    values.push_back (value);
    if ((sdbf.types[i] & SKIPFIELD) != SKIPFIELD) {
      if ((sdbf.types[i] & KNOWNFIELD) != KNOWNFIELD) {
        addValue (fieldValues, sdbf.names[i], string2real(unquote(value), 0));
      }
    }
  }
  // Now handle the standard fields.
  string srcName = getValue(values, sdbf.fieldNrs[NameNr]);
  string refType = getValue(values, sdbf.fieldNrs[RefTypeNr], "J2000");
  if (refType.empty()) {
    refType = "J2000";
  }
  checkRefType (refType);
  SourceInfo::Type srctype = string2type (getValue(values,
                                                   sdbf.fieldNrs[TypeNr]));
  double ra = string2pos (values, sdbf.fieldNrs[RaNr],
                          sdbf.fieldNrs[RahNr],
                          sdbf.fieldNrs[RadNr],
                          sdbf.fieldNrs[RamNr],
                          sdbf.fieldNrs[RasNr],
                          true);
  double dec = string2pos (values, sdbf.fieldNrs[DecNr],
                           sdbf.fieldNrs[DechNr],
                           sdbf.fieldNrs[DecdNr],
                           sdbf.fieldNrs[DecmNr],
                           sdbf.fieldNrs[DecsNr],
                           false);
  if (!(ra>-6.3 && ra<6.3 && dec>-1.6 && dec<1.6))
		throw std::runtime_error("RA " + to_string(ra)
             + " or DEC " + to_string(dec) + " radians is outside boundaries");
  int cat = string2int (values, sdbf.fieldNrs[CatNr], 2);
  double fluxI     = string2real (getValue (values, sdbf.fieldNrs[INr]), 1.);
  string fluxQ     = getValue (values, sdbf.fieldNrs[QNr]);
  string fluxU     = getValue (values, sdbf.fieldNrs[UNr]);
  string fluxV     = getValue (values, sdbf.fieldNrs[VNr]);
  string rm        = getValue (values, sdbf.fieldNrs[RotMeasNr]);
  string polFrac   = getValue (values, sdbf.fieldNrs[PolFracNr]);
  string polAng    = getValue (values, sdbf.fieldNrs[PolAngNr]);
  string refWavel  = getValue (values, sdbf.fieldNrs[RefWavelNr]);
  string shapeletI = getValue (values, sdbf.fieldNrs[IShapeletNr]);
  string shapeletQ = getValue (values, sdbf.fieldNrs[QShapeletNr]);
  string shapeletU = getValue (values, sdbf.fieldNrs[UShapeletNr]);
  string shapeletV = getValue (values, sdbf.fieldNrs[VShapeletNr]);

  vector<double> spinx(vector2real(string2vector(values,
                                                 sdbf.fieldNrs[SpInxNr],
                                                 vector<string>()),
                                   0.));
  double refFreq = string2real (values, sdbf.fieldNrs[RefFreqNr], 0);
  bool useLogSI = string2bool (values, sdbf.fieldNrs[LogSINr], true);
  bool useRM = false;
  double rmRefWavel = 0;
  if (rm.empty()) {
    if (!polFrac.empty() || !polAng.empty() || !refWavel.empty())
			throw std::runtime_error(
               "PolarizationAngle, PolarizedFraction, and ReferenceWavelength"
               " cannot be specified if RotationMeasure is not specified");
  } else {
    if (!fluxQ.empty() || !fluxU.empty()) {
      if (fluxQ.empty() || fluxU.empty() ||
                 !polFrac.empty() || !polAng.empty())
			throw std::runtime_error(
                 "PolarizationAngle/PolarizedFraction or Q/U must be "
                 "specified if RotationMeasure is specified");
      useRM = true;
      if (refWavel.empty()) {
        if (refFreq <= 0)
			throw std::runtime_error(
                   "For rotation measures the reference frequency or "
                   "wavelength must be given");
      }
      rmRefWavel = string2real (refWavel, C::c / refFreq);
    } else {
      if (polFrac.empty() || polAng.empty())
				throw std::runtime_error(
                 "PolarizationAngle/PolarizedFraction or Q/U must be "
                 "specified if RotationMeasure is specified");
      useRM = true;
      rmRefWavel = string2real (refWavel, 0);
      if (rmRefWavel != 0)
				 throw std::runtime_error("PolarizationAngle/PolarizedFraction can "
                 "only be given for ReferenceWavelength=0");
    }
  }               
  SourceInfo srcInfo(srcName, srctype, refType, useLogSI, spinx.size(), refFreq, useRM);
  if (srctype == SourceInfo::SHAPELET) {
    fillShapelet (srcInfo, shapeletI, shapeletQ, shapeletU, shapeletV);
  }
  add (fieldValues, INr, fluxI);
  double rmval = 0;
  if (!rm.empty()) {
    rmval = string2real(rm, 0.);
    add (fieldValues, RotMeasNr, rmval);
  }
  double fq = string2real(fluxQ, 0.);
  double fu = string2real(fluxU, 0.);
  if (useRM) {
    double pfrac = string2real(polFrac, 0.);
    double pang  = string2real(polAng, 0.);
    if (! fluxQ.empty()) {
      calcRMParam (pfrac, pang, fluxI, fq, fu,
                   spinx, rmval, refFreq, rmRefWavel);
    }
    add (fieldValues, PolFracNr, pfrac);
    add (fieldValues, PolAngNr, pang);
  } else {
    add (fieldValues, QNr, string2real(fluxQ, 0.));
    add (fieldValues, UNr, string2real(fluxU, 0.));
  }
  add (fieldValues, VNr, string2real(fluxV, 0.));
  addSpInx (fieldValues, spinx, refFreq);
  if (refFreq > 0) {
    add (fieldValues, RefFreqNr, refFreq);
  }
  string patch = getValue (values, sdbf.fieldNrs[PatchNr]);

  if (srctype == SourceInfo::GAUSSIAN) {
    add (fieldValues, MajorNr,
         string2real (values, sdbf.fieldNrs[MajorNr], 1));
    add (fieldValues, MinorNr,
         string2real (values, sdbf.fieldNrs[MinorNr], 1));
    add (fieldValues, OrientNr,
         string2real (values, sdbf.fieldNrs[OrientNr], 1));
  }
  // Add the source.
  // Do not check for duplicates yet.
  if (srcName.empty()) {
    if (patch.empty())
			throw std::runtime_error("Source and/or patch name must be filled in");
    if (matchSearchInfo (ra, dec, searchInfo)) {
      uint patchId = pdb.addPatch (patch, cat, fluxI, ra, dec, check);
      nrpatchfnd++;
      // Create an entry to collect the ra/dec/flux of the sources in the patch.
      patchSumInfo.insert (make_pair(patch, PatchSumInfo(patchId)));
    }
    nrpatch++;
  } else {
    if (matchSearchInfo (ra, dec, searchInfo)) {
      if (patch.empty()) {
        // Patch name is source name plus possible prefix and suffix.
        pdb.addSource (srcInfo, prefix + srcInfo.getName() + suffix,
                       cat, fluxI, fieldValues, ra, dec, check);
      } else {
        pdb.addSource (srcInfo, patch, fieldValues, ra, dec, check);
        // Add ra/dec/flux to patch sum info.
        map<string,PatchSumInfo>::iterator iter = patchSumInfo.find(patch);
        if(iter==patchSumInfo.end()) throw std::runtime_error("Patch name " + patch +
               " not defined before source using it");
        iter->second.add (ra, dec, fluxI);
      }
      nrsourcefnd++;
    }
    nrsource++;
  }
}

void make (const string& in, const string& out, const string& outType,
           const string& format, const string& prefix, const string& suffix,
           bool append, bool average, bool check, const SearchInfo& searchInfo)
{
  // Analyze the format string.
  SdbFormat sdbf = getFormat (format);
  // Create/open the sourcedb and lock it for write.
  ParmDBMeta ptm(outType, out);
  SourceDB pdb(ptm, !append);
  pdb.lock (true);
  int nrpatch     = 0;
  int nrsource    = 0;
  int nrpatchfnd  = 0;
  int nrsourcefnd = 0;
  map<string,PatchSumInfo> patchSumInfo;
  if (! in.empty()) {
    ifstream infile(in.c_str());
    if (!infile)
			throw std::runtime_error("File " + in + " could not be opened");
    casacore::Regex regexf("^[ \t]*[fF][oO][rR][mM][aA][tT][ \t]*=.*");
    string line;
    // Read first line.
    getInLine (infile, line);
    while (infile) {
      // Remove comment lines, empty lines, and possible format line.
      bool skip = true;
      for (uint i=0; i<line.size(); ++i) {
        if (line[i] == '#') {
          break;
        }
        if (line[i] != ' '  &&  line[i] != '\t') {
          if (line[i] == 'f'  ||  line[i] == 'F') {
            String sline(line);
            if (sline.matches (regexf)) {
              break;
            }
          }
          // Empty nor format line, thus use it.
          skip = false;
          break;
        }
      }
      if (!skip) {
        process (line, pdb, sdbf, prefix, suffix, check, nrpatch, nrsource,
                 nrpatchfnd, nrsourcefnd, patchSumInfo, searchInfo);
      }
      // Read next line
      getInLine (infile, line);
    }
  }
  // Write the calculated ra/dec/flux of the patches.
  if (average) {
    for (map<string,PatchSumInfo>::const_iterator iter=patchSumInfo.begin();
         iter!=patchSumInfo.end(); ++iter) {
      const PatchSumInfo& info = iter->second;
      if (info.getFlux() != 0) {
        pdb.updatePatch (info.getPatchId(), info.getFlux(),
                         info.getRa(), info.getDec());
      }
    }
  }
  cout << "Wrote " << nrpatchfnd << " patches (out of " << nrpatch << ") and "
       << nrsourcefnd << " sources (out of " << nrsource << ") into "
       << pdb.getParmDBMeta().getTableName() << endl;
  vector<string> dp(pdb.findDuplicatePatches());
  if (dp.size() > 0) {
    cerr << "Duplicate patches: " << dp << endl;
  }
  vector<string> ds(pdb.findDuplicateSources());
  if (ds.size() > 0) {
    cerr << "Duplicate sources: " << ds << endl;
  }
}

// Read the format from the file.
// It should be contained in a line like # format = .
string readFormat (string file, const string& catFile)
{
  // Use catalog itself if needed.
  if (file.empty()) {
    file = catFile;
  }
  if (file.empty()) {
    return string();
  }
  // Read file until format line is found or until non-comment is found.
  ifstream infile(file.c_str());
  if(!infile) throw std::runtime_error("File " + file
             + " containing format string could not be opened");
  string line;
  getInLine (infile, line);
  casacore::Regex regex("^[ \t]*#[ \t]*\\([ \t]*.*\\)[ \t]*=[ \t]*[fF][oO][rR][mM][aA][tT][ \t]*$");
  casacore::Regex regexs1("^[ \t]*#[ \t]*\\([ \t]*");
  casacore::Regex regexs2("\\)[ \t]*=[ \t]*[fF][oO][rR][mM][aA][tT][ \t]*$");
  while (infile) {
    uInt st=0;
    st = lskipws(line, st, line.size());   // skip whitespace
    if (st < line.size()) {
      if (line[st] != '#') {
        break;                             // data line
      }
      casacore::String sline(line);
      if (sline.matches (regex)) {
        sline.gsub (regexs1, String());
        sline.gsub (regexs2, String());
        return sline;
      }
    }
    getInLine (infile, line);
  }
  // See if a format line is given as "format=".
  casacore::Regex regexf("^[ \t]*[fF][oO][rR][mM][aA][tT][ \t]*=.*$");
  casacore::Regex regexf1("^[ \t]*[fF][oO][rR][mM][aA][tT][ \t]*=[ \t]*");
  casacore::String sline(line);
  if (sline.matches (regexf)) {
    sline.gsub (regexf1, String());
    return sline;
  }
  // No format line found, so use default.
  return string();
}

int main (int argc, char* argv[])
{
  try {
    // Get the inputs.
    Input inputs(1);
    inputs.version ("GvD 2013-May-16");
    inputs.create("in", "",
                  "Input file name", "string");
    inputs.create("out", "",
                  "Output sourcedb name", "string");
    inputs.create ("outtype", "casa",
                   "Output type (casa or blob)", "string");
    inputs.create("format", "<",
                  "Format of the input lines or name of file containing format",
                  "string");
    inputs.create("append", "true",
                  "Append to possibly existing sourcedb?", "bool");
    inputs.create("average", "true",
                  "Calculate average patch ra/dec and total flux?", "bool");
    inputs.create("patchprefix", "",
                  "Add this prefix to patch name if taken from source name",
                  "string");
    inputs.create("patchsuffix", "",
                  "Add this suffix to patch name if taken from source name",
                  "string");
    inputs.create("check", "false",
                  "Check immediately for duplicate entries?", "bool");
    inputs.create("center", "",
                  "Field center as ra,dec", "string");
    inputs.create("radius", "",
                  "Cone search radius", "string");
    inputs.create("width", "",
                  "Search box width as 1 or 2 values (e.g. 2deg,3deg)",
                  "string");
    ///    inputs.create("minflux", "",
    ///                  "Only select sources >= minflux Jy", "string");
    ///    inputs.create("maxflux", "",
    ///                  "Always select sources >= maxflux Jy", "string");
    ///    inputs.create("beammodel", "",
    ///                  "If given, apply beammodel to make fluxes apparent",
    ///                  "WSRT or LOFAR");
    inputs.readArguments(argc, argv);
    string in = inputs.getString("in");
    string out = inputs.getString("out");
    if (out.empty())
			throw std::runtime_error("no output sourcedb name given");
    string outType = boost::to_lower_copy(inputs.getString("outtype"));
    string format = inputs.getString("format");
    string prefix = inputs.getString("patchprefix");
    string suffix = inputs.getString("patchsuffix");
    bool append   = inputs.getBool("append");
    bool average  = inputs.getBool("average");
    bool check    = inputs.getBool("check");
    string center = inputs.getString ("center");
    string radius = inputs.getString ("radius");
    string width  = inputs.getString ("width");
///    double minFlux = inputs.getDouble ("minflux");
///    double maxFlux = inputs.getDouble ("maxflux");
///    string beamModel = inputs.getString ("beammodel");
    // Check if the format has to be read from a file.
    // It is if it starts with a <. The filename should follow it. An empty
    // filename means reading from the catalog file itself.
    if (! format.empty()  &&  format[0] == '<') {
      // Skip optional whitespace.
      uInt st=1;
      st = lskipws (format, st, format.size());
      // Read format from file.
      format = readFormat (format.substr(st), in);
    }
    // Use default if empty format.
    if (format.empty()) {
      cerr << "No format string found; using default format" << endl;
      format = "Name,Type,Ra,Dec,I,Q,U,V,MajorAxis,MinorAxis,Orientation";
    }
    make (in, out, outType, format, prefix, suffix, append, average, check,
          getSearchInfo (center, radius, width));
  } catch (AipsError& x) {
    cerr << "Caught AIPS error: " << x.what() << endl;
    return 1;
  } catch (std::exception& x) {
    cerr << "Caught exception: " << x.what() << endl;
    return 1;
  }
  
  return 0;
}
