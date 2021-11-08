// SkymodelToSourceDB.cc
//
// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SkymodelToSourceDB.h"

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>

using std::istream;

using dp3::parmdb::ParmDBMeta;
using dp3::parmdb::ParmMap;
using dp3::parmdb::ParmValue;
using dp3::parmdb::ParmValueSet;
using dp3::parmdb::PatchSumInfo;
using dp3::parmdb::SourceDB;
using dp3::parmdb::SourceInfo;

using casacore::IPosition;

// Define the sequence nrs of the various fields.
enum FieldNr {
  // First the standard fields.
  NameNr,
  TypeNr,
  RefTypeNr,
  RaNr,
  DecNr,
  INr,
  QNr,
  UNr,
  VNr,
  SpInxNr,
  LogSINr,
  RefFreqNr,
  MajorNr,
  MinorNr,
  OrientNr,
  RotMeasNr,
  PolFracNr,
  PolAngNr,
  RefWavelNr,
  IShapeletNr,
  QShapeletNr,
  UShapeletNr,
  VShapeletNr,
  NrKnownFields,
  // Now other fields
  CatNr = NrKnownFields,
  PatchNr,
  RahNr,
  RadNr,
  RamNr,
  RasNr,
  DechNr,
  DecdNr,
  DecmNr,
  DecsNr,
  // Nr of fields.
  NrFields
};

namespace dp3 {
namespace parmdb {
namespace skymodel_to_source_db {

std::vector<std::string> fillKnown() {
  // The order in which the names are pushed must match exactly with the
  // enum above.
  std::vector<std::string> names;
  names.reserve(NrFields);
  names.push_back("Name");
  names.push_back("Type");
  names.push_back("RefType");
  names.push_back("Ra");
  names.push_back("Dec");
  names.push_back("I");
  names.push_back("Q");
  names.push_back("U");
  names.push_back("V");
  names.push_back("SpectralIndex");
  names.push_back("LogarithmicSI");
  names.push_back("ReferenceFrequency");
  names.push_back("MajorAxis");
  names.push_back("MinorAxis");
  names.push_back("Orientation");
  names.push_back("RotationMeasure");
  names.push_back("PolarizedFraction");
  names.push_back("PolarizationAngle");
  names.push_back("ReferenceWavelength");
  names.push_back("IShapelet");
  names.push_back("QShapelet");
  names.push_back("UShapelet");
  names.push_back("VShapelet");
  names.push_back("Category");
  names.push_back("Patch");
  names.push_back("rah");
  names.push_back("rad");
  names.push_back("ram");
  names.push_back("ras");
  names.push_back("dech");
  names.push_back("decd");
  names.push_back("decm");
  names.push_back("decs");
  assert(names.size() == NrFields);
  return names;
}
// Define all field names.
std::vector<std::string> theFieldNames = fillKnown();

enum FieldType {
  KNOWNFIELD = 1,
  FIXEDVALUE = 2,
  SKIPFIELD = 4,
  DEFAULTVALUE = 8
};

struct SdbFormat {
  // Define fieldnr of the known fields.
  // -1 means not given.
  std::vector<int> fieldNrs;
  // Define the separators.
  std::vector<char> sep;
  // Define the field names.
  std::vector<std::string> names;
  // Define the field type.
  std::vector<int> types;
  // Define the fixed values.
  std::vector<std::string> values;
};

// Read a line and remove a possible carriage-return at the end.
void getInLine(istream& infile, string& line) {
  getline(infile, line);
  int sz = line.size();
  if (sz > 0 && line[sz - 1] == '\r') {
    line = line.substr(0, sz - 1);
  }
}

void checkRaDec(const SdbFormat& sdbf, int nr, int hnr, int dnr, int mnr,
                int snr) {
  int v = sdbf.fieldNrs[nr];
  int hv = sdbf.fieldNrs[hnr];
  int dv = sdbf.fieldNrs[dnr];
  int mv = sdbf.fieldNrs[mnr];
  int sv = sdbf.fieldNrs[snr];
  if (hv >= 0 && dv >= 0)
    throw std::runtime_error("rah and rad cannot be used both (same for dec)");
  if (v >= 0) {
    if (!(hv < 0 && dv < 0 && mv < 0 && sv < 0))
      throw std::runtime_error(
          "rah/d/m/s cannot be used if ra is given (same for dec)");
  }
  if (!(v >= 0 || hv >= 0 || dv >= 0 || mv >= 0 || sv >= 0))
    throw std::runtime_error("No Ra or Dec info given");
}

unsigned int ltrim(const string& value, unsigned int st, unsigned int end) {
  for (; st < end; ++st) {
    if (value[st] != ' ' && value[st] != '\t') {
      break;
    }
  }
  return st;
}

unsigned int rtrim(const string& value, unsigned int st, unsigned int end) {
  for (; end > st; --end) {
    if (value[end - 1] != ' ' && value[end - 1] != '\t') {
      break;
    }
  }
  return end;
}

// Get the next value by looking for the separator.
// The separator is ignored in parts enclosed in quotes or square brackets.
// Square brackets can be nested (they indicate arrays).
unsigned int nextValue(const string& str, char sep, unsigned int st,
                       unsigned int end) {
  unsigned int posbracket = 0;
  int nbracket = 0;
  while (st < end) {
    if (str[st] == '\'' || str[st] == '"') {
      // Ignore a quoted part.
      string::size_type pos = str.find(str[st], st + 1);
      if (pos == string::npos)
        throw std::runtime_error("Unbalanced quoted string at position " +
                                 std::to_string(st) + " in " + str);
      st = pos;
    } else if (str[st] == '[') {
      nbracket++;
      posbracket = st;
    } else if (str[st] == ']') {
      if (nbracket <= 0)
        throw std::runtime_error("Unbalanced square brackets at position " +
                                 std::to_string(st) + " in " + str);
      nbracket--;
    } else if (nbracket == 0) {
      if (str[st] == sep) {
        return st;
      }
    }
    ++st;
  }
  if (nbracket != 0)
    throw std::runtime_error("Unbalanced square brackets at position " +
                             std::to_string(posbracket) + " in " + str);
  return end;
}

SdbFormat getFormat(const string& format) {
  // Fill the map with known names.
  std::map<std::string, int> nameMap;
  for (unsigned int i = 0; i < theFieldNames.size(); ++i) {
    nameMap[boost::to_lower_copy(theFieldNames[i])] = i;
  }
  // Skip possible left and right whitespace.
  unsigned int end = format.size();
  unsigned int st = ltrim(format, 0, end);
  end = rtrim(format, st, end);
  // Initialize the format.
  SdbFormat sdbf;
  sdbf.fieldNrs.resize(NrFields);
  for (unsigned int i = 0; i < NrFields; ++i) {
    sdbf.fieldNrs[i] = -1;
  }
  // Use default if the string is empty.
  if (st >= end - 1) {
    char sep = ',';
    if (st == end - 1) sep = format[st];
    for (unsigned int i = 0; i < NrKnownFields; ++i) {
      sdbf.fieldNrs[i] = i;
      sdbf.sep.push_back(sep);
      sdbf.names.push_back(theFieldNames[i]);
      sdbf.types.push_back(KNOWNFIELD);
      sdbf.values.push_back("");
    }
    return sdbf;
  }
  // Parse the format string.
  unsigned int nr = 0;
  unsigned int i = st;
  while (i < end) {
    if ((format[i] >= 'a' && format[i] <= 'z') ||
        (format[i] >= 'A' && format[i] <= 'Z') ||
        (format[i] >= '0' && format[i] <= '9') || format[i] == '_' ||
        (format[i] == ':' && i > st)) {
      ++i;  // part of name
    } else {
      // End of name
      string name = format.substr(st, i - st);
      string lname = boost::to_lower_copy(name);
      int fieldType = 0;
      if (lname.empty() || lname == "dummy") {
        fieldType = SKIPFIELD;
      } else {
        /// Remove ASSERT branch once we're sure SpectralIndexDegree is not used
        /// anymore.
        if (lname == "spectralindexdegree")
          throw std::runtime_error(
              "Use SpectralIndex=[v1,v2,...] instead of SpectralIndexDegree "
              "and SpectralIndex:i");
        std::map<std::string, int>::const_iterator namepos =
            nameMap.find(lname);
        // Fill in fieldnr of a known field.
        if (namepos != nameMap.end()) {
          fieldType = KNOWNFIELD;
          sdbf.fieldNrs[namepos->second] = nr;
        }
      }
      // See if a default or fixed value is given.
      std::string fixedValue;
      i = ltrim(format, i, end);
      if (i < end && format[i] == '=') {
        i = ltrim(format, i + 1, end);
        // See if it is a fixed value.
        bool isDefault = true;
        if (i + 5 < end &&
            boost::to_lower_copy(format.substr(i, 5)) == "fixed") {
          isDefault = false;
          i += 5;
        }
        if (i >= end)
          throw std::runtime_error("No value given after " + name + '=');
        if (format[i] != '"' && format[i] != '\'')
          throw std::runtime_error("value after " + name +
                                   "= must start with a quote");
        // Find the ending quote.
        std::string::size_type pos = format.find(format[i], i + 1);
        if (pos == std::string::npos)
          throw std::runtime_error("No closing value quote given after " +
                                   name + '=');
        fixedValue = format.substr(i + 1, pos - i - 1);
        i = ltrim(format, pos + 1, end);
        if (isDefault) {
          fieldType |= DEFAULTVALUE;
        } else {
          fieldType |= FIXEDVALUE;
        }
      }
      // Now look for a separator.
      char sep = ' ';
      if (!((format[i] >= 'a' && format[i] <= 'z') ||
            (format[i] >= 'A' && format[i] <= 'Z') ||
            (format[i] >= '0' && format[i] <= '9') || format[i] == '_')) {
        sep = format[i];
        if (sep == '"' || sep == '\'')
          throw std::runtime_error(
              "A quote is found as separator; "
              "probably a quote around a value in the format string is "
              "missing");
        i = ltrim(format, i + 1, end);
      }
      sdbf.sep.push_back(sep);
      sdbf.names.push_back(name);
      sdbf.types.push_back(fieldType);
      sdbf.values.push_back(fixedValue);
      nr++;
      st = i;
    }
  }
  if (st < end) {
    // The last item was just a name which has to be processed.
    std::string name = format.substr(st, end - st);
    std::string lname = boost::to_lower_copy(name);
    int fieldType = 0;
    if (lname.empty() || lname == "dummy") {
      fieldType = SKIPFIELD;
    } else {
      std::map<std::string, int>::const_iterator namepos = nameMap.find(lname);
      // Fill in fieldnr of a known field.
      if (namepos != nameMap.end()) {
        sdbf.fieldNrs[namepos->second] = nr;
        fieldType = KNOWNFIELD;
      }
    }
    sdbf.sep.push_back(' ');
    sdbf.names.push_back(name);
    sdbf.types.push_back(fieldType);
    sdbf.values.push_back("");
  }
  // Make sure Ra and Dec are given correctly.
  checkRaDec(sdbf, RaNr, RahNr, RadNr, RamNr, RasNr);
  checkRaDec(sdbf, DecNr, DechNr, DecdNr, DecmNr, DecsNr);
  return sdbf;
}

SourceInfo::Type string2type(const std::string& str) {
  std::string s = boost::to_lower_copy(str);
  if (s == "point" || s.empty()) {
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

std::string unquote(const std::string& value) {
  std::string res(value);
  if (res.size() > 1) {
    int last = res.size() - 1;
    if (last >= 1 && ((res[0] == '"' && res[last] == '"') ||
                      (res[0] == '\'' && res[last] == '\''))) {
      res = res.substr(1, last - 1);
    }
  }
  return res;
}

std::string getValue(const std::vector<std::string>& values, int nr,
                     const std::string& defVal = std::string()) {
  if (nr < 0) {
    return defVal;
  }
  return unquote(values[nr]);
}

bool string2bool(const std::vector<std::string>& values, int nr, bool defVal) {
  std::string value = getValue(values, nr);
  if (value.empty()) {
    return defVal;
  }
  return dp3::common::strToBool(value);
}

int string2int(const std::vector<std::string>& values, int nr, int defVal) {
  std::string value = getValue(values, nr);
  if (value.empty()) {
    return defVal;
  }
  return dp3::common::strToInt(value);
}

double string2real(const std::string& value, double defVal) {
  if (value.empty()) {
    return defVal;
  }
  return dp3::common::strToDouble(value);
}

double string2real(const std::vector<std::string>& values, int nr,
                   double defVal) {
  return string2real(getValue(values, nr), defVal);
}

// Convert values in a possibly bracketed string to a vector of strings
// taking quoted or bracketed values into account.
// A comma isused as separator.
std::vector<std::string> string2vector(const std::string& value,
                                       const std::vector<std::string>& defVal) {
  std::vector<std::string> result;
  // Test if anything is given.
  if (value.empty()) {
    result = defVal;
  } else {
    unsigned int end = value.size();
    // If no brackets given, it is a single value.
    if (value.size() < 2 || value[0] != '[' || value[end - 1] != ']') {
      result.push_back(value);
    } else {
      // Skip opening and closing bracket and possible whitespace.
      unsigned int st = ltrim(value, 1, end - 1);
      end = rtrim(value, st, end - 1);
      while (st < end) {
        unsigned int pos = nextValue(value, ',', st, end);
        result.push_back(value.substr(st, rtrim(value, st, pos) - st));
        st = ltrim(value, pos + 1, end);
      }
    }
  }
  return result;
}

std::vector<std::string> string2vector(const std::vector<std::string>& values,
                                       int nr,
                                       const std::vector<std::string>& defVal) {
  return string2vector(getValue(values, nr), defVal);
}

std::vector<double> vector2real(const std::vector<std::string>& values,
                                double defVal) {
  std::vector<double> result;
  result.reserve(values.size());
  for (unsigned int i = 0; i < values.size(); ++i) {
    result.push_back(string2real(values[i], defVal));
  }
  return result;
}

double string2pos(const std::vector<std::string>& values, int pnr, int hnr,
                  int dnr, int mnr, int snr, bool canUseColon) {
  double deg = 0;
  bool fnd = false;
  if (pnr >= 0) {
    std::string value = getValue(values, pnr);
    if (!value.empty()) {
      if (!canUseColon) {
        if (value.find(':') != std::string::npos)
          throw std::runtime_error(
              "Colons cannot be used in declination value " + value);
      }
      casacore::Quantity q;
      if (!casacore::MVAngle::read(q, values[pnr]))
        throw std::runtime_error("Error in reading position " + values[pnr]);
      deg = q.getValue("deg");
      fnd = true;
    }
  } else {
    if (hnr >= 0) {
      std::string value = getValue(values, hnr);
      if (!value.empty()) {
        deg = string2real(values[hnr], 0);
        fnd = true;
      }
    } else if (dnr >= 0) {
      std::string value = getValue(values, dnr);
      if (!value.empty()) {
        deg = string2real(values[dnr], 0);
        fnd = true;
      }
    }
    double ms = 0;
    if (mnr >= 0) {
      std::string value = getValue(values, mnr);
      if (!value.empty()) {
        ms = string2real(values[mnr], 0);
        fnd = true;
      }
    }
    if (snr >= 0) {
      std::string value = getValue(values, snr);
      if (!value.empty()) {
        ms += string2real(values[snr], 0) / 60;
        fnd = true;
      }
    }
    if (deg < 0) {
      deg -= ms / 60;
    } else {
      deg += ms / 60;
    }
    if (hnr >= 0) {
      deg *= 15;
    }
  }
  if (fnd) {
    casacore::Quantity q(deg, "deg");
    return q.getValue("rad");
  }
  return 1e-9;
}

void checkRefType(const std::string& refType) {
  std::string type = boost::to_upper_copy(refType);
  if (type != "J2000" && type != "B1950" && type != "SUN" && type != "MOON" &&
      type != "VENUS" && type != "MARS" && type != "JUPITER" &&
      type != "SATURN" && type != "URANUS" && type != "NEPTUNE" &&
      type != "MERCURY") {
    throw std::runtime_error("Reference type " + refType + " is incorrect");
  }
}

// Get the search cone or box values.
SearchInfo GetSearchInfo(const std::string& center, const std::string& radius,
                         const std::string& width) {
  SearchInfo searchInfo;
  if (center.empty()) {
    searchInfo.search = false;
  } else {
    searchInfo.search = true;
    std::vector<std::string> pos;
    boost::algorithm::split(pos, center, boost::is_any_of(","));
    if (pos.size() != 2)
      throw std::runtime_error("center not specified as ra,dec");
    searchInfo.ra = string2pos(pos, 0, -1, -1, -1, -1, true);
    searchInfo.dec = string2pos(pos, 1, -1, -1, -1, -1, false);
    searchInfo.sinDec = sin(searchInfo.dec);
    searchInfo.cosDec = cos(searchInfo.dec);
    if (radius.empty() == width.empty())
      throw std::runtime_error(
          "radius OR width must be given if center is given (not both)");
    if (radius.empty()) {
      double raw, decw;
      searchInfo.asCone = false;
      pos.clear();
      boost::algorithm::split(pos, width, boost::is_any_of(","));
      if (pos.size() != 1 && pos.size() != 2)
        throw std::runtime_error("width should be specified as 1 or 2 values");
      raw = string2pos(pos, 0, -1, -1, -1, -1, true);
      if (pos.size() > 1) {
        decw = string2pos(pos, 1, -1, -1, -1, -1, false);
      } else {
        decw = raw;
      }
      searchInfo.raStart = searchInfo.ra - raw / 2;
      searchInfo.raEnd = searchInfo.ra + raw / 2;
      searchInfo.decStart = searchInfo.dec - decw / 2;
      searchInfo.decEnd = searchInfo.dec + decw / 2;
    } else {
      searchInfo.asCone = true;
      pos[0] = radius;
      searchInfo.cosRadius = cos(string2pos(pos, 0, -1, -1, -1, -1, false));
    }
  }
  return searchInfo;
}

bool matchSearchInfo(double ra, double dec, const SearchInfo& si) {
  if (!si.search) {
    return true;
  }
  bool match = false;
  if (si.asCone) {
    match = (si.cosRadius <=
             si.sinDec * sin(dec) + si.cosDec * cos(dec) * cos(si.ra - ra));
  } else {
    // Ra can be around 0 or 360 degrees, so make sure all cases are handled.
    ra -= casacore::C::_2pi;
    for (int i = 0; i < 4; ++i) {
      if (ra >= si.raStart && ra <= si.raEnd) {
        match = true;
        break;
      }
      ra += casacore::C::_2pi;
    }
    if (match) {
      match = (dec >= si.decStart && dec <= si.decEnd);
    }
  }
  return match;
}

void addValue(ParmMap& fieldValues, const std::string& name, double value) {
  fieldValues.define(name, ParmValueSet(ParmValue(value)));
}

void add(ParmMap& fieldValues, FieldNr field, double value) {
  fieldValues.define(theFieldNames[field], ParmValueSet(ParmValue(value)));
}

void addSpInx(ParmMap& fieldValues, const std::vector<double>& spinx,
              double refFreq) {
  if (spinx.size() > 0) {
    if (refFreq <= 0)
      throw std::runtime_error(
          "SpectralIndex given, but no ReferenceFrequency");
    /// Remove the following lines if not needed anymore for BBS.
    addValue(fieldValues, "SpectralIndexDegree", int(spinx.size()) - 1);
    for (unsigned int i = 0; i < spinx.size(); ++i) {
      std::ostringstream ostr;
      ostr << "SpectralIndex:" << i;
      addValue(fieldValues, ostr.str(), spinx[i]);
    }
  }
}

void readShapelet(const std::string& fileName, casacore::Array<double>& coeff,
                  double& scale) {
  std::ifstream file(fileName.c_str());
  if (!file)
    throw std::runtime_error("Shapelet file " + fileName +
                             " could not be opened");
  std::string line;
  getInLine(file, line);  // ra dec
  getInLine(file, line);  // order scale
  std::vector<std::string> parts;
  boost::algorithm::split(parts, line, boost::is_any_of(" "));
  if (parts.size() != 2)
    throw std::runtime_error("Expected 2 values in shapelet line " + line);
  int order = string2int(parts, 0, 0);
  scale = string2real(parts, 1, 0.);
  if (order <= 0)
    throw std::runtime_error("Invalid order in shapelet line " + line);
  coeff.resize(IPosition(2, order, order));
  double* coeffData = coeff.data();
  for (unsigned int i = 0; i < coeff.size(); ++i) {
    getInLine(file, line);  // index coeff
    std::vector<std::string> parts;
    boost::algorithm::split(parts, line, boost::is_any_of(" "));
    if (parts.size() != 2)
      throw std::runtime_error("Expected 2 values in shapelet line " + line);
    if (string2int(parts, 0, -1) != int(i))
      throw std::runtime_error("Expected shapelet line with index " +
                               std::to_string(i));
    *coeffData++ = string2real(parts, 1, 0.);
  }
}

void fillShapelet(SourceInfo& srcInfo, const std::string& shpI,
                  const std::string& shpQ, const std::string& shpU,
                  const std::string& shpV) {
  double scaleI = 0;
  double scaleQ = 0;
  double scaleU = 0;
  double scaleV = 0;
  casacore::Array<double> coeffI, coeffQ, coeffU, coeffV;
  readShapelet(shpI, coeffI, scaleI);
  if (shpQ.empty()) {
    coeffQ = coeffI;
    scaleQ = scaleI;
  } else {
    readShapelet(shpQ, coeffQ, scaleQ);
  }
  if (shpU.empty()) {
    coeffU = coeffI;
    scaleU = scaleI;
  } else {
    readShapelet(shpU, coeffU, scaleU);
  }
  if (shpV.empty()) {
    coeffV = coeffI;
    scaleV = scaleI;
  } else {
    readShapelet(shpV, coeffV, scaleV);
  }
  srcInfo.setShapeletCoeff(coeffI, coeffQ, coeffU, coeffV);
  srcInfo.setShapeletScale(scaleI, scaleQ, scaleU, scaleV);
}

// Calculate the polarization angle and polarized fraction given Q and U
// for a given reference wavelength.
// A spectral index can be used to calculate Stokes I.
void calcRMParam(double& polfrac, double& polang, double fluxi0, double fluxq,
                 double fluxu, const std::vector<double>& spinx, double rm,
                 double refFreq, double rmRefWavel) {
  // polfrac = sqrt(q^2 + u^2) / i
  // where i = i(0) * spinx
  // Compute spectral index for the RM reference wavelength as:
  // (v / v0) ^ (c0 + c1 * log10(v / v0) + c2 * log10(v / v0)^2 + ...)
  // Where v is the RM frequency and v0 is the spinx reference frequency.
  double si = 1;
  if (spinx.size() > 0) {
    if (rmRefWavel <= 0)
      throw std::runtime_error("No RM reference wavelength given");
    double rmFreq = casacore::C::c / rmRefWavel;
    double vv0 = rmFreq / refFreq;
    double factor = 1;
    double sum = 0;
    for (unsigned int i = 0; i < spinx.size(); ++i) {
      sum += factor * spinx[i];
      factor *= log10(vv0);
    }
    si = std::pow(vv0, sum);
  }
  double fluxi = fluxi0 * si;
  polfrac = sqrt(fluxq * fluxq + fluxu * fluxu) / fluxi;
  // Calculate polang(0) from Q and U given for the reference lambda.
  // polang = atan2(u,q) / 2
  // polang(lambda) = polang(0) + lambda^2 * rm
  // Scale between 0 and pi.
  double pa = 0.5 * atan2(fluxu, fluxq) - rmRefWavel * rmRefWavel * rm;
  polang = fmod(pa, casacore::C::pi);
  if (polang < 0) {
    polang += casacore::C::pi;
  }
}

void process(const std::string& line, dp3::parmdb::SourceDBBase& pdb,
             const SdbFormat& sdbf, const std::string& prefix,
             const std::string& suffix, bool check, int& nrpatch, int& nrsource,
             int& nrpatchfnd, int& nrsourcefnd,
             std::map<std::string, PatchSumInfo>& patchSumInfo,
             const SearchInfo& searchInfo) {
  //  cout << line << endl;
  // Hold the values.
  ParmMap fieldValues;
  std::vector<std::string> values;
  // Process the line.
  unsigned int end = line.size();
  unsigned int st = ltrim(line, 0, end);
  for (unsigned int i = 0; i < sdbf.names.size(); ++i) {
    std::string value;
    if ((sdbf.types[i] & FIXEDVALUE) == FIXEDVALUE) {
      value = sdbf.values[i];
    } else if (st < end) {
      unsigned int pos = nextValue(line, sdbf.sep[i], st, end);
      value = line.substr(st, rtrim(line, st, pos) - st);
      st = ltrim(line, pos + 1, end);
    }
    if (value.empty() && (sdbf.types[i] & DEFAULTVALUE) == DEFAULTVALUE) {
      value = sdbf.values[i];
    }
    values.push_back(value);
    if ((sdbf.types[i] & SKIPFIELD) != SKIPFIELD) {
      if ((sdbf.types[i] & KNOWNFIELD) != KNOWNFIELD) {
        addValue(fieldValues, sdbf.names[i], string2real(unquote(value), 0));
      }
    }
  }
  // Now handle the standard fields.
  std::string srcName = getValue(values, sdbf.fieldNrs[NameNr]);
  std::string refType = getValue(values, sdbf.fieldNrs[RefTypeNr], "J2000");
  if (refType.empty()) {
    refType = "J2000";
  }
  checkRefType(refType);
  SourceInfo::Type srctype =
      string2type(getValue(values, sdbf.fieldNrs[TypeNr]));
  double ra = string2pos(values, sdbf.fieldNrs[RaNr], sdbf.fieldNrs[RahNr],
                         sdbf.fieldNrs[RadNr], sdbf.fieldNrs[RamNr],
                         sdbf.fieldNrs[RasNr], true);
  double dec = string2pos(values, sdbf.fieldNrs[DecNr], sdbf.fieldNrs[DechNr],
                          sdbf.fieldNrs[DecdNr], sdbf.fieldNrs[DecmNr],
                          sdbf.fieldNrs[DecsNr], false);
  if (!(ra > -6.3 && ra < 6.3 && dec > -1.6 && dec < 1.6))
    throw std::runtime_error("RA " + std::to_string(ra) + " or DEC " +
                             std::to_string(dec) +
                             " radians is outside boundaries");
  int cat = string2int(values, sdbf.fieldNrs[CatNr], 2);
  double fluxI = string2real(getValue(values, sdbf.fieldNrs[INr]), 1.);
  std::string fluxQ = getValue(values, sdbf.fieldNrs[QNr]);
  std::string fluxU = getValue(values, sdbf.fieldNrs[UNr]);
  std::string fluxV = getValue(values, sdbf.fieldNrs[VNr]);
  std::string rm = getValue(values, sdbf.fieldNrs[RotMeasNr]);
  std::string polFrac = getValue(values, sdbf.fieldNrs[PolFracNr]);
  std::string polAng = getValue(values, sdbf.fieldNrs[PolAngNr]);
  std::string refWavel = getValue(values, sdbf.fieldNrs[RefWavelNr]);
  std::string shapeletI = getValue(values, sdbf.fieldNrs[IShapeletNr]);
  std::string shapeletQ = getValue(values, sdbf.fieldNrs[QShapeletNr]);
  std::string shapeletU = getValue(values, sdbf.fieldNrs[UShapeletNr]);
  std::string shapeletV = getValue(values, sdbf.fieldNrs[VShapeletNr]);

  std::vector<double> spinx(vector2real(
      string2vector(values, sdbf.fieldNrs[SpInxNr], std::vector<std::string>()),
      0.));
  double refFreq = string2real(values, sdbf.fieldNrs[RefFreqNr], 0);
  bool useLogSI = string2bool(values, sdbf.fieldNrs[LogSINr], true);
  bool useRM = false;
  double rmRefWavel = 0;
  if (rm.empty()) {
    if (!polFrac.empty() || !polAng.empty() || !refWavel.empty())
      throw std::runtime_error(
          "PolarizationAngle, PolarizedFraction, and ReferenceWavelength"
          " cannot be specified if RotationMeasure is not specified");
  } else {
    if (!fluxQ.empty() || !fluxU.empty()) {
      if (fluxQ.empty() || fluxU.empty() || !polFrac.empty() || !polAng.empty())
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
      rmRefWavel = string2real(refWavel, casacore::C::c / refFreq);
    } else {
      if (polFrac.empty() || polAng.empty())
        throw std::runtime_error(
            "PolarizationAngle/PolarizedFraction or Q/U must be "
            "specified if RotationMeasure is specified");
      useRM = true;
      rmRefWavel = string2real(refWavel, 0);
      if (rmRefWavel != 0)
        throw std::runtime_error(
            "PolarizationAngle/PolarizedFraction can "
            "only be given for ReferenceWavelength=0");
    }
  }
  SourceInfo srcInfo(srcName, srctype, refType, useLogSI, spinx.size(), refFreq,
                     useRM);
  if (srctype == SourceInfo::SHAPELET) {
    fillShapelet(srcInfo, shapeletI, shapeletQ, shapeletU, shapeletV);
  }
  add(fieldValues, INr, fluxI);
  double rmval = 0;
  if (!rm.empty()) {
    rmval = string2real(rm, 0.);
    add(fieldValues, RotMeasNr, rmval);
  }
  double fq = string2real(fluxQ, 0.);
  double fu = string2real(fluxU, 0.);
  if (useRM) {
    double pfrac = string2real(polFrac, 0.);
    double pang = string2real(polAng, 0.);
    if (!fluxQ.empty()) {
      calcRMParam(pfrac, pang, fluxI, fq, fu, spinx, rmval, refFreq,
                  rmRefWavel);
    }
    add(fieldValues, PolFracNr, pfrac);
    add(fieldValues, PolAngNr, pang);
  } else {
    add(fieldValues, QNr, string2real(fluxQ, 0.));
    add(fieldValues, UNr, string2real(fluxU, 0.));
  }
  add(fieldValues, VNr, string2real(fluxV, 0.));
  addSpInx(fieldValues, spinx, refFreq);
  if (refFreq > 0) {
    add(fieldValues, RefFreqNr, refFreq);
  }
  std::string patch = getValue(values, sdbf.fieldNrs[PatchNr]);

  if (srctype == SourceInfo::GAUSSIAN) {
    add(fieldValues, MajorNr, string2real(values, sdbf.fieldNrs[MajorNr], 1));
    add(fieldValues, MinorNr, string2real(values, sdbf.fieldNrs[MinorNr], 1));
    add(fieldValues, OrientNr, string2real(values, sdbf.fieldNrs[OrientNr], 1));
  }
  // Add the source.
  // Do not check for duplicates yet.
  if (srcName.empty()) {
    if (patch.empty())
      throw std::runtime_error("Source and/or patch name must be filled in");
    if (matchSearchInfo(ra, dec, searchInfo)) {
      unsigned int patchId = pdb.addPatch(patch, cat, fluxI, ra, dec, check);
      nrpatchfnd++;
      // Create an entry to collect the ra/dec/flux of the sources in the patch.
      patchSumInfo.insert(make_pair(patch, PatchSumInfo(patchId)));
    }
    nrpatch++;
  } else {
    if (matchSearchInfo(ra, dec, searchInfo)) {
      if (patch.empty()) {
        // Patch name is source name plus possible prefix and suffix.
        pdb.addSource(srcInfo, prefix + srcInfo.getName() + suffix, cat, fluxI,
                      fieldValues, ra, dec, check);
      } else {
        pdb.addSource(srcInfo, patch, fieldValues, ra, dec, check);
        // Add ra/dec/flux to patch sum info.
        std::map<std::string, PatchSumInfo>::iterator iter =
            patchSumInfo.find(patch);
        if (iter == patchSumInfo.end())
          throw std::runtime_error("Patch name " + patch +
                                   " not defined before source using it");
        iter->second.add(ra, dec, fluxI);
      }
      nrsourcefnd++;
    }
    nrsource++;
  }
}

static void ParseSkyModel(dp3::parmdb::SourceDBBase& pdb, std::istream& input,
                          const SdbFormat& sdbf, const std::string& prefix,
                          const std::string& suffix, bool check, int& nrpatch,
                          int& nrsource, int& nrpatchfnd, int& nrsourcefnd,
                          std::map<std::string, PatchSumInfo>& patchSumInfo,
                          const SearchInfo& searchInfo) {
  casacore::Regex regexf("^[ \t]*[fF][oO][rR][mM][aA][tT][ \t]*=.*");
  std::string line;
  // Read first line.
  getInLine(input, line);
  while (input) {
    // Remove comment lines, empty lines, and possible format line.
    bool skip = true;
    for (unsigned int i = 0; i < line.size(); ++i) {
      if (line[i] == '#') {
        break;
      }
      if (line[i] != ' ' && line[i] != '\t') {
        if (line[i] == 'f' || line[i] == 'F') {
          casacore::String sline(line);
          if (sline.matches(regexf)) {
            break;
          }
        }
        // Empty nor format line, thus use it.
        skip = false;
        break;
      }
    }
    if (!skip) {
      process(line, pdb, sdbf, prefix, suffix, check, nrpatch, nrsource,
              nrpatchfnd, nrsourcefnd, patchSumInfo, searchInfo);
    }
    // Read next line
    getInLine(input, line);
  }
}

SourceDB MakeSourceDb(const std::string& in, const std::string& out,
                      const std::string& outType, const std::string& format,
                      const std::string& prefix, const std::string& suffix,
                      bool append, bool average, bool check,
                      const SearchInfo& search_info) {
  SdbFormat sdbf = getFormat(format);
  // Create/open the sourcedb and lock it for write.
  ParmDBMeta ptm(outType, out);
  SourceDB pdb(ptm, false, !append);
  pdb.lock(true);
  int nrpatch = 0;
  int nrsource = 0;
  int nrpatchfnd = 0;
  int nrsourcefnd = 0;
  std::map<std::string, PatchSumInfo> patchSumInfo;
  if (!in.empty()) {
    std::ifstream infile(in.c_str());
    if (!infile)
      throw std::runtime_error("File " + in + " could not be opened");

    ParseSkyModel(pdb, infile, sdbf, prefix, suffix, check, nrpatch, nrsource,
                  nrpatchfnd, nrsourcefnd, patchSumInfo, search_info);
  }
  // Write the calculated ra/dec/flux of the patches.
  if (average) {
    for (std::map<std::string, PatchSumInfo>::const_iterator iter =
             patchSumInfo.begin();
         iter != patchSumInfo.end(); ++iter) {
      const PatchSumInfo& info = iter->second;
      if (info.getFlux() != 0) {
        pdb.updatePatch(info.getPatchId(), info.getFlux(), info.getRa(),
                        info.getDec());
      }
    }
  }
  std::cout << "Wrote " << nrpatchfnd << " patches (out of " << nrpatch
            << ") and " << nrsourcefnd << " sources (out of " << nrsource
            << ") into " << pdb.getParmDBMeta().getTableName() << '\n';
  casacore::Vector<std::string> dp(pdb.findDuplicatePatches());
  if (dp.size() > 0) {
    std::cerr << "Duplicate patches: " << dp << '\n';
  }
  casacore::Vector<std::string> ds(pdb.findDuplicateSources());
  if (ds.size() > 0) {
    std::cerr << "Duplicate sources: " << ds << '\n';
  }

  return pdb;
}

SourceDBSkymodel MakeSourceDBSkymodel(const std::string& filename,
                                      const std::string& format) {
  SdbFormat sdb_format = getFormat(format);
  SourceDBSkymodel source_db;
  if (filename.empty()) {
    return source_db;
  }
  int nrpatch = 0;
  int nrsource = 0;
  int nrpatchfnd = 0;
  int nrsourcefnd = 0;
  std::map<std::string, PatchSumInfo> patchSumInfo;

  std::ifstream file(filename);
  if (!file)
    throw std::runtime_error("File " + filename +
                             " could not be opened for reading.");

  ParseSkyModel(source_db, file, sdb_format, /*prefix*/ "", /*suffix*/ "",
                /*check*/ false, nrpatch, nrsource, nrpatchfnd, nrsourcefnd,
                patchSumInfo, GetSearchInfo("", "", ""));

  // Write the calculated ra/dec/flux of the patches.
  for (const auto& patch : patchSumInfo) {
    const auto& info = patch.second;
    if (info.getFlux() != 0) {
      source_db.updatePatch(info.getPatchId(), info.getFlux(), info.getRa(),
                            info.getDec());
    }
  }

  return source_db;
}

// Read the format from the file.
// It should be contained in a line like # format = .
std::string ReadFormat(std::string file, const std::string& cat_file) {
  // Use catalog itself if needed.
  if (file.empty()) {
    file = cat_file;
  }
  if (file.empty()) {
    return std::string();
  }
  // Read file until format line is found or until non-comment is found.
  std::ifstream infile(file.c_str());
  if (!infile)
    throw std::runtime_error("File " + file +
                             " containing format string could not be opened");
  std::string line;
  getInLine(infile, line);
  casacore::Regex regex(
      "^[ \t]*#[ \t]*\\([ \t]*.*\\)[ \t]*=[ \t]*[fF][oO][rR][mM][aA][tT][ "
      "\t]*$");
  casacore::Regex regexs1("^[ \t]*#[ \t]*\\([ \t]*");
  casacore::Regex regexs2("\\)[ \t]*=[ \t]*[fF][oO][rR][mM][aA][tT][ \t]*$");
  while (infile) {
    unsigned st = 0;
    st = dp3::common::lskipws(line, st, line.size());  // skip whitespace
    if (st < line.size()) {
      if (line[st] != '#') {
        break;  // data line
      }
      casacore::String sline(line);
      if (sline.matches(regex)) {
        sline.gsub(regexs1, casacore::String());
        sline.gsub(regexs2, casacore::String());
        return sline;
      }
    }
    getInLine(infile, line);
  }
  // See if a format line is given as "format=".
  casacore::Regex regexf("^[ \t]*[fF][oO][rR][mM][aA][tT][ \t]*=.*$");
  casacore::Regex regexf1("^[ \t]*[fF][oO][rR][mM][aA][tT][ \t]*=[ \t]*");
  casacore::String sline(line);
  if (sline.matches(regexf)) {
    sline.gsub(regexf1, casacore::String());
    return sline;
  }
  // No format line found, so use default.
  std::cerr << "No format string found; using default format\n";
  return "Name,Type,Ra,Dec,I,Q,U,V,MajorAxis,MinorAxis,Orientation";
}
}  // namespace skymodel_to_source_db
}  // namespace parmdb
}  // namespace dp3
