// makesourcedb.cc: Fill a SourceDB from an ASCII file
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

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

#include "SkymodelToSourceDB.h"

int main(int argc, char* argv[]) {
  try {
    // Get the inputs.
    casacore::Input inputs(1);
    inputs.version("GvD 2013-May-16");
    inputs.create("in", "", "Input file name", "string");
    inputs.create("out", "", "Output sourcedb name", "string");
    inputs.create("outtype", "casa", "Output type (casa or blob)", "string");
    inputs.create("format", "<",
                  "Format of the input lines or name of file containing format",
                  "string");
    inputs.create("append", "true", "Append to possibly existing sourcedb?",
                  "bool");
    inputs.create("average", "true",
                  "Calculate average patch ra/dec and total flux?", "bool");
    inputs.create("patchprefix", "",
                  "Add this prefix to patch name if taken from source name",
                  "string");
    inputs.create("patchsuffix", "",
                  "Add this suffix to patch name if taken from source name",
                  "string");
    inputs.create("check", "false", "Check immediately for duplicate entries?",
                  "bool");
    inputs.create("center", "", "Field center as ra,dec", "string");
    inputs.create("radius", "", "Cone search radius", "string");
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
    if (out.empty()) throw std::runtime_error("no output sourcedb name given");
    string outType = boost::to_lower_copy(inputs.getString("outtype"));
    string format = inputs.getString("format");
    string prefix = inputs.getString("patchprefix");
    string suffix = inputs.getString("patchsuffix");
    bool append = inputs.getBool("append");
    bool average = inputs.getBool("average");
    bool check = inputs.getBool("check");
    string center = inputs.getString("center");
    string radius = inputs.getString("radius");
    string width = inputs.getString("width");
    ///    double minFlux = inputs.getDouble ("minflux");
    ///    double maxFlux = inputs.getDouble ("maxflux");
    ///    string beamModel = inputs.getString ("beammodel");
    // Check if the format has to be read from a file.
    // It is if it starts with a <. The filename should follow it. An empty
    // filename means reading from the catalog file itself.

    if (!format.empty() && format[0] == '<') {
      // Skip optional whitespace.
      unsigned st = dp3::common::lskipws(format, 1, format.size());
      // Read format from file.
      format =
          dp3::parmdb::skymodel_to_source_db::ReadFormat(format.substr(st), in);
    }

    dp3::parmdb::skymodel_to_source_db::MakeSourceDb(
        in, out, outType, format, prefix, suffix, append, average, check,
        dp3::parmdb::skymodel_to_source_db::GetSearchInfo(center, radius,
                                                          width));
  } catch (std::exception& x) {
    std::cerr << "Caught exception: " << x.what() << '\n';
    return 1;
  }

  return 0;
}
