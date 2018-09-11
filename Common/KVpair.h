//# KVpair.h: Implements a KV pair as a pair<string, string>.
//#
//# Copyright (C) 2002-2004
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
//# $Id: KVpair.h 31468 2015-04-13 23:26:52Z amesfoort $

#ifndef LOFAR_COMMON_KVPAIR_H
#define LOFAR_COMMON_KVPAIR_H

// \file
// Implements a KV pair as a pair<string, string>.

//# Never #include <config.h> or #include <lofar_config.h> in a header file!
//# Includes
#include <ctime>
#include <utility>
#include <iosfwd>
#include <string>

namespace DP3 {
// \addtogroup Common
// @{

// Implements a KV pair as a pair<string, string>.
class KVpair : public std::pair<std::string, std::string>
{
public:
	// Note: while this class is not PVSS specific, it is mostly (only?)
	// used by the PVSSGateway, which uses valueType to map the enum values
	// below to PVSS types to query (write).
	// If you add types at all without PVSS support, document that below,
	// so that PVSS users can avoid them.
	KVpair(const std::string& aKey, const std::string& aValue, bool genTimestamp = false, bool timestampInKeyname = false);
	KVpair(const std::string& aKey, const char*   aValue, bool genTimestamp = false, bool timestampInKeyname = false);
	KVpair(const std::string& aKey, bool			 aValue, bool genTimestamp = false, bool timestampInKeyname = false);
	KVpair(const std::string& aKey, int			 aValue, bool genTimestamp = false, bool timestampInKeyname = false);
	KVpair(const std::string& aKey, double		 aValue, bool genTimestamp = false, bool timestampInKeyname = false);
	KVpair(const std::string& aKey, float		 aValue, bool genTimestamp = false, bool timestampInKeyname = false);
	KVpair(const std::string& aKey, time_t		 aValue, bool genTimestamp = false, bool timestampInKeyname = false);

	KVpair();
	~KVpair();

	// Copying is allowed
	KVpair(const KVpair&	that);
	KVpair& operator=(const KVpair& that);
	inline bool operator==(const KVpair& that) const { 
		return (first==that.first && second==that.second && timestamp==that.timestamp && valueType==that.valueType); 
	}

	// data-members
	double	timestamp;	// store also as double
	int16_t	valueType;	

	enum {
		VT_UNKNOWN = 0, VT_STRING, VT_BOOL, VT_INT, VT_DOUBLE, VT_FLOAT, VT_TIME_T
	};
};

// @} addgroup

std::ostream& operator<< (std::ostream& os, const DP3::KVpair& kv);


} // namespace LOFAR

#endif
