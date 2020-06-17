// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.

#include "FitsIOChecker.h"

#include <fitsio.h>

#include <sstream>
#include <stdexcept>

void FitsIOChecker::checkStatus(int status, const std::string& filename)
{
  if(status) {
    /* fits_get_errstatus returns at most 30 characters */
    char err_text[31];
    fits_get_errstatus(status, err_text);
    char err_msg[81];
    std::stringstream errMsg;
    errMsg << "CFITSIO reported error when performing IO on file '" << filename << "':" << err_text << " (";
    while(fits_read_errmsg(err_msg))
      errMsg << err_msg;
    errMsg << ')';
    throw std::runtime_error(errMsg.str());
  }
}

void FitsIOChecker::checkStatus(int status, const std::string& filename, const std::string& operation) 
{
  if(status) {
    /* fits_get_errstatus returns at most 30 characters */
    char err_text[31];
    fits_get_errstatus(status, err_text);
    char err_msg[81];
    std::stringstream errMsg;
    errMsg << "During operation " << operation << ", CFITSIO reported error when performing IO on file '" << filename << "': " << err_text << " (";
    while(fits_read_errmsg(err_msg))
      errMsg << err_msg;
    errMsg << ')';
    throw std::runtime_error(errMsg.str());
  }
}

