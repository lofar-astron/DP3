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

