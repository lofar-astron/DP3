#ifndef FITS_IO_CHECKER_H
#define FITS_IO_CHECKER_H

#include <string>

class FitsIOChecker
{
protected:
	static void checkStatus(int status, const std::string& filename);
	static void checkStatus(int status, const std::string& filename, const std::string& operation);
public:
	enum Unit {
		JanskyPerBeam,
		JanskyPerPixel,
		Jansky,
		Kelvin,
		MilliKelvin
	};
	static const char* UnitName(Unit unit) {
		switch(unit)
		{
			case JanskyPerBeam: return "Jansky/beam";
			case JanskyPerPixel: return "Jansky/pixel";
			case Jansky: return "Jansky";
			case Kelvin: return "Kelvin";
			case MilliKelvin: return "Milli-Kelvin";
		}
		return "";
	}
};

#endif
