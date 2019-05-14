#ifndef FITSWRITER_H
#define FITSWRITER_H

#include <fitsio.h>

#include <string>
#include <vector>
#include <map>
#include <cmath>

#include "Polarization.h"
#include "FitsIOChecker.h"

class FitsWriter : public FitsIOChecker
{
public:
	enum DimensionType { 
		FrequencyDimension, 
		PolarizationDimension, 
		AntennaDimension, 
		TimeDimension,
		MatrixDimension
	};
	
	FitsWriter() :
		_width(0), _height(0),
		_phaseCentreRA(0.0), _phaseCentreDec(0.0), _pixelSizeX(0.0), _pixelSizeY(0.0),
		_phaseCentreDL(0.0), _phaseCentreDM(0.0),
		_frequency(0.0), _bandwidth(0.0),
		_dateObs(0.0),
		_hasBeam(false),
		_beamMajorAxisRad(0.0), _beamMinorAxisRad(0.0), _beamPositionAngle(0.0),
		_polarization(Polarization::StokesI),
		_unit(JanskyPerBeam),
		_isUV(false),
		_telescopeName(), _observer(), _objectName(),
		_origin("AO/WSImager"), _originComment("Imager written by Andre Offringa"),
		_multiFPtr(nullptr)
	{
	}
	
	explicit FitsWriter(const class FitsReader& reader) :
		_width(0), _height(0),
		_phaseCentreRA(0.0), _phaseCentreDec(0.0), _pixelSizeX(0.0), _pixelSizeY(0.0),
		_phaseCentreDL(0.0), _phaseCentreDM(0.0),
		_frequency(0.0), _bandwidth(0.0),
		_dateObs(0.0),
		_hasBeam(false),
		_beamMajorAxisRad(0.0), _beamMinorAxisRad(0.0), _beamPositionAngle(0.0),
		_polarization(Polarization::StokesI),
		_unit(JanskyPerBeam),
		_isUV(false),
		_telescopeName(), _observer(), _objectName(),
		_origin("AO/WSImager"), _originComment("Imager written by Andre Offringa"),
		_multiFPtr(nullptr)
	{
		SetMetadata(reader);
	}
	
	~FitsWriter()
	{
		if(_multiFPtr != nullptr)
			FinishMulti();
	}
	
	template<typename NumType> void Write(const std::string& filename, const NumType* image) const;
	
	void WriteMask(const std::string& filename, const bool* mask) const;
	
	void StartMulti(const std::string& filename);
	
	template<typename NumType>
	void AddToMulti(const NumType* image)
	{
		if(_multiFPtr == 0)
			throw std::runtime_error("AddToMulti() called before StartMulti()");
		writeImage(_multiFPtr, _multiFilename, image, _currentPixel.data());
		size_t index = 2;
		_currentPixel[index]++;
		while(index < _currentPixel.size()-1 && _currentPixel[index] > long(_extraDimensions[index-2].size))
		{
			_currentPixel[index] = 1;
			++index;
			_currentPixel[index]++;
		}
	}
	
	void FinishMulti();
	
	void SetBeamInfo(double widthRad)
	{
		SetBeamInfo(widthRad, widthRad, 0.0);
	}
	void SetBeamInfo(double majorAxisRad, double minorAxisRad, double positionAngleRad)
	{
		_hasBeam = true;
		_beamMajorAxisRad = majorAxisRad;
		_beamMinorAxisRad = minorAxisRad;
		_beamPositionAngle = positionAngleRad;
	}
	void SetNoBeamInfo()
	{
		_hasBeam = false;
		_beamMajorAxisRad = 0.0;
		_beamMinorAxisRad = 0.0;
		_beamPositionAngle = 0.0;
	}
	void SetImageDimensions(size_t width, size_t height)
	{
		_width = width;
		_height = height;
	}
	void SetImageDimensions(size_t width, size_t height, double pixelSizeX, double pixelSizeY)
	{
		_width = width;
		_height = height;
		_pixelSizeX = pixelSizeX;
		_pixelSizeY = pixelSizeY;
	}
	void SetImageDimensions(size_t width, size_t height, double phaseCentreRA, double phaseCentreDec, double pixelSizeX, double pixelSizeY)
	{
		_width = width;
		_height = height;
		_phaseCentreRA = phaseCentreRA;
		_phaseCentreDec = phaseCentreDec;
		_pixelSizeX = pixelSizeX;
		_pixelSizeY = pixelSizeY;
	}
	void SetFrequency(double frequency, double bandwidth)
	{
		_frequency = frequency;
		_bandwidth = bandwidth;
	}
	void SetDate(double dateObs)
	{
		_dateObs = dateObs;
	}
	void SetPolarization(PolarizationEnum polarization)
	{
		_polarization = polarization;
	}
	Unit GetUnit() const { return _unit; }
	void SetUnit(Unit unit)
	{
		_unit = unit;
	}
	void SetIsUV(bool isUV)
	{
		_isUV = isUV;
	}
	void SetTelescopeName(const std::string& telescopeName)
	{
		_telescopeName = telescopeName;
	}
	void SetObserver(const std::string& observer)
	{
		_observer = observer;
	}
	void SetObjectName(const std::string& objectName)
	{
		_objectName = objectName;
	}
	void SetOrigin(const std::string& origin, const std::string& comment)
	{
		_origin = origin;
		_originComment = comment;
	}
	void SetHistory(const std::vector<std::string>& history)
	{
		_history = history;
	}
	void AddHistory(const std::string& historyLine)
	{
		_history.push_back(historyLine);
	}

	void SetMetadata(const class FitsReader& reader);
	
	double RA() const { return _phaseCentreRA; }
	double Dec() const { return _phaseCentreDec; }
	double Frequency() const { return _frequency; }
	double Bandwidth() const { return _bandwidth; }
	double BeamSizeMajorAxis() const { return _beamMajorAxisRad; }
	double BeamSizeMinorAxis() const { return _beamMinorAxisRad; }
	double BeamPositionAngle() const { return _beamPositionAngle; }
	
	void SetExtraKeyword(const std::string& name, const std::string& value)
	{
		if(_extraStringKeywords.count(name) != 0)
			_extraStringKeywords.erase(name);
		_extraStringKeywords.insert(std::make_pair(name, value));
	}
	void SetExtraKeyword(const std::string& name, double value)
	{
		if(_extraNumKeywords.count(name) != 0)
			_extraNumKeywords.erase(name);
		_extraNumKeywords.insert(std::make_pair(name, value));
	}
	void RemoveExtraKeyword(const std::string& name)
	{
		if(_extraNumKeywords.count(name) != 0)
			_extraNumKeywords.erase(name);
		if(_extraStringKeywords.count(name) != 0)
			_extraStringKeywords.erase(name);
	}
	void SetExtraStringKeywords(const std::map<std::string, std::string>& keywords)
	{
		_extraStringKeywords = keywords;
	}
	void SetExtraNumKeywords(const std::map<std::string, double>& keywords)
	{
		_extraNumKeywords = keywords;
	}
	void SetPhaseCentreShift(double dl, double dm)
	{
		_phaseCentreDL = dl;
		_phaseCentreDM = dm;
	}
	size_t Width() const { return _width; }
	size_t Height() const { return _height; }
	double PhaseCentreDL() const { return _phaseCentreDL; }
	double PhaseCentreDM() const { return _phaseCentreDM; }
	
	void CopyDoubleKeywordIfExists(class FitsReader& reader, const char* keywordName);
	void CopyStringKeywordIfExists(class FitsReader& reader, const char* keywordName);
	
	static void MJDToHMS(double mjd, int& hour, int& minutes, int& seconds, int& deciSec);
	
	void AddExtraDimension(enum DimensionType type, size_t size)
	{
		_extraDimensions.emplace_back(Dimension{type, size});
	}
	void SetTimeDirectionStart(double time) { _timeDirectionStart = time; }
	void SetTimeDirectionInc(double dTime) { _timeDirectionInc = dTime; }
private:
	struct Dimension
	{
		DimensionType type;
		size_t size;
	};
	
	template<typename T>
	static T setNotFiniteToZero(T num)
	{
		return std::isfinite(num) ? num : 0.0;
	}
	std::size_t _width, _height;
	double _phaseCentreRA, _phaseCentreDec, _pixelSizeX, _pixelSizeY;
	double _phaseCentreDL, _phaseCentreDM;
	double _frequency, _bandwidth;
	double _dateObs;
	bool _hasBeam;
	double _beamMajorAxisRad, _beamMinorAxisRad, _beamPositionAngle;
	PolarizationEnum _polarization;
	Unit _unit;
	bool _isUV;
	std::string _telescopeName, _observer, _objectName;
	std::string _origin, _originComment;
	std::vector<std::string> _history;
	std::vector<Dimension> _extraDimensions;
	std::map<std::string, std::string> _extraStringKeywords;
	std::map<std::string, double> _extraNumKeywords;
	double _timeDirectionStart, _timeDirectionInc;
	
	void julianDateToYMD(double jd, int &year, int &month, int &day) const;
	void writeHeaders(fitsfile*& fptr, const std::string& filename) const;
	void writeHeaders(fitsfile*& fptr, const std::string& filename, const std::vector<Dimension>& extraDimensions) const;
	void writeImage(fitsfile* fptr, const std::string& filename, const double* image, long* currentPixel) const;
	void writeImage(fitsfile* fptr, const std::string& filename, const float* image, long* currentPixel) const;
	template<typename NumType>
	void writeImage(fitsfile* fptr, const std::string& filename, const NumType* image, long* currentPixel) const;
	
	std::string _multiFilename;
	fitsfile *_multiFPtr;
	std::vector<long> _currentPixel;
};

#endif
