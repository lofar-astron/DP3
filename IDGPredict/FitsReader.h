#ifndef FITSREADER_H
#define FITSREADER_H

#include <string>
#include <vector>

#include <fitsio.h>

#include "Polarization.h"
#include "FitsIOChecker.h"

class FitsReader : public FitsIOChecker
{
	public:
		explicit FitsReader(const std::string &filename) 
		: FitsReader(filename, true, false)
		{ }
		explicit FitsReader(const std::string &filename, bool checkCType, bool allowMultipleImages=false) :
			_meta(filename, checkCType, allowMultipleImages)
		{
			initialize();
		}
		FitsReader(const FitsReader& source);
		FitsReader(FitsReader&& source);
		~FitsReader();
		
		FitsReader& operator=(const FitsReader& rhs);
		FitsReader& operator=(FitsReader&& rhs);
		
		template<typename NumType> void ReadIndex(NumType *image, size_t index);
		
		template<typename NumType> void Read(NumType *image)
		{
			ReadIndex(image, 0);
		}
		
		size_t ImageWidth() const { return _meta.imgWidth; }
		size_t ImageHeight() const { return _meta.imgHeight; }
		
		double PhaseCentreRA() const { return _meta.phaseCentreRA; }
		double PhaseCentreDec() const { return _meta.phaseCentreDec; }
		
		double PixelSizeX() const { return _meta.pixelSizeX; }
		double PixelSizeY() const { return _meta.pixelSizeY; }
		
		double PhaseCentreDL() const { return _meta.phaseCentreDL; }
		double PhaseCentreDM() const { return _meta.phaseCentreDM; }
		
		double Frequency() const { return _meta.frequency; }
		double Bandwidth() const { return _meta.bandwidth; }
		
		double DateObs() const { return _meta.dateObs; }
		PolarizationEnum Polarization() const { return _meta.polarization; }
		
		FitsIOChecker::Unit Unit() const { return _meta.unit; }
		
		bool HasBeam() const { return _meta.hasBeam; }
		double BeamMajorAxisRad() const { return _meta.beamMajorAxisRad; }
		double BeamMinorAxisRad() const { return _meta.beamMinorAxisRad; }
		double BeamPositionAngle() const { return _meta.beamPositionAngle; }
		
		const std::string& TelescopeName() const { return _meta.telescopeName; }
		const std::string& Observer() const { return _meta.observer; }
		const std::string& ObjectName() const { return _meta.objectName; }
		
		const std::string& Origin() const { return _meta.origin; }
		const std::string& OriginComment() const { return _meta.originComment; }
		
		const std::vector<std::string>& History() const { return _meta.history; }
		
		bool ReadDoubleKeyIfExists(const char* key, double& dest);
		bool ReadStringKeyIfExists(const char* key, std::string& dest) {
			std::string c;
			return ReadStringKeyIfExists(key, dest, c);
		}
		bool ReadStringKeyIfExists(const char* key, std::string& value, std::string& comment);
		bool ReadFloatKeyIfExists(const char* key, float& dest);
		
		static double ParseFitsDateToMJD(const char* valueStr);
		
		const std::string& Filename() const { return _meta.filename; }
		
		fitsfile* FitsHandle() const { return _fitsPtr; }
		
		size_t NMatrixElements() const { return _meta.nMatrixElements; }
		size_t NFrequencies() const { return _meta.nFrequencies; }
		size_t NAntennas() const { return _meta.nAntennas; }
		size_t NTimesteps() const { return _meta.nTimesteps; }
		
		double TimeDimensionStart() const { return _meta.timeDimensionStart; }
		double TimeDimensionIncr() const { return _meta.timeDimensionIncr; }
		
	private:
		double readDoubleKey(const char* key);
		std::string readStringKey(const char* key);
		void readHistory();
		bool readDateKeyIfExists(const char *key, double &dest);
		
		void initialize();
		
		fitsfile* _fitsPtr;
		
		struct MetaData {
			MetaData(const std::string &filename_, bool checkCType_, bool allowMultipleImages_) :
				filename(filename_),
				hasBeam(false),
				checkCType(checkCType_),
				allowMultipleImages(allowMultipleImages_)
			{ }
			std::string filename;
			size_t imgWidth, imgHeight;
			size_t nMatrixElements, nAntennas, nFrequencies, nTimesteps;
			double phaseCentreRA, phaseCentreDec;
			double pixelSizeX, pixelSizeY;
			double phaseCentreDL, phaseCentreDM;
			double frequency, bandwidth, dateObs;
			bool hasBeam;
			double beamMajorAxisRad, beamMinorAxisRad, beamPositionAngle;
			double timeDimensionStart, timeDimensionIncr;
			
			PolarizationEnum polarization;
			FitsIOChecker::Unit unit;
			std::string telescopeName, observer, objectName;
			std::string origin, originComment;
			std::vector<std::string> history;
			
			bool checkCType, allowMultipleImages;
		} _meta;
};

#endif
