#ifndef FACET_PREDICT_H
#define FACET_PREDICT_H

#ifdef HAVE_IDG

#include "DS9FacetFile.h"
#include "FacetImage.h"

#include <idg-api.h>

#include "IDGConfiguration.h"

#include "FitsReader.h"
#include "FitsWriter.h"

#include "../Common/UVector.h"

#include <algorithm>
#include <functional>
#include <string>
#include <vector>

class FacetPredict
{
public:
  FacetPredict(const std::string& fitsModelFile, const std::string& ds9RegionsFile) :
    _reader(fitsModelFile),
    _padding(1.0)
  {
    DS9FacetFile f(ds9RegionsFile);
    const size_t width = _reader.ImageWidth(), height = _reader.ImageHeight();
    ao::uvector<double> model(width * height);
    _reader.Read(model.data());
    _fullWidth = width;
    _fullHeight = height;
    _pixelSizeX = _reader.PixelSizeX();
    _pixelSizeY = _reader.PixelSizeY();
    

    FacetMap map;
    f.Read(map, _reader.PhaseCentreRA(), _reader.PhaseCentreDec(), _pixelSizeX, _pixelSizeY, _fullWidth, _fullHeight);
    std::cout << "Read " << map.NFacets() << " facet definitions.\n";
    
    bool makeSquare = true; // only necessary for IDG though
    size_t area = 0;
    for(size_t i=0; i!=map.NFacets(); ++i)
    {
      const Facet& facet = map[i];
      _directions.emplace_back(facet.RA(), facet.Dec());
      _images.emplace_back();
      FacetImage& image = _images.back();
      image.CopyFacetPart(facet, model.data(), _fullWidth, _fullHeight, _padding, makeSquare);
      area += image.Width() * image.Height();
    }
    std::cout << "Area covered: " << area/1024 << " Kpixels^2\n";
  }
  
  void SetMSInfo(std::vector<std::vector<double>>&& bands, size_t nr_stations)
  {
        _maxBaseline = 1.0/std::min(_pixelSizeX, _pixelSizeY);
    _maxW = _maxBaseline * 0.1;
    std::cout << "Predicting baselines up to " << _maxBaseline << " wavelengths.\n";
    _bands = std::move(bands);
    _nr_stations = nr_stations;
  }
  
  void SetMSInfo(double maxW, std::vector<std::vector<double>>&& bands, size_t nr_stations, double max_baseline)
  {
    _maxW = maxW;
    _bands = std::move(bands);
    _nr_stations = nr_stations;
          _maxBaseline = max_baseline;
  }
  
  bool IsStarted() const { return !_buffersets.empty(); }
  
  void StartIDG(bool saveFacets)
  {
    _buffersets.clear();
    idg::api::Type proxyType = idg::api::Type::CPU_OPTIMIZED;
    int buffersize = 256;
    idg::api::options_type options;
    IdgConfiguration::Read(proxyType, buffersize, options);
    ao::uvector<double> data;
    _metaData.clear();
    for(FacetImage& img : _images)
    {
      double
        dl = ( int(_reader.ImageWidth()/2) - (img.OffsetX()+int(img.Width()/2)) ) * _pixelSizeX,
        dm = (img.OffsetY()+int(img.Height()/2) - int(_reader.ImageHeight()/2)) * _pixelSizeY,
        dp = sqrt(1.0 - dl*dl - dm*dm) - 1.0;
      std::cout << "Initializing gridder " << _buffersets.size() << " (" << img.Width() << " x " << img.Height() << ", +" << img.OffsetX() << "," << img.OffsetY() << ", dl=" << dl*180.0/M_PI << " deg, dm=" << dm*180.0/M_PI << " deg)\n";
      
      data.assign(img.Width() * img.Height() * 4, 0.0);
      // TODO make full polarization
      std::copy(img.Data(), img.Data()+img.Width()*img.Height(), data.data());
      
      _buffersets.emplace_back(idg::api::BufferSet::create(proxyType));
      idg::api::BufferSet& bs = *_buffersets.back();
      options["padded_size"] = size_t(1.2*img.Width());
      //options["max_threads"] = int(1);
      bs.init(img.Width(), _pixelSizeX, _maxW+1.0, dl, dm, dp, options);
      bs.set_image(data.data());
      bs.init_buffers(buffersize, _bands, _nr_stations, _maxBaseline, options, idg::api::BufferSetType::degridding);
      
      if(saveFacets)
      {
        FitsWriter writer;
        writer.SetImageDimensions(img.Width(), img.Height(), _reader.PhaseCentreRA(), _reader.PhaseCentreDec(), _pixelSizeX, _pixelSizeY);
        writer.SetPhaseCentreShift(dl, dm);
        writer.Write("facet" + std::to_string(_metaData.size()) + ".fits", img.Data());
      }
      
      _metaData.emplace_back();
      FacetMetaData& m = _metaData.back();
      m.dl = dl;
      m.dm = dm;
      m.dp = dp;
      m.isInitialized = false;
      m.rowIdOffset = 0;
    }
  }
  
  void RequestPredict(size_t direction, size_t dataDescId, size_t rowId, size_t timeIndex, size_t antenna1, size_t antenna2, const double* uvw)
  {
    double uvwr2 = uvw[0]*uvw[0] + uvw[1]*uvw[1] + uvw[2]*uvw[2];
    if(uvw[2] > _maxW && uvwr2 <= _maxBaseline*_maxBaseline)
    {
      Flush();
      _maxW *= 1.5;
      std::cout << "Increasing maximum w to " << _maxW << '\n';
      StartIDG(false);
    }
    idg::api::BufferSet& bs = *_buffersets[direction];
    FacetMetaData& meta = _metaData[direction];
    if(!meta.isInitialized)
    {
      meta.rowIdOffset = rowId;
      meta.isInitialized = true;
    }
    size_t localRowId = rowId - meta.rowIdOffset;
    if(meta.uvws.size() <= localRowId*3)
      meta.uvws.resize((localRowId+1)*3);
    for(size_t i=0; i!=3; ++i)
      meta.uvws[localRowId*3 + i] = uvw[i];
    double uvwFlipped[3] = { uvw[0], -uvw[1], -uvw[2] }; // IDG uses a flipped coordinate system
    while (bs.get_degridder(dataDescId)->request_visibilities(rowId, timeIndex, antenna1, antenna2, uvwFlipped))
    {
      computePredictionBuffer(dataDescId, direction);
    }
  }
  
  std::function<void(size_t /*row*/, size_t /*direction*/, size_t /*dataDescId*/, const std::complex<float>* /*values*/)> PredictCallback;
  
  size_t NDirections() const { return _images.size(); }
  
  std::pair<double, double> Direction(size_t facet) const { return _directions[facet]; }
  
  void Flush()
    {
    for(size_t b=0; b!=_bands.size(); ++b)
    {
      for(size_t direction = 0; direction != _buffersets.size(); ++direction)
        computePredictionBuffer(b, direction);
    }
  }
	
private:
  void computePredictionBuffer(size_t dataDescId, size_t direction)
  {
    idg::api::BufferSet& bs = *_buffersets[direction];
    auto available_row_ids = bs.get_degridder(dataDescId)->compute();
    for(auto i : available_row_ids)
    {
      size_t row = i.first;
      std::complex<float>* values = i.second;
      size_t nChan = _bands[dataDescId].size();
      size_t localRow = row - _metaData[direction].rowIdOffset;
      double* uvw = &_metaData[direction].uvws[localRow*3];
      double
        dlFact = 2.0 * M_PI * _metaData[direction].dl,
        dmFact = 2.0 * M_PI * _metaData[direction].dm,
        dpFact = 2.0 * M_PI * _metaData[direction].dp;
      for(size_t ch=0; ch!=nChan; ++ch)
      {
        double
          lambda = wavelength(dataDescId, ch),
          u = uvw[0] / lambda,
          v = uvw[1] / lambda,
          w = uvw[2] / lambda;
        double angle = u * dlFact + v * dmFact + w * dpFact;
        float
          rotSin = sin(angle),
          rotCos = cos(angle);
        for(size_t p=0; p!=4; ++p) {
          std::complex<float> s = values[ch*4 + p];
          values[ch*4 + p] = std::complex<float>(
            s.real() * rotCos  -  s.imag() * rotSin,
            s.real() * rotSin  +  s.imag() * rotCos);
        }
      }
      PredictCallback(row, direction, dataDescId, values);
    }
    bs.get_degridder(dataDescId)->finished_reading();
    _metaData[direction].isInitialized = false;
  }
  
  constexpr static double c() { return 299792458.0L; }
  
  double wavelength(size_t dataDescId, size_t channel) const
  {
    return c() / _bands[dataDescId][channel];
  }
  
  std::vector<FacetImage> _images;
  std::vector<std::unique_ptr<idg::api::BufferSet>> _buffersets;
  struct FacetMetaData {
    double dl, dm, dp;
    bool isInitialized;
    size_t rowIdOffset;
    std::vector<double> uvws;
  };
  std::vector<FacetMetaData> _metaData;
  
  size_t _fullWidth, _fullHeight;
  double _pixelSizeX, _pixelSizeY;
  FitsReader _reader;
  double _padding;
  
  // MS info
  double _maxW;
  std::vector<std::vector<double>> _bands;
  size_t _nr_stations;
  double _maxBaseline;
  std::vector<std::pair<double, double>> _directions;
};

#else // HAVE_IDG

#include <algorithm>
#include <functional>
#include <string>
#include <vector>

class FacetPredict
{
public:
	FacetPredict(const std::string&, const std::string&)
	{ notCompiled(); }
	
	void SetMSInfo(std::vector<std::vector<double>>&& bands, size_t nr_stations)
	{ notCompiled(); }
	
	void SetMSInfo(double maxW, std::vector<std::vector<double>>&& bands, size_t nr_stations, double max_baseline)
	{ notCompiled(); }
	
	bool IsStarted() const
	{ notCompiled(); return false; }
	
	void StartIDG(bool)
	{ notCompiled(); }
	
	void RequestPredict(size_t, size_t, size_t, size_t, size_t, size_t, const double*)
	{ notCompiled(); }
	
	std::function<void(size_t, size_t, size_t, const std::complex<float>*)> PredictCallback;
	
	size_t NDirections() const
	{ notCompiled(); return 0; }
	
	std::pair<double, double> Direction(size_t) const
	{ notCompiled(); return std::pair<double,double>(); }
	
	void Flush()
	{ notCompiled(); }
	
private:
  void notCompiled() const {
    throw std::runtime_error("Facet prediction is not available, because DP3 was not compiled with IDG support");
  }
};

#endif // HAVE_IDG

#endif // header guard
