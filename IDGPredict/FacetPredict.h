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
  FacetPredict(const std::vector<std::string> fitsModelFiles, const std::string& ds9RegionsFile) :
    _padding(1.0),
    _bufferSize(0)
  {
    if(fitsModelFiles.empty())
      throw std::runtime_error("No fits files specified for IDG predict");
    _readers.reserve(fitsModelFiles.size());
    for(const std::string& file : fitsModelFiles)
      _readers.emplace_back(file);
    
    DS9FacetFile f(ds9RegionsFile);
    _fullWidth = _readers.front().ImageWidth();
    _fullHeight = _readers.front().ImageHeight();
    _refFrequency = _readers.front().Frequency();
    _pixelSizeX = _readers.front().PixelSizeX();
    _pixelSizeY = _readers.front().PixelSizeY();
    std::vector<ao::uvector<double>> models(_readers.size());
    for(size_t img=0; img!=_readers.size(); ++img)
    {
      if(_readers[img].ImageWidth() != _fullWidth || _readers[img].ImageHeight() != _fullHeight)
        throw std::runtime_error("Image for spectral term " + std::to_string(img) + " has inconsistent dimensions");
      if(_readers[img].PixelSizeX() != _pixelSizeX || _readers[img].PixelSizeY() != _pixelSizeY)
        throw std::runtime_error("Pixel size of spectral term " + std::to_string(img) + " is inconsistent with first spectral term");
      models[img].resize(_fullWidth * _fullHeight);
      _readers[img].Read(models[img].data());
    }

    FacetMap map;
    f.Read(map, _readers.front().PhaseCentreRA(), _readers.front().PhaseCentreDec(), _pixelSizeX, _pixelSizeY, _fullWidth, _fullHeight);
    std::cout << "Read " << map.NFacets() << " facet definitions.\n";
    
    bool makeSquare = true; // only necessary for IDG though
    size_t area = 0;
    for(size_t i=0; i!=map.NFacets(); ++i)
    {
      const Facet& facet = map[i];
      _directions.emplace_back(facet.RA(), facet.Dec());
      _images.emplace_back();
      FacetImage& image = _images.back();
      image.CopyFacetPart(facet, models, _fullWidth, _fullHeight, _padding, makeSquare);
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
    size_t nTerms = _readers.size();
    
    size_t maxChannels = 0;
    for(std::vector<double>& band : _bands)
      maxChannels = std::max(maxChannels, band.size());
    long int pageCount = sysconf(_SC_PHYS_PAGES), pageSize = sysconf(_SC_PAGE_SIZE);
    int64_t memory = (int64_t) pageCount * (int64_t) pageSize;
    uint64_t memPerTimestep = idg::api::BufferSet::get_memory_per_timestep(_nr_stations, maxChannels);
    memPerTimestep *= 2; // IDG uses two internal buffer
    // Allow the directions together to use 1/4th of the available memory for the vis buffers.
    size_t allocatableTimesteps = memory/4/_images.size()/nTerms / memPerTimestep;
    // TODO once a-terms are supported, this should include the size required for the a-terms.
    std::cout << "Allocatable timesteps per direction: " << allocatableTimesteps << '\n';
    
    int buffersize = std::max(allocatableTimesteps, size_t(1));
    if(_bufferSize != 0)
    {
      buffersize = _bufferSize;
      std::cout << "Buffer size manually set to " << buffersize << " timesteps\n";
    }
    idg::api::options_type options;
    IdgConfiguration::Read(proxyType, buffersize, options);
    std::vector<ao::uvector<double>> data(nTerms);
    _metaData.clear();
    FitsReader& reader = _readers.front();
    for(FacetImage& img : _images)
    {
      double
        dl = ( int(reader.ImageWidth()/2) - (img.OffsetX()+int(img.Width()/2)) ) * _pixelSizeX,
        dm = (img.OffsetY()+int(img.Height()/2) - int(reader.ImageHeight()/2)) * _pixelSizeY,
        dp = sqrt(1.0 - dl*dl - dm*dm) - 1.0;
      std::cout << "Initializing gridder " << _buffersets.size() << " (" << img.Width() << " x " << img.Height() << ", +" << img.OffsetX() << "," << img.OffsetY() << ", dl=" << dl*180.0/M_PI << " deg, dm=" << dm*180.0/M_PI << " deg)\n";
      
      // TODO make full polarization
      for(size_t term=0; term!=nTerms; ++term)
      {
        data[term].assign(img.Width() * img.Height() * 4, 0.0);
        std::copy(img.Data(term), img.Data(term)+img.Width()*img.Height(), data[term].data());
        
        _buffersets.emplace_back(idg::api::BufferSet::create(proxyType));
        idg::api::BufferSet& bs = *_buffersets.back();
        options["padded_size"] = size_t(1.2*img.Width());
        //options["max_threads"] = int(1);
        bs.init(img.Width(), _pixelSizeX, _maxW+1.0, dl, dm, dp, options);
        bs.set_image(data[term].data());
        bs.init_buffers(buffersize, _bands, _nr_stations, _maxBaseline, options, idg::api::BufferSetType::degridding);
      }
      
      if(saveFacets)
      {
        FitsWriter writer;
        writer.SetImageDimensions(img.Width(), img.Height(), reader.PhaseCentreRA(), reader.PhaseCentreDec(), _pixelSizeX, _pixelSizeY);
        writer.SetPhaseCentreShift(dl, dm);
        writer.Write("facet" + std::to_string(_metaData.size()) + ".fits", img.Data(0));
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
    size_t nTerms = _readers.size();
    double uvwr2 = uvw[0]*uvw[0] + uvw[1]*uvw[1] + uvw[2]*uvw[2];
    if(uvw[2] > _maxW && uvwr2 <= _maxBaseline*_maxBaseline)
    {
      Flush();
      _maxW *= 1.5;
      std::cout << "Increasing maximum w to " << _maxW << '\n';
      StartIDG(false);
    }
    for(size_t termIndex=0; termIndex!=nTerms; ++termIndex)
    {
      idg::api::BufferSet& bs = *_buffersets[direction*nTerms + termIndex];
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
  }
  
  std::function<void(size_t /*row*/, size_t /*direction*/, size_t /*dataDescId*/, const std::complex<float>* /*values*/)> PredictCallback;
  
  size_t NDirections() const { return _images.size(); }
  
  std::pair<double, double> Direction(size_t facet) const { return _directions[facet]; }
  
  void Flush()
  {
    for(size_t b=0; b!=_bands.size(); ++b)
    {
      for(size_t direction = 0; direction != _directions.size(); ++direction)
        computePredictionBuffer(b, direction);
    }
  }
  
  void SetBufferSize(size_t nTimesteps) { _bufferSize = nTimesteps; }
  
private:
  void computePredictionBuffer(size_t dataDescId, size_t direction)
  {
    size_t nTerms = _readers.size();
    typedef std::vector<std::pair<size_t, std::complex< float >*>> rowidlist_t;
    std::vector<rowidlist_t> available_row_ids(nTerms);
    for(size_t term=0; term!=nTerms; ++term)
    {
      idg::api::BufferSet& bs = *_buffersets[direction*nTerms + term];
      available_row_ids[term] = bs.get_degridder(dataDescId)->compute();
    }
    
    size_t nChan = _bands[dataDescId].size();
    double
      dlFact = 2.0 * M_PI * _metaData[direction].dl,
      dmFact = 2.0 * M_PI * _metaData[direction].dm,
      dpFact = 2.0 * M_PI * _metaData[direction].dp;
    for(size_t i=0; i!=available_row_ids[0].size(); ++i)
    {
      size_t row = available_row_ids[0][i].first;
      size_t localRow = row - _metaData[direction].rowIdOffset;
      const double* uvw = &_metaData[direction].uvws[localRow*3];
      
      // Correct the phase shift of the values for this facet
      for(size_t term=0; term!=nTerms; ++term)
      {
        std::complex<float>* values = available_row_ids[term][i].second;
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
      }
      
      // Apply polynomial-term corrections and add all to values of 'term 0'
      // The "polynomial spectrum" definition is used, equal to the one e.g. used by WSClean in component outputs
      // (see https://sourceforge.net/p/wsclean/wiki/ComponentList/ ) and in text files when 'logarithmic SI' is false. 
      // The definition is:
      //   S(nu) = term0 + term1 (nu/refnu - 1) + term2 (nu/refnu - 1)^2 + ...
      std::complex<float>* values0 = available_row_ids[0][i].second;
      for(size_t ch=0; ch!=nChan; ++ch)
      {
        double frequency = _bands[dataDescId][ch];
        double freqFactor = frequency / _refFrequency - 1.0;
        double polynomialFactor = 1.0;
        for(size_t term=1; term!=nTerms; ++term)
        {
          polynomialFactor *= freqFactor;
          const std::complex<float>* values = available_row_ids[term][i].second;
          for(size_t p=0; p!=4; ++p)
            values0[ch*4 + p] += values[ch*4 + p] * float(polynomialFactor);
        }
      }
      PredictCallback(row, direction, dataDescId, values0);
    }
    for(size_t term=0; term!=nTerms; ++term)
    {
      idg::api::BufferSet& bs = *_buffersets[direction*nTerms + term];
      bs.get_degridder(dataDescId)->finished_reading();
    }
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
  double _refFrequency;
  double _pixelSizeX, _pixelSizeY;
  std::vector<FitsReader> _readers;
  double _padding;
  size_t _bufferSize;
  
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
  FacetPredict(const std::vector<std::string>&, const std::string&)
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
  
  void SetBufferSize(size_t)
  { notCompiled(); }
  
private:
  void notCompiled() const {
    throw std::runtime_error("Facet prediction is not available, because DP3 was not compiled with IDG support");
  }
};

#endif // HAVE_IDG

#endif // header guard
