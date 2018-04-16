#ifndef DPPP_STMANPARSETKEYS_H
#define DPPP_STMANPARSETKEYS_H

#include <Common/ParameterSet.h>

#include <string>

#include <casacore/casa/Containers/Record.h>

namespace LOFAR {
  namespace DPPP {
    struct StManParsetKeys
    {
      casacore::String    stManName;
      uint            dyscoDataBitRate;     //# Bits per data float, or 0 if data column is not compressed
      uint            dyscoWeightBitRate;   //# Bits per weight float, or 0 if weight column is not compressed
      std::string     dyscoDistribution;    //# Distribution assumed for compression; e.g. "Uniform" or "TruncatedGaussian"
      double          dyscoDistTruncation;  //# For truncated distributions, the truncation point (e.g. 3 for 3 sigma in TruncGaus).
      std::string     dyscoNormalization;   //# Kind of normalization; "AF", "RF" or "Row".     
      
      void Set(const ParameterSet& parset, const std::string& prefix)
      {
        stManName = toLower(parset.getString(prefix+"storagemanager",
                                             parset.getString(prefix+"storagemanager.name",
                                                              string())));
        if(stManName == "dysco") {
          dyscoDataBitRate    = parset.getInt(
                                    prefix+"storagemanager.databitrate", 10);
          dyscoWeightBitRate  = parset.getInt(
                                    prefix+"storagemanager.weightbitrate", 12);
          dyscoDistribution   = parset.getString(
                                    prefix+"storagemanager.distribution", 
                                    "TruncatedGaussian");
          dyscoDistTruncation = parset.getDouble(
                                    prefix+"storagemanager.disttruncation", 2.5);
          dyscoNormalization  = parset.getString(
                                    prefix+"storagemanager.normalization", "AF");
        }
      }
      
      casacore::Record GetDyscoSpec() const
      {
        casacore::Record dyscoSpec;
        dyscoSpec.define ("distribution", dyscoDistribution);
        dyscoSpec.define ("normalization", dyscoNormalization);
        dyscoSpec.define ("distributionTruncation", dyscoDistTruncation);
        // DPPP uses bitrate of 0 to disable compression of the data/weight column.
        // However, Dysco does not allow the data or weight bitrates to be set to 0,
        // so we set the values to something different. The values are not actually used.
        uint dataBitRate = dyscoDataBitRate;
        if(dataBitRate == 0)
          dataBitRate = 16;
        dyscoSpec.define ("dataBitCount", dataBitRate);
        uint weightBitRate = dyscoWeightBitRate;
        if(weightBitRate == 0)
          weightBitRate = 16;
        dyscoSpec.define ("weightBitCount", weightBitRate);
        return dyscoSpec;
      }
    };
  }
}
#endif
