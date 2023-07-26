// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Interface to sagecal model prediction

#ifndef DP3_STEPS_SAGECALPREDICT_H_
#define DP3_STEPS_SAGECALPREDICT_H_

#include <dp3/steps/Step.h>
#include <dp3/base/DP3.h>
#include <dp3/base/DPBuffer.h>
#include "../base/ModelComponent.h"
#include "../base/Patch.h"
#include "../base/PredictBuffer.h"
#include "../base/SourceDBUtil.h"
#include "../common/ParameterSet.h"
#include "ApplyCal.h"
#include "ResultStep.h"

#include "../base/ComponentInfo.h"

#if defined(HAVE_LIBDIRAC) || defined(HAVE_LIBDIRAC_CUDA)
#include <Dirac_radio.h>
// Remove 'complex' def here as we do not need it afterwards
#undef complex
#endif

namespace dp3 {
namespace steps {

class SagecalPredict : public ModelDataStep {
#if defined(HAVE_LIBDIRAC) || defined(HAVE_LIBDIRAC_CUDA)
  class IOData {
   public:
    IOData(){};
    std::vector<double> data_;
    std::vector<double> u_;
    std::vector<double> v_;
    std::vector<double> w_;
    std::vector<double> freqs_;
    std::vector<clus_source_t> cluster_arr_;
    std::vector<baseline_t> baseline_arr_;

    size_t n_baselines, n_stations, n_channels;
    double f0;
    double fdelta;

   private:
  };
  class BeamData {
   public:
    BeamData() { sources_prec = false; };
    size_t n_stations;            /* N */
    std::vector<double> time_utc; /* time coord UTC (s), size tileszx1,
                       convert from MJD (s) to JD (days) */
    std::vector<int> n_elem;      /* no of elements in each station, size Nx1 */
    /* position (ITRF) of stations (m)
     later changed to logitude,latitude,height (rad,rad,m) */
    std::vector<double> sx; /* x: size Nx1 */
    std::vector<double> sy; /* y: ... */
    std::vector<double> sz; /* z: ... */
    /* x,y,z coords of elements, projected, converted to ITRF (m) */
    std::vector<double*> xx; /* x coord pointer, size Nx1, each *x: x coord of
                    station, size Nelem[]x1 */
    std::vector<double*> yy; /* y ... */
    std::vector<double*> zz; /* z ... */

    /* LOFAR HBA tile beam pointing center */
    double b_ra0;
    double b_dec0;

    /* pointing center of beams (only one) (could be different from phase
     * center) */
    double p_ra0;
    double p_dec0;

    /* flag to indicate this beam is only a dipole */
    bool isDipole;

    /* beamformer type STAT_NONE, STAT_SINGLE or STAT_TILE */
    int beamformer_type;

    elementcoeff ecoeff;
    bool sources_prec;

   private:
  };
#endif /* HAVE_LIBDIRAC || HAVE_LIBDIRAC_CUDA */

 public:
  SagecalPredict(const common::ParameterSet&, const std::string& prefix,
                 MsType input_type = MsType::kRegular);
  // The following constructor is to be used within DDECal
  SagecalPredict(const common::ParameterSet&, const std::string& prefix,
                 const std::vector<std::string>& given_source_patterns,
                 MsType input_type = MsType::kRegular);

  ~SagecalPredict() override;

  common::Fields getRequiredFields() const override {
    // We need UVW and data
    common::Fields fields = kUvwField;
    if ((operation_ == Operation::kAdd) ||
        (operation_ == Operation::kSubtract)) {
      fields |= kDataField;
    }
    return fields;
  }

  common::Fields getProvidedFields() const override {
    // We only modify the output data
    common::Fields fields = kDataField;
    return fields;
  }

  ///
  /// buffer, copies the buffer and change its data.
  bool process(std::unique_ptr<dp3::base::DPBuffer> buffer) override;

  /// Update the general info.
  void updateInfo(const base::DPInfo&) override;

  void finish() override;

  void show(std::ostream&) const override;

  void showTimings(std::ostream& os, double duration) const override;

  base::Direction GetFirstDirection() const override;

  void setNThreads(const size_t n_threads) { n_threads_ = n_threads; }

 private:
  enum class Operation { kReplace, kAdd, kSubtract };

  const MsType ms_type_;

  void SetOperation(const std::string& operation);

  unsigned int nPol(const std::string& parmName);

  void setCorrectType(std::vector<std::string>& solTabs);

  void updateFromH5(const double startTime);

  /// The actual constructor
  void init(const common::ParameterSet&, const std::string& prefix,
            const std::vector<std::string>& source_patterns);

  std::string name_;     ///< The name of the step (or prefix)
  Operation operation_;  ///< Operation to use on the DATA column

  std::string h5_name_;  ///< HDF5 parameter file
  std::vector<std::string>
      directions_list_;  ///< list of directions (one or more patches)
  std::vector<std::shared_ptr<base::Patch>> patch_list_;
  std::string source_db_name_;
  std::vector<std::pair<std::shared_ptr<base::ModelComponent>,
                        std::shared_ptr<base::Patch>>>
      source_list_;
  bool any_orientation_is_absolute_{false};  ///< Any of the Gaussian sources
                                             ///< has absolute orientation

  // H5 solutions handling
  schaapcommon::h5parm::H5Parm h5_parm_;
  std::string solset_name_;
  std::string soltab_name_;
  bool invert_{false};
  bool parm_on_disk_{false};
  bool use_amp_phase_{false};
  double sigma_mmse_;
  JonesParameters::CorrectType corr_type_;
  JonesParameters::InterpolationType interp_type_;
  JonesParameters::MissingAntennaBehavior missing_ant_behavior_;
  schaapcommon::h5parm::SolTab sol_tab_;
  schaapcommon::h5parm::SolTab sol_tab2_;
  std::vector<casacore::String> parm_expressions_;
  unsigned int timeslots_per_parmupdate_;
  unsigned int timestep_;
  double time_interval_;
  double time_last_;
  std::vector<std::string> soltab_names_;

  std::vector<double> params_;

  base::Direction phase_ref_;

  common::NSTimer timer_;

  // This should be set in updateInfo()
  size_t n_threads_{0};

#if defined(HAVE_LIBDIRAC) || defined(HAVE_LIBDIRAC_CUDA)
  IOData iodata_;
  BeamData beam_data_;
  void loadData(std::unique_ptr<dp3::base::DPBuffer>& buffer);
  void writeData(std::unique_ptr<dp3::base::DPBuffer>& buffer);
  void readAuxData(const base::DPInfo&);
  int beam_mode_{DOBEAM_NONE};
  std::vector<int> ignore_list_;
#endif /* HAVE_LIBDIRAC || HAVE_LIBDIRAC_CUDA */
};

}  // namespace steps
}  // namespace dp3

#endif  // header guard
