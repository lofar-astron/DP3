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

  // The following class is unique to each instance of SagecalPredict
  // and will be updated during each process() call.
  class BeamData {
   public:
    BeamData() { sources_prec_ = false; }
    // The following are for beam calculation, but will be updated
    // at each process() call, hence they belong here
    /* time coord UTC (s), size tileszx1,
     convert from MJD (s) to JD (days) */
    std::vector<double> time_utc_;
    bool sources_prec_;
    /* LOFAR HBA tile beam pointing center */
    double b_ra0_;
    double b_dec0_;
    /* pointing center of beams (only one) (could be different from phase
     * center) */
    double p_ra0_;
    double p_dec0_;

   private:
  };

  // The following class is a singleton.
  // Only updated during updateInfo(), and reused (read only) by all
  // SagecalPredict instances.
  class BeamDataSingle {
   private:
    explicit BeamDataSingle() : is_updated_(false) {}
    std::mutex mutex_;
    bool is_updated_;

   public:
    BeamDataSingle(BeamDataSingle const&) = delete;
    BeamDataSingle& operator=(BeamDataSingle const&) = delete;
    ~BeamDataSingle();

    static BeamDataSingle* get_instance() {
      static BeamDataSingle instance;
      return &instance;
    }

    void update_metadata(const base::DPInfo& _info, const double freq_f0,
                         const size_t n_channels,
                         std::vector<double>& freq_channel,
                         const int beam_mode);

    size_t n_stations_; /* N */
    /* number of elements in each station, size Nx1 */
    std::vector<int> n_elem_;
    /* position (ITRF) of stations (m)
     later changed to longitude,latitude,height (rad,rad,m) */
    std::vector<double> sx_; /* x: size Nx1 */
    std::vector<double> sy_; /* y: ... */
    std::vector<double> sz_; /* z: ... */
    /* x,y,z coords of elements, projected, converted to ITRF (m) */
    std::vector<double*> xx_; /* x coord pointer, size Nx1, each *x: x coord of
                    station, size Nelem[]x1 */
    std::vector<double*> yy_; /* y ... */
    std::vector<double*> zz_; /* z ... */

    /* flag to indicate this beam is only a dipole */
    bool is_dipole_;

    /* beamformer type STAT_NONE, STAT_SINGLE or STAT_TILE */
    int beamformer_type_;

    /* LOFAR HBA tile beam pointing center */
    double b_ra0_;
    double b_dec0_;

    /* pointing center of beams (only one) (could be different from phase
     * center) */
    double p_ra0_;
    double p_dec0_;

    int beam_mode_{DOBEAM_NONE};
    elementcoeff ecoeff;
  };
#endif /* HAVE_LIBDIRAC || HAVE_LIBDIRAC_CUDA */

  // Singleton to open and read H5 file only by one thread
  // But H5Parm itself is assumed thread-safe
  class H5ParmSingle {
   private:
    explicit H5ParmSingle() : is_opened_(false) {}
    std::mutex mutex_;
    bool is_opened_;

   public:
    H5ParmSingle(H5ParmSingle const&) = delete;
    H5ParmSingle& operator=(H5ParmSingle const&) = delete;
    ~H5ParmSingle(){};

    static H5ParmSingle* get_instance() {
      static H5ParmSingle instance;
      return &instance;
    }

    schaapcommon::h5parm::H5Parm& open_file(const std::string& h5_name);

    schaapcommon::h5parm::H5Parm h5_parm_;
  };

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
  H5ParmSingle* h5_parm_reference_{nullptr};
  std::string solset_name_;
  std::string soltab_name_;
  bool invert_{false};
  bool parm_on_disk_{false};
  bool use_amp_phase_{false};
  double sigma_mmse_;
  GainType gain_type_;
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

#if defined(HAVE_LIBDIRAC) || defined(HAVE_LIBDIRAC_CUDA)
  IOData iodata_;
  BeamData runtime_beam_data_;
  BeamDataSingle* beam_data_{nullptr};
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
