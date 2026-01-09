// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Interface to sagecal model prediction

#ifndef DP3_STEPS_SAGECALPREDICT_H_
#define DP3_STEPS_SAGECALPREDICT_H_

#include <Dirac_radio.h>
// Remove 'complex' def here as we do not need it afterwards
#undef complex

#include "base/DP3.h"
#include "base/DPBuffer.h"
#include "base/ComponentInfo.h"
#include "base/ModelComponent.h"
#include "base/PredictBuffer.h"
#include "common/ParameterSet.h"
#include "model/Patch.h"
#include "model/SourceDBUtil.h"
#include "steps/Step.h"
#include "ApplyCal.h"
#include "ResultStep.h"

namespace dp3 {
namespace steps {

class SagecalPredict : public ModelDataStep {
  struct IOData {
    std::vector<double> data;
    std::vector<double> u;
    std::vector<double> v;
    std::vector<double> w;
    std::vector<double> frequencies;
    std::vector<clus_source_t> cluster_arr;
    std::vector<baseline_t> baseline_arr;

    size_t n_baselines;
    size_t n_stations;
    size_t n_channels;
    double f0;
    double fdelta;
  };

  // The following struct is unique to each instance of SagecalPredict
  // and will be updated during each process() call.
  struct BeamData {
    // The following are for beam calculation, but will be updated
    // at each process() call, hence they belong here
    // time coord UTC (s), convert from MJD (s) to JD (days)
    std::vector<double> time_utc;
    bool sources_are_precessed = false;
    // LOFAR HBA tile beam pointing center
    double tile_ra;
    double tile_dec;
    // pointing center of beams (only one) (can be different from phase center)
    double reference_ra;
    double reference_dec;
  };

  // The following struct is a singleton.
  // Only updated during updateInfo(), and reused (read only) by all
  // SagecalPredict instances.
  struct BeamDataSingle {
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

    /* number of elements in each station, size Nx1 */
    std::vector<int> n_elem;
    /* position (ITRF) of stations (m)
     later changed to longitude,latitude,height (rad,rad,m) */
    std::vector<double> sx; /* x: size Nx1 */
    std::vector<double> sy; /* y: ... */
    /* x,y,z coords of elements, projected, converted to ITRF (m) */
    std::vector<std::vector<double>> xx;
    std::vector<std::vector<double>> yy;
    std::vector<std::vector<double>> zz;

    /* pointers the xx, yy and zz data, since libdirac has double** arguments */
    std::vector<double*> xx_ptr;
    std::vector<double*> yy_ptr;
    std::vector<double*> zz_ptr;

    /* beamformer type STAT_NONE, STAT_SINGLE or STAT_TILE */
    int beamformer_type;

    /* LOFAR HBA tile beam pointing center */
    double tile_ra;
    double tile_dec;

    /* pointing center of beams (only one) (could be different from phase
     * center) */
    double reference_ra;
    double reference_dec;

    int beam_mode{DOBEAM_NONE};
    elementcoeff ecoeff;

   private:
    explicit BeamDataSingle() : is_updated_(false) {}
    std::mutex mutex_;
    bool is_updated_;
  };

  // Singleton to open and read H5 file only by one thread
  // But H5Parm itself is assumed thread-safe
  class H5ParmSingle {
   public:
    H5ParmSingle(H5ParmSingle const&) = delete;
    H5ParmSingle& operator=(H5ParmSingle const&) = delete;
    ~H5ParmSingle(){};

    static H5ParmSingle& GetInstance() {
      static H5ParmSingle instance;
      return instance;
    }

    schaapcommon::h5parm::H5Parm& OpenFile(const std::string& h5_name);

   private:
    explicit H5ParmSingle() {}

    /// Protects concurrent access to h5_parm_.
    std::mutex mutex_;

    // Do not use a plain H5Parm object, since its implicitly declared
    // assignment operator of H5Parm is deprecated,
    // Using a (unique) pointer does allow re-assigning to this variable.
    std::unique_ptr<schaapcommon::h5parm::H5Parm> h5_parm_;
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

  void SetCorrectType(schaapcommon::h5parm::H5Parm& h5_parm,
                      const std::vector<std::string>& sol_tabs);

  void updateFromH5(const double startTime);

  /// The actual constructor
  void init(const common::ParameterSet&, const std::string& prefix,
            const std::vector<std::string>& source_patterns);

  void LoadData(const dp3::base::DPBuffer& buffer);
  void WriteData(dp3::base::DPBuffer& buffer) const;
  void ReadAuxData(const base::DPInfo&);

  std::string name_;     ///< The name of the step (or prefix)
  Operation operation_;  ///< Operation to use on the DATA column

  std::string h5_name_;  ///< HDF5 parameter file
  std::vector<std::string>
      directions_list_;  ///< list of directions (one or more patches)
  std::vector<std::shared_ptr<model::Patch>> patch_list_;
  std::string source_db_name_;
  std::vector<std::pair<std::shared_ptr<base::ModelComponent>,
                        std::shared_ptr<model::Patch>>>
      source_list_;
  bool any_orientation_is_absolute_{false};  ///< Any of the Gaussian sources
                                             ///< has absolute orientation

  // H5 solutions handling
  std::string solset_name_;
  std::string soltab_name_;
  bool invert_{false};
  bool parm_on_disk_{false};
  bool use_amp_phase_{false};
  schaapcommon::h5parm::GainType gain_type_;
  schaapcommon::h5parm::JonesParameters::InterpolationType interp_type_;
  schaapcommon::h5parm::JonesParameters::MissingAntennaBehavior
      missing_ant_behavior_;
  /// Non-owning pointer to a solution table in the H5ParmSingle object.
  schaapcommon::h5parm::SolTab* sol_tab_;
  /// Non-owning pointer to a solution table in the H5ParmSingle object.
  schaapcommon::h5parm::SolTab* sol_tab2_;
  std::vector<casacore::String> parm_expressions_;
  unsigned int timeslots_per_parmupdate_;
  unsigned int timestep_;
  double time_interval_;
  double time_last_;
  std::vector<std::string> soltab_names_;

  std::vector<double> params_;

  base::Direction phase_ref_;

  common::NSTimer timer_;

  IOData iodata_;
  BeamData runtime_beam_data_;
  BeamDataSingle* beam_data_{nullptr};
  int beam_mode{DOBEAM_NONE};
  std::vector<int> ignore_list_;
};

}  // namespace steps
}  // namespace dp3

#endif  // header guard
