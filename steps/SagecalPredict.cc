// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <stddef.h>
#include <string>
#include <sstream>
#include <utility>
#include <vector>
#include <iostream>

#include "../common/ParameterSet.h"
#include "../common/Timer.h"
#include "../common/StreamUtil.h"

#include "../parmdb/ParmDBMeta.h"
#include "../parmdb/PatchInfo.h"
#include "../parmdb/SkymodelToSourceDB.h"
#include "../parmdb/SourceDB.h"

#include <dp3/base/DPInfo.h>
#include <dp3/base/DPBuffer.h>
#include "../base/FlagCounter.h"
#include "../base/Simulate.h"
#include "../base/Simulator.h"
#include "../base/Stokes.h"
#include "../base/PointSource.h"
#include "../base/GaussianSource.h"
#include "../base/Telescope.h"
#include "../base/ComponentInfo.h"
#include "../base/SkyModelCache.h"

#include <casacore/casa/Arrays/Cube.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/Precession.h>
#include <casacore/measures/Measures/Nutation.h>
#include <casacore/casa/Quanta/Quantum.h>

#include <schaapcommon/h5parm/h5parm.h>

#include <aocommon/threadpool.h>

#include "SagecalPredict.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::operator<<;

using schaapcommon::h5parm::H5Parm;
using schaapcommon::h5parm::SolTab;

namespace dp3 {
namespace steps {

SagecalPredict::SagecalPredict(const common::ParameterSet& parset,
                               const std::string& prefix, MsType input_type)
    : ms_type_(input_type),
      name_(prefix),
      operation_(Operation::kReplace),
      h5_name_(parset.getString(prefix + "applycal.parmdb", "")),
      directions_list_(parset.getStringVector(prefix + "directions",
                                              std::vector<std::string>())),
      solset_name_(parset.getString(prefix + "applycal.solset", "")),
      soltab_name_(parset.getString(prefix + "applycal.correction", "")),
      invert_(false),
      parm_on_disk_(!h5_name_.empty()),
      use_amp_phase_(false),
      sigma_mmse_(0.0),
      interp_type_(JonesParameters::InterpolationType::NEAREST),
      timeslots_per_parmupdate_(0),
      timestep_(0),
      time_interval_(-1),
      time_last_(-1) {
  init(parset, prefix, std::vector<std::string>());
}

SagecalPredict::SagecalPredict(const common::ParameterSet& parset,
                               const std::string& prefix,
                               const std::vector<std::string>& source_patterns,
                               MsType input_type)
    : ms_type_(input_type),
      name_(prefix),
      operation_(Operation::kReplace),
      h5_name_(parset.getString(prefix + "applycal.parmdb", "")),
      directions_list_(std::vector<std::string>()),
      solset_name_(parset.getString(prefix + "applycal.solset", "")),
      soltab_name_(parset.getString(prefix + "applycal.correction", "")),
      invert_(false),
      parm_on_disk_(!h5_name_.empty()),
      use_amp_phase_(false),
      sigma_mmse_(0.0),
      interp_type_(JonesParameters::InterpolationType::NEAREST),
      timeslots_per_parmupdate_(0),
      timestep_(0),
      time_interval_(-1),
      time_last_(-1) {
  init(parset, prefix, source_patterns);
}

SagecalPredict::~SagecalPredict() {
#if defined(HAVE_LIBDIRAC) || defined(HAVE_LIBDIRAC_CUDA)

  for (size_t patch_index = 0; patch_index < patch_list_.size();
       ++patch_index) {
    delete[] iodata_.cluster_arr_[patch_index].ll;
    delete[] iodata_.cluster_arr_[patch_index].mm;
    delete[] iodata_.cluster_arr_[patch_index].nn;
    delete[] iodata_.cluster_arr_[patch_index].sI;
    delete[] iodata_.cluster_arr_[patch_index].sQ;
    delete[] iodata_.cluster_arr_[patch_index].sU;
    delete[] iodata_.cluster_arr_[patch_index].sV;
    // Extra data
    for (int ci = 0; ci < iodata_.cluster_arr_[patch_index].N; ci++) {
      if (iodata_.cluster_arr_[patch_index].stype[ci] == STYPE_GAUSSIAN) {
        exinfo_gaussian* exg = static_cast<exinfo_gaussian*>(
            iodata_.cluster_arr_[patch_index].ex[ci]);
        delete exg;
      }
    }
    delete[] iodata_.cluster_arr_[patch_index].stype;
    delete[] iodata_.cluster_arr_[patch_index].ex;
    delete[] iodata_.cluster_arr_[patch_index].sI0;
    delete[] iodata_.cluster_arr_[patch_index].sQ0;
    delete[] iodata_.cluster_arr_[patch_index].sU0;
    delete[] iodata_.cluster_arr_[patch_index].sV0;
    delete[] iodata_.cluster_arr_[patch_index].f0;
    delete[] iodata_.cluster_arr_[patch_index].spec_idx;
    delete[] iodata_.cluster_arr_[patch_index].spec_idx1;
    delete[] iodata_.cluster_arr_[patch_index].spec_idx2;
    delete[] iodata_.cluster_arr_[patch_index].ra;
    delete[] iodata_.cluster_arr_[patch_index].dec;
    delete[] iodata_.cluster_arr_[patch_index].p;
  }

#endif /* HAVE_LIBDIRAC || HAVE_LIBDIRAC_CUDA */
}

#if defined(HAVE_LIBDIRAC) || defined(HAVE_LIBDIRAC_CUDA)
void SagecalPredict::BeamDataSingle::update_metadata(
    const DPInfo& _info, const double freq_f0, const size_t n_channels,
    std::vector<double>& freq_channel, const int beam_mode) {
  std::lock_guard<std::mutex> lock(mutex_);
  if (is_updated_) {
    return;
  }
  beam_mode_ = beam_mode;
  n_stations_ = _info.nantenna();
  casacore::Table ms(_info.msName());
  is_dipole_ = false;
  try {
    n_elem_.resize(n_stations_);
    sx_.resize(n_stations_);
    sy_.resize(n_stations_);
    sz_.resize(n_stations_);
    xx_.resize(n_stations_);
    yy_.resize(n_stations_);
    zz_.resize(n_stations_);
  } catch (const std::bad_alloc& e) {
    throw std::runtime_error("Allocating memory failure");
  }

  if (!ms.keywordSet().isDefined("LOFAR_ANTENNA_FIELD"))
    throw std::runtime_error("subtable LOFAR_ANTENNA_FIELD is missing");
  casacore::Table antfield(ms.keywordSet().asTable("LOFAR_ANTENNA_FIELD"));
  casacore::ROArrayColumn<double> position(antfield, "POSITION");
  casacore::ROArrayColumn<double> offset(antfield, "ELEMENT_OFFSET");
  casacore::ROArrayColumn<double> coord(antfield, "COORDINATE_AXES");
  casacore::ROArrayColumn<bool> eflag(antfield, "ELEMENT_FLAG");
  casacore::ROArrayColumn<double> tileoffset(antfield, "TILE_ELEMENT_OFFSET");
  casacore::Table _field(ms.keywordSet().asTable("FIELD"));
  // check if TILE_ELEMENT_OFFSET has any rows, if no rows present,
  //      we know this is LBA
  bool isHBA = tileoffset.hasContent(0);
  /* read positions, also setup memory for element coords */
  for (size_t ci = 0; ci < n_stations_; ci++) {
    casacore::Array<double> _pos = position(ci);
    double* tx = _pos.data();
    sz_[ci] = tx[2];

    casacore::MPosition stnpos(casacore::MVPosition(tx[0], tx[1], tx[2]),
                               casacore::MPosition::ITRF);
    casacore::Array<double> _radpos = stnpos.getAngle("rad").getValue();
    tx = _radpos.data();

    sx_[ci] = tx[0];
    sy_[ci] = tx[1];
    n_elem_[ci] = offset.shape(ci)[1];
  }

  if (isHBA) {
    beamformer_type_ = STAT_TILE;
    double tempT[3 * HBA_TILE_SIZE];
    /* now read in element offsets, also transform them to local coordinates */
    for (size_t ci = 0; ci < n_stations_; ci++) {
      casacore::Array<double> _off = offset(ci);
      double* off = _off.data();
      casacore::Array<double> _coord = coord(ci);
      double* coordmat = _coord.data();
      casacore::Array<bool> _eflag = eflag(ci);
      bool* ef = _eflag.data();
      casacore::Array<double> _toff = tileoffset(ci);
      double* toff = _toff.data();

      double* tempC = new double[3 * n_elem_[ci]];
      my_dgemm('T', 'N', n_elem_[ci], 3, 3, 1.0, off, 3, coordmat, 3, 0.0,
               tempC, n_elem_[ci]);
      my_dgemm('T', 'N', HBA_TILE_SIZE, 3, 3, 1.0, toff, 3, coordmat, 3, 0.0,
               tempT, HBA_TILE_SIZE);

      /* now inspect the element flag table to see if any of the dipoles are
       * flagged */
      size_t fcount = 0;
      for (int cj = 0; cj < n_elem_[ci]; cj++) {
        if (ef[2 * cj] == 1 || ef[2 * cj + 1] == 1) {
          fcount++;
        }
      }

      xx_[ci] = new double[HBA_TILE_SIZE + (n_elem_[ci] - fcount)];
      yy_[ci] = new double[HBA_TILE_SIZE + (n_elem_[ci] - fcount)];
      zz_[ci] = new double[HBA_TILE_SIZE + (n_elem_[ci] - fcount)];
      my_dcopy(HBA_TILE_SIZE, &tempT[0], 1, &(xx_[ci][0]), 1);
      my_dcopy(HBA_TILE_SIZE, &tempT[HBA_TILE_SIZE], 1, &(yy_[ci][0]), 1);
      my_dcopy(HBA_TILE_SIZE, &tempT[2 * HBA_TILE_SIZE], 1, &(zz_[ci][0]), 1);
      /* copy unflagged tile centroids */
      fcount = 0;
      for (int cj = 0; cj < n_elem_[ci]; cj++) {
        if (!(ef[2 * cj] == 1 || ef[2 * cj + 1] == 1)) {
          xx_[ci][HBA_TILE_SIZE + fcount] = tempC[cj];
          yy_[ci][HBA_TILE_SIZE + fcount] = tempC[cj + n_elem_[ci]];
          zz_[ci][HBA_TILE_SIZE + fcount] = tempC[cj + 2 * n_elem_[ci]];
          fcount++;
        }
      }
      n_elem_[ci] = fcount;

      delete[] tempC;
    }
  } else { /* LBA */
    beamformer_type_ = STAT_SINGLE;
    /* read in element offsets, also transform them to local coordinates */
    for (size_t ci = 0; ci < n_stations_; ci++) {
      casacore::Array<double> _off = offset(ci);
      double* off = _off.data();
      casacore::Array<double> _coord = coord(ci);
      double* coordmat = _coord.data();
      casacore::Array<bool> _eflag = eflag(ci);
      bool* ef = _eflag.data();

      double* tempC = new double[3 * n_elem_[ci]];
      my_dgemm('T', 'N', n_elem_[ci], 3, 3, 1.0, off, 3, coordmat, 3, 0.0,
               tempC, n_elem_[ci]);

      /* now inspect the element flag table to see if any of the dipoles are
       * flagged */
      size_t fcount = 0;
      for (int cj = 0; cj < n_elem_[ci]; cj++) {
        if (ef[2 * cj] == 1 || ef[2 * cj + 1] == 1) {
          fcount++;
        }
      }

      xx_[ci] = new double[(n_elem_[ci] - fcount)];
      yy_[ci] = new double[(n_elem_[ci] - fcount)];
      zz_[ci] = new double[(n_elem_[ci] - fcount)];
      /* copy unflagged coords for each dipole */
      fcount = 0;
      for (int cj = 0; cj < n_elem_[ci]; cj++) {
        if (!(ef[2 * cj] == 1 || ef[2 * cj + 1] == 1)) {
          xx_[ci][fcount] = tempC[cj];
          yy_[ci][fcount] = tempC[cj + n_elem_[ci]];
          zz_[ci][fcount] = tempC[cj + 2 * n_elem_[ci]];
          fcount++;
        }
      }
      n_elem_[ci] = fcount;
      delete[] tempC;
    }
  }
  /* read beam pointing direction */
  casacore::ROArrayColumn<double> point_dir(_field, "REFERENCE_DIR");
  casacore::Array<double> pdir = point_dir(0);
  double* pc = pdir.data();
  p_ra0_ = pc[0];
  p_dec0_ = pc[1];
  /* read tile beam pointing direction */
  casacore::ROArrayColumn<double> tile_dir(_field, "LOFAR_TILE_BEAM_DIR");
  casacore::Array<double> tdir = tile_dir(0);
  double* tc = tdir.data();
  b_ra0_ = tc[0];
  b_dec0_ = tc[1];

  if (beam_mode_ == DOBEAM_FULL || beam_mode_ == DOBEAM_ELEMENT) {
    set_elementcoeffs((freq_f0 < 100e6 ? ELEM_LBA : ELEM_HBA), freq_f0,
                      &ecoeff);
  } else if (beam_mode_ == DOBEAM_FULL_WB || beam_mode_ == DOBEAM_ELEMENT_WB) {
    set_elementcoeffs_wb((freq_f0 < 100e6 ? ELEM_LBA : ELEM_HBA),
                         freq_channel.data(), n_channels, &ecoeff);
  }

  is_updated_ = true;
}

SagecalPredict::BeamDataSingle::~BeamDataSingle() {
  for (size_t ci = 0; ci < n_stations_; ci++) {
    delete[] xx_[ci];
    delete[] yy_[ci];
    delete[] zz_[ci];
  }
  if (beam_mode_ == DOBEAM_FULL || beam_mode_ == DOBEAM_ELEMENT ||
      beam_mode_ == DOBEAM_FULL_WB || beam_mode_ == DOBEAM_ELEMENT_WB) {
    free_elementcoeffs(ecoeff);
  }
}
#endif /* HAVE_LIBDIRAC || HAVE_LIBDIRAC_CUDA */

schaapcommon::h5parm::H5Parm& SagecalPredict::H5ParmSingle::open_file(
    const std::string& h5_name) {
  std::lock_guard<std::mutex> lock(mutex_);
  if (is_opened_) {
    return h5_parm_;
  }

  h5_parm_ = H5Parm(h5_name, false);

  is_opened_ = true;

  return h5_parm_;
}

void SagecalPredict::SetOperation(const std::string& operation) {
  if (operation == "replace") {
    operation_ = Operation::kReplace;
  } else if (operation == "add") {
    operation_ = Operation::kAdd;
  } else if (operation == "subtract") {
    operation_ = Operation::kSubtract;
  } else {
    throw std::invalid_argument("Invalid operation " + operation);
  }
}

unsigned int SagecalPredict::nPol(const std::string& parmName) {
  if (!sol_tab_.HasAxis("pol")) {
    return 1;
  } else {
    return sol_tab_.GetAxis("pol").size;
  }
}

void SagecalPredict::setCorrectType(std::vector<std::string>& solTabs) {
  if (soltab_name_ == "fulljones") {
    if (solTabs.size() != 2) {
      throw std::runtime_error(
          "The soltab parameter requires two soltabs for fulljones "
          "correction (amplitude and phase)");
    }
    sol_tab_ = h5_parm_.GetSolTab(solTabs[0]);
    sol_tab2_ = h5_parm_.GetSolTab(solTabs[1]);
    soltab_name_ = solTabs[0] + "," + solTabs[1];
    gain_type_ = GainType::kFullJones;
  } else if (soltab_name_ == "gain") {
    sol_tab_ = h5_parm_.GetSolTab(solTabs[0]);
    if (solTabs.size() == 2) {
      sol_tab2_ = h5_parm_.GetSolTab(solTabs[1]);
      soltab_name_ = solTabs[0] + "," + solTabs[1];
      gain_type_ = GainType::kDiagonalComplex;
    } else {
      soltab_name_ = solTabs[0];
      gain_type_ =
          JonesParameters::H5ParmTypeStringToGainType(sol_tab_.GetType());
    }
  } else {
    sol_tab_ = h5_parm_.GetSolTab(soltab_name_);
    gain_type_ =
        JonesParameters::H5ParmTypeStringToGainType(sol_tab_.GetType());
  }
  if (gain_type_ == GainType::kDiagonalPhase && nPol("") == 1) {
    gain_type_ = GainType::kScalarPhase;
  }
  if (gain_type_ == GainType::kDiagonalAmplitude && nPol("") == 1) {
    gain_type_ = GainType::kScalarAmplitude;
  }
}

void SagecalPredict::init(
    const common::ParameterSet& parset, const std::string& prefix,
    const std::vector<std::string>& given_source_patterns) {
  if (ms_type_ == MsType::kBda) {
    throw std::runtime_error("Does not support BDA data");
  }
  // Override 'applycal.correction' if substeps like 'applycal.steps=[step1,..]'
  // and 'applycal.step1.correction=.. , applycal.step2.correction=..
  // are defined.
  // Try to find 'applycal.stepname.correction'
  // where 'stepname' given as 'applycal.steps=[*]'
  if (parset.isDefined(prefix + "applycal.steps")) {
    common::ParameterValue step_names(
        parset.getString(prefix + "applycal.steps"));
    std::vector<std::string> substep_names;
    if (step_names.isVector()) {
      substep_names = step_names.getStringVector();
    } else {
      substep_names.push_back(step_names.getString());
    }
    // We can only handle one or two substeps (amplitude and/or phase)
    if (substep_names.size() > 2) {
      throw std::runtime_error(
          "Only one or two types of gain application can be done at a time"
          "applycal.steps should specify only one or two steps");
    }
    // update soltab_names_ from substeps
    if (!substep_names.empty()) {
      for (auto substep_name : substep_names) {
        if (!substep_name.empty()) {
          auto substep_soltab_name = parset.getString(
              prefix + "applycal." + substep_name + ".correction");
          soltab_names_.push_back(substep_soltab_name);
        }
      }
      // update soltab_name_
      soltab_name_ = "gain";
    }
  }

  std::vector<std::string> source_patterns;
  // Order of precedence for getting the sources list:
  // 1) parset.directions key, 2) h5parm directions, 3) parset.sources key
  // withing DDECal, 3) will be invoked because source_patterns in the
  // constructor
  if (!given_source_patterns.empty()) {
    source_patterns = given_source_patterns;
  } else {
    source_patterns =
        parset.getStringVector(prefix + "sources", std::vector<std::string>());
  }
  if (!source_patterns.empty() && !directions_list_.empty()) {
    throw std::runtime_error(
        "Both sources and directions specified, use only one of them");
  }

  // Try to get HDF5 file if any
  if (!parm_on_disk_ && !directions_list_.empty()) {
    throw std::runtime_error(
        "H5parm name not specified but directions specified, give both of "
        "them");
  }
  std::vector<std::string> h5directions;
  if (parm_on_disk_) {
    // Create a singleton to only open and read the H5 file once
    h5_parm_reference_ = SagecalPredict::H5ParmSingle::get_instance();
    h5_parm_ = h5_parm_reference_->open_file(h5_name_);
    // Check to see soltab is initialized at the constructor
    if (soltab_name_.empty()) {
      soltab_name_ = parset.getString(prefix + "applycal.correction");
    }
    // soltab_names_ will be already updated if 'substeps' are defined
    if (soltab_names_.empty()) {
      soltab_names_ = parset.getStringVector(
          prefix + "applycal.soltab",
          std::vector<std::string>{"amplitude000", "phase000"});
    }
    setCorrectType(soltab_names_);
    h5directions = sol_tab_.GetStringAxis("dir");
    if (h5directions.empty())
      throw std::runtime_error("H5Parm has empty dir axis");
    // Also do sanity check for time,freq axes
    if (sol_tab_.HasAxis("freq") && sol_tab_.HasAxis("time") &&
        sol_tab_.GetAxisIndex("freq") < sol_tab_.GetAxisIndex("time"))
      throw std::runtime_error("H5Parm fastest varying axis should be freq");
  }
  missing_ant_behavior_ = JonesParameters::StringToMissingAntennaBehavior(
      parset.getString(prefix + "applycal.missingantennabehavior", "error"));

  if (!h5directions.empty() && directions_list_.empty()) {
    directions_list_ = h5directions;
  }

  if (!directions_list_.empty() && source_patterns.empty()) {
    // Build source patterns from the given directions list
    // only if source_patterns is empty
    for (size_t ci = 0; ci < directions_list_.size(); ci++) {
      std::string dir = directions_list_[ci];
      std::vector<std::string> dir_vec;
      if (dir.size() <= 2 || dir[0] != '[' || dir[dir.size() - 1] != ']')
        throw std::runtime_error(
            "Invalid direction string: expecting array between square "
            "brackets, "
            "e.g. [a, b, c, ...]");
      dir_vec =
          common::stringtools::tokenize(dir.substr(1, dir.size() - 2), ",");
      for (auto source : dir_vec) {
        // only include source if it is not already in the source_patterns
        if (find(source_patterns.begin(), source_patterns.end(), source) ==
            source_patterns.end()) {
          source_patterns.push_back(source);
        }
      }
    }
  }

  // Reconcile source patterns and h5parm directions, if h5parm is given
  if (parm_on_disk_) {
    for (auto dir : source_patterns) {
      if (find(h5directions.begin(), h5directions.end(), "[" + dir + "]") ==
          h5directions.end()) {
        throw std::runtime_error("Direction " + dir + " not found in " +
                                 h5_name_);
      }
    }
  }

  name_ = prefix;
  SetOperation(parset.getString(prefix + "operation", "replace"));
  source_db_name_ = parset.getString(prefix + "sourcedb");

  patch_list_.clear();

  base::SourceDBWrapper source_db =
      base::SkyModelCache::GetInstance().GetSkyModel(source_db_name_);
  source_db.Filter(source_patterns,
                   base::SourceDBWrapper::FilterMode::kPattern);

  try {
    patch_list_ = source_db.MakePatchList();
    if (patch_list_.empty()) {
      std::stringstream ss;
      ss << source_patterns;
      throw std::runtime_error("Couldn't find patch for directions " +
                               ss.str());
    }
  } catch (std::exception& exception) {
    throw std::runtime_error(std::string("Something went wrong while reading "
                                         "the source model. The error was: ") +
                             exception.what());
  }

  source_list_ = makeSourceList(patch_list_);
  any_orientation_is_absolute_ = source_db.CheckAnyOrientationIsAbsolute();
#if defined(HAVE_LIBDIRAC) || defined(HAVE_LIBDIRAC_CUDA)
  const bool apply_beam = parset.getBool(prefix + "usebeammodel", false);
  const bool use_channel_freq = parset.getBool(prefix + "usechannelfreq", true);
  beam_mode_ = DOBEAM_NONE;
  if (apply_beam) {
    const std::string& beam_mode_string =
        parset.getString(prefix + "beammode", "full");
    if (beam_mode_string == "full" || beam_mode_string == "default") {
      beam_mode_ = use_channel_freq ? DOBEAM_FULL_WB : DOBEAM_FULL;
    } else if (beam_mode_string == "arrayfactor" ||
               beam_mode_string == "array_factor") {
      beam_mode_ = use_channel_freq ? DOBEAM_ARRAY_WB : DOBEAM_ARRAY;
    } else if (beam_mode_string == "element") {
      beam_mode_ = use_channel_freq ? DOBEAM_ELEMENT_WB : DOBEAM_ELEMENT;
    }
  }
  if (parm_on_disk_) {
    ignore_list_.resize(patch_list_.size());
    memset(ignore_list_.data(), 0, sizeof(int) * patch_list_.size());
  }
  // Create singleton if apply_beam is set
  beam_data_ = SagecalPredict::BeamDataSingle::get_instance();
#endif /* HAVE_LIBDIRAC || HAVE_LIBDIRAC_CUDA */
}

void SagecalPredict::updateFromH5(const double startTime) {
  time_last_ = startTime + time_interval_ * (-0.5 + timeslots_per_parmupdate_);

  const double lastMSTime =
      info().startTime() + info().ntime() * time_interval_;
  if (time_last_ > lastMSTime &&
      !casacore::nearAbs(time_last_, lastMSTime, 1e-3)) {
    time_last_ = lastMSTime;
  }

  std::vector<double> times(info().ntime());
#pragma GCC ivdep
  for (size_t t = 0; t < times.size(); ++t) {
    times[t] = info().startTime() + (t + 0.5) * info().timeInterval();
  }

  const size_t n_chan = info().chanFreqs().size();
  const size_t n_dir = directions_list_.size();

  size_t dir_index = 0;
  for (auto dir : directions_list_) {
    hsize_t direction_index = sol_tab_.GetDirIndex(dir);
    std::unique_ptr<JonesParameters> Jones_params_ =
        std::make_unique<JonesParameters>(
            info().chanFreqs(), times, info().antennaNames(), gain_type_,
            interp_type_, direction_index, &sol_tab_, &sol_tab2_, invert_,
            sigma_mmse_, parm_expressions_.size(), missing_ant_behavior_);

    // shape : ncorr x n_stations x ntime*nfreq (ncorr:2,4)
    const casacore::Cube<casacore::Complex>& gains = Jones_params_->GetParms();
    const size_t n_corr = gains.shape()[0];
    const size_t n_stat = gains.shape()[1];
    const size_t n_timefreq = gains.shape()[2];
    const size_t n_time = n_timefreq / n_chan;

    if (!(n_corr == 2 || n_corr == 4)) {
      throw std::runtime_error("H5 solutions correlations has to be 2 or 4.");
    }

    // If not allocated, setup memory
    if (params_.size() == 0) {
      // storage: 8*n_dir*n_stat*n_time*n_chan
      params_.resize(8 * n_dir * n_stat * n_time * n_chan);
      memset(params_.data(), 0,
             sizeof(double) * 8 * n_dir * n_stat * n_time * n_chan);
    } else if (!dir_index) {  // if already allocated, reset at the first dir
      memset(params_.data(), 0,
             sizeof(double) * 8 * n_dir * n_stat * n_time * n_chan);
    }

    // Copy to params
    for (size_t chan = 0; chan < n_chan; chan++) {
      for (size_t t = 0; t < n_time; t++) {
        const size_t time_freq_offset = t * n_chan + chan;
        if (n_corr == 2) {
#pragma GCC ivdep
          for (size_t ant = 0; ant < n_stat; ant++) {
            const casacore::Complex* xx = &gains(0, ant, time_freq_offset);
            const casacore::Complex* yy = &gains(1, ant, time_freq_offset);
            params_[chan * 8 * n_dir * n_stat * n_time +
                    t * 8 * n_dir * n_stat + dir_index * 8 * n_stat + 8 * ant] =
                xx->real();
            params_[chan * 8 * n_dir * n_stat * n_time +
                    t * 8 * n_dir * n_stat + dir_index * 8 * n_stat + 8 * ant +
                    1] = xx->imag();
            params_[chan * 8 * n_dir * n_stat * n_time +
                    t * 8 * n_dir * n_stat + dir_index * 8 * n_stat + 8 * ant +
                    6] = yy->real();
            params_[chan * 8 * n_dir * n_stat * n_time +
                    t * 8 * n_dir * n_stat + dir_index * 8 * n_stat + 8 * ant +
                    7] = yy->imag();
          }
        } else if (n_corr == 4) {
#pragma GCC ivdep
          for (size_t ant = 0; ant < n_stat; ant++) {
            const casacore::Complex* xx = &gains(0, ant, time_freq_offset);
            const casacore::Complex* xy = &gains(1, ant, time_freq_offset);
            const casacore::Complex* yx = &gains(2, ant, time_freq_offset);
            const casacore::Complex* yy = &gains(3, ant, time_freq_offset);
            params_[chan * 8 * n_dir * n_stat * n_time +
                    t * 8 * n_dir * n_stat + dir_index * 8 * n_stat + 8 * ant] =
                xx->real();
            params_[chan * 8 * n_dir * n_stat * n_time +
                    t * 8 * n_dir * n_stat + dir_index * 8 * n_stat + 8 * ant +
                    1] = xx->imag();
            params_[chan * 8 * n_dir * n_stat * n_time +
                    t * 8 * n_dir * n_stat + dir_index * 8 * n_stat + 8 * ant +
                    2] = xy->real();
            params_[chan * 8 * n_dir * n_stat * n_time +
                    t * 8 * n_dir * n_stat + dir_index * 8 * n_stat + 8 * ant +
                    3] = xy->imag();
            params_[chan * 8 * n_dir * n_stat * n_time +
                    t * 8 * n_dir * n_stat + dir_index * 8 * n_stat + 8 * ant +
                    4] = yx->real();
            params_[chan * 8 * n_dir * n_stat * n_time +
                    t * 8 * n_dir * n_stat + dir_index * 8 * n_stat + 8 * ant +
                    5] = yx->imag();
            params_[chan * 8 * n_dir * n_stat * n_time +
                    t * 8 * n_dir * n_stat + dir_index * 8 * n_stat + 8 * ant +
                    6] = yy->real();
            params_[chan * 8 * n_dir * n_stat * n_time +
                    t * 8 * n_dir * n_stat + dir_index * 8 * n_stat + 8 * ant +
                    7] = yy->imag();
          }
        }
      }
    }
    dir_index++;
  }
}

bool SagecalPredict::process(std::unique_ptr<DPBuffer> buffer) {
  timer_.start();
#if defined(HAVE_LIBDIRAC) || defined(HAVE_LIBDIRAC_CUDA)
  // Determine the various sizes.
  const size_t nDr = patch_list_.size();
  const size_t nBl = info().nbaselines();
  const size_t nCh = info().nchan();
  const size_t nCr = info().ncorr();

  buffer->GetData().resize({nBl, nCh, nCr});

  if (parm_on_disk_ && (buffer->GetTime() > time_last_)) {
    timestep_ = 0;
    // Update solutions (also aux info such as time_last_ ...)
    updateFromH5(buffer->GetTime());
  } else {
    timestep_++;
  }
  loadData(buffer);

  // update time
  if (beam_mode_ != DOBEAM_NONE) {
    runtime_beam_data_.time_utc_[0] = (buffer->GetTime() / 86400.0 + 2400000.5);
    // precess source positions
    if (!runtime_beam_data_.sources_prec_) {
      casacore::Precession prec(casacore::Precession::IAU2000);
      casacore::RotMatrix rotat_prec(
          prec(runtime_beam_data_.time_utc_[0] - 2400000.5));
      casacore::Nutation nut(casacore::Nutation::IAU2000);
      casacore::RotMatrix rotat_nut(
          nut(runtime_beam_data_.time_utc_[0] - 2400000.5));

      casacore::RotMatrix rotat = rotat_prec * rotat_nut;
      rotat.transpose();
      for (size_t cl = 0; cl < nDr; cl++) {
#pragma GCC ivdep
        for (int ci = 0; ci < iodata_.cluster_arr_[cl].N; ci++) {
          casacore::MVDirection pos(
              casacore::Quantity(iodata_.cluster_arr_[cl].ra[ci], "rad"),
              casacore::Quantity(iodata_.cluster_arr_[cl].dec[ci], "rad"));
          casacore::MVDirection newdir = rotat * pos;
          iodata_.cluster_arr_[cl].ra[ci] = newdir.get()[0];
          iodata_.cluster_arr_[cl].dec[ci] = newdir.get()[1];
        }
      }
      casacore::MVDirection pos(
          casacore::Quantity(runtime_beam_data_.p_ra0_, "rad"),
          casacore::Quantity(runtime_beam_data_.p_dec0_, "rad"));
      casacore::MVDirection newdir = rotat * pos;
      runtime_beam_data_.p_ra0_ = newdir.get()[0];
      runtime_beam_data_.p_dec0_ = newdir.get()[1];
      casacore::MVDirection pos_tile(
          casacore::Quantity(runtime_beam_data_.b_ra0_, "rad"),
          casacore::Quantity(runtime_beam_data_.b_dec0_, "rad"));
      casacore::MVDirection newdir_tile = rotat * pos_tile;
      runtime_beam_data_.b_ra0_ = newdir_tile.get()[0];
      runtime_beam_data_.b_dec0_ = newdir_tile.get()[1];
      runtime_beam_data_.sources_prec_ = true;
    }
  }

  int operation = SIMUL_ONLY;
  if (operation_ == Operation::kAdd) {
    operation = SIMUL_ADD;
  } else if (operation_ == Operation::kSubtract) {
    operation = SIMUL_SUB;
  }
  const int tile_size = 1;
  const double time_smear_factor = 1.0;
  const size_t n_threads = aocommon::ThreadPool::GetInstance().NThreads();
#ifdef HAVE_LIBDIRAC /* mutually exclusive with HAVE_LIBDIRAC_CUDA */
  if (!parm_on_disk_) {
    if (beam_mode_ == DOBEAM_NONE) {
      predict_visibilities_multifreq(
          iodata_.u_.data(), iodata_.v_.data(), iodata_.w_.data(),
          iodata_.data_.data(), iodata_.n_stations, iodata_.n_baselines,
          tile_size, iodata_.baseline_arr_.data(), iodata_.cluster_arr_.data(),
          nDr, iodata_.freqs_.data(), iodata_.n_channels, iodata_.fdelta,
          time_smear_factor, phase_ref_.dec, n_threads, operation);
    } else {
      predict_visibilities_multifreq_withbeam(
          iodata_.u_.data(), iodata_.v_.data(), iodata_.w_.data(),
          iodata_.data_.data(), iodata_.n_stations, iodata_.n_baselines,
          tile_size, iodata_.baseline_arr_.data(), iodata_.cluster_arr_.data(),
          nDr, iodata_.freqs_.data(), iodata_.n_channels, iodata_.fdelta,
          time_smear_factor, phase_ref_.dec, beam_data_->beamformer_type_,
          runtime_beam_data_.b_ra0_, runtime_beam_data_.b_dec0_,
          runtime_beam_data_.p_ra0_, runtime_beam_data_.p_dec0_, iodata_.f0,
          beam_data_->sx_.data(), beam_data_->sy_.data(),
          runtime_beam_data_.time_utc_.data(), beam_data_->n_elem_.data(),
          beam_data_->xx_.data(), beam_data_->yy_.data(),
          beam_data_->zz_.data(), &beam_data_->ecoeff, beam_mode_, n_threads,
          operation);
    }
  } else {
    if (beam_mode_ == DOBEAM_NONE) {
      predict_visibilities_multifreq_withsol(
          iodata_.u_.data(), iodata_.v_.data(), iodata_.w_.data(),
          &params_[timestep_ * 8 * nDr * nCh * info().nantenna()],
          iodata_.data_.data(), ignore_list_.data(), iodata_.n_stations,
          iodata_.n_baselines, tile_size, iodata_.baseline_arr_.data(),
          iodata_.cluster_arr_.data(), nDr, iodata_.freqs_.data(),
          iodata_.n_channels, iodata_.fdelta, time_smear_factor, phase_ref_.dec,
          n_threads, operation, -1, 0.0, false);

    } else {
      predict_visibilities_multifreq_withsol_withbeam(
          iodata_.u_.data(), iodata_.v_.data(), iodata_.w_.data(),
          &params_[timestep_ * 8 * nDr * nCh * info().nantenna()],
          iodata_.data_.data(), ignore_list_.data(), iodata_.n_stations,
          iodata_.n_baselines, tile_size, iodata_.baseline_arr_.data(),
          iodata_.cluster_arr_.data(), nDr, iodata_.freqs_.data(),
          iodata_.n_channels, iodata_.fdelta, time_smear_factor, phase_ref_.dec,
          beam_data_->beamformer_type_, runtime_beam_data_.b_ra0_,
          runtime_beam_data_.b_dec0_, runtime_beam_data_.p_ra0_,
          runtime_beam_data_.p_dec0_, iodata_.f0, beam_data_->sx_.data(),
          beam_data_->sy_.data(), runtime_beam_data_.time_utc_.data(),
          beam_data_->n_elem_.data(), beam_data_->xx_.data(),
          beam_data_->yy_.data(), beam_data_->zz_.data(), &beam_data_->ecoeff,
          beam_mode_, n_threads, operation, -1, 0.0, false);
    }
  }
#endif                    /* HAVE_LIBDIRAC */
#ifdef HAVE_LIBDIRAC_CUDA /* mutually exclusive with HAVE_LIBDIRAC */
  if (!parm_on_disk_) {
    predict_visibilities_multifreq_withbeam_gpu(
        iodata_.u_.data(), iodata_.v_.data(), iodata_.w_.data(),
        iodata_.data_.data(), iodata_.n_stations, iodata_.n_baselines,
        tile_size, iodata_.baseline_arr_.data(), iodata_.cluster_arr_.data(),
        nDr, iodata_.freqs_.data(), iodata_.n_channels, iodata_.fdelta,
        time_smear_factor, phase_ref_.dec, beam_data_->beamformer_type_,
        runtime_beam_data_.b_ra0_, runtime_beam_data_.b_dec0_,
        runtime_beam_data_.p_ra0_, runtime_beam_data_.p_dec0_, iodata_.f0,
        beam_data_->sx_.data(), beam_data_->sy_.data(),
        runtime_beam_data_.time_utc_.data(), beam_data_->n_elem_.data(),
        beam_data_->xx_.data(), beam_data_->yy_.data(), beam_data_->zz_.data(),
        &beam_data_->ecoeff, beam_mode_, n_threads, operation);
  } else {
    predict_visibilities_withsol_withbeam_gpu(
        iodata_.u_.data(), iodata_.v_.data(), iodata_.w_.data(),
        &params_[timestep_ * 8 * nDr * nCh * info().nantenna()],
        iodata_.data_.data(), ignore_list_.data(), iodata_.n_stations,
        iodata_.n_baselines, tile_size, iodata_.baseline_arr_.data(),
        iodata_.cluster_arr_.data(), nDr, iodata_.freqs_.data(),
        iodata_.n_channels, iodata_.fdelta, time_smear_factor, phase_ref_.dec,
        beam_data_->beamformer_type_, runtime_beam_data_.b_ra0_,
        runtime_beam_data_.b_dec0_, runtime_beam_data_.p_ra0_,
        runtime_beam_data_.p_dec0_, iodata_.f0, beam_data_->sx_.data(),
        beam_data_->sy_.data(), runtime_beam_data_.time_utc_.data(),
        beam_data_->n_elem_.data(), beam_data_->xx_.data(),
        beam_data_->yy_.data(), beam_data_->zz_.data(), &beam_data_->ecoeff,
        beam_mode_, n_threads, operation, -1, 0.0, false);
  }
#endif /* HAVE_LIBDIRAC_CUDA */
  writeData(buffer);
#endif /* HAVE_LIBDIRAC || HAVE_LIBDIRAC_CUDA */

  timer_.stop();
  getNextStep()->process(std::move(buffer));

  return false;
}

void SagecalPredict::updateInfo(const DPInfo& _info) {
  const size_t n_correlations = _info.ncorr();

  if (n_correlations != 4) {
    throw std::invalid_argument("Can only predict with all 4 correlations.");
  }

  Step::updateInfo(_info);
  try {
    casacore::MDirection dirJ2000(casacore::MDirection::Convert(
        _info.phaseCenter(), casacore::MDirection::J2000)());
    casacore::Quantum<casacore::Vector<double>> angles = dirJ2000.getAngle();
    // casacore::Quantum<casacore::Vector<double>> angles =
    // _info.phaseCenter().getAngle();
    phase_ref_ =
        base::Direction(angles.getBaseValue()[0], angles.getBaseValue()[1]);
  } catch (casacore::AipsError&) {
    throw std::runtime_error("Casacore error.");
  }

#if defined(HAVE_LIBDIRAC) || defined(HAVE_LIBDIRAC_CUDA)
  const size_t n_directions = patch_list_.size();
  const size_t n_channels = _info.nchan();
  const size_t n_stations = _info.nantenna();
  const size_t n_given_baselines = _info.nbaselines();

  // Some data will not have autocorrelations, in that case
  // n_baselines = n_stations*(n_stations-1)/2, otherwise
  // n_baselines = n_stations*(n_stations-1)/2+n_stations
  // we only need baselines without autocorrelations
  size_t n_baselines = n_given_baselines;
  if (n_given_baselines == n_stations * (n_stations - 1) / 2 + n_stations) {
    n_baselines = n_stations * (n_stations - 1) / 2;
  }

  try {
    iodata_.cluster_arr_.resize(n_directions);
    iodata_.baseline_arr_.resize(n_baselines);
    iodata_.data_.resize(n_baselines * n_correlations * n_channels * 2);
    iodata_.u_.resize(n_baselines);
    iodata_.v_.resize(n_baselines);
    iodata_.w_.resize(n_baselines);
    iodata_.freqs_.resize(n_channels);
    iodata_.n_baselines = n_baselines;
    iodata_.n_stations = n_stations;
    iodata_.n_channels = n_channels;
  } catch (const std::bad_alloc& e) {
    throw std::runtime_error("Allocating memory failure");
  }

  for (size_t patch_index = 0; patch_index < n_directions; ++patch_index) {
    // Count sources of this patch
    size_t n_sources =
        std::count_if(source_list_.begin(), source_list_.end(),
                      [&](const std::pair<std::shared_ptr<base::ModelComponent>,
                                          std::shared_ptr<base::Patch>>& item) {
                        return item.second == patch_list_[patch_index];
                      });

    iodata_.cluster_arr_[patch_index].id = patch_index;
    iodata_.cluster_arr_[patch_index].nchunk = 1;
    iodata_.cluster_arr_[patch_index].N = n_sources;
    try {
      iodata_.cluster_arr_[patch_index].ll = new double[n_sources];
      iodata_.cluster_arr_[patch_index].mm = new double[n_sources];
      iodata_.cluster_arr_[patch_index].nn = new double[n_sources];
      iodata_.cluster_arr_[patch_index].sI = new double[n_sources];
      iodata_.cluster_arr_[patch_index].sQ = new double[n_sources];
      iodata_.cluster_arr_[patch_index].sU = new double[n_sources];
      iodata_.cluster_arr_[patch_index].sV = new double[n_sources];
      iodata_.cluster_arr_[patch_index].stype = new unsigned char[n_sources];
      // Extra data (stored as a generic pointer)
      iodata_.cluster_arr_[patch_index].ex = new void*[n_sources];
      iodata_.cluster_arr_[patch_index].sI0 = new double[n_sources];
      iodata_.cluster_arr_[patch_index].sQ0 = new double[n_sources];
      iodata_.cluster_arr_[patch_index].sU0 = new double[n_sources];
      iodata_.cluster_arr_[patch_index].sV0 = new double[n_sources];
      iodata_.cluster_arr_[patch_index].f0 = new double[n_sources];
      iodata_.cluster_arr_[patch_index].spec_idx = new double[n_sources];
      iodata_.cluster_arr_[patch_index].spec_idx1 = new double[n_sources];
      iodata_.cluster_arr_[patch_index].spec_idx2 = new double[n_sources];
      iodata_.cluster_arr_[patch_index].ra = new double[n_sources];
      iodata_.cluster_arr_[patch_index].dec = new double[n_sources];

      iodata_.cluster_arr_[patch_index].p =
          new int[iodata_.cluster_arr_[patch_index].nchunk];
    } catch (const std::bad_alloc& e) {
      throw std::runtime_error("Allocating memory failure");
    }

    dp3::base::ComponentInfo component_info;
    size_t src_index = 0;
    for (size_t scan_index = 0; scan_index < source_list_.size();
         ++scan_index) {
      if (source_list_[scan_index].second == patch_list_[patch_index]) {
        auto source = source_list_[scan_index].first;
        component_info.inspect(source);
        // Insert this source
        iodata_.cluster_arr_[patch_index].ra[src_index] = component_info.ra_;
        iodata_.cluster_arr_[patch_index].dec[src_index] = component_info.dec_;
        iodata_.cluster_arr_[patch_index].sI[src_index] = component_info.sI_;
        iodata_.cluster_arr_[patch_index].sQ[src_index] = component_info.sQ_;
        iodata_.cluster_arr_[patch_index].sU[src_index] = component_info.sU_;
        iodata_.cluster_arr_[patch_index].sV[src_index] = component_info.sV_;
        iodata_.cluster_arr_[patch_index].sI0[src_index] = component_info.sI_;
        iodata_.cluster_arr_[patch_index].sQ0[src_index] = component_info.sQ_;
        iodata_.cluster_arr_[patch_index].sU0[src_index] = component_info.sU_;
        iodata_.cluster_arr_[patch_index].sV0[src_index] = component_info.sV_;
        iodata_.cluster_arr_[patch_index].stype[src_index] =
            (component_info.source_type_ == dp3::base::ComponentInfo::kPoint
                 ? STYPE_POINT
                 : STYPE_GAUSSIAN);

        iodata_.cluster_arr_[patch_index].spec_idx[src_index] =
            component_info.spectrum_[0];
        iodata_.cluster_arr_[patch_index].spec_idx1[src_index] =
            component_info.spectrum_[1] / log(10.0);
        iodata_.cluster_arr_[patch_index].spec_idx2[src_index] =
            component_info.spectrum_[2] / (log(10.0) * log(10.0));
        iodata_.cluster_arr_[patch_index].f0[src_index] = component_info.f0_;
        if (iodata_.cluster_arr_[patch_index].stype[src_index] ==
            STYPE_GAUSSIAN) {
          exinfo_gaussian* exg = new exinfo_gaussian;
          // The following will be updated later
          exg->eX = component_info.g_major_;
          exg->eY = component_info.g_minor_;
          exg->eP = component_info.g_pa_;
          iodata_.cluster_arr_[patch_index].ex[src_index] =
              static_cast<void*>(exg);
        } else {
          iodata_.cluster_arr_[patch_index].ex[src_index] = nullptr;
        }
        src_index++;
      }
    }
    assert(src_index == n_sources);
  }
  const std::vector<double>& freqs = info().chanFreqs();
  const std::vector<double>& widths = info().chanWidths();
  iodata_.f0 = 0.0;
  iodata_.fdelta = 0.0;
  for (size_t ch = 0; ch < n_channels; ch++) {
    iodata_.freqs_[ch] = freqs[ch];
    iodata_.f0 += freqs[ch];
    iodata_.fdelta += widths[ch];
  }
  iodata_.f0 /= (double)n_channels;

  int param_offset = 0;
  for (size_t patch_index = 0; patch_index < n_directions; ++patch_index) {
    for (int chunk = 0; chunk < iodata_.cluster_arr_[patch_index].nchunk;
         chunk++) {
      iodata_.cluster_arr_[patch_index].p[chunk] =
          8 * param_offset * iodata_.n_stations;
      param_offset++;
    }
  }
#endif /* HAVE_LIBDIRAC || HAVE_LIBDIRAC_CUDA */

#if defined(HAVE_LIBDIRAC) || defined(HAVE_LIBDIRAC_CUDA)
  // Update source positions
  for (size_t cl = 0; cl < n_directions; cl++) {
#pragma GCC ivdep
    for (int ci = 0; ci < iodata_.cluster_arr_[cl].N; ci++) {
      const double c_dec = cos(iodata_.cluster_arr_[cl].dec[ci]);
      const double s_dec = sin(iodata_.cluster_arr_[cl].dec[ci]);
      const double c_dec0 = cos(phase_ref_.dec);
      const double s_dec0 = sin(phase_ref_.dec);
      const double c_radiff =
          cos(iodata_.cluster_arr_[cl].ra[ci] - phase_ref_.ra);
      const double s_radiff =
          sin(iodata_.cluster_arr_[cl].ra[ci] - phase_ref_.ra);
      iodata_.cluster_arr_[cl].ll[ci] = c_dec * s_radiff;
      iodata_.cluster_arr_[cl].mm[ci] =
          s_dec * c_dec0 - c_dec * s_dec0 * c_radiff;
      iodata_.cluster_arr_[cl].nn[ci] =
          s_dec * s_dec0 + c_dec * c_dec0 * c_radiff - 1.0;
      if (iodata_.cluster_arr_[cl].stype[ci] == STYPE_GAUSSIAN) {
        exinfo_gaussian* exg =
            static_cast<exinfo_gaussian*>(iodata_.cluster_arr_[cl].ex[ci]);
        exg->eX /= (2.0 * (sqrt(2.0 * log(2.0))));
        exg->eY /= (2.0 * (sqrt(2.0 * log(2.0))));
        exg->eP += M_PI_2;
        if (any_orientation_is_absolute_) {
          const double phi = acos(iodata_.cluster_arr_[cl].nn[ci] + 1.0);
          const double xi = atan2(-iodata_.cluster_arr_[cl].ll[ci],
                                  iodata_.cluster_arr_[cl].mm[ci]);
          /* negate angles */
          exg->cxi = cos(xi);
          exg->sxi = sin(-xi);
          exg->cphi = cos(phi);
          exg->sphi = sin(-phi);
          if (iodata_.cluster_arr_[cl].nn[ci] + 1.0 < PROJ_CUT) {
            /* only then consider projection */
            exg->use_projection = 1;
          } else {
            exg->use_projection = 0;
          }
        } else {
          exg->use_projection = 0;
        }
      }
    }
  }

  // If beam calculation is enabled
  if (beam_mode_ != DOBEAM_NONE) {
    readAuxData(_info);
  }
#endif /* HAVE_LIBDIRAC || HAVE_LIBDIRAC_CUDA */

  // Setup reading the H5 solutions
  if (parm_on_disk_) {
    time_interval_ = _info.timeInterval();
    // Read in solutions for all timeslots
    timeslots_per_parmupdate_ = info().ntime();
    if (gain_type_ == GainType::kDiagonalComplex ||
        gain_type_ == GainType::kFullJones) {
      use_amp_phase_ = true;
    }
    if (gain_type_ == GainType::kDiagonalComplex) {
      if (use_amp_phase_) {
        parm_expressions_.push_back("Gain:0:0:Ampl");
        parm_expressions_.push_back("Gain:0:0:Phase");
        parm_expressions_.push_back("Gain:1:1:Ampl");
        parm_expressions_.push_back("Gain:1:1:Phase");
      } else {
        parm_expressions_.push_back("Gain:0:0:Real");
        parm_expressions_.push_back("Gain:0:0:Imag");
        parm_expressions_.push_back("Gain:1:1:Real");
        parm_expressions_.push_back("Gain:1:1:Imag");
      }
    } else if (gain_type_ == GainType::kFullJones) {
      if (use_amp_phase_) {
        parm_expressions_.push_back("Gain:0:0:Ampl");
        parm_expressions_.push_back("Gain:0:0:Phase");
        parm_expressions_.push_back("Gain:0:1:Ampl");
        parm_expressions_.push_back("Gain:0:1:Phase");
        parm_expressions_.push_back("Gain:1:0:Ampl");
        parm_expressions_.push_back("Gain:1:0:Phase");
        parm_expressions_.push_back("Gain:1:1:Ampl");
        parm_expressions_.push_back("Gain:1:1:Phase");
      } else {
        parm_expressions_.push_back("Gain:0:0:Real");
        parm_expressions_.push_back("Gain:0:0:Imag");
        parm_expressions_.push_back("Gain:0:1:Real");
        parm_expressions_.push_back("Gain:0:1:Imag");
        parm_expressions_.push_back("Gain:1:0:Real");
        parm_expressions_.push_back("Gain:1:0:Imag");
        parm_expressions_.push_back("Gain:1:1:Real");
        parm_expressions_.push_back("Gain:1:1:Imag");
      }
    } else if (gain_type_ == GainType::kTec) {
      if (nPol("TEC") == 1) {
        parm_expressions_.push_back("TEC");
      } else {
        parm_expressions_.push_back("TEC:0");
        parm_expressions_.push_back("TEC:1");
      }
    } else if (gain_type_ == GainType::kClock) {
      if (nPol("Clock") == 1) {
        parm_expressions_.push_back("Clock");
      } else {
        parm_expressions_.push_back("Clock:0");
        parm_expressions_.push_back("Clock:1");
      }
    } else if (gain_type_ == GainType::kRotationAngle) {
      parm_expressions_.push_back("{Common,}RotationAngle");
    } else if (gain_type_ == GainType::kScalarPhase) {
      parm_expressions_.push_back("{Common,}ScalarPhase");
    } else if (gain_type_ == GainType::kRotationMeasure) {
      parm_expressions_.push_back("RotationMeasure");
    } else if (gain_type_ == GainType::kScalarAmplitude) {
      parm_expressions_.push_back("{Common,}ScalarAmplitude");
    } else if (gain_type_ == GainType::kDiagonalPhase) {
      parm_expressions_.push_back("Phase:0");
      parm_expressions_.push_back("Phase:1");
    } else if (gain_type_ == GainType::kDiagonalAmplitude) {
      parm_expressions_.push_back("Amplitude:0");
      parm_expressions_.push_back("Amplitude:1");
    } else {
      throw std::runtime_error(
          "Correction type " +
          JonesParameters::GainTypeToHumanReadableString(gain_type_) +
          " unknown");
    }
  }
}

void SagecalPredict::finish() { getNextStep()->finish(); }

void SagecalPredict::show(std::ostream& os) const {
#ifdef HAVE_LIBDIRAC
  os << "SagecalPredict " << name_ << '\n';
#endif
#ifdef HAVE_LIBDIRAC_CUDA
  os << "SagecalPredict (GPU) " << name_ << '\n';
#endif
#if defined(HAVE_LIBDIRAC) || defined(HAVE_LIBDIRAC_CUDA)
  os << "  sourcedb:                " << source_db_name_ << '\n';
  os << "   number of directions:      " << patch_list_.size() << '\n';
  os << "   number of components:   " << source_list_.size() << '\n';
  os << "   correct freq smearing:  " << std::boolalpha << true << '\n';
  os << "  apply beam:              ";
  switch (beam_mode_) {
    case DOBEAM_NONE:
      os << std::boolalpha << false << '\n';
      break;
    case DOBEAM_FULL:
      os << "mode: full beam" << '\n';
      break;
    case DOBEAM_FULL_WB:
      os << "mode: full beam" << '\n';
      os << "  use channelfreq: true" << '\n';
      break;
    case DOBEAM_ARRAY:
      os << "mode: array beam" << '\n';
      break;
    case DOBEAM_ARRAY_WB:
      os << "mode: array beam" << '\n';
      os << "  use channelfreq: true" << '\n';
      break;
    case DOBEAM_ELEMENT:
      os << "mode: element beam" << '\n';
      break;
    case DOBEAM_ELEMENT_WB:
      os << "mode: element beam" << '\n';
      os << "  use channelfreq: true" << '\n';
      break;
  }
  os << "  operation:   ";
  switch (operation_) {
    case Operation::kReplace:
      os << "replace";
      break;
    case Operation::kAdd:
      os << "add";
      break;
    case Operation::kSubtract:
      os << "subtract";
      break;
  }
  os << '\n';
#endif /* HAVE_LIBDIRAC || HAVE_LIBDIRAC_CUDA */
  if (parm_on_disk_) {
    os << "H5 name " << h5_name_ << "\n";
    os << " SolSet " << h5_parm_.GetSolSetName() << "\n";
    os << " SolTab " << soltab_name_ << "\n";
    os << "Time interval " << time_interval_ << "\n";
    os << "Timeslots per parmupdate " << timeslots_per_parmupdate_ << "\n";
  }
}

void SagecalPredict::showTimings(std::ostream& os,
                                 [[maybe_unused]] double duration) const {
  os << "  ";
  base::FlagCounter::showPerc1(os, timer_.getElapsed(), duration);
  os << " SagecalPredict " << name_ << '\n';
}

base::Direction SagecalPredict::GetFirstDirection() const {
  return patch_list_.front()->direction();
}

#if defined(HAVE_LIBDIRAC) || defined(HAVE_LIBDIRAC_CUDA)
void SagecalPredict::loadData(std::unique_ptr<dp3::base::DPBuffer>& buffer) {
  const size_t nBl = info().nbaselines();
  const size_t nCh = info().nchan();
  [[maybe_unused]] const size_t nCr = info().ncorr();

  assert(iodata_.n_baselines >= nBl);
  assert(iodata_.n_stations == getInfo().nantenna());
  assert(iodata_.n_channels == nCh);
  assert(4 == nCr);

  const DPBuffer::UvwType& uvw = buffer->GetUvw();
  size_t row0 = 0;
  // load data, skipping autocorrelations
  for (size_t bl = 0; bl < nBl; bl++) {
    uint ant1 = info().getAnt1()[bl];
    uint ant2 = info().getAnt2()[bl];
    if (ant1 != ant2) {
      // shape nBl x 3
      iodata_.u_[row0] = uvw(bl, 0);
      iodata_.v_[row0] = uvw(bl, 1);
      iodata_.w_[row0] = uvw(bl, 2);

      iodata_.baseline_arr_[row0].sta1 = ant1;
      iodata_.baseline_arr_[row0].sta2 = ant2;
      iodata_.baseline_arr_[row0].flag = 0;

#pragma GCC ivdep
      for (uint ch = 0; ch < nCh; ch++) {
        const std::complex<float>* data_pointer = &buffer->GetData()(bl, ch, 0);
        iodata_.data_[iodata_.n_baselines * 8 * ch + row0 * 8 + 0] =
            data_pointer[0].real();
        iodata_.data_[iodata_.n_baselines * 8 * ch + row0 * 8 + 1] =
            data_pointer[0].imag();
        iodata_.data_[iodata_.n_baselines * 8 * ch + row0 * 8 + 2] =
            data_pointer[1].real();
        iodata_.data_[iodata_.n_baselines * 8 * ch + row0 * 8 + 3] =
            data_pointer[1].imag();
        iodata_.data_[iodata_.n_baselines * 8 * ch + row0 * 8 + 4] =
            data_pointer[2].real();
        iodata_.data_[iodata_.n_baselines * 8 * ch + row0 * 8 + 5] =
            data_pointer[2].imag();
        iodata_.data_[iodata_.n_baselines * 8 * ch + row0 * 8 + 6] =
            data_pointer[3].real();
        iodata_.data_[iodata_.n_baselines * 8 * ch + row0 * 8 + 7] =
            data_pointer[3].imag();
      }
      row0++;
    }
  }

  /* rescale u,v,w by 1/c */
  const double inv_c = 1.0 / casacore::C::c;
  my_dscal(iodata_.n_baselines, inv_c, iodata_.u_.data());
  my_dscal(iodata_.n_baselines, inv_c, iodata_.v_.data());
  my_dscal(iodata_.n_baselines, inv_c, iodata_.w_.data());
}

void SagecalPredict::writeData(std::unique_ptr<DPBuffer>& buffer) {
  const size_t nBl = info().nbaselines();
  const size_t nCh = info().nchan();

  assert(iodata_.n_baselines >= nBl);
  assert(iodata_.n_channels == nCh);

  DPBuffer::DataType& data = buffer->GetData();
  size_t row0 = 0;
  // load data, skipping autocorrelations
  for (size_t bl = 0; bl < nBl; bl++) {
    uint ant1 = info().getAnt1()[bl];
    uint ant2 = info().getAnt2()[bl];
    if (ant1 != ant2) {
#pragma GCC ivdep
      for (uint ch = 0; ch < nCh; ch++) {
        data(bl, ch, 0) = std::complex<float>(
            iodata_.data_[iodata_.n_baselines * 8 * ch + row0 * 8 + 0],
            iodata_.data_[iodata_.n_baselines * 8 * ch + row0 * 8 + 1]);
        data(bl, ch, 1) = std::complex<float>(
            iodata_.data_[iodata_.n_baselines * 8 * ch + row0 * 8 + 2],
            iodata_.data_[iodata_.n_baselines * 8 * ch + row0 * 8 + 3]);
        data(bl, ch, 2) = std::complex<float>(
            iodata_.data_[iodata_.n_baselines * 8 * ch + row0 * 8 + 4],
            iodata_.data_[iodata_.n_baselines * 8 * ch + row0 * 8 + 5]);
        data(bl, ch, 3) = std::complex<float>(
            iodata_.data_[iodata_.n_baselines * 8 * ch + row0 * 8 + 6],
            iodata_.data_[iodata_.n_baselines * 8 * ch + row0 * 8 + 7]);
      }
      row0++;
    }
  }
}

void SagecalPredict::readAuxData(const DPInfo& _info) {
  beam_data_->update_metadata(_info, iodata_.f0, iodata_.n_channels,
                              iodata_.freqs_, beam_mode_);
  runtime_beam_data_.p_ra0_ = beam_data_->p_ra0_;
  runtime_beam_data_.p_dec0_ = beam_data_->p_dec0_;
  runtime_beam_data_.b_ra0_ = beam_data_->b_ra0_;
  runtime_beam_data_.b_dec0_ = beam_data_->b_dec0_;
  const size_t tile_size = 1;
  try {
    runtime_beam_data_.time_utc_.resize(tile_size);
  } catch (const std::bad_alloc& e) {
    throw std::runtime_error("Allocating memory failure");
  }
}
#endif /* HAVE_LIBDIRAC || HAVE_LIBDIRAC_CUDA */

}  // namespace steps
}  // namespace dp3
