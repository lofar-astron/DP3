// MSWriter.h: DPPP step writing to an MS
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step writing to an MS
/// @author Ger van Diepen

#ifndef DPPP_MSWRITER_H
#define DPPP_MSWRITER_H

#include "Step.h"

#include "../base/StManParsetKeys.h"

#include <casacore/tables/Tables/ColumnDesc.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ArrayColumn.h>

namespace casacore {
class Table;
}

namespace dp3 {
namespace common {
class ParameterSet;
}

namespace steps {
class InputStep;

/// @brief DPPP step writing to an MS

/// This class is a Step creating a new MeasurementSet and writing
/// all data in it.
/// Most meta information (subtables and meta columns in main table) is
/// copied from the input MeasurementSet given by the MSReader object.
/// <br>
/// In principle the new MS uses the same storage managers as used in the
/// input MS, but in case of an MS stored with LofarStMan it will use the
/// optimal storage managers (ISM for slowly varying meta data, TSM for
/// bulk data, SSM for others).
///
/// The SPECTRAL_WINDOW table will be changed to reflect the channels
/// being used or averaged.
/// The OBSERVATION table will be updated for the correct start and end time.
/// The HISTORY table gets an entry containing the parset values and the
/// DPPP version.

class MSWriter : public Step {
 public:
  explicit MSWriter(InputStep& reader, const std::string& out_name,
                    const common::ParameterSet&, const std::string& prefix);

  ~MSWriter() override;

  /// Process the next data chunk.
  /// It returns false when at the end.
  bool process(const base::DPBuffer&) override;

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Show the step parameters.
  void show(std::ostream&) const override;

  /// Add some data to the MeasurementSet written/updated.
  /// Calls addToMS from the previous step, with the current output msname.
  void addToMS(const string&) override;

  /// Update the general info.
  void updateInfo(const base::DPInfo&) override;

  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override;

  /// Write the parset info into the HISTORY table of the MS.
  static void WriteHistory(casacore::Table& ms,
                           const common::ParameterSet& parset);

  static void UpdateBeam(const std::string& out_name,
                         const std::string& out_col_name,
                         const base::DPInfo& info);

 private:
  /// Create an array column description and add to table with given
  /// stoage manager (if given).
  void MakeArrayColumn(casacore::ColumnDesc desc,
                       const casacore::IPosition& shape,
                       casacore::DataManager* dm, casacore::Table& table,
                       bool make_direct_column = false);

  /// Create the MS by cloning all subtables from the input MS.
  /// All output columns in the main table are using normal storage managers.
  /// The SPECTRAL_WINDOW table is adapted as needed.
  void CreateMs(const std::string& out_name, const base::DPInfo& info,
                unsigned int tile_size, unsigned int tile_n_chan);

  /// Update the SPECTRAL_WINDOW table for averaged channels.
  void UpdateSpw(const string& out_name, const base::DPInfo& info);

  /// Update the OBSERVATION table with the correct start and end time.
  void UpdateObs(const string& out_name);

  /// Update the FIELD table with the new phase center.
  void UpdatePhaseCentre(const string& out_name, const base::DPInfo& info);

  /// Update @ref buffer_ with the provided @a buffer.
  void UpdateInternalBuffer(const base::DPBuffer& buffer);

  /// Process the data in @ref buffer.
  ///
  /// This function does not access @ref buffer_.
  void ProcessBuffer(base::DPBuffer& buffer);

  /// Write the data, flags, etc.
  ///
  /// This function does not access @ref buffer_.
  void WriteData(casacore::Table& out, const base::DPBuffer& buf);

  /// Write the full resolution flags (flags before any averaging).
  ///
  /// This function does not access @ref buffer_.
  void WriteFullResFlags(casacore::Table& out, const base::DPBuffer& buf);

  /// Write all meta data columns for a time slot (ANTENNA1, etc.)
  ///
  /// This function does not access @ref buffer_.
  void WriteMeta(casacore::Table& out, const base::DPBuffer& buf);

  /// Copy meta data columns for a time slot (ANTENNA1, etc.)
  /// It also copies all time info if possible.
  void CopyMeta(const casacore::Table& in, casacore::Table& out,
                bool copy_time_info);

  /// Copy the contents of a scalar column.
  template <typename T>
  void FillSca(const T& value, casacore::Table& out,
               const casacore::String& column_name) {
    casacore::ScalarColumn<T> out_col(out, column_name);
    out_col.fillColumn(value);
  }

  /// Copy the contents of an array column.
  template <typename T>
  void FillArr(const casacore::Array<T>& value, casacore::Table& out,
               const casacore::String& column_name) {
    casacore::ArrayColumn<T> out_col(out, column_name);
    out_col.fillColumn(value);
  }

  /// Copy the contents of a scalar column.
  template <typename T>
  void CopySca(const casacore::Table& in, casacore::Table& out,
               const casacore::String& column_name) {
    casacore::ROScalarColumn<T> in_col(in, column_name);
    casacore::ScalarColumn<T> out_col(out, column_name);
    out_col.putColumn(in_col.getColumn());
  }

  /// Copy the contents of an array column.
  template <typename T>
  void CopyArr(const casacore::Table& in, casacore::Table& out,
               const casacore::String& column_name) {
    casacore::ROArrayColumn<T> in_col(in, column_name);
    casacore::ArrayColumn<T> out_col(out, column_name);
    out_col.putColumn(in_col.getColumn());
  }

  InputStep& reader_;
  string name_;
  string out_name_;
  base::DPBuffer buffer_;
  casacore::Table ms_;
  common::ParameterSet parset_;  ///< parset for writing history
  casacore::String data_col_name_;
  casacore::String flag_col_name_;
  casacore::String weight_col_name_;
  double interval_;
  bool overwrite_;  ///< Overwrite an existing output MS?
  bool copy_corr_data_;
  bool copy_model_data_;
  bool write_full_res_flags_;
  unsigned int tile_size_;
  unsigned int tile_n_chan_;
  unsigned int nr_corr_;
  unsigned int nr_chan_;
  unsigned int nr_bl_;
  unsigned int nr_times_;
  unsigned int n_chan_avg_;      ///< nr of channels in input averaged to 1
  unsigned int n_time_avg_;      ///< nr of times in input averaged to 1
  unsigned int nr_times_flush_;  ///< flush every N time slots (0=no flush)
  unsigned int nr_done_;         ///< nr of time slots written
  std::string vds_dir_;          ///< directory where to put VDS file
  std::string cluster_desc_;     ///< name of clusterdesc file
  base::StManParsetKeys st_man_keys_;

  /// The total time spent in the writer.
  common::NSTimer timer_;

  /// The time spent updating the buffer.
  ///
  /// This update process will spend time reading from the input.
  common::NSTimer update_buffer_timer_;

  /// The time spent writing data to the output MS.
  common::NSTimer writer_timer_;
};

}  // namespace steps
}  // namespace dp3

#endif
