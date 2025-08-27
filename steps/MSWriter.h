// MSWriter.h: DP3 step writing to an MS
// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_MSWRITER_H_
#define DP3_STEPS_MSWRITER_H_

#include "OutputStep.h"

#include <memory>
#include <thread>

#include <aocommon/lane.h>

#include <casacore/tables/Tables/ColumnDesc.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ArrayColumn.h>

#include "../base/StManParsetKeys.h"
#include "../common/Timer.h"

namespace casacore {
class Table;
}

namespace dp3 {
namespace common {
class ParameterSet;
}

namespace steps {
class InputStep;

std::unique_ptr<casacore::DataManager> MakeStMan(
    const std::string& type_name, const std::string& instance_name,
    const casacore::Record& record = casacore::Record());

/// @brief DP3 step writing to an MS

/// This class is a Step creating a new MeasurementSet and writing
/// all data in it.
/// Most meta information (subtables and meta columns in main table) is
/// copied from the input MeasurementSet given by the MsReader object.
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

class MSWriter : public OutputStep {
 public:
  explicit MSWriter(const std::string& out_name, const common::ParameterSet&,
                    const std::string& prefix);

  ~MSWriter() override;

  common::Fields getRequiredFields() const override {
    return kDataField | kFlagsField | kWeightsField | kUvwField;
  }

  /// Process the next data chunk.
  /// It returns false when at the end.
  bool process(std::unique_ptr<base::DPBuffer> buffer) override;

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Show the step parameters.
  void show(std::ostream&) const override;

  /// Update the general info.
  void updateInfo(const base::DPInfo&) override;

  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override;

  /// Write the parset info into the HISTORY table of the MS.
  static void WriteHistory(casacore::Table& ms,
                           const common::ParameterSet& parset);

  static void UpdateBeam(casacore::Table& main_table,
                         const std::string& out_col_name,
                         const base::DPInfo& info);

  /// Update the FIELD table with the new phase center.
  static void UpdatePhaseCentre(const std::string& out_name,
                                const casacore::MDirection& new_phase_dir);

  /// Update the SPECTRAL_WINDOW table for averaged channels.
  static void UpdateSpw(const std::string& out_name, const base::DPInfo& info);

  /// Update the OBSERVATION table with the correct start and end time.
  static void UpdateObs(const std::string& out_name, const base::DPInfo& info);

  static std::string InsertNumberInFilename(const std::string& name,
                                            size_t number);

 private:
  void StartNewMs();
  void FinishMs();

  /// Create the MS by cloning all subtables from the input MS.
  /// All output columns in the main table are using normal storage managers.
  /// The SPECTRAL_WINDOW table is adapted as needed.
  void CreateMs(const std::string& out_name, unsigned int tile_size,
                unsigned int tile_n_chan);

  /// Copy subtables from 'original_table` to the output MS.
  void CopySubTables(casacore::Table& original_table);

  /// Process the data in @ref buffer.
  ///
  /// This function does not access @ref internal_buffer_.
  void ProcessBuffer(base::DPBuffer& buffer);

  /// Write the data, flags, etc.
  ///
  /// This function does not access @ref internal_buffer_.
  void WriteData(casacore::Table& out, base::DPBuffer& buf);

  /// Write all meta data columns for a time slot (ANTENNA1, etc.)
  ///
  /// This function does not access @ref internal_buffer_.
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

  /// Name of the step
  std::string name_;
  /// Requested filename of the output ms (before adding chunk number)
  std::string out_name_;
  /// The filename of the currently being written chunk
  std::string chunk_name_;
  casacore::Table ms_;
  common::ParameterSet parset_;  ///< parset for writing history
  casacore::String data_col_name_;
  casacore::String flag_col_name_;
  casacore::String weight_col_name_;
  bool overwrite_;  ///< Overwrite an existing output MS?
  bool copy_corr_data_;
  bool copy_model_data_;
  unsigned int tile_size_;
  unsigned int tile_n_chan_;
  unsigned int nr_times_flush_;  ///< flush every N time slots (0=no flush)
  unsigned int nr_done_;         ///< nr of time slots written
  /// In seconds, split measurement set in chunks with a given maximum size.
  /// Zero means no chunking.
  double chunk_duration_ = 0.0;
  /// Start time of last started chunk, to determine length
  double chunk_start_time_ = 0.0;
  size_t current_chunk_index_ = 0;

  std::string vds_dir_;       ///< directory where to put VDS file
  std::string cluster_desc_;  ///< name of clusterdesc file
  base::StManParsetKeys st_man_keys_;
  // For now, metadata compression is turned on only when explicitly turned on
  // in compilation. Once we're confident enough about this feature, we can
  // enable it by default.
  bool scalar_flags_ = METADATA_COMPRESSION_DEFAULT;
  bool uvw_compression_ = METADATA_COMPRESSION_DEFAULT;
  bool antenna_compression_ = METADATA_COMPRESSION_DEFAULT;

  /// The total time spent in the writer.
  common::NSTimer timer_;

  /// The time spent writing data to the output MS.
  ///
  /// This timer can either run in the main thread or the dedicated write
  /// thread. When it is running in the write thread the time may exceed the
  /// time of \ref timer_.
  common::NSTimer writer_timer_;

  /// The time spent creating a task for the write queue.
  ///
  /// Creating a task requires an expensive deep copy. When this measures the
  /// time the main thread spends in that operation.
  common::NSTimer create_task_timer_;

  /// The size of the write buffer.
  ///
  /// On machines with "fast" I/O the writing is usually done when the next
  /// write task arrives. On machines with "slow" I/O the writing is never
  /// finished when the next task arrives. The buffer is optimized for machines
  /// with faster I/O. When an unrelated process causes temporary I/O the buffer
  /// can keep the processing to stall on a full buffer.
  ///
  /// On the other hand the buffer requires additional memory, so making the
  /// buffer large takes more memory while the "slow" I/O will just delay the
  /// final part of the processing. Based on experiments locally and on DAS6 the
  /// selected size seems a nice trade-off.
  aocommon::Lane<std::unique_ptr<base::DPBuffer>> write_queue_{3};

  /// Creates task for the \ref write_queue_.
  void CreateTask(std::unique_ptr<base::DPBuffer> buffer);

  /// The thread used to process \ref write_queue_.
  ///
  /// Note since writing with multiple threads to the same MeasurementSet will
  /// cause issues in Casacore only one thread is used.
  std::thread write_queue_thread_;

  /// The task running in \ref write_queue_thread_.
  void WriteQueueProcess();

  /// Does the writer use a separate thread?
  ///
  /// A seperate thread can only be used when the write step in the chain.
  bool use_write_thread_{false};

  /// Is the write thread running?
  ///
  /// There are two places where the write thread can be stopped:
  /// - finish, there it needs to stop so the finish process can write to the
  ///   MS.
  /// - destructor, destroying a non-detached running thread is a contract
  ///   violation.
  bool is_write_thread_active_{false};

  /// Stops the write thread when it is active.
  void StopWriteThread();
};

}  // namespace steps
}  // namespace dp3

#endif
