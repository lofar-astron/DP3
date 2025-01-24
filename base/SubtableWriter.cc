#include "SubtableWriter.h"

#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Containers/Record.h>
#include <casacore/measures/Measures/MFrequency.h>
#include <casacore/measures/TableMeasures/TableMeasDesc.h>
#include <casacore/tables/Tables/ArrColDesc.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScaColDesc.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/SetupNewTab.h>
#include <casacore/measures/Measures/Stokes.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <cassert>

namespace dp3::base {
using namespace casacore;

SubtableWriter::SubtableWriter(std::string path, const int nr_channels)
    : path_{path} {
  TableDesc desc = MS::requiredTableDesc();
  SetupNewTable main_table(path_, desc, Table::New);

  const int kNCorrelations = 4;
  casacore::IPosition dataShape(2, kNCorrelations, nr_channels);
  desc.rwColumnDesc("SIGMA").setShape(IPosition(1, kNCorrelations));
  desc.rwColumnDesc("SIGMA").setOptions(ColumnDesc::Option::FixedShape |
                                        ColumnDesc::Option::Direct);
  desc.rwColumnDesc("WEIGHT").setShape(IPosition(1, kNCorrelations));
  desc.rwColumnDesc("WEIGHT").setOptions(ColumnDesc::Option::FixedShape |
                                         ColumnDesc::Option::Direct);
  desc.rwColumnDesc("FLAG").setShape(dataShape);
  desc.rwColumnDesc("FLAG").setOptions(ColumnDesc::Option::FixedShape |
                                       ColumnDesc::Option::Direct);
  desc.rwColumnDesc("UVW").setOptions(ColumnDesc::Option::FixedShape |
                                      ColumnDesc::Option::Direct);
  desc.rwColumnDesc("FLAG_CATEGORY").setShape(IPosition(3, 0, 0, 0));

  IncrementalStMan storage_manager;
  main_table.bindColumn("TIME", storage_manager);
  main_table.bindColumn("TIME_CENTROID", storage_manager);
  main_table.bindColumn("ANTENNA1", storage_manager);
  main_table.bindColumn("DATA_DESC_ID", storage_manager);
  main_table.bindColumn("INTERVAL", storage_manager);
  main_table.bindColumn("EXPOSURE", storage_manager);
  main_table.bindColumn("PROCESSOR_ID", storage_manager);
  main_table.bindColumn("SCAN_NUMBER", storage_manager);
  main_table.bindColumn("STATE_ID", storage_manager);
  main_table.bindColumn("SIGMA", storage_manager);

  main_table.bindColumn("FLAG_CATEGORY", storage_manager);
  main_table.bindColumn("ARRAY_ID", storage_manager);
  main_table.bindColumn("FEED1", storage_manager);
  main_table.bindColumn("FEED2", storage_manager);
  main_table.bindColumn("FIELD_ID", storage_manager);
  main_table.bindColumn("FLAG_ROW", storage_manager);
  main_table.bindColumn("OBSERVATION_ID", storage_manager);

  ms_ = MeasurementSet(main_table);
  ms_.createDefaultSubtables(Table::New);

  ArrayColumnDesc<std::complex<float>> data_column_desc =
      ArrayColumnDesc<std::complex<float>>(
          MS::columnName(casacore::MSMainEnums::DATA));
  data_column_desc.setShape(dataShape);
  data_column_desc.setOptions(ColumnDesc::Direct | ColumnDesc::FixedShape);
  ms_.addColumn(data_column_desc);

  ArrayColumnDesc<float> weight_spectrum_column_desc = ArrayColumnDesc<float>(
      MS::columnName(casacore::MSMainEnums::WEIGHT_SPECTRUM));
  weight_spectrum_column_desc.setShape(dataShape);
  weight_spectrum_column_desc.setOptions(ColumnDesc::Direct |
                                         ColumnDesc::FixedShape);
  ms_.addColumn(weight_spectrum_column_desc);
};

void SubtableWriter::WriteDataDescEntry(size_t spectralWindowId,
                                        size_t polarizationId, bool flagRow) {
  MSDataDescription data_description_table = ms_.dataDescription();
  ScalarColumn<int> spectral_window_id_column(
      data_description_table, data_description_table.columnName(
                                  MSDataDescriptionEnums::SPECTRAL_WINDOW_ID));
  ScalarColumn<int> polarization_id_column(
      data_description_table, data_description_table.columnName(
                                  MSDataDescriptionEnums::POLARIZATION_ID));
  ScalarColumn<bool> flag_row_column(
      data_description_table,
      data_description_table.columnName(MSDataDescriptionEnums::FLAG_ROW));

  const size_t index = data_description_table.nrow();
  data_description_table.addRow();

  spectral_window_id_column.put(index, spectralWindowId);
  polarization_id_column.put(index, polarizationId);
  flag_row_column.put(index, flagRow);
}

void SubtableWriter::WriteFeedEntries(const std::vector<AntennaInfo> &antennas,
                                      double time) {
  // Open feed table
  MSFeed feed_table = ms_.feed();
  // Define feed table columns
  ScalarColumn<int> antenna_id_col(
      feed_table, feed_table.columnName(MSFeedEnums::ANTENNA_ID));
  ScalarColumn<int> feed_id_col(feed_table,
                                feed_table.columnName(MSFeedEnums::FEED_ID));
  ScalarColumn<int> spectral_window_id_col(
      feed_table, feed_table.columnName(MSFeedEnums::SPECTRAL_WINDOW_ID));
  ScalarColumn<double> time_col(feed_table,
                                feed_table.columnName(MSFeedEnums::TIME));
  ScalarColumn<int> num_receptors_col(
      feed_table, feed_table.columnName(MSFeedEnums::NUM_RECEPTORS));
  ScalarColumn<int> beam_id_col(feed_table,
                                feed_table.columnName(MSFeedEnums::BEAM_ID));
  ArrayColumn<double> beam_offset_col(
      feed_table, feed_table.columnName(MSFeedEnums::BEAM_OFFSET));
  ArrayColumn<casacore::String> polarization_type_col(
      feed_table, feed_table.columnName(MSFeedEnums::POLARIZATION_TYPE));
  ArrayColumn<std::complex<float>> pol_response_col(
      feed_table, feed_table.columnName(MSFeedEnums::POL_RESPONSE));
  ArrayColumn<double> position_col(
      feed_table, feed_table.columnName(MSFeedEnums::POSITION));
  ArrayColumn<double> receptor_angle_col(
      feed_table, feed_table.columnName(MSFeedEnums::RECEPTOR_ANGLE));

  // Get last index in feed table
  // this is because it is appending in the feed the antennas that have been
  // specified as argument
  size_t row_index = feed_table.nrow();
  feed_table.addRow(antennas.size());

  for (size_t antIndex = 0; antIndex != antennas.size(); ++antIndex) {
    antenna_id_col.put(row_index, antIndex);
    feed_id_col.put(row_index, 0);
    spectral_window_id_col.put(row_index, -1);
    time_col.put(row_index, time);
    num_receptors_col.put(row_index, 2);
    beam_id_col.put(row_index, -1);

    // Define the been offset
    casacore::Array<double> beam_offset(IPosition(2, 2, 2), 0.0);
    beam_offset_col.put(row_index, beam_offset);

    // Define the polarization names
    casacore::Vector<casacore::String> polType(2);
    polType[0] = 'X';
    polType[1] = 'Y';
    polarization_type_col.put(row_index, polType);

    // Define polarization response
    casacore::Array<std::complex<float>> pol_response(IPosition(2, 2, 2));
    pol_response(IPosition(2, 0, 0)) = 1.;
    pol_response(IPosition(2, 1, 0)) = 0.;
    pol_response(IPosition(2, 0, 1)) = 0.;
    pol_response(IPosition(2, 1, 1)) = 1.;
    pol_response_col.put(row_index, pol_response);

    // Define feed position
    casacore::Vector<double> position(3);
    position[0] = 0.0;
    position[1] = 0.0;
    position[2] = 0.0;
    position_col.put(row_index, position);

    // Define feed receptor angle
    casacore::Vector<double> receptorAngle(2);
    receptorAngle[0] = 0.0;
    receptorAngle[1] = M_PI * 0.5;
    receptor_angle_col.put(row_index, receptorAngle);

    ++row_index;
  }
}

void SubtableWriter::WriteBandInfo(const std::string &name,
                                   const std::vector<ChannelInfo> &channels,
                                   double reference_frequency,
                                   double total_bandwidth, bool flag_row) {
  MSSpectralWindow spw_table = ms_.spectralWindow();

  ScalarColumn<int> num_chan_col = ScalarColumn<int>(
      spw_table, spw_table.columnName(MSSpectralWindowEnums::NUM_CHAN));
  ScalarColumn<casacore::String> name_col = ScalarColumn<casacore::String>(
      spw_table, spw_table.columnName(MSSpectralWindowEnums::NAME));
  ScalarColumn<double> ref_freq_col = ScalarColumn<double>(
      spw_table, spw_table.columnName(MSSpectralWindowEnums::REF_FREQUENCY));

  ArrayColumn<double> chan_freq_col = ArrayColumn<double>(
      spw_table, spw_table.columnName(MSSpectralWindowEnums::CHAN_FREQ));
  ArrayColumn<double> chan_width_col = ArrayColumn<double>(
      spw_table, spw_table.columnName(MSSpectralWindowEnums::CHAN_WIDTH));
  ScalarColumn<int> meas_freq_ref_col = ScalarColumn<int>(
      spw_table, spw_table.columnName(MSSpectralWindowEnums::MEAS_FREQ_REF));
  ArrayColumn<double> effective_bw_col = ArrayColumn<double>(
      spw_table, spw_table.columnName(MSSpectralWindowEnums::EFFECTIVE_BW));
  ArrayColumn<double> resolution_col = ArrayColumn<double>(
      spw_table, spw_table.columnName(MSSpectralWindowEnums::RESOLUTION));
  ScalarColumn<double> total_bw_col = ScalarColumn<double>(
      spw_table, spw_table.columnName(MSSpectralWindowEnums::TOTAL_BANDWIDTH));
  ScalarColumn<bool> flag_row_col = ScalarColumn<bool>(
      spw_table, spw_table.columnName(MSSpectralWindowEnums::FLAG_ROW));

  const size_t n_channels = channels.size();
  size_t row_index = spw_table.nrow();
  spw_table.addRow();
  num_chan_col.put(row_index, n_channels);
  name_col.put(row_index, name);
  ref_freq_col.put(row_index, reference_frequency);

  casacore::Vector<double> chan_freq_vec(n_channels);
  casacore::Vector<double> chan_width_vec(n_channels);
  casacore::Vector<double> effective_bw_vec(n_channels);
  casacore::Vector<double> resolution_vec(n_channels);

  for (size_t ch = 0; ch != channels.size(); ++ch) {
    chan_freq_vec[ch] = channels[ch].channel_frequency;
    chan_width_vec[ch] = channels[ch].channel_width;
    effective_bw_vec[ch] = channels[ch].effective_bandwidth;
    resolution_vec[ch] = channels[ch].resolution;
  }
  chan_freq_col.put(row_index, chan_freq_vec);
  chan_width_col.put(row_index, chan_width_vec);
  meas_freq_ref_col.put(row_index, casacore::MFrequency::Types::TOPO);
  effective_bw_col.put(row_index, effective_bw_vec);
  resolution_col.put(row_index, resolution_vec);

  total_bw_col.put(row_index, total_bandwidth);
  flag_row_col.put(row_index, flag_row);

  WriteDataDescEntry(row_index, 0, false);
}

void SubtableWriter::WriteAntennas(const std::vector<AntennaInfo> &antennas,
                                   const std::array<double, 9> &coordinate_axes,
                                   double time) {
  MSAntenna antenna_table = ms_.antenna();
  ScalarColumn<casacore::String> name_col = ScalarColumn<casacore::String>(
      antenna_table, antenna_table.columnName(MSAntennaEnums::NAME));
  ScalarColumn<casacore::String> station_col = ScalarColumn<casacore::String>(
      antenna_table, antenna_table.columnName(MSAntennaEnums::STATION));
  ScalarColumn<casacore::String> type_col = ScalarColumn<casacore::String>(
      antenna_table, antenna_table.columnName(MSAntennaEnums::TYPE));
  ScalarColumn<casacore::String> mount_col = ScalarColumn<casacore::String>(
      antenna_table, antenna_table.columnName(MSAntennaEnums::MOUNT));
  ArrayColumn<double> position_col = ArrayColumn<double>(
      antenna_table, antenna_table.columnName(MSAntennaEnums::POSITION));
  ScalarColumn<double> dish_diameter_col = ScalarColumn<double>(
      antenna_table, antenna_table.columnName(MSAntennaEnums::DISH_DIAMETER));
  ScalarColumn<casacore::Bool> flag_row = ScalarColumn<casacore::Bool>(
      antenna_table, antenna_table.columnName(MSAntennaEnums::FLAG_ROW));

  size_t row_index = antenna_table.nrow();
  antenna_table.addRow(antennas.size());

  for (std::vector<AntennaInfo>::const_iterator antenna_iterator =
           antennas.begin();
       antenna_iterator != antennas.end(); ++antenna_iterator) {
    name_col.put(row_index, antenna_iterator->name);
    station_col.put(row_index, antenna_iterator->station);
    type_col.put(row_index, antenna_iterator->type);
    mount_col.put(row_index, antenna_iterator->mount);
    dish_diameter_col.put(row_index, antenna_iterator->diameter);
    flag_row.put(row_index, antenna_iterator->flag);

    Vector<double> position_array(3);
    position_array[0] = antenna_iterator->x;
    position_array[1] = antenna_iterator->y;
    position_array[2] = antenna_iterator->z;
    position_col.put(row_index, position_array);

    ++row_index;
  }

  casacore::Matrix<double> coordinateAxesVec(3, 3);
  std::copy(coordinate_axes.begin(), coordinate_axes.end(),
            coordinateAxesVec.begin());

  MSObservation obs_table = ms_.observation();
  auto telescope_name =
      obs_table.columnName(MSObservationEnums::TELESCOPE_NAME);

  antenna_table.rwKeywordSet().define(telescope_name + "_COORDINATE_AXES",
                                      coordinateAxesVec);

  WriteFeedEntries(antennas, time);
}

void SubtableWriter::WriteLinearPolarizations(bool flag_row, const int n_pol) {
  assert(n_pol == 2 or n_pol == 4);
  MSPolarization pol_table = ms_.polarization();
  ScalarColumn<int> num_corr_col = ScalarColumn<int>(
      pol_table, pol_table.columnName(MSPolarizationEnums::NUM_CORR));
  ArrayColumn<int> corr_type_col = ArrayColumn<int>(
      pol_table, pol_table.columnName(MSPolarizationEnums::CORR_TYPE));
  ArrayColumn<int> corr_product_col = ArrayColumn<int>(
      pol_table, pol_table.columnName(MSPolarizationEnums::CORR_PRODUCT));
  ScalarColumn<bool> flag_row_col = ScalarColumn<bool>(
      pol_table, pol_table.columnName(MSPolarizationEnums::FLAG_ROW));

  size_t row_index = pol_table.nrow();
  pol_table.addRow(1);
  num_corr_col.put(row_index, n_pol);

  casacore::Vector<int> c_type_vec(n_pol);
  if (n_pol == 2) {
    c_type_vec[0] = Stokes::XX;
    c_type_vec[1] = Stokes::YY;
  } else if (n_pol == 4) {
    c_type_vec[0] = Stokes::XX;
    c_type_vec[1] = Stokes::XY;
    c_type_vec[2] = Stokes::YX;
    c_type_vec[3] = Stokes::YY;
  }
  corr_type_col.put(row_index, c_type_vec);

  casacore::Array<int> c_prod_arr(IPosition(2, 2, n_pol));
  if (n_pol == 2) {
    c_prod_arr(IPosition(2, 0, 0)) = 0;
    c_prod_arr(IPosition(2, 1, 0)) = 0;
    c_prod_arr(IPosition(2, 0, 1)) = 1;
    c_prod_arr(IPosition(2, 1, 1)) = 1;
  } else if (n_pol == 4) {
    c_prod_arr(IPosition(2, 0, 0)) = 0;
    c_prod_arr(IPosition(2, 1, 0)) = 0;
    c_prod_arr(IPosition(2, 0, 1)) = 0;
    c_prod_arr(IPosition(2, 1, 1)) = 1;
    c_prod_arr(IPosition(2, 0, 2)) = 1;
    c_prod_arr(IPosition(2, 1, 2)) = 0;
    c_prod_arr(IPosition(2, 0, 3)) = 1;
    c_prod_arr(IPosition(2, 1, 3)) = 1;
  }

  corr_product_col.put(row_index, c_prod_arr);

  flag_row_col.put(row_index, flag_row);
}

void SubtableWriter::WriteSource(const SourceInfo &source) {
  // The source table is not a required table so it has to be created
  // from scratch

  // Get default table description for source
  TableDesc source_table_desc = MSSource::requiredTableDesc();
  // Add rest frequency column to the default
  casacore::ArrayColumnDesc<double> restFrequencyColumnDesc =
      ArrayColumnDesc<double>(
          MSSource::columnName(MSSourceEnums::REST_FREQUENCY));
  source_table_desc.addColumn(restFrequencyColumnDesc);

  // Add frequency information to source table
  TableMeasRefDesc meas_ref(MFrequency::DEFAULT);
  TableMeasValueDesc meas_val(
      source_table_desc, MSSource::columnName(MSSourceEnums::REST_FREQUENCY));
  TableMeasDesc<MFrequency> rest_freq_col_meas(meas_val, meas_ref);
  // write makes the Measure column persistent.
  rest_freq_col_meas.write(source_table_desc);

  // Creating the table using the helper class
  SetupNewTable source_table_setup(ms_.sourceTableName(), source_table_desc,
                                   Table::New);
  MSSource source_table(source_table_setup);

  // Link table to the main measurement set
  ms_.rwKeywordSet().defineTable(MS::keywordName(casacore::MSMainEnums::SOURCE),
                                 source_table);
  // Open columns to write data into
  ScalarColumn<int> source_id_col = ScalarColumn<int>(
      source_table, MSSource::columnName(MSSourceEnums::SOURCE_ID));
  ScalarColumn<double> time_col = ScalarColumn<double>(
      source_table, MSSource::columnName(MSSourceEnums::TIME));
  ScalarColumn<double> interval_col = ScalarColumn<double>(
      source_table, MSSource::columnName(MSSourceEnums::INTERVAL));
  ScalarColumn<int> spectral_window_id_col = ScalarColumn<int>(
      source_table, MSSource::columnName(MSSourceEnums::SPECTRAL_WINDOW_ID));
  ScalarColumn<int> num_lines_col = ScalarColumn<int>(
      source_table, MSSource::columnName(MSSourceEnums::NUM_LINES));
  ScalarColumn<casacore::String> name_col = ScalarColumn<casacore::String>(
      source_table, MSSource::columnName(MSSourceEnums::NAME));
  ScalarColumn<int> calibration_group_col = ScalarColumn<int>(
      source_table, MSSource::columnName(MSSourceEnums::CALIBRATION_GROUP));
  ScalarColumn<casacore::String> code_col = ScalarColumn<casacore::String>(
      source_table, MSSource::columnName(MSSourceEnums::CODE));
  ArrayColumn<double> direction_col = ArrayColumn<double>(
      source_table, MSSource::columnName(MSSourceEnums::DIRECTION));
  ArrayColumn<double> proper_motion_col = ArrayColumn<double>(
      source_table, MSSource::columnName(MSSourceEnums::PROPER_MOTION));

  // Write by appending to the end one row
  size_t row_index = source_table.nrow();
  source_table.addRow();

  source_id_col.put(row_index, source.source_id);
  time_col.put(row_index, source.time);
  interval_col.put(row_index, source.interval);
  spectral_window_id_col.put(row_index, source.spectral_window_id);
  num_lines_col.put(row_index, source.num_lines);
  name_col.put(row_index, source.name);
  calibration_group_col.put(row_index, source.calibration_group);
  code_col.put(row_index, source.code);

  casacore::Vector<double> direction(2);
  direction[0] = source.ra;
  direction[1] = source.dec;
  direction_col.put(row_index, direction);

  casacore::Vector<double> proper_motion(2);
  proper_motion[0] = source.proper_motion[0];
  proper_motion[1] = source.proper_motion[1];
  proper_motion_col.put(row_index, proper_motion);
}

void SubtableWriter::WriteField(const FieldInfo &field) {
  MSField field_table = ms_.field();
  ScalarColumn<casacore::String> name_col = ScalarColumn<casacore::String>(
      field_table, field_table.columnName(MSFieldEnums::NAME));
  ScalarColumn<casacore::String> code_col = ScalarColumn<casacore::String>(
      field_table, field_table.columnName(MSFieldEnums::CODE));
  ScalarColumn<double> time_col = ScalarColumn<double>(
      field_table, field_table.columnName(MSFieldEnums::TIME));
  ScalarColumn<int> num_poly_col = ScalarColumn<int>(
      field_table, field_table.columnName(MSFieldEnums::NUM_POLY));
  ArrayColumn<double> delay_dir_col = ArrayColumn<double>(
      field_table, field_table.columnName(MSFieldEnums::DELAY_DIR));
  ArrayColumn<double> phase_dir_col = ArrayColumn<double>(
      field_table, field_table.columnName(MSFieldEnums::PHASE_DIR));
  ArrayColumn<double> ref_dir_col = ArrayColumn<double>(
      field_table, field_table.columnName(MSFieldEnums::REFERENCE_DIR));
  ScalarColumn<int> source_id_col = ScalarColumn<int>(
      field_table, field_table.columnName(MSFieldEnums::SOURCE_ID));
  ScalarColumn<bool> flag_row_col = ScalarColumn<bool>(
      field_table, field_table.columnName(MSFieldEnums::FLAG_ROW));

  size_t index = field_table.nrow();
  field_table.addRow();

  name_col.put(index, field.name);
  code_col.put(index, field.code);
  time_col.put(index, field.time);
  num_poly_col.put(index, field.num_poly);

  casacore::Array<double> arr(IPosition(2, 2, 1));
  arr(IPosition(2, 0, 0)) = field.delay_direction_ra;
  arr(IPosition(2, 1, 0)) = field.delay_direction_dec;
  delay_dir_col.put(index, arr);

  arr(IPosition(2, 0, 0)) = field.phase_direction_ra;
  arr(IPosition(2, 1, 0)) = field.phase_direction_dec;
  phase_dir_col.put(index, arr);

  arr(IPosition(2, 0, 0)) = field.reference_direction_ra;
  arr(IPosition(2, 1, 0)) = field.reference_direction_dec;
  ref_dir_col.put(index, arr);

  source_id_col.put(index, field.source_id);
  flag_row_col.put(index, field.flag_row);
}

void SubtableWriter::WriteObservation(const ObservationInfo &observation) {
  MSObservation obs_table = ms_.observation();

  ScalarColumn<casacore::String> telescope_name_col(
      obs_table, obs_table.columnName(MSObservationEnums::TELESCOPE_NAME));
  ArrayColumn<double> time_range_col(
      obs_table, obs_table.columnName(MSObservationEnums::TIME_RANGE));
  ScalarColumn<casacore::String> observer_col(
      obs_table, obs_table.columnName(MSObservationEnums::OBSERVER));
  ScalarColumn<casacore::String> schedule_type_col(
      obs_table, obs_table.columnName(MSObservationEnums::SCHEDULE_TYPE));
  ScalarColumn<casacore::String> project_col(
      obs_table, obs_table.columnName(MSObservationEnums::PROJECT));
  ScalarColumn<double> release_date_col(
      obs_table, obs_table.columnName(MSObservationEnums::RELEASE_DATE));
  ScalarColumn<bool> flag_row_col(
      obs_table, obs_table.columnName(MSObservationEnums::FLAG_ROW));

  size_t rowIndex = obs_table.nrow();
  obs_table.addRow();

  casacore::Vector<double> timeRange(2);
  timeRange[0] = observation.start_time;
  timeRange[1] = observation.end_time;

  telescope_name_col.put(rowIndex, observation.telescope_name);
  time_range_col.put(rowIndex, timeRange);
  observer_col.put(rowIndex, observation.observer);
  schedule_type_col.put(rowIndex, observation.schedule_type);
  project_col.put(rowIndex, observation.project);
  release_date_col.put(rowIndex, observation.release_date);
  flag_row_col.put(rowIndex, observation.flag_row);

  // Following are specialized for AARTFAAC,...
  if (observation.telescope_name == "AARTFAAC") {
    obs_table.addColumn((ScalarColumnDesc<casacore::String>(
        observation.telescope_name + "_ANTENNA_TYPE")));
    obs_table.addColumn((ScalarColumnDesc<casacore::Int>(
        observation.telescope_name + "_RCU_MODE")));
    obs_table.addColumn((ScalarColumnDesc<casacore::Int>(
        observation.telescope_name + "_FLAG_WINDOW_SIZE")));
    ScalarColumn<casacore::String> antenna_type_col =
        ScalarColumn<casacore::String>(
            obs_table, observation.telescope_name + "_ANTENNA_TYPE");
    ScalarColumn<int> rcu_mode_col =
        ScalarColumn<int>(obs_table, observation.telescope_name + "_RCU_MODE");
    ScalarColumn<int> flag_window_size_col = ScalarColumn<int>(
        obs_table, observation.telescope_name + "_FLAG_WINDOW_SIZE");
    antenna_type_col.put(rowIndex, observation.antenna_type);
    rcu_mode_col.put(rowIndex, observation.rcu_mode);
    flag_window_size_col.put(rowIndex, observation.flag_window_size);
  }
}

void SubtableWriter::WriteHistoryItem(const std::string &commandLine,
                                      const std::string &application,
                                      const std::vector<std::string> &params) {
  MSHistory history_table = ms_.history();

  casacore::ScalarColumn<double> time_col(
      history_table, MSHistory::columnName(MSHistoryEnums::TIME));
  casacore::ScalarColumn<int> obs_id_col(
      history_table, MSHistory::columnName(MSHistoryEnums::OBSERVATION_ID));
  casacore::ScalarColumn<casacore::String> message_col(
      history_table, MSHistory::columnName(MSHistoryEnums::MESSAGE));
  casacore::ScalarColumn<casacore::String> application_col(
      history_table, MSHistory::columnName(MSHistoryEnums::APPLICATION));
  casacore::ScalarColumn<casacore::String> priority_col(
      history_table, MSHistory::columnName(MSHistoryEnums::PRIORITY));
  casacore::ScalarColumn<casacore::String> origin_col(
      history_table, MSHistory::columnName(MSHistoryEnums::ORIGIN));
  casacore::ArrayColumn<casacore::String> parms_col(
      history_table, MSHistory::columnName(MSHistoryEnums::APP_PARAMS));
  casacore::ArrayColumn<casacore::String> cli_col(
      history_table, MSHistory::columnName(MSHistoryEnums::CLI_COMMAND));

  casacore::Vector<casacore::String> app_params(params.size());
  for (size_t i = 0; i != params.size(); ++i) {
    app_params[i] = params[i];
  }

  casacore::Vector<casacore::String> cliVec(1);
  cliVec[0] = commandLine;

  size_t rowIndex = history_table.nrow();
  const float HOURS_IN_DAY = 24.;
  const float SECONDS_IN_HOURS = 3600.;
  const float DAY_IN_SECONDS = HOURS_IN_DAY * SECONDS_IN_HOURS;

  history_table.addRow();
  time_col.put(rowIndex, casacore::Time().modifiedJulianDay() * DAY_IN_SECONDS);
  obs_id_col.put(rowIndex, 0);
  message_col.put(rowIndex, "Preprocessed");
  application_col.put(rowIndex, application);
  priority_col.put(rowIndex, "NORMAL");
  origin_col.put(rowIndex, "standalone");
  parms_col.put(rowIndex, app_params);
  cli_col.put(rowIndex, cliVec);
};

}  // namespace dp3::base
