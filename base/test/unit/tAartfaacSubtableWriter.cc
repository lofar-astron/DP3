#include <base/AartfaacSubtableWriter.h>

#include <filesystem>

#include <boost/filesystem.hpp>  // for the unique_path generation
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/ScalarColumn.h>

namespace {
template <typename T>
void TestScalarValueFromTable(casacore::Table &table, const std::string &field,
                              const T &expected_value, int row_id = 0) {
  casacore::ScalarColumn<T> table_column(table, field);
  BOOST_CHECK_EQUAL(table_column.get(row_id), expected_value);
}
template <typename T>
std::vector<T> GetArrayColumnFromTable(casacore::Table &table,
                                       const std::string &field,
                                       int row_id = 0) {
  return casacore::ArrayColumn<T>(table, field).get(row_id).tovector();
}

}  // namespace

BOOST_AUTO_TEST_SUITE(aartfaacsubtablewriter)

struct CleanUpFixture {
  CleanUpFixture() {
    const std::filesystem::path tmp_path =
        std::filesystem::temp_directory_path() /
        boost::filesystem::unique_path("tmp%%%%%%%.ms").string();
    path = tmp_path.string();
  };
  ~CleanUpFixture() { std::filesystem::remove_all(path); };

  std::string path;
};

struct WriterFixture : CleanUpFixture {
  WriterFixture() : writer(path){};
  dp3::base::AartfaacSubtableWriter writer;
};

const size_t kDirectionDimensions = 2u;
const size_t kRangeDimensions = 2u;
const size_t kPositionDimensions = 3u;

BOOST_FIXTURE_TEST_CASE(measurementset_create, CleanUpFixture) {
  // Use CleanUpFixture instead of WriterFixture to make
  // it clear we are testing here the constructor of the
  // writer class
  dp3::base::AartfaacSubtableWriter writer(path);

  BOOST_CHECK(std::filesystem::exists(path));
}

BOOST_FIXTURE_TEST_CASE(observation_table, WriterFixture) {
  dp3::base::AartfaacSubtableWriter::ObservationInfo observation;
  observation.telescope_name = "MyArray";
  observation.start_time = 1.0;
  observation.end_time = 3.0;
  observation.observer = "ME";
  observation.flag_row = false;
  observation.flag_window_size = 10;
  observation.project = "LC_TEST";
  observation.release_date = 51232.0;
  observation.schedule_type = "ATYPICAL";

  writer.WriteObservation(observation);

  casacore::MeasurementSet ms(path);

  casacore::Table obs_table = ms.observation();
  TestScalarValueFromTable(obs_table, "TELESCOPE_NAME",
                           casacore::String(observation.telescope_name));
  TestScalarValueFromTable(obs_table, "OBSERVER",
                           casacore::String(observation.observer));
  TestScalarValueFromTable(obs_table, "FLAG_ROW",
                           casacore::Bool(observation.flag_row));
  TestScalarValueFromTable(obs_table, "AARTFAAC_FLAG_WINDOW_SIZE",
                           casacore::Int(observation.flag_window_size));
  TestScalarValueFromTable(obs_table, "PROJECT",
                           casacore::String(observation.project));
  TestScalarValueFromTable(obs_table, "RELEASE_DATE",
                           casacore::Double(observation.release_date));
  TestScalarValueFromTable(obs_table, "SCHEDULE_TYPE",
                           casacore::String(observation.schedule_type));

  casacore::ArrayColumn<casacore::Double> time_range(obs_table, "TIME_RANGE");

  std::vector<double> stored_time_range =
      GetArrayColumnFromTable<casacore::Double>(obs_table, "TIME_RANGE");
  BOOST_REQUIRE_EQUAL(stored_time_range.size(), kRangeDimensions);
  BOOST_CHECK_EQUAL(stored_time_range[0], observation.start_time);
  BOOST_CHECK_EQUAL(stored_time_range[1], observation.end_time);
}

BOOST_FIXTURE_TEST_CASE(antenna_table, WriterFixture) {
  const size_t kNAntennas = 2;

  std::vector<dp3::base::AartfaacSubtableWriter::AntennaInfo> antennas(
      kNAntennas);

  for (size_t idx = 0u; idx < kNAntennas; ++idx) {
    antennas[idx].flag = idx % 2;
    antennas[idx].mount = "ALTAX";
    antennas[idx].name = "MY_ANTN" + std::to_string(idx);
    antennas[idx].type = "LBLA";
    antennas[idx].x = 1.0 * idx;
    antennas[idx].y = 2.0 * idx;
    antennas[idx].z = 3.0 * idx;
    antennas[idx].diameter = 1.23;
  }

  std::array<double, 9> axes = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  double time = 123.45;
  writer.WriteAntennae(antennas, axes, time);

  casacore::MeasurementSet ms(writer.GetPath());
  casacore::Table antenna_table = ms.antenna();
  BOOST_CHECK_EQUAL(antenna_table.nrow(), kNAntennas);

  for (size_t idx = 0u; idx < kNAntennas; ++idx) {
    TestScalarValueFromTable(antenna_table, "MOUNT",
                             casacore::String(antennas[idx].mount), idx);
    TestScalarValueFromTable(antenna_table, "NAME",
                             casacore::String(antennas[idx].name), idx);
    TestScalarValueFromTable(antenna_table, "DISH_DIAMETER",
                             casacore::Double(antennas[idx].diameter), idx);
    TestScalarValueFromTable(antenna_table, "FLAG_ROW",
                             casacore::Bool(antennas[idx].flag), idx);

    std::vector<casacore::Double> stored_positions =
        GetArrayColumnFromTable<casacore::Double>(antenna_table, "POSITION",
                                                  idx);
    BOOST_REQUIRE_EQUAL(stored_positions.size(), kPositionDimensions);
    BOOST_CHECK_EQUAL(stored_positions[0], antennas[idx].x);
    BOOST_CHECK_EQUAL(stored_positions[1], antennas[idx].y);
    BOOST_CHECK_EQUAL(stored_positions[2], antennas[idx].z);
  }
  casacore::Array<casacore::Double> coordinates_axes;

  antenna_table.keywordSet().get("AARTFAAC_COORDINATE_AXES", coordinates_axes);
  std::vector<double> stored_coordinate_axes = coordinates_axes.tovector();

  BOOST_CHECK_EQUAL_COLLECTIONS(stored_coordinate_axes.begin(),
                                stored_coordinate_axes.end(), axes.begin(),
                                axes.end());
}

BOOST_FIXTURE_TEST_CASE(write_band_info, WriterFixture) {
  const int kNChannels = 5;
  const double kReferenceFrequency = 1.0e6;
  const double kBandwidth = 5.0e6;
  const bool kFlag = true;
  const std::string kBandName = "test_band";

  std::vector<dp3::base::AartfaacSubtableWriter::ChannelInfo> channels(
      kNChannels);

  for (int idx = 0; idx < kNChannels; ++idx) {
    channels[idx].channel_frequency = 40.0e6 + 500.0e5 * idx;
    channels[idx].channel_width = 1.0e6;
    channels[idx].effective_bandwidth = 1.0e6;
    channels[idx].resolution = 1.0e6;
  }

  writer.WriteBandInfo("test_band", channels, kReferenceFrequency, kBandwidth,
                       kFlag);

  casacore::MeasurementSet ms(writer.GetPath());
  casacore::Table window_table = ms.spectralWindow();
  BOOST_CHECK_EQUAL(window_table.nrow(), 1);

  TestScalarValueFromTable(window_table, "NUM_CHAN", casacore::Int(kNChannels));
  TestScalarValueFromTable(window_table, "NAME", casacore::String(kBandName));
  TestScalarValueFromTable(window_table, "REF_FREQUENCY",
                           casacore::Double(kReferenceFrequency));
  TestScalarValueFromTable(window_table, "TOTAL_BANDWIDTH",
                           casacore::Double(kBandwidth));
  TestScalarValueFromTable(window_table, "FLAG_ROW", casacore::Bool(kFlag));

  std::vector<casacore::Double> channels_frequencies =
      GetArrayColumnFromTable<casacore::Double>(window_table, "CHAN_FREQ");
  std::vector<casacore::Double> channels_width =
      GetArrayColumnFromTable<casacore::Double>(window_table, "CHAN_WIDTH");
  std::vector<casacore::Double> channels_effective_bw =
      GetArrayColumnFromTable<casacore::Double>(window_table, "EFFECTIVE_BW");
  std::vector<casacore::Double> channels_resolution =
      GetArrayColumnFromTable<casacore::Double>(window_table, "RESOLUTION");

  BOOST_REQUIRE_EQUAL(channels_frequencies.size(), kNChannels);
  BOOST_REQUIRE_EQUAL(channels_width.size(), kNChannels);
  BOOST_REQUIRE_EQUAL(channels_effective_bw.size(), kNChannels);
  BOOST_REQUIRE_EQUAL(channels_resolution.size(), kNChannels);

  for (int idx = 0; idx < kNChannels; ++idx) {
    BOOST_CHECK_EQUAL(channels[idx].channel_frequency,
                      channels_frequencies[idx]);
    BOOST_CHECK_EQUAL(channels[idx].channel_width, channels_width[idx]);
    BOOST_CHECK_EQUAL(channels[idx].effective_bandwidth,
                      channels_effective_bw[idx]);
    BOOST_CHECK_EQUAL(channels[idx].resolution, channels_resolution[idx]);
  }
}

BOOST_FIXTURE_TEST_CASE(polarization, WriterFixture) {
  const bool kFlag = true;
  const int kNCorrelations = 4;
  writer.WriteLinearPolarizations(kFlag);

  casacore::MeasurementSet ms(writer.GetPath());
  casacore::Table pol_table = ms.polarization();
  TestScalarValueFromTable(pol_table, "FLAG_ROW", casacore::Bool(kFlag));
  TestScalarValueFromTable(pol_table, "NUM_CORR",
                           casacore::Int(kNCorrelations));
}

BOOST_FIXTURE_TEST_CASE(source, WriterFixture) {
  dp3::base::AartfaacSubtableWriter::SourceInfo sinfo;
  sinfo.source_id = -1;
  sinfo.time = 54321;
  sinfo.interval = 1;
  sinfo.spectral_window_id = 1;
  sinfo.num_lines = 1;
  sinfo.name = "CasA";
  sinfo.code = "A1";
  sinfo.calibration_group = 1;
  sinfo.ra = 1.23;
  sinfo.dec = 2.3;
  sinfo.proper_motion[0] = 0.1;
  sinfo.proper_motion[0] = 0.3;

  writer.WriteSource(sinfo);

  casacore::MeasurementSet ms(writer.GetPath());
  casacore::Table source_table = ms.source();
  TestScalarValueFromTable(source_table, "SOURCE_ID",
                           casacore::Int(sinfo.source_id));
  TestScalarValueFromTable(source_table, "TIME", casacore::Double(sinfo.time));
  TestScalarValueFromTable(source_table, "INTERVAL",
                           casacore::Double(sinfo.interval));
  TestScalarValueFromTable(source_table, "SPECTRAL_WINDOW_ID",
                           casacore::Int(sinfo.spectral_window_id));
  TestScalarValueFromTable(source_table, "NUM_LINES",
                           casacore::Int(sinfo.num_lines));
  TestScalarValueFromTable(source_table, "NAME", casacore::String(sinfo.name));
  TestScalarValueFromTable(source_table, "CODE", casacore::String(sinfo.code));
  TestScalarValueFromTable(source_table, "CALIBRATION_GROUP",
                           casacore::Int(sinfo.calibration_group));
  std::vector<casacore::Double> stored_direction =
      GetArrayColumnFromTable<casacore::Double>(source_table, "DIRECTION");
  std::vector<casacore::Double> stored_motion =
      GetArrayColumnFromTable<casacore::Double>(source_table, "PROPER_MOTION");

  BOOST_REQUIRE_EQUAL(stored_direction.size(), kDirectionDimensions);
  BOOST_REQUIRE_EQUAL(stored_motion.size(), kDirectionDimensions);
  BOOST_CHECK_EQUAL(stored_direction[0], sinfo.ra);
  BOOST_CHECK_EQUAL(stored_direction[1], sinfo.dec);
  BOOST_CHECK_EQUAL(stored_motion[0], sinfo.proper_motion[0]);
  BOOST_CHECK_EQUAL(stored_motion[1], sinfo.proper_motion[1]);
}

BOOST_FIXTURE_TEST_CASE(field, WriterFixture) {
  dp3::base::AartfaacSubtableWriter::FieldInfo finfo;
  finfo.name = "test-field-name";
  finfo.code = "field_01";
  finfo.time = 54321;
  finfo.num_poly = 1;
  finfo.delay_direction_ra = 1.23;
  finfo.delay_direction_dec = 2.23;
  finfo.phase_direction_ra = 3.12;
  finfo.phase_direction_dec = 1.23;
  finfo.reference_direction_ra = 0.23;
  finfo.reference_direction_dec = 0.32;
  finfo.source_id = 1;
  finfo.flag_row = true;

  writer.WriteField(finfo);

  casacore::MeasurementSet ms(writer.GetPath());
  casacore::Table field_table = ms.field();

  TestScalarValueFromTable(field_table, "NAME", casacore::String(finfo.name));
  TestScalarValueFromTable(field_table, "CODE", casacore::String(finfo.code));
  TestScalarValueFromTable(field_table, "TIME", casacore::Double(finfo.time));
  TestScalarValueFromTable(field_table, "NUM_POLY",
                           casacore::Int(finfo.num_poly));
  TestScalarValueFromTable(field_table, "SOURCE_ID",
                           casacore::Int(finfo.source_id));
  TestScalarValueFromTable(field_table, "FLAG_ROW",
                           casacore::Bool(finfo.flag_row));

  std::vector<casacore::Double> delay_pointing =
      GetArrayColumnFromTable<casacore::Double>(field_table, "DELAY_DIR");

  BOOST_REQUIRE_EQUAL(delay_pointing.size(), kDirectionDimensions);
  BOOST_CHECK_EQUAL(delay_pointing[0], finfo.delay_direction_ra);
  BOOST_CHECK_EQUAL(delay_pointing[1], finfo.delay_direction_dec);

  std::vector<casacore::Double> phase_pointing =
      GetArrayColumnFromTable<casacore::Double>(field_table, "PHASE_DIR");
  BOOST_REQUIRE_EQUAL(phase_pointing.size(), kDirectionDimensions);
  BOOST_CHECK_EQUAL(phase_pointing[0], finfo.phase_direction_ra);
  BOOST_CHECK_EQUAL(phase_pointing[1], finfo.phase_direction_dec);

  std::vector<casacore::Double> reference_pointing =
      GetArrayColumnFromTable<casacore::Double>(field_table, "REFERENCE_DIR");
  BOOST_REQUIRE_EQUAL(reference_pointing.size(), kDirectionDimensions);
  BOOST_CHECK_EQUAL(reference_pointing[0], finfo.reference_direction_ra);
  BOOST_CHECK_EQUAL(reference_pointing[1], finfo.reference_direction_dec);
}

BOOST_FIXTURE_TEST_CASE(history, WriterFixture) {
  const std::string kCommandLineString = "DP3 one two --test-optional";
  const std::string kApplication = "DP3";
  const std::vector<std::string> kParams = {"one", "two", "--test-optional"};
  const size_t kCliCommandSize = 1;
  writer.WriteHistoryItem(kCommandLineString, kApplication, kParams);

  casacore::MeasurementSet ms(writer.GetPath());
  casacore::Table history_table = ms.history();

  std::vector<casacore::String> params =
      GetArrayColumnFromTable<casacore::String>(history_table, "APP_PARAMS");
  std::vector<casacore::String> cli_command =
      GetArrayColumnFromTable<casacore::String>(history_table, "CLI_COMMAND");

  BOOST_REQUIRE_EQUAL(cli_command.size(), kCliCommandSize);
  BOOST_CHECK_EQUAL(cli_command.back(), kCommandLineString);
  TestScalarValueFromTable<casacore::String>(history_table, "APPLICATION",
                                             kApplication);
  BOOST_CHECK_EQUAL_COLLECTIONS(kParams.begin(), kParams.end(), params.begin(),
                                params.end());
}

BOOST_AUTO_TEST_SUITE_END()
