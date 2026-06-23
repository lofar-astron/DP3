// Copyright (C) 2026 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include <EveryBeam/telescope/telescope.h>
#include <EveryBeam/griddedresponse/griddedresponse.h>
#include <EveryBeam/pointresponse/pointresponse.h>

namespace dp3::test {

class MockTelescope;

/// The class returned by MockTelescope::GetPointResponse().
class MockPointResponse : public everybeam::pointresponse::PointResponse {
 public:
  MockPointResponse(const MockTelescope& telescope);

  ~MockPointResponse() override;

  /// Response calls fail, since DP3 only calls the other overload.
  void Response(everybeam::BeamMode beam_mode,
                std::complex<float>* response_matrix, double ra, double dec,
                double freq, size_t station_id, size_t field_id) override {
    BOOST_TEST(false, "Unexpected MockPointResponse::Response call.");
  }

  /// This overload checks input values and returns a fixed response value.
  aocommon::MC2x2 Response(everybeam::BeamMode beam_mode, size_t station_idx,
                           double freq, const everybeam::vector3r_t& direction,
                           std::mutex* mutex = nullptr) override;

 private:
  std::size_t response_count_ = 0;
};

class MockTelescope : public everybeam::telescope::Telescope {
 public:
  /**
   * Constructor.
   * @param expected_frequencies List of expected frequencies for
   * MockPointResponse::Reponse() calls: The first call should receive the first
   * frequency etc. The MockPointResponse destructor checks the call count.
   * @param expected_direction Expected direction for
   * MockPointResponse::Response() calls. MockPointResponse currently assumes
   * all calls receive this direction.
   */
  MockTelescope(const std::vector<double>& expected_frequencies,
                const everybeam::vector3r_t& expected_direction)
      : everybeam::telescope::Telescope(0, {}),
        expected_frequencies_(expected_frequencies),
        expected_direction_(expected_direction) {}

  std::unique_ptr<everybeam::griddedresponse::GriddedResponse>
  GetGriddedResponse(const aocommon::CoordinateSystem&) const override {
    BOOST_TEST(false, "Unexpected MockTelescope::GetGriddedResponse call.");
    return nullptr;
  }

  std::unique_ptr<everybeam::pointresponse::PointResponse> GetPointResponse(
      double) const override {
    return std::make_unique<MockPointResponse>(*this);
  }

  // Emulate a homogeneous telescope, so SelectStationIndices() works.
  bool IsHomogeneous() const override { return true; }

  const std::vector<double>& ExpectedFrequencies() const {
    return expected_frequencies_;
  }

  const everybeam::vector3r_t& ExpectedDirection() const {
    return expected_direction_;
  }

 private:
  const std::vector<double>& expected_frequencies_;
  const everybeam::vector3r_t& expected_direction_;
};

}  // namespace dp3::test
