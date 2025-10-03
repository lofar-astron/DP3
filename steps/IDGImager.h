// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_IDGIMAGER_H_
#define DP3_STEPS_IDGIMAGER_H_

#ifdef HAVE_IDG
#include <idg-api.h>
#endif

#include <dp3/steps/Step.h>

#include "../common/ParameterSet.h"
#include "../common/Timer.h"
namespace dp3::steps {

class IDGImager : public Step {
 public:
  IDGImager(const common::ParameterSet&, const std::string& prefix);

  common::Fields getRequiredFields() const override {
    return kWeightsField | kUvwField | kDataField | kFlagsField;
  }

  common::Fields getProvidedFields() const override { return kDataField; }
  /// Process the data. The dummy step forwards the data to its next step.
  bool process(std::unique_ptr<base::DPBuffer>) override;

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Update the general info.
  void updateInfo(const base::DPInfo&) override;

  /// Show the step parameters.
  void show(std::ostream&) const override;

  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override;

  size_t GetImageSize() { return grid_size_; }

 private:
  void ConfigureIDG();
  void InitBuffers(const double max_w, const double max_baseline);

#ifdef HAVE_IDG
  static idg::api::Type ProxyTypeFromString(const std::string& idg_type);
  static void CopyVisibilitiesToBuffer(
      idg::api::GridderBuffer& buffer, const std::vector<int>& antenna1,
      const std::vector<int>& antenna2, const base::DPBuffer::UvwType& uvw,
      const base::DPBuffer::WeightsType& weights,
      const base::DPBuffer::DataType& visibilities);
#endif

  static void PrepareVisibilities(base::DPBuffer::DataType& visibilities,
                                  base::DPBuffer::WeightsType& weights,
                                  const base::DPBuffer::FlagsType& flags);

  static xt::xtensor<float, 2> ComputeWeightMap(
      size_t grid_size_, const base::DPBuffer::WeightsType& weights,
      const xt::xtensor<size_t, 3>& uv_pixels);

  static void ComputeUVPixels(const base::DPBuffer::UvwType& uvw,
                              const float uv_max,
                              const std::vector<double>& frequencies,
                              xt::xtensor<size_t, 3>& uv_pixels,
                              size_t grid_size);

  static void ApplyWeights(const xt::xtensor<float, 2>& weight_map,
                           const xt::xtensor<size_t, 3>& uv_pixels,
                           base::DPBuffer::DataType& visibilities);

  static void WeightVisibilities(const base::DPBuffer::UvwType& uvw,
                                 const base::DPBuffer::WeightsType& weights,
                                 const std::vector<double> frequencies,
                                 const size_t grid_size,
                                 base::DPBuffer::DataType& visibilities);

  static void WriteFITS(const std::string file_path, const base::DPInfo& info,
                        double time, double dl, double dm, double pixel_scale,
                        xt::xtensor<double, 3> image);

  static inline double DegToRad(double angle);
  static std::string FormatName(std::string name, int time);

  std::string name_;
  common::NSTimer timer_;
  size_t grid_size_;
  float pixel_scale_;
  double dl_;
  double dm_;
  std::string image_name_;
  int buffer_size_;
  int current_image_idx_ = 0;
  bool init_buffers_ = false;

#ifdef HAVE_IDG
  idg::api::Type proxy_;
  std::unique_ptr<idg::api::BufferSet> buffer_set_;
  idg::api::options_type options_;
#endif
};

}  // namespace dp3::steps

#endif
