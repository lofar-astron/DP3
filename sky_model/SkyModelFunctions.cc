#include "ReadSkyModel.h"
#include "SkyModelFunctions.h"
#include "SkyModelSelection.h"

#include "base/PointSource.h"
#include "base/GaussianSource.h"

#include "common/ParameterValue.h"
#include "common/ProximityClustering.h"

#include <cassert>
#include <cmath>
#include <mutex>
#include <numeric>
#include <set>
#include <vector>

namespace dp3::sky_model {

using base::Direction;
using base::GaussianSource;
using base::ModelComponent;
using base::PointSource;
using base::Stokes;

using sky_model::Source;
using sky_model::SourceInfo;

static PointSource::Ptr MakePointSource(const Source& src) {
  // Fetch direction.
  if (src.getInfo().getRefType() != "J2000")
    throw std::runtime_error("Reference type should be J2000");
  const Direction direction(src.getRa(), src.getDec());

  // Fetch stokes vector.
  Stokes stokes;
  stokes.I = src.getI();
  stokes.V = src.getV();
  if (!src.getInfo().getUseRotationMeasure()) {
    stokes.Q = src.getQ();
    stokes.U = src.getU();
  }

  PointSource::Ptr source;
  switch (src.getInfo().getType()) {
    case SourceInfo::POINT: {
      source = PointSource::Ptr(new PointSource(direction, stokes));
    } break;

    case SourceInfo::GAUSSIAN: {
      GaussianSource::Ptr gauss(new GaussianSource(direction, stokes));

      const double deg2rad = (M_PI / 180.0);
      gauss->setPositionAngle(src.getOrientation() * deg2rad);
      gauss->setPositionAngleIsAbsolute(
          src.getInfo().getPositionAngleIsAbsolute());

      const double arcsec2rad = (M_PI / 3600.0) / 180.0;
      gauss->setMajorAxis(src.getMajorAxis() * arcsec2rad);
      gauss->setMinorAxis(src.getMinorAxis() * arcsec2rad);
      source = gauss;
    } break;

    default: {
      throw std::runtime_error(
          "Only point sources and Gaussian sources are"
          " supported at this time.");
    }
  }

  // Fetch spectral index attributes (if applicable).
  bool isLogarithmic = src.getInfo().getHasLogarithmicSI();
  if (src.getSpectralTerms().size() > 0) {
    source->setSpectralTerms(src.getInfo().getSpectralTermsRefFreq(),
                             isLogarithmic, src.getSpectralTerms().begin(),
                             src.getSpectralTerms().end());
  }

  // Fetch rotation measure attributes (if applicable).
  if (src.getInfo().getUseRotationMeasure()) {
    source->setRotationMeasure(src.getPolarizedFraction(),
                               src.getPolarizationAngle(),
                               src.getRotationMeasure());
  }

  return source;
}

std::vector<std::shared_ptr<Patch>> MakePatches(
    const sky_model::SkyModel& source_db,
    const std::vector<std::string>& patch_names) {
  // Create a component list for each patch name.
  std::vector<std::vector<std::shared_ptr<ModelComponent>>> componentsList(
      patch_names.size());

  for (const auto& source : source_db.GetSources()) {
    // Use the source if its patch matches a patch name.
    for (size_t i = 0; i < patch_names.size(); ++i) {
      if (source.getPatchName() == patch_names[i]) {
        componentsList[i].push_back(MakePointSource(source));
        break;
      }
    }
  }

  std::vector<std::shared_ptr<Patch>> patchList;
  patchList.reserve(componentsList.size());
  for (unsigned int i = 0; i < componentsList.size(); ++i) {
    if (componentsList[i].empty())
      throw std::runtime_error("No sources found for patch " + patch_names[i]);
    auto ppatch = std::make_shared<Patch>(
        patch_names[i], componentsList[i].begin(), componentsList[i].end());

    const sky_model::PatchInfo& info = source_db.GetPatch(patch_names[i]);
    // Set the direction and apparent flux of the patch.
    const Direction patchDirection(info.getRa(), info.getDec());
    ppatch->SetDirection(patchDirection);
    ppatch->SetBrightness(info.apparentBrightness());
    ///    ppatch->computeDirection();
    patchList.push_back(std::move(ppatch));
  }
  return patchList;
}

std::vector<std::pair<std::shared_ptr<ModelComponent>, std::shared_ptr<Patch>>>
makeSourceList(std::vector<std::shared_ptr<Patch>>& patchList) {
  const size_t nSources =
      std::accumulate(patchList.begin(), patchList.end(), 0,
                      [](size_t left, const std::shared_ptr<Patch>& right) {
                        return left + right->NComponents();
                      });

  std::vector<
      std::pair<std::shared_ptr<ModelComponent>, std::shared_ptr<Patch>>>
      sourceList;
  sourceList.reserve(nSources);

  for (const std::shared_ptr<Patch>& patch : patchList) {
    for (const std::shared_ptr<ModelComponent>& component : *patch) {
      sourceList.emplace_back(component, patch);
    }
  }

  return sourceList;
}

std::vector<std::shared_ptr<Patch>> makeOnePatchPerComponent(
    const std::vector<std::shared_ptr<Patch>>& patchList) {
  size_t numComponents = 0;

  for (const auto& patch : patchList) {
    numComponents += patch->NComponents();
  }

  std::vector<std::shared_ptr<Patch>> largePatchList;
  largePatchList.reserve(numComponents);

  for (const std::shared_ptr<Patch>& patch : patchList) {
    size_t compNum = 0;
    for (auto compIt = patch->begin(); compIt != patch->end(); ++compIt) {
      std::shared_ptr<Patch> ppatch = std::make_shared<Patch>(
          patch->Name() + "_" + std::to_string(compNum), compIt, compIt + 1);
      ppatch->SetDirection((*compIt)->direction());
      largePatchList.push_back(std::move(ppatch));
      ++compNum;
    }
  }

  return largePatchList;
}

std::vector<std::shared_ptr<Patch>> clusterProximateSources(
    const std::vector<std::shared_ptr<Patch>>& patch_list,
    double proximity_limit) {
  // Create a list of all source positions and a lookup table to go from index
  // back to patch & comp
  std::vector<std::pair<double, double>> sources;
  std::vector<std::pair<size_t, size_t>> lookup_table;
  size_t patchIndex = 0;
  for (const std::shared_ptr<Patch>& patch : patch_list) {
    size_t compIndex = 0;
    for (const auto& comp : *patch) {
      sources.emplace_back(comp->direction().ra, comp->direction().dec);
      lookup_table.emplace_back(patchIndex, compIndex);
      ++compIndex;
    }
    ++patchIndex;
  }

  // Call the clustering algorithm
  common::ProximityClustering clustering(sources);
  std::vector<std::vector<size_t>> clusters = clustering.Group(proximity_limit);

  // Make a new patch list from the results
  std::vector<std::shared_ptr<Patch>> clusteredPatchList;
  for (const auto& groups : clusters) {
    std::vector<std::shared_ptr<ModelComponent>> componentsInPatch;
    double alpha = 0.0;
    double delta = 0.0;
    for (const auto& component : groups) {
      auto lookupIndices = lookup_table[component];
      const auto& patch = patch_list[lookupIndices.first];
      const std::shared_ptr<ModelComponent>& comp =
          *(patch->begin() + lookupIndices.second);
      alpha += comp->direction().ra;
      delta += comp->direction().dec;
      componentsInPatch.emplace_back(comp);
    }
    const std::shared_ptr<Patch>& firstPatch =
        patch_list[lookup_table[groups.front()].first];
    std::shared_ptr<Patch> ppatch = std::make_shared<Patch>(
        firstPatch->Name() + "_" + std::to_string(clusteredPatchList.size()),
        componentsInPatch.begin(), componentsInPatch.end());
    ppatch->SetDirection(
        Direction(alpha / groups.size(), delta / groups.size()));
    clusteredPatchList.push_back(std::move(ppatch));
  }

  return clusteredPatchList;
}

std::vector<std::string> MakePatchList(
    const sky_model::SkyModel& source_db,
    const std::vector<std::string>& patterns) {
  if (patterns.empty()) return source_db.GetPatchNames();

  std::set<std::string> patches;
  for (const auto& pattern : patterns) {
    if (pattern.empty()) continue;

    if (pattern[0] == '@')
      patches.insert(pattern);
    else {
      const std::vector<std::string> match = source_db.FindPatches(pattern);
      patches.insert(match.begin(), match.end());
    }
  }
  return std::vector<std::string>(patches.begin(), patches.end());
}

std::vector<std::vector<std::string>> MakeDirectionList(
    const std::vector<std::string>& packed_directions,
    const std::string& sky_model_filename) {
  std::vector<std::vector<std::string>> directions;

  if (packed_directions.empty() && !sky_model_filename.empty()) {
    // Use all patches from the skymodel if the user did not give directions.
    SkyModelSelection sky_model(ReadSkyModel(sky_model_filename));
    sky_model.SelectAllPatches();
    const std::vector<std::shared_ptr<Patch>> patches =
        std::move(sky_model).MakePatchList();

    for (const auto& patch : patches) {
      directions.emplace_back(1, patch->Name());
    }
  } else {
    for (const std::string& direction : packed_directions) {
      directions.push_back(common::ParameterValue(direction).getStringVector());
    }
  }

  return directions;
}

bool CheckPolarized(const sky_model::SkyModel& source_db,
                    const std::vector<std::string>& patch_names) {
  for (const auto& source : source_db.GetSources()) {
    const std::string& source_patch_name = source.getPatchName();
    for (const auto& patch_name : patch_names)
      if (patch_name == source_patch_name)
        if (source.getV() != 0.0 || source.getQ() != 0.0 ||
            source.getU() != 0.0)
          return true;
  }
  return false;
}

bool CheckAnyOrientationIsAbsolute(
    const sky_model::SkyModel& source_db,
    const std::vector<std::string>& patch_names) {
  for (const auto& source : source_db.GetSources()) {
    const std::string& source_patch_name = source.getPatchName();
    for (const auto& patch_name : patch_names)
      if (patch_name == source_patch_name)
        if (source.getInfo().getPositionAngleIsAbsolute()) {
          return true;
        }
  }
  return false;
}

void SetPatchIndices(std::vector<std::shared_ptr<Patch>>& patch_list) {
  for (size_t index = 0; index != patch_list.size(); ++index)
    patch_list[index]->SetIndex(index);
}

}  // namespace dp3::sky_model
