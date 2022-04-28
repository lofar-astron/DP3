// SourceDBUtil.cc: Helper functions to extract patch and source information
// from a SourceDB.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SourceDBUtil.h"

#include "Exceptions.h"
#include "PointSource.h"
#include "GaussianSource.h"

#include "../parmdb/SourceDB.h"
#include "../parmdb/SkymodelToSourceDB.h"

#include "../common/ParameterValue.h"
#include "../common/ProximityClustering.h"

#include <cassert>
#include <mutex>
#include <numeric>
#include <set>
#include <string_view>
#include <vector>

namespace dp3 {
namespace base {
using parmdb::SourceData;
using parmdb::SourceInfo;

static PointSource::Ptr MakePointSource(const SourceData& src) {
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

      const double deg2rad = (casacore::C::pi / 180.0);
      gauss->setPositionAngle(src.getOrientation() * deg2rad);
      gauss->setPositionAngleIsAbsolute(
          src.getInfo().getPositionAngleIsAbsolute());

      const double arcsec2rad = (casacore::C::pi / 3600.0) / 180.0;
      gauss->setMajorAxis(src.getMajorAxis() * arcsec2rad);
      gauss->setMinorAxis(src.getMinorAxis() * arcsec2rad);
      source = gauss;
    } break;

    default: {
      throw Exception(
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

std::vector<std::shared_ptr<Patch>> makePatches(
    parmdb::SourceDB& sourceDB, const std::vector<string>& patchNames,
    unsigned int nModel) {
  // Create a component list for each patch name.
  std::vector<std::vector<std::shared_ptr<ModelComponent>>> componentsList(
      nModel);

  // Loop over all sources.
  sourceDB.lock();
  sourceDB.rewind();
  SourceData src;
  while (!sourceDB.atEnd()) {
    sourceDB.getNextSource(src);
    // Use the source if its patch matches a patch name.
    for (unsigned int i = 0; i < nModel; ++i) {
      if (src.getPatchName() == patchNames[i]) {
        componentsList[i].push_back(MakePointSource(src));
        break;
      }
    }
  }
  sourceDB.unlock();

  std::vector<std::shared_ptr<Patch>> patchList;
  patchList.reserve(componentsList.size());
  for (unsigned int i = 0; i < componentsList.size(); ++i) {
    if (componentsList[i].empty())
      throw Exception("No sources found for patch " + patchNames[i]);
    auto ppatch = std::make_shared<Patch>(
        patchNames[i], componentsList[i].begin(), componentsList[i].end());
    std::vector<parmdb::PatchInfo> patchInfo(
        sourceDB.getPatchInfo(-1, patchNames[i]));
    if (patchInfo.size() != 1)
      throw std::runtime_error("Patch \"" + patchNames[i] +
                               "\" defined more than once in SourceDB");
    // Set the direction and apparent flux of the patch.
    const Direction patchDirection(patchInfo[0].getRa(), patchInfo[0].getDec());
    ppatch->setDirection(patchDirection);
    ppatch->setBrightness(patchInfo[0].apparentBrightness());
    ///    ppatch->computeDirection();
    patchList.push_back(std::move(ppatch));
  }
  return patchList;
}

std::vector<std::shared_ptr<Patch>> MakePatches(
    const parmdb::SourceDBSkymodel& source_db,
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
      throw Exception("No sources found for patch " + patch_names[i]);
    auto ppatch = std::make_shared<Patch>(
        patch_names[i], componentsList[i].begin(), componentsList[i].end());

    const parmdb::PatchInfo& info = source_db.GetPatch(patch_names[i]);
    // Set the direction and apparent flux of the patch.
    const Direction patchDirection(info.getRa(), info.getDec());
    ppatch->setDirection(patchDirection);
    ppatch->setBrightness(info.apparentBrightness());
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
                        return left + right->nComponents();
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
    numComponents += patch->nComponents();
  }

  std::vector<std::shared_ptr<Patch>> largePatchList;
  largePatchList.reserve(numComponents);

  for (const std::shared_ptr<Patch>& patch : patchList) {
    size_t compNum = 0;
    for (auto compIt = patch->begin(); compIt != patch->end(); ++compIt) {
      std::shared_ptr<Patch> ppatch(new Patch(
          patch->name() + "_" + std::to_string(compNum), compIt, compIt + 1));
      ppatch->setDirection((*compIt)->direction());
      largePatchList.push_back(std::move(ppatch));
      compNum++;
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

  // Make a new patch list according from the results
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
    std::shared_ptr<Patch> ppatch(new Patch(
        firstPatch->name() + "_" + std::to_string(clusteredPatchList.size()),
        componentsInPatch.begin(), componentsInPatch.end()));
    ppatch->setDirection(
        Direction(alpha / groups.size(), delta / groups.size()));
    clusteredPatchList.push_back(std::move(ppatch));
  }

  return clusteredPatchList;
}

std::vector<string> makePatchList(parmdb::SourceDB& sourceDB,
                                  std::vector<string> patterns) {
  if (patterns.empty()) {
    patterns.push_back("*");
  }

  std::set<string> patches;
  std::vector<string>::iterator it = patterns.begin();
  while (it != patterns.end()) {
    if (!it->empty() && (*it)[0] == '@') {
      patches.insert(*it);
      it = patterns.erase(it);
    } else {
      std::vector<string> match(sourceDB.getPatches(-1, *it));
      patches.insert(match.begin(), match.end());
      ++it;
    }
  }

  return std::vector<string>(patches.begin(), patches.end());
}

std::vector<std::string> MakePatchList(
    const parmdb::SourceDBSkymodel& source_db,
    const std::vector<std::string>& patterns) {
  if (patterns.empty()) return source_db.FindPatches("*");

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
    const std::string& source_db_filename) {
  std::vector<std::vector<std::string>> directions;

  if (packed_directions.empty() && !source_db_filename.empty()) {
    // Use all patches from the SourceDB if the user did not give directions.
    const std::vector<std::shared_ptr<Patch>> patches =
        SourceDB(source_db_filename, std::vector<std::string>{},
                 SourceDB::FilterMode::kPattern)
            .MakePatchList();

    for (const auto& patch : patches) {
      directions.emplace_back(1, patch->name());
    }
  } else {
    for (const std::string& direction : packed_directions) {
      directions.push_back(common::ParameterValue(direction).getStringVector());
    }
  }

  return directions;
}

bool checkPolarized(parmdb::SourceDB& sourceDB,
                    const std::vector<string>& patchNames,
                    unsigned int nModel) {
  bool polarized = false;

  // Loop over all sources.
  const std::lock_guard<parmdb::SourceDB> lock{sourceDB};
  sourceDB.rewind();
  SourceData src;
  while (!sourceDB.atEnd()) {
    sourceDB.getNextSource(src);
    // Use the source if its patch matches a patch name.
    for (unsigned int i = 0; i < nModel; ++i) {
      if (src.getPatchName() == patchNames[i]) {
        // Determine whether source is unpolarized.
        if (src.getV() > 0.0 || src.getQ() > 0.0 || src.getU() > 0.0) {
          polarized = true;
          break;
        }
      }
    }
    if (polarized) {
      break;
    }
  }
  return polarized;
}

bool CheckPolarized(const parmdb::SourceDBSkymodel& source_db,
                    const std::vector<std::string>& patch_names) {
  for (const auto& source : source_db.GetSources()) {
    const std::string& source_patch_name = source.getPatchName();
    for (const auto& patch_name : patch_names)
      if (patch_name == source_patch_name)
        if (source.getV() > 0.0 || source.getQ() > 0.0 || source.getU() > 0.0)
          return true;
  }
  return false;
}

bool CheckAnyOrientationIsAbsolute(parmdb::SourceDB& source_db,
                                   const std::vector<string>& patch_names,
                                   unsigned int n_model) {
  bool anyOrientationIsAbsolute = false;

  // Loop over all sources.
  const std::lock_guard<parmdb::SourceDB> lock{source_db};
  source_db.rewind();
  SourceData src;
  while (!source_db.atEnd()) {
    source_db.getNextSource(src);
    // Use the source if its patch matches a patch name.
    for (unsigned int i = 0; i < n_model; ++i) {
      if (src.getPatchName() == patch_names[i]) {
        if (src.getInfo().getPositionAngleIsAbsolute()) {
          return true;
        }
      }
    }
  }
  return false;
}

bool CheckAnyOrientationIsAbsolute(
    const parmdb::SourceDBSkymodel& source_db,
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

static bool HasSkymodelExtension(const std::string& source_db_name) {
  static const std::string_view kSymodelExtension = ".skymodel";
  static const std::string_view kTxtExtension = ".txt";
  return (source_db_name.size() >= kSymodelExtension.size() &&
          std::equal(kSymodelExtension.rbegin(), kSymodelExtension.rend(),
                     source_db_name.rbegin())) ||
         (source_db_name.size() >= kTxtExtension.size() &&
          std::equal(kTxtExtension.rbegin(), kTxtExtension.rend(),
                     source_db_name.rbegin()));
}

SourceDB::SourceDB(const std::string& source_db_name,
                   const std::vector<std::string>& filter,
                   FilterMode filter_mode) {
  if (std::find(filter.cbegin(), filter.cend(), "") != filter.end()) {
    throw std::runtime_error("Empty source pattern not allowed");
  }

  if (HasSkymodelExtension(source_db_name)) {
    InitialiseUsingSkymodel(source_db_name, filter, filter_mode);
  } else {
    InitialiseUsingSourceDb(source_db_name, filter, filter_mode);
  }
}

std::vector<std::shared_ptr<Patch>> SourceDB::MakePatchList() {
  assert(!HoldsAlternative<std::monostate>() &&
         "The constructor should have properly initialized the source_db_");
  if (HoldsAlternative<parmdb::SourceDBSkymodel>())
    return base::MakePatches(Get<parmdb::SourceDBSkymodel>(), patch_names_);

  return base::makePatches(Get<parmdb::SourceDB>(), patch_names_,
                           patch_names_.size());
}

bool SourceDB::CheckPolarized() {
  assert(!HoldsAlternative<std::monostate>() &&
         "The constructor should have properly initialized the source_db_");

  if (HoldsAlternative<parmdb::SourceDBSkymodel>())
    return base::CheckPolarized(Get<parmdb::SourceDBSkymodel>(), patch_names_);

  return base::checkPolarized(Get<parmdb::SourceDB>(), patch_names_,
                              patch_names_.size());
}

bool SourceDB::CheckAnyOrientationIsAbsolute() {
  assert(!HoldsAlternative<std::monostate>() &&
         "The constructor should have properly initialized the source_db_");

  if (HoldsAlternative<parmdb::SourceDBSkymodel>())
    return base::CheckAnyOrientationIsAbsolute(Get<parmdb::SourceDBSkymodel>(),
                                               patch_names_);

  return base::CheckAnyOrientationIsAbsolute(Get<parmdb::SourceDB>(),
                                             patch_names_, patch_names_.size());
}

void SourceDB::InitialiseUsingSkymodel(const std::string& source_db_name,
                                       const std::vector<std::string>& filter,
                                       FilterMode filter_mode) {
  assert(HasSkymodelExtension(source_db_name) &&
         "Use InitialiseUsingSourceDb instead.");

  source_db_ = parmdb::skymodel_to_source_db::MakeSourceDBSkymodel(
      source_db_name,
      parmdb::skymodel_to_source_db::ReadFormat("", source_db_name));
  switch (filter_mode) {
    case FilterMode::kPattern:
      patch_names_ =
          base::MakePatchList(Get<parmdb::SourceDBSkymodel>(), filter);
      break;
    case FilterMode::kValue:
      patch_names_ = filter;
      break;
  }
}

void SourceDB::InitialiseUsingSourceDb(const std::string& source_db_name,
                                       const std::vector<std::string>& filter,
                                       FilterMode filter_mode) {
  assert(!HasSkymodelExtension(source_db_name) &&
         "Use InitialiseUsingSourceSkymodel instead.");

  source_db_ =
      parmdb::SourceDB(parmdb::ParmDBMeta("", source_db_name), true, false);
  switch (filter_mode) {
    case FilterMode::kPattern:
      patch_names_ = base::makePatchList(Get<parmdb::SourceDB>(), filter);
      break;
    case FilterMode::kValue:
      patch_names_ = filter;
      break;
  }
}

}  // namespace base
}  // namespace dp3
