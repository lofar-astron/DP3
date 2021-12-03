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

#include "../common/ProximityClustering.h"

#include <boost/utility/string_view.hpp>

#include <cassert>
#include <set>
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

std::vector<Patch::ConstPtr> makePatches(parmdb::SourceDB& sourceDB,
                                         const std::vector<string>& patchNames,
                                         unsigned int nModel) {
  // Create a component list for each patch name.
  std::vector<std::vector<ModelComponent::Ptr>> componentsList(nModel);

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

  std::vector<Patch::ConstPtr> patchList;
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

std::vector<Patch::ConstPtr> MakePatches(
    const parmdb::SourceDBSkymodel& source_db,
    const std::vector<std::string>& patch_names) {
  // Create a component list for each patch name.
  std::vector<std::vector<ModelComponent::Ptr>> componentsList(
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

  std::vector<Patch::ConstPtr> patchList;
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

std::vector<std::pair<ModelComponent::ConstPtr, Patch::ConstPtr>>
makeSourceList(const std::vector<Patch::ConstPtr>& patchList) {
  std::vector<Patch::ConstPtr>::const_iterator pIter = patchList.begin();
  std::vector<Patch::ConstPtr>::const_iterator pEnd = patchList.end();

  unsigned int nSources = 0;
  for (; pIter != pEnd; ++pIter) {
    nSources += (*pIter)->nComponents();
  }

  std::vector<std::pair<ModelComponent::ConstPtr, Patch::ConstPtr>> sourceList;
  sourceList.reserve(nSources);

  pIter = patchList.begin();

  for (; pIter != pEnd; ++pIter) {
    Patch::const_iterator sIter = (*pIter)->begin();
    Patch::const_iterator sEnd = (*pIter)->end();
    for (; sIter != sEnd; ++sIter) {
      sourceList.push_back(make_pair(*sIter, *pIter));
    }
  }

  return sourceList;
}

std::vector<Patch::ConstPtr> makeOnePatchPerComponent(
    const std::vector<Patch::ConstPtr>& patchList) {
  size_t numComponents = 0;

  for (const auto& patch : patchList) {
    numComponents += patch->nComponents();
  }

  std::vector<Patch::ConstPtr> largePatchList;
  largePatchList.reserve(numComponents);

  for (const auto& patch : patchList) {
    size_t compNum = 0;
    for (auto compIt = patch->begin(); compIt != patch->end(); ++compIt) {
      Patch::Ptr ppatch(new Patch(patch->name() + "_" + std::to_string(compNum),
                                  compIt, compIt + 1));
      ppatch->setDirection((*compIt)->direction());
      largePatchList.push_back(std::move(ppatch));
      compNum++;
    }
  }

  return largePatchList;
}

std::vector<Patch::ConstPtr> clusterProximateSources(
    const std::vector<Patch::ConstPtr>& patch_list, double proximity_limit) {
  // Create a list of all source positions and a lookup table to go from index
  // back to patch & comp
  std::vector<std::pair<double, double>> sources;
  std::vector<std::pair<size_t, size_t>> lookup_table;
  size_t patchIndex = 0;
  for (const auto& patch : patch_list) {
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
  std::vector<Patch::ConstPtr> clusteredPatchList;
  for (const auto& groups : clusters) {
    std::vector<ModelComponent::ConstPtr> componentsInPatch;
    double alpha = 0.0;
    double delta = 0.0;
    for (const auto& component : groups) {
      auto lookupIndices = lookup_table[component];
      const auto& patch = patch_list[lookupIndices.first];
      const ModelComponent::ConstPtr& comp =
          *(patch->begin() + lookupIndices.second);
      alpha += comp->direction().ra;
      delta += comp->direction().dec;
      componentsInPatch.emplace_back(comp);
    }
    const auto& firstPatch = patch_list[lookup_table[groups.front()].first];
    Patch::Ptr ppatch(new Patch(
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

bool checkPolarized(parmdb::SourceDB& sourceDB,
                    const std::vector<string>& patchNames,
                    unsigned int nModel) {
  bool polarized = false;

  // Loop over all sources.
  sourceDB.lock();
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
  sourceDB.unlock();
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

static bool HasSkymodelExtension(const std::string& source_db_name) {
  static const boost::string_view kSymodelExtension = ".skymodel";
  static const boost::string_view kTxtExtension = ".txt";
  return (source_db_name.size() >= kSymodelExtension.size() &&
          std::equal(kSymodelExtension.rbegin(), kSymodelExtension.rend(),
                     source_db_name.rbegin())) ||
         (source_db_name.size() >= kTxtExtension.size() &&
          std::equal(kTxtExtension.rbegin(), kTxtExtension.rend(),
                     source_db_name.rbegin()));
}

SourceDB::SourceDB(const std::string& source_db_name,
                   const std::vector<std::string>& source_patterns) {
  if (std::find(source_patterns.cbegin(), source_patterns.cend(), "") !=
      source_patterns.end()) {
    throw std::runtime_error("Empty source pattern not allowed");
  }

  if (HasSkymodelExtension(source_db_name)) {
    source_db_ = parmdb::skymodel_to_source_db::MakeSourceDBSkymodel(
        source_db_name,
        parmdb::skymodel_to_source_db::ReadFormat("", source_db_name));
    patch_names_ =
        base::MakePatchList(Get<parmdb::SourceDBSkymodel>(), source_patterns);
  } else {
    source_db_ =
        parmdb::SourceDB(parmdb::ParmDBMeta("", source_db_name), true, false);
    patch_names_ =
        base::makePatchList(Get<parmdb::SourceDB>(), source_patterns);
  }
}

std::vector<Patch::ConstPtr> SourceDB::MakePatchList() {
  assert(!HoldsAlternative<monostate>() &&
         "The constructor should have properly initialized the source_db_");
  if (HoldsAlternative<parmdb::SourceDBSkymodel>())
    return base::MakePatches(Get<parmdb::SourceDBSkymodel>(), patch_names_);

  return base::makePatches(Get<parmdb::SourceDB>(), patch_names_,
                           patch_names_.size());
}

bool SourceDB::CheckPolarized() {
  assert(!HoldsAlternative<monostate>() &&
         "The constructor should have properly initialized the source_db_");

  if (HoldsAlternative<parmdb::SourceDBSkymodel>())
    return base::CheckPolarized(Get<parmdb::SourceDBSkymodel>(), patch_names_);

  return base::checkPolarized(Get<parmdb::SourceDB>(), patch_names_,
                              patch_names_.size());
}

}  // namespace base
}  // namespace dp3
