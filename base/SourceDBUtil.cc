// SourceDBUtil.cc: Helper functions to extract patch and source information
// from a SourceDB.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// $Id$

#include "SourceDBUtil.h"

#include "Exceptions.h"
#include "PointSource.h"
#include "GaussianSource.h"

#include "../parmdb/SourceDB.h"

#include "../common/ProximityClustering.h"

#include <set>
#include <vector>

namespace dp3 {
namespace base {
using parmdb::SourceData;
using parmdb::SourceDB;
using parmdb::SourceInfo;

std::vector<Patch::ConstPtr> makePatches(SourceDB& sourceDB,
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
          source->setSpectralTerms(
              src.getInfo().getSpectralTermsRefFreq(), isLogarithmic,
              src.getSpectralTerms().begin(), src.getSpectralTerms().end());
        }

        // Fetch rotation measure attributes (if applicable).
        if (src.getInfo().getUseRotationMeasure()) {
          source->setRotationMeasure(src.getPolarizedFraction(),
                                     src.getPolarizationAngle(),
                                     src.getRotationMeasure());
        }

        componentsList[i].push_back(source);
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

std::vector<string> makePatchList(SourceDB& sourceDB,
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

bool checkPolarized(SourceDB& sourceDB, const std::vector<string>& patchNames,
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

}  // namespace base
}  // namespace dp3
