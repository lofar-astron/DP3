//# EstimateNew.cc: Estimate Jones matrices for several directions and stations
//#
//# Copyright (C) 2012
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id$

#include <lofar_config.h>
#include <DPPP/EstimateNew.h>
#include <Common/LofarLogger.h>
#include <Common/StreamUtil.h> ///

#include <scimath/Fitting/LSQFit.h>


namespace LOFAR {
  namespace DPPP {

    EstimateNew::EstimateNew()
    {}

    void EstimateNew::update (size_t maxndir, size_t nBaseline, size_t nStation,
                              size_t nChannel, size_t maxIter,
                              bool propagateSolution)
    {
      itsNrBaselines = nBaseline;
      itsNrStations  = nStation;
      itsNrChannels  = nChannel;
      itsMaxIter     = maxIter;
      itsNrDir       = maxndir;
      itsPropagateSolution = propagateSolution;
      itsSolveStation.resize (nStation);
      itsUnknowns.resize (maxndir * nStation * 4 * 2);
      itsSolution.resize (itsUnknowns.size());
      std::fill (itsSolution.begin(), itsSolution.end(), 0);
      itsDerivIndex.resize (maxndir*4*2*4);
      itsM.resize  (maxndir*4);
      itsdM.resize (maxndir*4*4);
      itsdR.resize (maxndir*8);
      itsdI.resize (maxndir*8);
    }

    // Initialize the solution to the defaultGain for sources/stations not to solve.
    // Set to 0 and diagonal to defaultGaib for solvable ones if no propagation.
    void EstimateNew::initSolution (const vector<vector<int> >& unknownsIndex,
                                    const vector<uint>& srcSet,
				    double defaultGain)
    {
      uint dr=0;
      double* solution = &(itsSolution[0]);
      for (size_t i=0; i<itsNrDir; ++i) {
        if (dr < srcSet.size()  &&  srcSet[dr] == i) {
          // This source is to be solved.
          for (size_t st=0; st<itsNrStations; ++st) {
            if (unknownsIndex[dr][st] >= 0) {
              if (! itsPropagateSolution) {
                std::fill (solution, solution+8, 0);
              }
              // Solvable; set diagonal to 1 if it is 0.
              if (solution[0] == 0) solution[0] = defaultGain;
              if (solution[6] == 0) solution[6] = defaultGain;
            } else {
              // Set non-solvable station-source to 0.
              std::fill (solution, solution+8, 0);
              solution[0] = solution[6] = defaultGain;
            }
            solution += 8;
          }
          dr++;
        } else {
          // Set entire non-solvable source to 0.
          std::fill (solution, solution + 8*itsNrStations, 0.);
          solution += 8*itsNrStations;
        }
      }
    }

    // Clear the solution for unsolvable stations
    // (essentially changing defaultGain to 0).
    void EstimateNew::clearNonSolvable (const vector<vector<int> >& unknownsIndex,
                                        const vector<uint>& srcSet)
    {
      uint dr=0;
      double* solution = &(itsSolution[0]);
      for (size_t i=0; i<itsNrDir; ++i) {
        if (dr < srcSet.size()  &&  srcSet[dr] == i) {
          // This source is to be solved.
          for (size_t st=0; st<itsNrStations; ++st) {
            if (unknownsIndex[dr][st] < 0) {
              // Set non-solvable station-source to 0.
              std::fill (solution, solution+8, 0);
            }
            solution += 8;
          }
          dr++;
        } else {
          // Skip entire source.
          solution += 8*itsNrStations;
        }
      }
    }

    void EstimateNew::fillSolution (const vector<vector<int> >& unknownsIndex,
                                    const vector<uint>& srcSet)
    {
      // Copy the solution for the direction-stations to solve.
      // Note that the solution vector contains all possible sources/stations.
      const double* unknowns = &(itsUnknowns[0]);
      for (size_t dr=0; dr<srcSet.size(); ++dr) {
        size_t inx = srcSet[dr] * itsNrStations * 8;
        double* solution = &(itsSolution[inx]);
        for (size_t st=0; st<itsNrStations; ++st) {
          if (unknownsIndex[dr][st] >= 0) {
            for (size_t k=0; k<8; ++k) {
              *solution++ = *unknowns++;
            }
          } else {
            solution += 8;
          }
        }
      }
    }

    // Form the partial derivative index for a particular baseline.
    // Partial derivatives are only needed if direction-station is solved for.
    // Each visibility provides information about two (complex) unknowns per
    // station per direction. A visibility is measured by a specific
    // interferometer, which is the combination of two stations. Thus, in total
    // each visibility provides information about (no. of directions) x 2 x 2
    // x 2 (scalar) unknowns = (no. of directions) x 8. For each of these
    // unknowns the value of the partial derivative of the model with respect
    // to the unknown has to be computed.
    uint EstimateNew::fillDerivIndex (size_t ndir,
                                      const vector<vector<int> >& unknownsIndex,
                                      const Baseline& baseline)
    {
      // Per direction a baseline has 32 equations with information about
      // 16 unknowns: real and imag part of p00,p01,p10,p11,q00,q01,q10,q11
      // where p and q are the stations forming the baseline.
      // However, only fill if a station has to be solved.
      size_t n = 0;
      for (size_t cr=0; cr<4; ++cr) {
        for (size_t dr=0; dr<ndir; ++dr) {
          if (unknownsIndex[dr][baseline.first] >= 0) {
            size_t idx0 = unknownsIndex[dr][baseline.first] + (cr/2)*4;
            itsDerivIndex[n++] = idx0;
            itsDerivIndex[n++] = idx0 + 1;
            itsDerivIndex[n++] = idx0 + 2;
            itsDerivIndex[n++] = idx0 + 3;
          }
          if (unknownsIndex[dr][baseline.second] >= 0) {
            size_t idx1 = unknownsIndex[dr][baseline.second] + (cr%2)*4;
            itsDerivIndex[n++] = idx1;
            itsDerivIndex[n++] = idx1 + 1;
            itsDerivIndex[n++] = idx1 + 2;
            itsDerivIndex[n++] = idx1 + 3;
          }
        }
      }
      // Return nr of partial derivatives per correlation.
      return n/4;
    }

    // Note that the cursors are passed by value, so a copy is made.
    // In this way no explicit reset of the cursor is needed on a next call.
    bool EstimateNew::estimate (const vector<vector<int> >& unknownsIndex,
                                const vector<uint>& srcSet,
                                const_cursor<Baseline> baselines,
                                vector<const_cursor<fcomplex> > data,
                                vector<const_cursor<dcomplex> > model,
                                const_cursor<bool> flag,
                                const_cursor<float> weight,
                                const_cursor<dcomplex> mix,
				double defaultGain,
                                bool solveBoth,
                                uint verbose)
    {
      initSolution (unknownsIndex, srcSet, defaultGain);
      // Determine if a station has to be solved for any source.
      itsSolveStation = false;
      size_t nUnknowns = 0;
      const size_t nDirection = srcSet.size();
      for (size_t dr=0; dr<nDirection; ++dr) {
        uint drOrig = srcSet[dr];
        const double* solution = &(itsSolution[drOrig*itsNrStations*8]);
        for (size_t st=0; st<itsNrStations; ++st) {
          if (unknownsIndex[dr][st] >= 0) {
            itsSolveStation[st] = true;
            std::copy (solution, solution+8, itsUnknowns.begin()+nUnknowns);
            nUnknowns += 8;
          }
          solution += 8;
        }
      }
      if (verbose > 12) {
        cout<<"unkindex="<<unknownsIndex<<endl;
      }
      // Initialize LSQ solver.
      casa::LSQFit solver(nUnknowns);
      // Iterate until convergence.
      itsNrIter = 0;
      while (!solver.isReady()  &&  itsNrIter < itsMaxIter) {
        if (verbose > 12) {
          cout<<endl<<"iteration " << itsNrIter << endl;
        }
        for (size_t bl=0; bl<itsNrBaselines; ++bl) {
          const size_t p = baselines->first;
          const size_t q = baselines->second;
          // Only compute if no autocorr and if stations need to be solved.
          if (p != q  &&  ((itsSolveStation[p] || itsSolveStation[q])  &&
                           (!solveBoth ||
                            (itsSolveStation[p] && itsSolveStation[q])))) {
            // Create partial derivative index for current baseline.
            size_t nPartial = fillDerivIndex (srcSet.size(), unknownsIndex,
                                              *baselines);
            if (verbose > 13) {
              cout<<"derinx="<<itsDerivIndex<<endl;
            }
            // Generate equations for each channel.
            for (size_t ch=0; ch<itsNrChannels; ++ch) {
              for (size_t dr=0; dr<nDirection; ++dr) {
                uint drOrig = srcSet[dr];
                // Jones matrix for station P.
                const double *Jp = &(itsSolution[(drOrig*itsNrStations + p)*8]);
                const dcomplex Jp_00(Jp[0], Jp[1]);
                const dcomplex Jp_01(Jp[2], Jp[3]);
                const dcomplex Jp_10(Jp[4], Jp[5]);
                const dcomplex Jp_11(Jp[6], Jp[7]);

                // Jones matrix for station Q, conjugated.
                const double *Jq = &(itsSolution[(drOrig*itsNrStations + q)*8]);
                const dcomplex Jq_00(Jq[0], -Jq[1]);
                const dcomplex Jq_01(Jq[2], -Jq[3]);

                const dcomplex Jq_10(Jq[4], -Jq[5]);
                const dcomplex Jq_11(Jq[6], -Jq[7]);

                // Fetch model visibilities for the current direction.
                const dcomplex xx = model[dr][0];
                const dcomplex xy = model[dr][1];
                const dcomplex yx = model[dr][2];
                const dcomplex yy = model[dr][3];

                // Precompute terms involving conj(Jq) and the model
                // visibilities.
                const dcomplex Jq_00xx_01xy = Jq_00 * xx + Jq_01 * xy;
                const dcomplex Jq_00yx_01yy = Jq_00 * yx + Jq_01 * yy;
                const dcomplex Jq_10xx_11xy = Jq_10 * xx + Jq_11 * xy;
                const dcomplex Jq_10yx_11yy = Jq_10 * yx + Jq_11 * yy;

                // Precompute (Jp x conj(Jq)) * vec(data), where 'x'
                // denotes the Kronecker product. This is the model
                // visibility for the current direction, with the
                // current Jones matrix estimates applied. This is
                // stored in M.
                // Also, precompute the partial derivatives of M with
                // respect to all 16 parameters (i.e. 2 Jones matrices
                // Jp and Jq, 4 complex scalars per Jones matrix, 2 real
                // scalars per complex scalar, 2 * 4 * 2 = 16). These
                // partial derivatives are stored in dM.
                // Note that conj(Jq) is used and that q01 and q10 are swapped.

                itsM[dr * 4] = Jp_00 * Jq_00xx_01xy + Jp_01 * Jq_00yx_01yy;
                // Derivatives of M00 wrt p00, p01, q00, q01
                itsdM[dr * 16] = Jq_00xx_01xy;
                itsdM[dr * 16 + 1] = Jq_00yx_01yy;
                itsdM[dr * 16 + 2] = Jp_00 * xx + Jp_01 * yx;
                itsdM[dr * 16 + 3] = Jp_00 * xy + Jp_01 * yy;

                itsM[dr * 4 + 1] = Jp_00 * Jq_10xx_11xy + Jp_01 * Jq_10yx_11yy;
                // Derivatives of M01 wrt p00, p01, q10, q11
                itsdM[dr * 16 + 4] = Jq_10xx_11xy;
                itsdM[dr * 16 + 5] = Jq_10yx_11yy;
                itsdM[dr * 16 + 6] = itsdM[dr * 16 + 2];
                itsdM[dr * 16 + 7] = itsdM[dr * 16 + 3];

                itsM[dr * 4 + 2] = Jp_10 * Jq_00xx_01xy + Jp_11 * Jq_00yx_01yy;
                // Derivatives of M10 wrt p10, p11, q00, q01
                itsdM[dr * 16 + 8] = itsdM[dr * 16];
                itsdM[dr * 16 + 9] = itsdM[dr * 16 + 1];
                itsdM[dr * 16 + 10] = Jp_10 * xx + Jp_11 * yx;
                itsdM[dr * 16 + 11] = Jp_10 * xy + Jp_11 * yy;

                itsM[dr * 4 + 3] = Jp_10 * Jq_10xx_11xy + Jp_11 * Jq_10yx_11yy;
                // Derivatives of M11 wrt p10, p11, q10, q11
                itsdM[dr * 16 + 12] = itsdM[dr * 16 + 4];
                itsdM[dr * 16 + 13] = itsdM[dr * 16 + 5];
                itsdM[dr * 16 + 14] = itsdM[dr * 16 + 10];
                itsdM[dr * 16 + 15] = itsdM[dr * 16 + 11];
              }
              if (verbose > 14) {
                cout<<"M="<<itsM<<endl;
                cout<<"dM="<<itsdM<<endl;
              }

              // Now compute the equations (per pol) for D*M=A where
              //  D is the NxN demixing weight matrix
              //  M is the model visibilities vector for the N directions
              //  A is the shifted observed visibilities vector for N directions
              // Note that each element in the vectors is a 2x2 matrix
              // (xx,xy,yx,yy) of complex values.
              // A complex multiplication of (a,b) and (c,d) gives (ac-bd,ad+bc)
              // Thus real partial derivatives wrt a,b,c,d are c,-d,a,-b.
              // Imaginary partial derivatives wrt a,b,c,d are d,c,b,a
              for (size_t cr=0; cr<4; ++cr) {
                // Only use visibility if not flagged.
                if (!flag[cr]) {
                  // For each direction a set of equations is generated.
                  for (size_t tg=0; tg<nDirection; ++tg) {
                    dcomplex visibility(0.0, 0.0);
                    // Each direction is dependent on all directions.
                    size_t off = 0;
                    for (size_t dr=0; dr<nDirection; ++dr) {
                      bool do1 = unknownsIndex[dr][p] >= 0;
                      bool do2 = unknownsIndex[dr][q] >= 0;
                      // Only generate equations if a station has to be solved
                      // for this direction.
                      if ((do1 && do2)  ||  (!solveBoth && (do1 || do2))) {
                        // Look-up mixing weight.
                        const dcomplex mix_weight = *mix;
                        // Sum weighted model visibilities.
                        visibility += mix_weight * itsM[dr * 4 + cr];

                        // Compute weighted partial derivatives.
                        if (do1) {
                          dcomplex der(mix_weight * itsdM[dr * 16 + cr * 4]);
                          itsdR[off]     = real(der);
                          itsdI[off]     = imag(der);
                          itsdR[off + 1] = -imag(der);
                          itsdI[off + 1] = real(der);
                          off += 2;
                        }
                        if (do2) {
                          dcomplex der(mix_weight * itsdM[dr * 16 + cr * 4 + 1]);
                          itsdR[off]     = real(der);
                          itsdI[off]     = imag(der);
                          itsdR[off + 1] = -imag(der);
                          itsdI[off + 1] = real(der);
                          off += 2;
                        }
                        if (do1) {
                          dcomplex der(mix_weight * itsdM[dr * 16 + cr * 4 + 2]);
                          itsdR[off]     = real(der);
                          itsdI[off]     = imag(der);
                          itsdR[off + 1] = imag(der);  // conjugate
                          itsdI[off + 1] = -real(der);
                          off += 2;
                        }
                        if (do2) {
                          dcomplex der(mix_weight * itsdM[dr * 16 + cr * 4 + 3]);
                          itsdR[off]     = real(der);
                          itsdI[off]     = imag(der);
                          itsdR[off + 1] = imag(der);
                          itsdI[off + 1] = -real(der);
                          off += 2;
                        }
                      }
                      // Move to next source direction.
                      mix.forward(1);
                    } // Source directions.

                    // Compute the residual.
                    dcomplex residual(data[tg][cr]);
                    residual -= visibility;

                    // Update the normal equations.
                    solver.makeNorm(nPartial,
                                    &(itsDerivIndex[cr * nPartial]), &(itsdR[0]),
                                    static_cast<double>(weight[cr]),
                                    real(residual));
                    solver.makeNorm(nPartial,
                                    &(itsDerivIndex[cr * nPartial]), &(itsdI[0]),
                                    static_cast<double>(weight[cr]),
                                    imag(residual));
                    if (verbose > 14) {
                      cout<<"makeres "<<real(residual)<<' '<<weight[cr]
                          <<' '<<nPartial;
                      for (uint i=0; i<nPartial; ++i) {
                        cout << ' '<<itsDerivIndex[cr*nPartial+i]<<' '<<itsdR[i];
                      }
                      cout<<endl;
                    }

                    // Move to next target direction.
                    mix.backward(1, nDirection);
                    mix.forward(0);
                  } // Target directions.

                  // Reset cursor to the start of the correlation.
                  mix.backward(0, nDirection);
                }

                // Move to the next correlation.
                mix.forward(2);
              } // Correlations.

              // Move to the next channel.
              mix.backward(2, 4);
              mix.forward(3);

              for (size_t dr=0; dr<nDirection; ++dr) {
                model[dr].forward(1);
                data[dr].forward(1);
              }
              flag.forward(1);
              weight.forward(1);
            } // Channels.

            // Reset cursors to the start of the baseline.
            for (size_t dr=0; dr<nDirection; ++dr) {
              model[dr].backward(1, itsNrChannels);
              data[dr].backward(1, itsNrChannels);
            }
            flag.backward(1, itsNrChannels);
            weight.backward(1, itsNrChannels);
            mix.backward(3, itsNrChannels);
          }

          // Move cursors to the next baseline.
          for (size_t dr=0; dr<nDirection; ++dr) {
            model[dr].forward(2);
            data[dr].forward(2);
          }
          flag.forward(2);
          weight.forward(2);
          mix.forward(4);
          ++baselines;
        } // Baselines.

        // Reset all cursors for the next iteration.
        for (size_t dr=0; dr<nDirection; ++dr) {
          model[dr].backward(2, itsNrBaselines);
          data[dr].backward(2, itsNrBaselines);
        }
        flag.backward(2, itsNrBaselines);
        weight.backward(2, itsNrBaselines);
        mix.backward(4, itsNrBaselines);
        baselines -= itsNrBaselines;

        // Perform LSQ iteration.
        casa::uInt rank;
        bool status = solver.solveLoop(rank, &(itsUnknowns[0]), true);
        ASSERT(status);
        // Copy the unknowns to the full solution.
        fillSolution (unknownsIndex, srcSet);
        if (verbose > 13) {
          cout<<"unknowns="<<nUnknowns<<' '<<itsUnknowns<<endl;
          cout<<"solution="<<itsSolution<<endl;
        }
        // Update iteration count.
        itsNrIter++;
      }
      bool converged = (solver.isReady() == casa::LSQFit::SOLINCREMENT  ||
                        solver.isReady() == casa::LSQFit::DERIVLEVEL);
      ///      clearNonSolvable (unknownsIndex, srcSet);
      return converged;
    }

  } //# namespace DPPP
} //# namespace LOFAR
