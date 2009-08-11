#  $Id$
#
#  Copyright (C) 2008-2009
#  ASTRON (Netherlands Foundation for Research in Astronomy)
#  P.O.Box 2, 7990 AA Dwingeloo, The Netherlands, seg@astron.nl
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

include(LofarPackageVersion)

lofar_add_library(dppp
  Package__Version.cc
  FlaggerStatistics.cc
  MsInfo.cc
  MsFile.cc
  ComplexMedianFlagger.cc
  FrequencyFlagger.cc
  ComplexMedianFlagger2.cc
  MADFlagger.cc
  BandpassCorrector.cc
  DataSquasher.cc
  DataBuffer.cc
  TimeBuffer.cc
  RunDetails.cc
  Pipeline.cc)

lofar_add_bin_program(CS1_IDPPP IDPPP.cc PipelineProcessControl.cc)
lofar_add_bin_program(versiondppp versiondppp.cc)
