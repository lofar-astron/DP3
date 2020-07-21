// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.

#ifndef IDG_CONFIGURATION_H
#define IDG_CONFIGURATION_H

#include <idg-api.h>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include <fstream>

class IdgConfiguration {
public:
  static void Read(idg::api::Type& proxyType, int& bufferSize, idg::api::options_type& options)
  {
    namespace po = boost::program_options; 
    po::options_description desc("Options"); 
    desc.add_options() 
    ("proxy", "idg proxy")
    ("max_nr_w_layers", po::value<int>(), "")
    ("buffersize", po::value<int>(), ""); 

    po::variables_map vm;
    std::ifstream ifs("idg.conf");
    if (ifs.fail())
    {
      // Logger::Debug << "could not open config file\n";
    }
    else
    {
      try 
      { 
        po::store(po::parse_config_file(ifs, desc), vm);
      }
      catch(po::error&)
      { }
    }
    if(vm.count("proxy")) 
    {
      std::string proxy(vm["proxy"].as<std::string>());
      boost::to_lower(proxy);
      if (proxy == "cpu-optimized") proxyType = idg::api::Type::CPU_OPTIMIZED;
      if (proxy == "cpu-reference") proxyType = idg::api::Type::CPU_REFERENCE;
      if (proxy == "cuda-generic") proxyType = idg::api::Type::CUDA_GENERIC;
      if (proxy == "hybrid-cuda-cpu-optimized") proxyType = idg::api::Type::HYBRID_CUDA_CPU_OPTIMIZED;
    }
    if (vm.count("buffersize")) 
    {
      bufferSize = vm["buffersize"].as<int>();
    }
    if (vm.count("max_nr_w_layers")) 
    {
      options["max_nr_w_layers"] = vm["max_nr_w_layers"].as<int>();
    }
  }
};

#endif
