/***************************************************************************
 *   Copyright (C) 2007 by Adriaan Renting, ASTRON                         *
 *   renting@astron.nl                                                     *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <lofar_config.h>
#include <libgen.h>
#include <PLC/ACCmain.h>
#include <casa/Exceptions.h>
#include <CS1_SPWCombine/CombinerProcessControl.h>

int main(int argc, char *argv[])
{
  try
  {
    INIT_LOGGER(basename(argv[0]));
    LOFAR::CS1::CombinerProcessControl myProcess;
    return LOFAR::ACC::PLC::ACCmain(argc, argv, &myProcess);
  } //try
  catch(casa::AipsError& err)
  {
    std::cerr << "Aips++ error detected: " << err.getMesg() << std::endl;
    return -2;
  }
  catch(...)
  {
    std::cerr << "** PROBLEM **: Unhandled exception caught." << std::endl;
    return -3;
  }
}
