/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/
/**
 * @file 3dCorrectedNormalCurrent.cpp
 * @ingroup surfaceTools
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2017/06/01
 *
 *
 * This file is part of the DGtalTools.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

// Shape constructors
#include "EstimatorHelpers.h"

// // Integral Invariant includes
// #include "DGtal/geometry/surfaces/estimation/IIGeometricFunctors.h"
// #include "DGtal/geometry/surfaces/estimation/IntegralInvariantVolumeEstimator.h"
// #include "DGtal/geometry/surfaces/estimation/IntegralInvariantCovarianceEstimator.h"

// #ifdef WITH_VISU3D_QGLVIEWER
// #include "DGtal/io/viewers/Viewer3D.h"
// #endif

using namespace DGtal;
using namespace functors;

/**
 @page Doc3DCorrectedNormalCurrent 3dCorrectedNormalCurrent

 @brief  Computes and visualizes the 3d corrected normal current of digital surfaces.

 @b Usage:  3dCorrectedNormalCurrent -i file.vol

 @b Allowed @b options @b are:

 @code
  -h [ --help ]                         display this message
  -i [ --input ] arg                    .vol file
 @endcode

 @b Example:
 @see
 @ref 3dCorrectedNormalCurrent.cpp
 @ref Doc3DCorrectedNormalCurrent
 */

///////////////////////////////////////////////////////////////////////////////

/**
 * Missing parameter error message.
 *
 * @param param
 */
void missingParam( std::string param )
{
  trace.error() << " Parameter: " << param << " is required.";
  trace.info() << std::endl;
}

namespace po = boost::program_options;

int main( int argc, char** argv )
{
  typedef Z3i::KSpace                KSpace;
  typedef EstimatorHelpers< KSpace > EH;
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are");
  general_opt.add_options() ("help,h", "display this message");
  EH::optionsImplicitShape ( general_opt );
  EH::optionsDigitizedShape( general_opt );

  po::variables_map vm;
  bool parseOK = EH::args2vm( general_opt, argc, argv, vm );
  // try
  //   {
  //     po::store( po::parse_command_line( argc, argv, general_opt ), vm );
  //   }
  // catch( const std::exception & ex )
  //   {
  //     parseOK = false;
  //     trace.error() << " Error checking program options: " << ex.what() << std::endl;
  //   }
  bool neededArgsGiven=true;

  if (parseOK && !(vm.count("polynomial"))){
    missingParam("--polynomial");
    neededArgsGiven=false;
  }
  
  if(!neededArgsGiven || !parseOK || vm.count("help") || argc <= 1 )
    {
      trace.info()<< "Builds the 3d corrected normal currents" <<std::endl
                  << general_opt << "\n"
                  << "Basic usage: "<<std::endl
                  << "\t 3dCorrectedNormalCurrent -p \"3x^2+5y^2+7z^2-1\" "<<std::endl
                  << std::endl;
      return 0;
    }

  // Digital space.
  KSpace K;
  
  trace.beginBlock( "Make implicit shape and digitized shape" );
  auto shape  = EH::makeImplicitShape ( vm );
  auto dshape = EH::makeDigitizedShape( vm, shape, K );
  trace.endBlock();


  return 0;
}

///////////////////////////////////////////////////////////////////////////////
