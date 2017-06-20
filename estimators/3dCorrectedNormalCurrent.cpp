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

#include "EstimatorHelpers.h"
#include "CorrectedNormalCurrent.h"

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
  typedef Z3i::KSpace                     KSpace;
  typedef EstimatorHelpers< KSpace >      EH;
  typedef EH::Surface                     Surface;
  typedef CorrectedNormalCurrent<Surface> Current;
  // parse command line ----------------------------------------------
  po::options_description general_opt( "Allowed options are" );
  general_opt.add_options()
    ( "help,h", "display this message" )
    ( "m-coef", po::value<double>()->default_value( 1.0 ), "the coefficient k that defines the radius of the ball used in measures, that is r := k h^b" )
    ( "m-pow", po::value<double>()->default_value( 1.0 ), "the coefficient b that defines the radius of the ball used in measures, that is r := k h^b" )
    ;
  EH::optionsImplicitShape   ( general_opt );
  EH::optionsDigitizedShape  ( general_opt );
  EH::optionsNoisyImage      ( general_opt );
  EH::optionsNormalEstimators( general_opt );

  
  po::variables_map vm;
  bool parseOK = EH::args2vm( general_opt, argc, argv, vm );
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
  unsigned int nb = 0;
  auto shape   = EH::makeImplicitShape ( vm );
  auto dshape  = EH::makeDigitizedShape( vm, shape, K );
  auto size    = dshape->getDomain().upperBound() - dshape->getDomain().lowerBound();
  trace.info() << "- Domain size is " << ( size[ 0 ] + 1 ) << " x " << ( size[ 1 ] + 1 )
	       << " x " << ( size[ 2 ] + 1 ) << std::endl;
  auto bimage  = EH::makeNoisyOrNotImage( vm, dshape );
  std::for_each( bimage->cbegin(), bimage->cend(), [&nb] ( bool v ) { nb += v ? 1 : 0; } );
  trace.info() << "- digital shape has " << nb << " voxels." << std::endl;
  auto surface = EH::makeDigitalSurface( K, bimage );
  trace.info() << "- surface component has " << surface->size()<< " surfels." << std::endl;
  trace.endBlock();

  trace.beginBlock( "Compute true normal estimations" );
  auto h        = vm[ "gridstep" ].as<double>();
  auto surfels  = EH::computeDepthFirstSurfelRange( surface );
  auto tnormals = EH::computeTrueNormals( K, shape, surfels );
  trace.endBlock();
  trace.beginBlock( "Compute VCM normal estimations" );
  auto vnormals = EH::computeVCMNormals( vm, surface, surfels );
  auto vstat    = EH::measureAngleDeviation( tnormals, vnormals );
  trace.endBlock();
  trace.beginBlock( "Compute II normal estimations" );
  auto inormals = EH::computeIINormals( vm, K, bimage, surfels );
  EH::orientVectors( tnormals, inormals ); // necessary for II
  auto istat    = EH::measureAngleDeviation( tnormals, inormals );
  trace.endBlock();
  trace.info() << "- VCM h=" << h << " L1=" << vstat.mean() // L1
	       << " L2=" << sqrt( vstat.unbiasedVariance()
				  + vstat.mean()*vstat.mean() ) // L2
	       << " Loo=" << vstat.max() // Loo
	       << std::endl;
  trace.info() << "- II  h=" << h << " L1=" << istat.mean() // L1
	       << " L2=" << sqrt( istat.unbiasedVariance()
				  + istat.mean()*istat.mean() ) // L2
	       << " Loo=" << istat.max() // Loo
	       << std::endl;

  trace.beginBlock( "Computing corrected normal current" );
  Current C( surface, h );
  C.setCorrectedNormals( surfels.begin(), surfels.end(), vnormals.begin() );
  const double mcoef = vm["m-coef"].as<double>();
  const double mpow  = vm["m-pow"].as<double>();
  const double r     = mcoef*pow( h, mpow );
  trace.info() << C << " m-ball-r = " << r << std::endl;
  double            area = 0.0;
  Statistic<double> meanCurv;
  for ( auto v : C )
    {
      area   += C.mu0( v );
      auto m0 = C.mu0Ball( v, r );
      auto m1 = C.mu1Ball( v, r );
      meanCurv.addValue( (m1/(2.0*h*m0)) );
      // std::cout << v
      // 		<< " mu0 = " << m0 << " mu1 = " << m1
      // 		<< " mu1/(r*mu0) = " << (m1/(2.0*h*m0)) << std::endl;
    }
  meanCurv.terminate();
  trace.info() << "- area = " << area << std::endl;
  trace.info() << "- mean curv: avg = " << meanCurv.mean() << std::endl;
  trace.info() << "- mean curv: min = " << meanCurv.min() << std::endl;
  trace.info() << "- mean curv: max = " << meanCurv.max() << std::endl;
  trace.endBlock();

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
