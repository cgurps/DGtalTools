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
 * @file 3dFastCorrectedNormalCurrent.cpp
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

#include "LaplaceOperatorHelpers.h"
#include "FastCorrectedNormalCurrent.h"

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
 @page Doc3DCorrectedNormalCurrent 3dFastCorrectedNormalCurrent

 @brief  Computes and visualizes the 3d corrected normal current of digital surfaces.

 @b Usage:  3dFastCorrectedNormalCurrent -i file.vol

 @b Allowed @b options @b are:

 @code
  -h [ --help ]                         display this message
  -i [ --input ] arg                    .vol file
 @endcode

 @b Example:
 @see
 @ref 3dFastCorrectedNormalCurrent.cpp
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

//namespace po = DGtal::po;

int main( int argc, char** argv )
{
#ifdef WITH_VISU3D_QGLVIEWER
  QApplication application( argc, argv );
#endif

  typedef Z3i::Space                      Space;
  typedef Z3i::KSpace                     KSpace;
  typedef EstimatorHelpers< KSpace >      EH;
  typedef EH::Surface                     Surface;
  typedef EH::RealVector                  RealVector;
  typedef EH::BinaryImage                 BinaryImage;
  typedef EH::ImplicitShape               ImplicitShape;
  typedef EH::SurfaceContainer            SurfaceContainer;
  typedef FastCorrectedNormalCurrent<SurfaceContainer> Current;
  typedef Current::Vertex                 Vertex;
  typedef double                          Scalar;

  // parse command line ----------------------------------------------
  po::options_description general_opt( "Allowed options are" );
  general_opt.add_options()
    ( "help,h", "display this message" )
    ( "m-coef", po::value<Scalar>()->default_value( 3.0 ), "the coefficient k that defines the radius of the ball used in measures, that is r := k h^b" )
    ( "m-pow", po::value<Scalar>()->default_value( 0.33, "0.33" ), "the coefficient b that defines the radius of the ball used in measures, that is r := k h^b" );

  EH::optionsImplicitShape   ( general_opt );
  EH::optionsDigitizedShape  ( general_opt );
  EH::optionsVolFile         ( general_opt );
  EH::optionsNoisyImage      ( general_opt );
  EH::optionsNormalEstimators( general_opt );

  po::variables_map vm;
  bool parseOK = EH::args2vm( general_opt, argc, argv, vm );
  bool neededArgsGiven=true;

  if ( vm.count( "polynomial-list" ) )
    {
      trace.info() << "List of predefined polynomials:" << std::endl;
      auto L = EH::getPolynomialList();
      for ( auto p : L ) {
	trace.info() << "  " << p.first << " -> " << p.second << std::endl;
      }
    }

  if (parseOK && ( ! vm.count("polynomial") ) && ( ! vm.count( "input" ) ) ) {
    missingParam("--polynomial or --input");
    neededArgsGiven=false;
  }

  if(!neededArgsGiven || !parseOK || vm.count("help") || argc <= 1 )
    {
      trace.info()<< "Builds the 3d corrected normal currents" <<std::endl
                  << general_opt << "\n"
                  << "Basic usage: "<<std::endl
		  << "\t 3dFastCorrectedNormalCurrent -e II -t 3 -r 3 -Q H --tics 1 --minValue -0.3 --maxValue 0.3 --polynomial-list -p goursat -g 0.5 -V Measure --m-coef 3 -N 0.2" << std::endl
                  << "\t 3dFastCorrectedNormalCurrent -p \"3x^2+5y^2+7z^2-1\" "<<std::endl
                  << "\t 3dFastCorrectedNormalCurrent -i \"file.vol\" "<<std::endl
		  << "\t 3dFastCorrectedNormalCurrent -e II -t 3 -r 3 -Q H --tics 1 --minValue -0.3 --maxValue 0.3 -V Measure --m-coef 3  -i ~/Images/3d/vol/fandisk-128.vol -m 0 -M 255 -N 0.4 -g 0.25" << std::endl
                  << "\t 3dFastCorrectedNormalCurrent --polynomial-list  // to get the list of predefined polynomials. "<<std::endl
                  << std::endl;
      return 0;
    }

  //while( vm["gridstep"].as<Scalar>() > 0.01 )
  {
    KSpace                        K;
    unsigned int                 nb = 0;
    CountedPtr<BinaryImage>  bimage( nullptr );
    CountedPtr<ImplicitShape> shape( nullptr );

    trace.beginBlock( "Make Shape" );
    if ( vm.count( "polynomial" ) )
      {
        shape        = EH::makeImplicitShape ( vm );
        auto dshape  = EH::makeImplicitDigitalShapeFromImplicitShape( vm, shape, K );
        auto size    = dshape->getDomain().upperBound() - dshape->getDomain().lowerBound();
        trace.info() << "- Domain size is " << ( size[ 0 ] + 1 )
	  	   << " x " << ( size[ 1 ] + 1 )
	  	   << " x " << ( size[ 2 ] + 1 ) << std::endl;
        bimage       = EH::makeNoisyOrNotBinaryImageFromImplicitDigitalShape( vm, dshape );
      }
    else if ( vm.count( "input" ) )
      {
        bimage       = EH::makeNoisyOrNotBinaryImageFromVolFile( vm, K );
      }
    std::for_each( bimage->cbegin(), bimage->cend(),
	  	 [&nb] ( bool v ) { nb += v ? 1 : 0; } );
    trace.info() << "- digital shape has " << nb << " voxels." << std::endl;
    auto surface = EH::makeDigitalSurfaceFromBinaryImage( K, bimage );
    if ( surface == 0 ) {
        trace.info() << "- surface is empty (either empty or full volume). "
	  	   << std::endl;
        return 1;
    }
    trace.info() << "- surface component has " << surface->size()<< " surfels." << std::endl;
    trace.endBlock();

    trace.beginBlock( "Compute surfels" );
    auto surfels  = EH::computeDepthFirstSurfelRange( surface );
    trace.endBlock();

    const Scalar h  = vm["gridstep"].as<Scalar>();
    std::vector<RealVector> normals;
    if ( vm.count( "polynomial" ) )
    {
      trace.beginBlock( "Compute estimated normals" );
      normals = EH::computeIINormals   ( vm, K, bimage, surfels );
      trace.endBlock();
    }

    typedef DGtal::LaplaceOperatorHelpers<Scalar, KSpace> LaplaceOperatorHelpers;
    typedef LaplaceOperatorHelpers::DenseVector DenseVector;
    LaplaceOperatorHelpers laplaceOperatorHelpers;
    laplaceOperatorHelpers.generateIndexedDigitalSurface( surface );
    laplaceOperatorHelpers.generateTriangulatedSurface( surface );
    laplaceOperatorHelpers.setNormals( surfels.begin(), surfels.end(), normals.begin() );
    laplaceOperatorHelpers.setVariableMap( vm );

    trace.info() << "There are " << std::distance( surfels.begin(), surfels.end() ) << " surfels." << std::endl;
    trace.info() << "Indexed Digital Surface has " << laplaceOperatorHelpers.myIndexedDigitalSurface.nbVertices() << " faces." << std::endl;
    trace.info() << "Triangulation has " << laplaceOperatorHelpers.myTriangulatedSurface.nbVertices() << " vertices." << std::endl;

    // DIGITAL PART
    {
      LaplaceOperatorHelpers::SparseMatrix combiLaplace             = laplaceOperatorHelpers.computeCombinatorialLaplace();

      trace.beginBlock("Laplace Heat");
      DenseVector heatCurvature  = laplaceOperatorHelpers.computeHeatCurvatureFromFunctor();
      trace.endBlock();
      trace.beginBlock("Combi");
      DenseVector combiCurvature = 0.5 * ( combiLaplace * laplaceOperatorHelpers.computeDualDigitalEmbedding() ).rowwise().lpNorm<2>() ;
      trace.endBlock();

      DenseVector realCurvature = Eigen::Map<DenseVector>( EH::computeMeanCurvatures( K, shape, h, surfels ).data(), surfels.size() ).cwiseAbs();
      trace.beginBlock("II");
      DenseVector IICurvature   = Eigen::Map<DenseVector>( EH::computeIIMeanCurvatures( vm, K, bimage, surfels ).data(), surfels.size() ).cwiseAbs();
      trace.endBlock();
      trace.beginBlock("JET");
      DenseVector jetCurvature  = Eigen::Map<DenseVector>( EH::computeJetMeanCurvatures( vm, surface, surfels ).data(), surfels.size() ).cwiseAbs();
      trace.endBlock();

      DenseVector realReorD = laplaceOperatorHelpers.reorder( realCurvature, surfels.begin(), false ).cwiseAbs();

      laplaceOperatorHelpers.exportErrors( realReorD, heatCurvature, "heat_errors.dat" );
      laplaceOperatorHelpers.exportErrors( realReorD, combiCurvature, "combi_errors.dat" );
      laplaceOperatorHelpers.exportErrors( realCurvature, IICurvature, "II_errors.dat" );
      laplaceOperatorHelpers.exportErrors( realCurvature, jetCurvature, "jet_errors.dat" );

      laplaceOperatorHelpers.exportDualQuantityToOBJ( realCurvature, "real_curvature.obj", surfels.begin() );
      laplaceOperatorHelpers.exportDualQuantityToOBJ( IICurvature  , "II_curvature.obj"  , surfels.begin() );
      laplaceOperatorHelpers.exportDualQuantityToOBJ( jetCurvature , "jet_curvature.obj" , surfels.begin() );
      laplaceOperatorHelpers.exportDualQuantityToOBJ( heatCurvature, "heat_curvature.obj" );
      laplaceOperatorHelpers.exportDualQuantityToOBJ( combiCurvature, "combi_curvature.obj" );
    }

    // TRIANGULATION PART
    {
      LaplaceOperatorHelpers::SparseMatrix cotanLaplace  = laplaceOperatorHelpers.computeCotanLaplace();
      LaplaceOperatorHelpers::SparseMatrix rLocalLaplace = laplaceOperatorHelpers.computeRLocalLaplace();

      trace.beginBlock("Laplace Heat Triangulation");
      DenseVector heatCurvature  = laplaceOperatorHelpers.computeHeatCurvatureFromTriangulationFunctor();
      trace.endBlock();
      trace.beginBlock("Cotan");
      DenseVector cotanCurvature  = 0.5 * ( cotanLaplace * laplaceOperatorHelpers.computeDualTriangulationEmbedding() ).rowwise().lpNorm<2>() ;
      trace.endBlock();
      trace.beginBlock("rLocal");
      DenseVector rLocalCurvature = 0.5 * ( rLocalLaplace * laplaceOperatorHelpers.computeDualTriangulationEmbedding() ).rowwise().lpNorm<2>() ;
      trace.endBlock();
      DenseVector realCurvature = laplaceOperatorHelpers.computeRealCurvatureTriangulation( shape ).cwiseAbs();

      laplaceOperatorHelpers.exportErrors( realCurvature, cotanCurvature, "cotan_errors.dat" );
      laplaceOperatorHelpers.exportErrors( realCurvature, rLocalCurvature, "rLocal_errors.dat" );
      laplaceOperatorHelpers.exportErrors( realCurvature, heatCurvature, "heat_triangulation_errors.dat" );

      laplaceOperatorHelpers.exportTriangulationQuantityToOBJ( cotanCurvature, "cotan_curvature.obj" );
      laplaceOperatorHelpers.exportTriangulationQuantityToOBJ( rLocalCurvature, "rLocal_curvature.obj" );
      laplaceOperatorHelpers.exportTriangulationQuantityToOBJ( heatCurvature, "heat_triangulation_curvature.obj" );
    }


    vm.at("gridstep").value() = vm["gridstep"].as<Scalar>() / 1.3;
  }

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
