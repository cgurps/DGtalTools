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

  // Digital space.
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
  LaplaceOperatorHelpers laplaceOperatorHelpers;
  laplaceOperatorHelpers.generateIndexedDigitalSurface( surface );
  laplaceOperatorHelpers.generateTriangulatedSurface( surface );
  laplaceOperatorHelpers.setNormals( surfels.begin(), surfels.end(), normals.begin() );
  LaplaceOperatorHelpers::SparseMatrix heatLaplace              = laplaceOperatorHelpers.computeHeatLaplace( vm );
  LaplaceOperatorHelpers::SparseMatrix heatLaplaceTriangulation = laplaceOperatorHelpers.computeHeatLaplaceTriangulation( vm );
  LaplaceOperatorHelpers::SparseMatrix combiLaplace             = laplaceOperatorHelpers.computeCombinatorialLaplace();
  LaplaceOperatorHelpers::SparseMatrix cotanLaplace             = laplaceOperatorHelpers.computeCotanLaplace( vm );
  LaplaceOperatorHelpers::SparseMatrix rLocalLaplace            = laplaceOperatorHelpers.computeRLocalLaplace( vm );

  LaplaceOperatorHelpers::DenseVector heatCurvature              = 0.5 * ( heatLaplace
      * laplaceOperatorHelpers.computeDualDigitalEmbedding( vm ) ).rowwise().lpNorm<2>() ;
  LaplaceOperatorHelpers::DenseVector heatCurvatureTriangulation = 0.5 * ( heatLaplaceTriangulation
      * laplaceOperatorHelpers.computeDualTriangulationEmbedding( vm ) ).rowwise().lpNorm<2>() ;
  LaplaceOperatorHelpers::DenseVector combiCurvature             = 0.5 * ( combiLaplace
      * laplaceOperatorHelpers.computeDualDigitalEmbedding( vm ) ).rowwise().lpNorm<2>() ;
  LaplaceOperatorHelpers::DenseVector cotanCurvature             = 0.5 * ( cotanLaplace
      * laplaceOperatorHelpers.computeDualTriangulationEmbedding( vm ) ).rowwise().lpNorm<2>() ;
  LaplaceOperatorHelpers::DenseVector rLocalCurvature            = 0.5 * ( rLocalLaplace
      * laplaceOperatorHelpers.computeDualTriangulationEmbedding( vm ) ).rowwise().lpNorm<2>() ;
  LaplaceOperatorHelpers::DenseVector realCurvature              = Eigen::Map<LaplaceOperatorHelpers::DenseVector>( EH::computeMeanCurvatures( K, shape, h, surfels ).data(), surfels.size() );

  laplaceOperatorHelpers.exportDualQuantityToOBJ( realCurvature,
      "real_curvature.obj", surfels.begin() );
  laplaceOperatorHelpers.exportDualQuantityToOBJ( heatCurvature,
      "heat_curvature.obj"   );
  laplaceOperatorHelpers.exportDualQuantityToOBJ( heatCurvatureTriangulation,
      "heat_curvature_triangulation.obj", true );
  laplaceOperatorHelpers.exportDualQuantityToOBJ( combiCurvature,
      "combi_curvature.obj"  );
  laplaceOperatorHelpers.exportDualQuantityToOBJ( cotanCurvature,
      "cotan_curvature.obj", true  );
  laplaceOperatorHelpers.exportDualQuantityToOBJ( rLocalCurvature,
      "rLocal_curvature.obj", true );

//  if ( ( quantity == "HII" ) || ( quantity == "GII" ) )
//    {
//      trace.beginBlock( "Compute II curvature estimations" );
//      normals         = EH::computeIINormals( vm, K, bimage, surfels );
//      measured_values = ( ( quantity == "H" ) || ( quantity == "Mu1" )
//			  || ( quantity == "HII" ) )
//	? EH::computeIIMeanCurvatures    ( vm, K, bimage, surfels )
//	: EH::computeIIGaussianCurvatures( vm, K, bimage, surfels );
//      for ( unsigned int i = 0; i < measured_values.size(); ++i )
//	measured_curv.addValue( measured_values[ i ] );
//      measured_curv.terminate();
//      trace.info() << "- II curv: avg = " << measured_curv.mean() << std::endl;
//      trace.info() << "- II curv: min = " << measured_curv.min() << std::endl;
//      trace.info() << "- II curv: max = " << measured_curv.max() << std::endl;
//      trace.endBlock();
//    }
//  else // if ( view == "Measure" || view == "Error" )
//    {
//      trace.beginBlock( "Compute normal estimations" );
//      normals  = EH::computeNormals( vm, K, shape, surface, bimage, surfels );
//      trace.endBlock();
//
//      trace.beginBlock( "Computing corrected normal current" );
//      Current C( fastDS, h, vm.count( "crisp" ) );
//      C.setCorrectedNormals( surfels.begin(), surfels.end(), normals.begin() );
//      trace.info() << C << " m-ball-r = " << mr << "(continuous)"
//		   << " " << (mr/h) << " (discrete)" << std::endl;
//      double              area = 0.0;
//      double              intG = 0.0;
//      std::vector<double> mu0( surfels.size() );
//      std::vector<double> mu1( surfels.size() );
//      std::vector<double> mu2( surfels.size() );
//      std::vector<double> muOmega( surfels.size() );
//      bool       mu0_needed = true;
//      bool       mu1_needed = false;
//      bool       mu2_needed = false;
//      bool   muOmega_needed = false;
//      if ( quantity == "Mu0" )     mu0_needed = true;
//      if ( quantity == "Mu1" )     mu1_needed = true;
//      if ( quantity == "Mu2" )     mu2_needed = true;
//      if ( quantity == "MuOmega" ) muOmega_needed = true;
//      if ( quantity == "H" )       mu0_needed = mu1_needed = true;
//      if ( quantity == "G" )       mu0_needed = mu2_needed = true;
//      if ( quantity == "Omega" )   mu0_needed = muOmega_needed = true;
//      trace.info() << "computeAllMu0" << std::endl;
//      if ( mu0_needed ) C.computeAllMu0();
//      trace.info() << "computeAllMu1" << std::endl;
//      if ( mu1_needed ) C.computeAllMu1();
//      trace.info() << "computeAllMu2" << std::endl;
//      if ( mu2_needed ) C.computeAllMu2();
//      trace.info() << "computeAllMuOmega" << std::endl;
//      if ( muOmega_needed ) C.computeAllMuOmega();
//      //#pragma omp parallel for schedule(dynamic)
//      trace.info() << "compute measures" << std::endl;
//      Vertex              i = 0;
//      Vertex              j = surfels.size();
//      for ( auto aSurfel : surfels )
//	{
//	  // std::cout << i << " / " << j << std::endl;
//	  trace.progressBar( i, j );
//	  Vertex v = C.getVertex( aSurfel );
//	  area    += C.mu0( v );
//	  if ( mu0_needed ) mu0[ i ] = C.mu0Ball( v, mr );
//	  if ( mu1_needed ) mu1[ i ] = C.mu1Ball( v, mr );
//	  if ( mu2_needed ) mu2[ i ] = C.mu2Ball( v, mr );
//	  if ( muOmega_needed ) muOmega[ i ] = C.muOmegaBall( v, mr );
//	  ++i;
//	}
//      // Computing total Gauss curvature.
//      if ( mu2_needed )
//	{
//	  trace.info() << "compute total Gauss curvature" << std::endl;
//	  for ( auto f : fastDS.allFaces() ) intG += C.mu2( f );
//	}
//      if ( quantity == "Mu0" )          measured_values = mu0;
//      else if ( quantity == "Mu1" )     measured_values = mu1;
//      else if ( quantity == "Mu2" )     measured_values = mu2;
//      else if ( quantity == "MuOmega" ) measured_values = muOmega;
//      else if ( quantity == "H" )
//	{
//	  measured_values.resize( surfels.size() );
//	  std::transform( mu0.cbegin(), mu0.cend(),
//			  mu1.cbegin(), measured_values.begin(),
//			  [] ( double m0, double m1 ) { return m1 / (2.0*m0); } );
//	}
//      else if ( quantity == "G" )
//	{
//	  measured_values.resize( surfels.size() );
//	  std::transform( mu0.cbegin(), mu0.cend(),
//			  mu2.cbegin(), measured_values.begin(),
//			  [] ( double m0, double m2 ) { return m2 / m0; } );
//	}
//      else if ( quantity == "Omega" )
//	{
//	  measured_values.resize( surfels.size() );
//	  std::transform( mu0.cbegin(), mu0.cend(),
//			  muOmega.cbegin(), measured_values.begin(),
//			  [] ( double m0, double m2 ) { return m2 / sqrt( m0 ); } );
//	}
//      for ( i = 0; i < j; ++i ) measured_curv.addValue( measured_values[ i ] );
//      measured_curv.terminate();
//      trace.info() << "- CNC area      = " << area << std::endl;
//      trace.info() << "- CNC total G   = " << intG << std::endl;
//      trace.info() << "- CNC curv: avg = " << measured_curv.mean() << std::endl;
//      trace.info() << "- CNC curv: min = " << measured_curv.min() << std::endl;
//      trace.info() << "- CNC curv: max = " << measured_curv.max() << std::endl;
//      trace.endBlock();
//    }
//
//#ifdef WITH_VISU3D_QGLVIEWER
//  if ( view != "None" ) {
//    typedef Viewer3D<Space,KSpace> MyViewever3D;
//    typedef Display3DFactory<Space,KSpace> MyDisplay3DFactory;
//
//    trace.beginBlock( "View measure" );
//    trace.info() << "view mode is " << view << std::endl;
//    trace.info() << "#mvalues=" << measured_values.size() << std::endl;
//    trace.info() << "#evalues=" << expected_values.size() << std::endl;
//    MyViewever3D viewer( K );
//    viewer.show();
//    if ( view == "Measure" )
//      displayed_values = measured_values;
//    else if ( view == "Truth" )
//      displayed_values = expected_values;
//    else if ( view == "Error" )
//      displayed_values = EH::absoluteDifference( measured_values, expected_values );
//    trace.info() << "#surfels=" << surfels.size() << std::endl;
//    trace.info() << "#dvalues=" << displayed_values.size() << std::endl;
//    EH::viewSurfelValues( viewer, vm, surfels, displayed_values, normals );
//    EH::viewSurfaceIsolines( viewer, vm, surface, surfels, displayed_values );
//    viewer << MyViewever3D::updateDisplay;
//    application.exec();
//    trace.endBlock();
//  }
//#endif
//  if ( ! vm.count( "polynomial" ) ) return 0;
//
//  trace.beginBlock( "Output statistics" );
//  auto error_fname = vm[ "error" ].as<std::string>();
//  std::ofstream ferr;
//  ferr.open ( error_fname.c_str(),
//	      std::ofstream::out | std::ofstream::app );
//  if ( ! ferr.good() )
//    trace.warning() << "Unable to open file " << error_fname << std::endl;
//  ferr << "######################################################################"
//       << std::endl;
//  ferr << "# ";
//  for ( int i = 0; i < argc; ++i ) ferr << " " << argv[ i ];
//  ferr << std::endl;
//  ferr << "#---------------------------------------------------------------------"
//       << std::endl;
//  ferr << "# Q=" << quantity
//       << " P="  << vm[ "polynomial" ].as<std::string>()
//       << " mr= " << mr << "(continuous)"
//       << " mrd=" << (mr/h) << " (discrete)" << std::endl;
//  ferr << "# h size l1 l2 loo"
//       << " m_mean m_dev m_min m_max"
//       << " exp_mean exp_dev exp_min exp_max " << std::endl;
//  ferr << h << " " << measured_values.size()
//       << " " << EH::normL1 ( measured_values, expected_values )
//       << " " << EH::normL2 ( measured_values, expected_values )
//       << " " << EH::normLoo( measured_values, expected_values )
//       << " " << measured_curv.mean()
//       << " " << sqrt( measured_curv.variance() )
//       << " " << measured_curv.min()
//       << " " << measured_curv.max()
//       << " " << expected_curv.mean()
//       << " " << sqrt( expected_curv.variance() )
//       << " " << expected_curv.min()
//       << " " << expected_curv.max()
//       << std::endl;
//  ferr << "#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
//       << std::endl;
//  ferr.close();
//  trace.endBlock();

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
