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

#pragma once

/**
 * @file EstimatorHelpers.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 * @author Marion Foare (\c marion.foare@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2017/06/17
 *
 * Header file for module EstimatorHelpers.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(EstimatorHelpers_RECURSES)
#error Recursive header files inclusion detected in EstimatorHelpers.h
#else // defined(EstimatorHelpers_RECURSES)
/** Prevents recursive inclusion of headers. */
#define EstimatorHelpers_RECURSES

#if !defined EstimatorHelpers_h
/** Prevents repeated inclusion of headers. */
#define EstimatorHelpers_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <sstream>
#include <string>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
// always include EigenSupport.h before any other Eigen headers
// #include "DGtal/math/linalg/EigenSupport.h"
#include "DGtal/base/Common.h"
#include "DGtal/base/Clone.h"
#include "DGtal/math/MPolynomial.h"
#include "DGtal/math/Statistic.h"
#include "DGtal/io/readers/MPolynomialReader.h"
#include "DGtal/shapes/implicit/ImplicitPolynomial3Shape.h"
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/shapes/ShapeGeometricFunctors.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/SurfelAdjacency.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"
#include "DGtal/geometry/volumes/KanungoNoise.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/geometry/surfaces/estimation/TrueDigitalSurfaceLocalEstimator.h"
#include "DGtal/geometry/surfaces/estimation/VoronoiCovarianceMeasureOnDigitalSurface.h"
#include "DGtal/geometry/surfaces/estimation/VCMDigitalSurfaceLocalEstimator.h"
#include "DGtal/geometry/surfaces/estimation/IIGeometricFunctors.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantCovarianceEstimator.h"


//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{
  namespace po = boost::program_options;

  /////////////////////////////////////////////////////////////////////////////
  // template class EstimatorHelpers
  template < typename TKSpace >
  struct EstimatorHelpers
  {
    typedef TKSpace                                  KSpace;
    typedef typename KSpace::Space                   Space;
    typedef typename Space::Point                    Point;
    typedef typename Space::Vector                   Vector;
    typedef typename Space::RealVector               RealVector;
    typedef typename Space::RealPoint                RealPoint;
    typedef typename RealVector::Component           Scalar;
    typedef MPolynomial< Space::dimension, Scalar >  PolynomialN;
    typedef ImplicitPolynomial3Shape<Space>          ImplicitShape;
    typedef GaussDigitizer< Space, ImplicitShape >   ImplicitDigitalShape;
    typedef typename ImplicitDigitalShape::Domain    Domain;
    typedef ImageContainerBySTLVector<Domain, bool>  BinaryImage;
    typedef LightImplicitDigitalSurface< KSpace, BinaryImage > SurfaceContainer;
    typedef DigitalSurface< SurfaceContainer >       Surface;
    typedef typename Surface::Surfel                 Surfel;
    typedef functors::ShapeGeometricFunctors::ShapeNormalVectorFunctor<ImplicitShape> NormalFunctor;
    typedef TrueDigitalSurfaceLocalEstimator<KSpace, ImplicitShape, NormalFunctor> TrueNormalEstimator;
    typedef DGtal::Statistic<Scalar>                 AngleDevStatistic;
    
    // ------------------- parsing related functions -----------------------------
    
    /// Parses command line given as \a argc, \a argv according to \a
    /// opts and fills in storage map \a vm.
    ///
    /// @param[in] opts the allowed options.
    /// @param[in] argc the size of array \a argv.
    ///
    /// @param[in] argv an array of C-style strings containing the
    /// desired command-line options.
    ///
    /// @param[inout] vm the map options to variables that is updated.
    static bool args2vm( const po::options_description& opts,
			 int argc, char** argv,
			 po::variables_map& vm )
    {
      bool parseOK = true;
      try
	{
	  po::store( po::parse_command_line( argc, argv, opts ), vm );
	}
      catch( const std::exception & ex )
	{
	  parseOK = false;
	  trace.error() << "[EstimatorHelpers::args2vm]"
			<< " Error checking program options: " << ex.what() << std::endl;
	}
      return parseOK;
    }
    
    /// Parses command line given as a vector of strings \a args
    /// according to \a opts and fills in storage map \a vm.
    ///
    /// @param[in] opts the allowed options.
    /// @param[in] args a vector of strings containing the
    /// desired command-line options cut by words.
    ///
    /// @param[inout] vm the map options to variables that is updated.
    static bool vectorString2vm( const po::options_description& opts,
				 std::vector< std::string > args,
				 po::variables_map& vm )
    {
      bool parseOK = true;
      try
	{
	  po::store( po::command_line_parser( args ).options( opts ).run(), vm );
	}
      catch( const std::exception & ex )
	{
	  parseOK = false;
	  trace.error() << "[EstimatorHelpers::vectorString2vm]"
			<< " Error checking program options: " << ex.what() << std::endl;
	}
      return parseOK;
    }
    
    /// Parses a string, cuts it into words, analyzes it
    /// according to \a opts and fills in storage map \a vm.
    ///
    /// @param[in] opts the allowed options.
    /// @param[in] s a string containing the
    /// desired command-line options separated by spaces.
    ///
    /// @param[inout] vm the map options to variables that is updated.
    static bool string2vm( const po::options_description& opts,
			   const std::string& s,
			   po::variables_map& vm )
    {
      // construct a stream from the string
      std::stringstream sstr( s );
      // use stream iterators to copy the stream to the vector as
      // whitespace separated strings
      std::istream_iterator<std::string> it(sstr);
      std::istream_iterator<std::string> end;
      std::vector<std::string> results(it, end);
      return vectorString2vm( opts, results, vm );
    }

    
    // ------------------- setting options related functions ---------------------------    
    /// Add options for implicit shape given by a polynomial.
    static void optionsImplicitShape( po::options_description& desc )
    {
      desc.add_options()
	( "polynomial,p", po::value<std::string>(), "the implicit polynomial whose zero-level defines the shape of interest." );
    }

    /// Add options for implicit shape digitization.
    static void optionsDigitizedShape( po::options_description& desc )
    {
      desc.add_options()
	("minAABB,a",  po::value<Scalar>()->default_value( -10.0 ), "the min value of the AABB bounding box (domain)" )
	("maxAABB,A",  po::value<Scalar>()->default_value( 10.0 ), "the max value of the AABB bounding box (domain)" )
	("gridstep,g", po::value< Scalar >()->default_value( 1.0 ), "the gridstep that defines the digitization (often called h). " );
    }

    /// Add options for implicit shape digitization.
    static void optionsNoisyImage( po::options_description& desc )
    {
      desc.add_options()
	("noise,N", po::value<Scalar>()->default_value( 0.0 ), "the Kanungo noise level l=arg, with l^d the probability that a point at distance d is flipped inside/outside." );
    }

    /// Add options for implicit shape digitization.
    static void optionsNormalEstimators( po::options_description& desc )
    {
      desc.add_options()
	("estimator,e", po::value<std::string>()->default_value( "True" ), "the chosen normal estimator: True | VCM | II | Trivial" )
	("R-radius,R", po::value<Scalar>()->default_value( 5 ), "the constant for parameter R in R(h)=R h^alpha (VCM)." )
	("r-radius,r", po::value<Scalar>()->default_value( 3 ), "the constant for parameter r in r(h)=r h^alpha (VCM,II,Trivial)." )
	("kernel,k", po::value<std::string>()->default_value( "hat" ), "the function chi_r, either hat or ball." )
	("alpha", po::value<Scalar>()->default_value( 0.0 ), "the parameter alpha in r(h)=r h^alpha (VCM)." )
	("trivial-radius,t", po::value<Scalar>()->default_value( 3 ), "the parameter t defining the radius for the Trivial estimator. Also used for reorienting the VCM." )
	("embedding,E", po::value<int>()->default_value( 0 ), "the surfel -> point embedding for VCM estimator: 0: Pointels, 1: InnerSpel, 2: OuterSpel." );
    }

    // ------------------- Shapes related functions --------------------------------
    
    /// Builds a 3D implicit shape from argument "-polynomial".
    ///
    /// @param[in] vm the options sets in the variable map (arguments
    /// given to the program). Recognized parameters are given in \ref
    /// optionsImplicitShape.
    ///
    /// @return a smart pointer on the created implicit shape.
    static CountedPtr<ImplicitShape>
    makeImplicitShape( const po::variables_map& vm )
    {
      typedef MPolynomialReader< Space::dimension, Scalar> Polynomial3Reader;
      std::string poly_str = vm[ "polynomial" ].as<std::string>();
      PolynomialN poly;
      Polynomial3Reader reader;
      std::string::const_iterator iter
	= reader.read( poly, poly_str.begin(), poly_str.end() );
      if ( iter != poly_str.end() )
	{
	  trace.error() << "[EstimatorHelpers::makeImplicitShape]"
			<< " ERROR reading polynomial: I read only <"
			<< poly_str.substr( 0, iter - poly_str.begin() )
			<< ">, and I built P=" << poly << std::endl;
	}
      return CountedPtr<ImplicitShape>( new ImplicitShape( poly ) );
    }

    /// Makes the Gauss digitization of the given implicit shape
    /// according to parameters given as arguments to the program.
    ///
    /// @param[in] vm the options sets in the variable map (arguments
    /// given to the program). Recognized parameters are given in \ref
    /// optionsDigitizedShape.
    ///
    /// @param[in] shape the implicit shape.
    /// @param[out] K the Khalimsky space whose domain encompasses the digital shape.
    /// @return a smart pointer on the created implicit digital shape.
    static CountedPtr<ImplicitDigitalShape>
    makeDigitizedShape( const po::variables_map& vm,
			CountedPtr<ImplicitShape> shape,
			KSpace& K )
    {
      Scalar min_x = vm[ "minAABB" ].as<Scalar>();
      Scalar max_x = vm[ "maxAABB" ].as<Scalar>();
      Scalar h     = vm[ "gridstep" ].as<Scalar>();
      RealPoint p1( min_x, min_x, min_x );
      RealPoint p2( max_x, max_x, max_x );
      CountedPtr<ImplicitDigitalShape> dshape( new ImplicitDigitalShape() );
      dshape->attach( *shape );
      dshape->init( p1, p2, h );
      Domain domain = dshape->getDomain();
      if ( ! K.init( domain.lowerBound(), domain.upperBound(), true ) )
	trace.error() << "[EstimatorHelpers::makeImplicitDigitalShape]"
		      << " Error building Khalimsky space K=" << K << std::endl;
      return dshape;
    }

    /// Vectorizes an implicitly defined digital shape into a binary image.
    ///
    /// @param[in] dshape a smart pointer on an implicit digital shape.
    /// @return a smart pointer on a binary image that samples the digital shape.
    static CountedPtr<BinaryImage>
    makeImage( CountedPtr<ImplicitDigitalShape> dshape )
    {
      const Domain shapeDomain    = dshape->getDomain();
      CountedPtr<BinaryImage> img ( new BinaryImage( shapeDomain ) );
      std::transform( shapeDomain.begin(), shapeDomain.end(),
		      img->begin(),
		      [dshape] ( const Point& p ) { return (*dshape)(p); } );
      return img;
    }

    /// Vectorizes an implicitly defined digital shape into a binary
    /// image with added Kanungo noise.
    ///
    /// @param[in] dshape a smart pointer on an implicit digital shape.
    ///
    /// @param[in] noiseLevel the parameter that tunes the amount of
    /// Kanungo noise (0.0 is none, 0.5 is quite a lot, 1.0 means
    /// random image).
    ///
    /// @return a smart pointer on a binary image that samples the
    /// digital shape with added noise.
    static CountedPtr<BinaryImage>
    makeNoisyImage( CountedPtr<ImplicitDigitalShape> dshape, Scalar noiseLevel )
    {
      typedef KanungoNoise< ImplicitDigitalShape, Domain > KanungoPredicate;
      const Domain shapeDomain    = dshape->getDomain();
      CountedPtr<BinaryImage> img ( new BinaryImage( shapeDomain ) );
      KanungoPredicate noisy_dshape( *dshape, shapeDomain, noiseLevel );
      std::transform( shapeDomain.begin(), shapeDomain.end(),
		      img->begin(),
		      [&noisy_dshape] ( const Point& p ) { return noisy_dshape(p); } );
      return img;
    }

    /// Vectorizes an implicitly defined digital shape into a binary
    /// image, and possibly add Kanungo noise to the result depending
    /// on parameters given in \a vm.
    ///
    /// @param[in] vm the options sets in the variable map (arguments
    /// given to the program). Recognized parameters are given in \ref
    /// optionsNoisyImage.
    ///
    /// @param[in] dshape a smart pointer on an implicit digital shape.
    ///
    /// @return a smart pointer on a binary image that samples the
    /// digital shape with added noise.
    static CountedPtr<BinaryImage>
    makeNoisyOrNotImage( const po::variables_map& vm,
			 CountedPtr<ImplicitDigitalShape> dshape )
    {
      Scalar noiseLevel = vm[ "noise" ].as<Scalar>();
      if ( noiseLevel == 0.0 )
	return makeImage( dshape );
      else
	return makeNoisyImage( dshape, noiseLevel );
    }

    /// Builds a digital surface from a space \a K and a binary image \a bimage.
    ///
    /// @param[in] K the Khalimsky space whose domain encompasses the digital shape.
    /// @param[in] bimage a binary image representing the characteristic function of a digital shape.
    ///
    /// @return a smart pointer on a (light) digital surface that
    /// represents the boundary of the digital shape (at least a big
    /// component).
    static CountedPtr<Surface>
    makeDigitalSurface( const KSpace& K, CountedPtr<BinaryImage> bimage )
    {
      SurfelAdjacency< KSpace::dimension > surfAdj( true );

      // We have to search for a surfel that belong to a big connected component.
      CountedPtr<Surface> ptrSurface;
      Surfel              bel;
      Scalar              minsize    = bimage->extent().norm();
      unsigned int        nb_surfels = 0;
      unsigned int        tries      = 0;
      do {
        try { // Search initial bel
          bel = Surfaces<KSpace>::findABel( K, *bimage, 10000 );
        } catch (DGtal::InputException e) {
          trace.error() << "[EstimatorHelpers::makeDigitalSurface]"
			<< " ERROR Unable to find bel." << std::endl;
          return ptrSurface;
        }
	// this pointer will be acquired by the surface.
        SurfaceContainer* surfContainer = new SurfaceContainer( K, *bimage, surfAdj, bel );
        ptrSurface = CountedPtr<Surface>( new Surface( surfContainer ) ); // acquired
        nb_surfels = ptrSurface->size();
      } while ( ( nb_surfels < 2 * minsize ) && ( tries++ < 150 ) );
      if( tries >= 150 )
        {
          trace.error() << "[EstimatorHelpers::makeDigitalSurface]"
			<< "ERROR cannot find a proper bel in a big enough component."
			<< std::endl;
          return ptrSurface;
        }
      return ptrSurface;
    }

    /// Given a digital surface, returns a vector of surfels in a depth-first
    /// traversal order (best for Integral Invariant).
    ///
    /// @param[in] surface a smort pointer on a digital surface.
    /// @return a range of surfels as a vector.
    static std::vector< Surfel >
    computeDepthFirstSurfelRange( CountedPtr<Surface> surface )
    {
      typedef DepthFirstVisitor< Surface > Visitor;
      typedef GraphVisitorRange< Visitor > VisitorRange;
      std::vector< Surfel > result;
      VisitorRange range( new Visitor( *surface, *( surface->begin() ) ) );
      std::for_each( range.begin(), range.end(),
		     [&result] ( Surfel s ) { result.push_back( s ); } );
      return result;
    }
       

    /// Given a space \a K, an implicit \a shape, a sequence of \a
    /// surfels, and a gridstep \a h, returns the true normals at the
    /// specified surfels, in the same order.
    ///
    /// @param[in] K the Khalimsky space whose domain encompasses the digital shape.
    /// @param[in] shape the implicit shape.
    /// @param[in] surfels the sequence of surfels at which we compute the normals
    ///
    /// @return the vector containing the true normals, in the same
    /// order as \a surfels.
    static std::vector< RealVector >
    computeTrueNormals( const KSpace&             K,
			CountedPtr<ImplicitShape> shape,
			std::vector< Surfel >     surfels )
    {
      std::vector< RealVector > n_true_estimations;
      TrueNormalEstimator       true_estimator;
      Scalar                    h = 1.0; // useless here.
      true_estimator.attach( *shape );
      true_estimator.setParams( K, NormalFunctor(), 20, 0.1, 0.01 );
      true_estimator.init( h, surfels.begin(), surfels.end() );
      true_estimator.eval( surfels.begin(), surfels.end(),
			   std::back_inserter( n_true_estimations ) );
      return n_true_estimations;
    }

    /// Given a digital surface \a surface, a sequence of \a surfels,
    /// and some parameters \a vm, returns the normal Voronoi
    /// Covariance Measure (VCM) estimation at the specified surfels,
    /// in the same order.
    ///
    /// @param[in] vm the options sets in the variable map (arguments
    /// given to the program). Recognized parameters are given in \ref
    /// optionsNoisyImage.
    /// @param[in] surface the digital surface
    /// @param[in] surfels the sequence of surfels at which we compute the normals
    ///
    /// @return the vector containing the estimated normals, in the
    /// same order as \a surfels.
    static std::vector< RealVector >
    computeVCMNormals( const po::variables_map& vm,
		       CountedPtr<Surface>      surface,
		       std::vector< Surfel >    surfels )
    {
      typedef ExactPredicateLpSeparableMetric<Space,2> Metric;
      std::vector< RealVector > n_estimations;
      std::string kernel = vm[ "kernel" ].as<std::string>();
      Scalar      h      = vm[ "gridstep" ].as<Scalar>();
      Scalar      R      = vm[ "R-radius" ].as<Scalar>();
      Scalar      r      = vm[ "r-radius" ].as<Scalar>();
      Scalar      t      = vm[ "trivial-radius" ].as<Scalar>();
      Scalar      alpha  = vm[ "alpha" ].as<Scalar>();
      int      embedding = vm[ "embedding" ].as<int>();
      // Adjust parameters according to gridstep if specified.
      if ( alpha != 0.0 ) R *= pow( h, alpha-1.0 );
      if ( alpha != 0.0 ) r *= pow( h, alpha-1.0 );
      Surfel2PointEmbedding embType = embedding == 0 ? Pointels :
                                      embedding == 1 ? InnerSpel : OuterSpel;
      trace.info() << "- R=" << R << " r=" << r << " t=" << t << std::endl;
      if ( kernel == "hat" ) {
	typedef functors::HatPointFunction<Point,Scalar>             KernelFunction;
	typedef VoronoiCovarianceMeasureOnDigitalSurface
	  < SurfaceContainer, Metric, KernelFunction >               VCMOnSurface;
	typedef functors::VCMNormalVectorFunctor<VCMOnSurface>       NormalFunctor;
	typedef VCMDigitalSurfaceLocalEstimator
	  < SurfaceContainer, Metric, KernelFunction, NormalFunctor> VCMNormalEstimator;
	KernelFunction chi_r( 1.0, r );
	VCMNormalEstimator estimator;
	estimator.attach( *surface );
	estimator.setParams( embType, R, r, chi_r, t, Metric(), true );
	estimator.init( h, surfels.begin(), surfels.end() );
	estimator.eval( surfels.begin(), surfels.end(),
			std::back_inserter( n_estimations ) );
      } else if ( kernel == "ball" ) {
	typedef functors::BallConstantPointFunction<Point,Scalar>    KernelFunction;
	typedef VoronoiCovarianceMeasureOnDigitalSurface
	  < SurfaceContainer, Metric, KernelFunction >               VCMOnSurface;
	typedef functors::VCMNormalVectorFunctor<VCMOnSurface>       NormalFunctor;
	typedef VCMDigitalSurfaceLocalEstimator
	  < SurfaceContainer, Metric, KernelFunction, NormalFunctor> VCMNormalEstimator;
	KernelFunction chi_r( 1.0, r );
	VCMNormalEstimator estimator;
	estimator.attach( *surface );
	estimator.setParams( embType, R, r, chi_r, t, Metric(), true );
	estimator.init( h, surfels.begin(), surfels.end() );
	estimator.eval( surfels.begin(), surfels.end(),
			std::back_inserter( n_estimations ) );
      }
      return n_estimations;
    }

    /// Given a digital surface \a surface, a sequence of \a surfels,
    /// and some parameters \a vm, returns the normal Integral
    /// Invariant (VCM) estimation at the specified surfels, in the
    /// same order.
    ///
    /// @param[in] vm the options sets in the variable map (arguments
    /// given to the program). Recognized parameters are given in \ref
    /// optionsNoisyImage.
    /// @param[in] K the Khalimsky space whose domain encompasses the digital shape.
    /// @param[in] surfels the sequence of surfels at which we compute the normals
    ///
    /// @return the vector containing the estimated normals, in the
    /// same order as \a surfels.
    ///
    /// @note It is better to have surfels in a specific order, as
    /// given for instance by computeDepthFirstSurfelRange.
    static std::vector< RealVector >
    computeIINormals( const po::variables_map& vm,
		      const KSpace&            K,
		      CountedPtr<BinaryImage>  bimage,
		      std::vector< Surfel >    surfels )
    {
      typedef functors::IINormalDirectionFunctor<Space> IINormalFunctor;
      typedef IntegralInvariantCovarianceEstimator
	<KSpace, BinaryImage, IINormalFunctor>          IINormalEstimator;
      std::vector< RealVector > n_estimations;
      Scalar h     = vm[ "gridstep" ].as<Scalar>();
      Scalar r     = vm[ "r-radius" ].as<Scalar>();
      Scalar alpha = vm[ "alpha"    ].as<Scalar>();
      if ( alpha != 0.0 ) r *= pow( h, alpha-1.0 );
      trace.info() << " r=" << r << std::endl;
      IINormalEstimator ii_estimator( K, *bimage );
      ii_estimator.setParams( r );
      ii_estimator.init( h, surfels.begin(), surfels.end() );
      ii_estimator.eval( surfels.begin(), surfels.end(),
			 std::back_inserter( n_estimations ) );
      return n_estimations;
    }
    
    /// Orient \a v so that it points in the same direction as \a
    /// ref_v (scalar product is then non-negative afterwards).
    ///
    /// @param[in]    ref_v the vectors having the reference orientation.
    /// @param[inout] v the vectors to reorient.
    static void
    orientVectors( const std::vector< RealVector > & ref_v,
		   std::vector< RealVector > &       v )
    {
      std::transform( ref_v.cbegin(), ref_v.cend(), v.cbegin(), v.begin(), 
		      [] ( RealVector rw, RealVector w )
		      { return rw.dot( w ) >= 0.0 ? w : -w; } );
    }

    /// Computes the statistic that measures the angle differences
    /// between the two arrays of unit vectors.
    ///
    /// @param[in] v1 the first array of unit vectors (normals)
    /// @param[in] v2 the second array of unit vectors (normals)
    /// @return their angle difference as a statistic.
    static AngleDevStatistic
    measureAngleDeviation( const std::vector< RealVector > & v1,
			   const std::vector< RealVector > & v2 )
    {
      AngleDevStatistic stat;
      for ( auto it1 = v1.cbegin(), it2 = v2.cbegin(), itE1 = v1.end();
	    it1 != itE1; ++it1, ++it2 )
	{
          Scalar angle_error = acos( (*it1).dot( *it2 ) );
          stat.addValue( angle_error );
	}
      stat.terminate();
      return stat;
    }

  }; // END of class EstimatorHelpers

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined EstimatorHelpers_h

#undef EstimatorHelpers_RECURSES
#endif // else defined(EstimatorHelpers_RECURSES)
