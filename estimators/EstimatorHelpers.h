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
#include "DGtal/io/readers/MPolynomialReader.h"
#include "DGtal/shapes/implicit/ImplicitPolynomial3Shape.h"
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/topology/DigitalSurface.h"


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
	("minAABB,a",  po::value<double>()->default_value( -10.0 ), "the min value of the AABB bounding box (domain)" )
	("maxAABB,A",  po::value<double>()->default_value( 10.0 ), "the max value of the AABB bounding box (domain)" )
	("gridstep,g", po::value< double >()->default_value( 1.0 ), "the gridstep that defines the digitization (often called h). " );
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

    static CountedPtr<BinaryImage>
    makeImage( CountedPtr<ImplicitDigitalShape> dshape )
    {
      const Domain shapeDomain    = dshape->getDomain();
      CountedPtr<BinaryImage> img ( new BinaryImage( shapeDomain ) );
      std::transform( shapeDomain.begin(), shapeDomain.end(),
		      img->begin(),
		      [dshape] ( const Point& p ) { return (*dshape)(p); } );
      // KanungoPredicate* noisified_dshape = new KanungoPredicate( dshape, shapeDomain, noiseLevel );
      return img;
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
