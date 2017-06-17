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
#include <string>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
// always include EigenSupport.h before any other Eigen headers
#include "DGtal/math/linalg/EigenSupport.h"
#include "DGtal/base/Common.h"
#include "DGtal/base/Clone.h"
#include "DGtal/math/MPolynomial.h"
#include "DGtal/io/readers/MPolynomialReader.h"
#include "DGtal/shapes/implicit/ImplicitPolynomial3Shape.h"


//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{
  namespace po = boost::program_options;

  /////////////////////////////////////////////////////////////////////////////
  // template class EstimatorHelpers
  template < typename TKSpace >
  struct EstimatorHelpers
  {
    typedef TKSpace                             KSpace;
    typedef typename KSpace::Space              Space;
    typedef typename Space::Point               Point;
    typedef typename Space::Vector              Vector;
    typedef typename Space::RealVector          RealVector;
    typedef typename Space::RealPoint           RealPoint;
    typedef typename RealVector::Component      Scalar;
    typedef MPolynomial< Space::dimension, Scalar > PolynomialN;
    typedef ImplicitPolynomial3Shape<Space>     ImplicitShape;

    /// Add options for implicit shape given by a polynomial.
    static void optionsImplicitShape( po::options_description& desc )
    {
      desc.add_options()
	( "polynomial,p", po::value<std::string>(), "the implicit polynomial whose zero-level defines the shape of interest." );
    }
    
    /// Builds a 3D implicit shape from argument "-polynomial".
    static CountedPtr<ImplicitShape> makeImplicitShape( const po::variables_map& vm )
    {
      typedef MPolynomialReader< Space::dimension, Scalar> Polynomial3Reader;
      std::string poly_str = vm[ "polynomial" ].as<std::string>();
      PolynomialN poly;
      Polynomial3Reader reader;
      std::string::const_iterator iter
	= reader.read( poly, poly_str.begin(), poly_str.end() );
      if ( iter != poly_str.end() )
	{
	  trace.error() << "ERROR reading polynomial: I read only <"
			<< poly_str.substr( 0, iter - poly_str.begin() )
			<< ">, and I built P=" << poly << std::endl;
	}
      return CountedPtr<ImplicitShape>( new ImplicitShape( poly ) );
    }
    
    
  }; // end of class EstimatorHelpers

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined EstimatorHelpers_h

#undef EstimatorHelpers_RECURSES
#endif // else defined(EstimatorHelpers_RECURSES)
