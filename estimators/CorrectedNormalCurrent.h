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
 * @file CorrectedNormalCurrent.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 *
 * @date 2017/06/19
 *
 * Header file for module CorrectedNormalCurrent.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(CorrectedNormalCurrent_RECURSES)
#error Recursive header files inclusion detected in CorrectedNormalCurrent.h
#else // defined(CorrectedNormalCurrent_RECURSES)
/** Prevents recursive inclusion of headers. */
#define CorrectedNormalCurrent_RECURSES

#if !defined CorrectedNormalCurrent_h
/** Prevents repeated inclusion of headers. */
#define CorrectedNormalCurrent_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <map>
#include "DGtal/base/Common.h"
#include "DGtal/topology/CDigitalSurfaceContainer.h"

//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // class CorrectedNormalCurrent
  /**
     Description of class 'CorrectedNormalCurrent' <p> \brief Aim:
     Represent a current over a digital surface, whose tangent plane
     is described by an additional normal vector field. It is useful
     to define geometric measures onto digital surface, which provide
     area, mean and gaussian curvature information.
     
     @tparam TSurface any type of digital surface.
  */
  template <typename TSurface>
  class CorrectedNormalCurrent
  {
  public:
    typedef TSurface                                  Surface;
    typedef typename Surface::DigitalSurfaceContainer DigitalSurfaceContainer;
    BOOST_CONCEPT_ASSERT(( concepts::CDigitalSurfaceContainer<DigitalSurfaceContainer> ));

    typedef typename Surface::KSpace                  KSpace;
    typedef typename Surface::Cell                    Cell;
    typedef typename Surface::SCell                   SCell;
    typedef typename Surface::Surfel                  Surfel;
    typedef typename Surface::Vertex                  Vertex;
    typedef typename Surface::Arc                     Arc;
    typedef typename Surface::Face                    Face;
    typedef typename Surface::ConstIterator           ConstIterator;
    typedef typename KSpace::Space                    Space;
    typedef typename Space::Point                     Point;
    typedef typename Space::Vector                    Vector;
    typedef typename Space::RealPoint                 RealPoint;
    typedef typename Space::RealVector                RealVector;
    typedef typename RealVector::Component            Scalar;
    typedef std::map< Surfel, RealVector >            NormalVectorField;

    // Checks that dimension is 3.
    BOOST_STATIC_ASSERT(( Kspace::dimension == 3 ));
    
    // ----------------------- Standard services ------------------------------
  public:
  
    /**
     * Destructor.
     */
    ~CorrectedNormalCurrent() {}

    /// Default constructor. The object is invalid.
    CorrectedNormalCurrent()
      : theSurface( nullptr ) {}

    /**
       Constructor from surface. The surface is shared by the
       current. Computes also its natural tangent bundle
       (i.e. using the plane defined by each surfel).
 
       @param surface the digital surface which defines the current.
       @param h the digitization gridstep of the surface.
    */
    CorrectedNormalCurrent( Alias<Surface> surface, Scalar h = 1.0 )
      : theSurface( surface )
    {
      setParams( h );
      computeTrivialNormals();
      setCorrectedNormals( myTrivialNormals );
    }
  
    /**
     * Copy constructor.
     * @param other the object to clone.
     */
    CorrectedNormalCurrent ( const CorrectedNormalCurrent & other ) = default;
  
    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
    CorrectedNormalCurrent & operator= ( const CorrectedNormalCurrent & other ) = default;

    /// Sets the digitization grid step of the digital surface.
    void setParams( Scalar h )
    {
      myH = h;
      myEmbedder.init( myH );
    }
    
    /// Attaches the given surface to this current. Computes also its
    /// natural tangent bundle.
    /// @param surface the digital surface which defines the current.
    void attach( Alias<Surface> surface )
    {
      theSurface = surface;
      myTrivialNormals.clear();
      computeTrivialNormals();
      setCorrectedNormals( myTrivialNormals );
    }
    
    /// @return the Khalimsky space associated with this current.
    const KSpace& space() const
    {
      ASSERT( theSurface != 0 );
      return theSurface->container().space();
    }
    
    /// @return the surface associated with this current.
    Surface& surface()
    {
      ASSERT( theSurface != 0 );
      return *theSurface;
    }

    /// @return an iterator on the first surfel of the surface.
    ConstIterator begin() const
    {
      ASSERT( theSurface != 0 );
      return theSurface->begin();
    }

    /**
       @return a ConstIterator after the last surfel of the surface.
    */
    ConstIterator end() const
    {
      ASSERT( theSurface != 0 );
      return theSurface->end();
    }

    // ----------------------- Normal services --------------------------------------
  public:

    /// Computes the trivial normals of the attached surface.
    void computeTrivialNormals()
    {
      const KSpace& K = space();
      for ( auto s : *this )
	{
	  Dimension  k = K.sOrthDir( s );
	  bool  direct = K.sDirect( s, k );
	  RealVector t = RealVector::zero;
	  t[ k ]       = direct ? -1.0 : 1.0;
	  myTrivialNormals[ s ] = t ;
	}
    }

    /// Sets the corrected normal vector field of the current.
    /// @param nvf a map Surfel -> RealVector
    void setCorrectedNormals( const NormalVectorField& nvf )
    {
      myCorrectedNormals = nvf;
    }

    /// Sets the corrected normal vector field of the current.
    template <typename SurfelIterator, typename RealVectorIterator>
    void setCorrectedNormals( SurfelIterator itS, SurfelIterator itSEnd,
			      RealVectorIterator itRV )
    {
      myCorrectedNormals.clear();
      for ( ; itS != itSEnd; ++itS )
	myCorrectedNormals[ *itS ] = *itRV++;
    }

    // ----------------------- Measure services -------------------------------
  public:
    /// @param[in] c any cell.
    /// @return the centroid of the cell c in real space;
    RealPoint uCentroid( Cell c ) const
    {
      Point    kp = space().uKCoords( c );
      Point    lo = Point( kp[ 0 ] / 2, kp[ 1 ] / 2, kp[ 2 ] / 2 );
      Point    hi = Point( (kp[ 0 ]+1) / 2, (kp[ 1 ]+1) / 2, (kp[ 2 ]+1) / 2 );
      return   0.5 * ( myEmbedder( lo ) + myEmbedder( hi ) );
    }
    
    /// @param[in] c any signed cell.
    /// @return the centroid of the signed cell c in real space;
    RealPoint sCentroid( SCell c ) const
    {
      Point    kp = space().sKCoords( c );
      Point    lo = Point( kp[ 0 ] / 2, kp[ 1 ] / 2, kp[ 2 ] / 2 );
      Point    hi = Point( (kp[ 0 ]+1) / 2, (kp[ 1 ]+1) / 2, (kp[ 2 ]+1) / 2 );
      return   0.5 * ( myEmbedder( lo ) + myEmbedder( hi ) );
    }
    
    /// Computes the d-dimensional Hausdorff measure of a d-dimensional cell.
    Scalar Hmeasure( Cell c ) const
    {
      const KSpace & K = space();
      Scalar H = 1.0;
      for ( Dimension k = K.udim( c ); k != 0; ++k ) H *= myH;
      return H;
    }

    /// Computes an approximation of the Hausdorff measure of
    /// intersection of the ball of radius \a r and center \a p with
    /// the given cell. Cells are embedded naturally in the grid of
    /// step h.
    Scalar Hmeasure( Cell c, RealPoint p, Scalar r ) const
    {
      const KSpace & K = space();
      auto       faces = K.uFaces( c );
      faces.push_back( c );
      bool       first = true;
      Scalar     d_max = 0.0;
      Scalar     d_min = 0.0;
      for ( Cell f : faces )
	{
	  RealPoint x = uCentroid( f );
	  Scalar    d = ( x - p ).norm();
	  if ( first ) { d_min = d; first = false; }
	  d_max = std::max( d_max, d );
	  d_min = std::min( d_min, d );
	}
      Scalar H = Hmeasure( c );
      if      ( d_max <= r     ) return H;
      else if ( r     <= d_min ) return 0.0;
      return H * ( r - d_min ) / ( d_max - d_min );
    }

    /// \f$ \mu_0 \f$ Lipschitz-Killing measure. It corresponds to a
    /// corrected area measure, and is non null only on 2-cells (or Vertex).
    ///
    /// @param[in] c any 2-dimensional cell (or a Vertex in 3D digital
    /// surfaces).
    ///
    /// @return the corrected area measure \f$ \mu_0 := \cos \alpha
    /// d\mathcal{H}^2 \f$.
    Scalar mu0( Cell c ) const
    {
      const KSpace & K = space();
      if ( K.udim( c ) != 2 ) return 0.0;
      return Hmeasure( c ) * myTrivialNormals[ c ].dot( myCorrectedNormals[ c ] );
    }

    /// \f$ \mu_1 \f$ Lipschitz-Killing measure. It corresponds to a
    /// measure of mean curvature, since normal vectors may turn
    /// around an edge, and is non null only on 1-cells (or Arc). The
    /// measure is oriented, but gives the same result of an arc and
    /// its opposite.
    ///
    /// @param[in] a any arc (ie a 1-cell in-between two 2-cells on 3D digital
    /// surfaces).
    ///
    /// @return the corrected mean curvature measure \f$ \mu_1 := \Psi
    /// \vec{e} \cdot \vec{e}_1 d\mathcal{H}^1 \f$.
    Scalar mu1( Arc a ) const
    {
      const KSpace & K = space();
      Arc            b = K.opposite( a );
      auto      face_a = K.facesAroundArc( a );
      if ( face_a.size() != 1 ) return 0.0; // boundary face.
      auto      face_b = K.facesAroundArc( b );
      if ( face_b.size() != 1 ) return 0.0; // boundary face.
      Surfel    s_plus = K.tail( a );
      Surfel   s_minus = K.head( a );
      SCell      linel = K.separator( a ); // oriented 1-cell
      SCell        pta = K.pivot( face_a[ 0 ] ); 
      SCell        ptb = K.pivot( face_b[ 0 ] );
      RealPoint      a = sCentroid( pta );
      RealVector     e = sCentroid( ptb ) - a;
      RealPoint     s0 = sCentroid( s_plus );
      RealPoint     s1 = sCentroid( s_minus );
      // s_plus must be to the left of e.
      if ( e.crossProduct( s0 - a ) < 0.0 ) std::swap( s_plus, s_minus );
      // Computes u_+ and u_-, then their cross product.
      RealVector     u_p = myCorrectedNormals[ s_plus ];
      RealVector     u_m = myCorrectedNormals[ s_minus ];
      RealVector  psi_e1 = u_p.crossProduct( u_m );
      return  e.dot( psi_e1 );
      // RealVector u_cross = u_p.crossProduct( u_m );
      // Scalar         psi = asin( fabs( u_cross ) );
      // RealVector      e1 = u_cross.getNormalized();
      // return  psi * e.dot( e1 );
    }
    
    // ----------------------- Interface --------------------------------------
  public:

    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const;

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const;

    // ------------------------- Private Datas --------------------------------
  private:

    /// A smart securable pointer onto the digital surface.
    CountedPtrOrPtr<Surface> theSurface;
    /// The digitization grid step.
    Scalar                   myH;
    /// The standard embedding with gridstep h.
    RegularPointEmbedder<Space> myEmbedder;
    /// The natural normal vector field.
    NormalVectorField        myTrivialNormals;
    /// The corrected normal vector field.
    NormalVectorField        myCorrectedNormals;
    
    // ------------------------- Hidden services ------------------------------
  protected:


    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class CorrectedNormalCurrent


  /**
   * Overloads 'operator<<' for displaying objects of class 'CorrectedNormalCurrent'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'CorrectedNormalCurrent' to write.
   * @return the output stream after the writing.
   */
  template <typename Surface>
  std::ostream&
  operator<< ( std::ostream & out, const CorrectedNormalCurrent<Surface> & object );


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "CorrectedNormalCurrent.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined CorrectedNormalCurrent_h

#undef CorrectedNormalCurrent_RECURSES
#endif // else defined(CorrectedNormalCurrent_RECURSES)
