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
#include "DGtal/graph/DistanceBreadthFirstVisitor.h"

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
    typedef typename Surface::VertexRange             VertexRange;
    typedef typename Surface::ArcRange                ArcRange;
    typedef typename KSpace::Space                    Space;
    typedef typename Space::Point                     Point;
    typedef typename Space::Vector                    Vector;
    typedef typename Space::RealPoint                 RealPoint;
    typedef typename Space::RealVector                RealVector;
    typedef typename RealVector::Component            Scalar;
    typedef std::map< Surfel, RealVector >            NormalVectorField;
    typedef std::map< Cell,   RealPoint >             CentroidMap;
    typedef std::map< Vertex, Scalar >                Measure0Map;
    typedef std::map< Arc,    Scalar >                Measure1Map;

    // Checks that dimension is 3.
    BOOST_STATIC_ASSERT(( KSpace::dimension == 3 ));
    
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
      computeCellCentroids();
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
      computeCellCentroids();
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

    /// Computes the trivial normals of the attached surface.
    void computeCellCentroids()
    {
      const KSpace&   K = space();
      for ( auto s : *this )
	{
	  Cell  cell2 = K.unsigns( s );
	  Dimension k = K.uOrthDir( cell2 );
	  Dimension i = (k+1)%3;
	  Dimension j = (k+2)%3;
	  Cell cell1a = K.uIncident( cell2, i, false );
	  Cell cell1b = K.uIncident( cell2, i, true );
	  Cell cell1c = K.uIncident( cell2, j, false );
	  Cell cell1d = K.uIncident( cell2, j, true );
	  Cell cell0a = K.uIncident( cell1a, j, false );
	  Cell cell0b = K.uIncident( cell1a, j, true );
	  Cell cell0c = K.uIncident( cell1b, j, false );
	  Cell cell0d = K.uIncident( cell1b, j, true );
	  myCellCentroids[ cell2 ]  = computeCentroid( cell2 );
	  myCellCentroids[ cell1a ] = computeCentroid( cell1a );
	  myCellCentroids[ cell1b ] = computeCentroid( cell1b );
	  myCellCentroids[ cell1c ] = computeCentroid( cell1c );
	  myCellCentroids[ cell1d ] = computeCentroid( cell1d );
	  myCellCentroids[ cell0a ] = computeCentroid( cell0a );
	  myCellCentroids[ cell0b ] = computeCentroid( cell0b );
	  myCellCentroids[ cell0c ] = computeCentroid( cell0c );
	  myCellCentroids[ cell0d ] = computeCentroid( cell0d );
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
    RealPoint uCentroid( Cell c )
    {
      return myCellCentroids[ c ];
    }

    /// @param[in] c any signed cell.
    /// @return the centroid of the signed cell c in real space;
    RealPoint sCentroid( SCell c )
    {
      return myCellCentroids[ space().unsigns( c ) ];
    }
    
    /// @param[in] c any cell.
    /// @return the centroid of the cell c in real space;
    RealPoint computeCentroid( Cell c ) const
    {
      Point    kp = space().uKCoords( c );
      return   0.5 * myEmbedder( kp );
      // Point    lo = Point( (kp[ 0 ]-0) / 2, (kp[ 1 ]-0) / 2, (kp[ 2 ]-0) / 2 );
      // Point    hi = Point( (kp[ 0 ]+1) / 2, (kp[ 1 ]+1) / 2, (kp[ 2 ]+1) / 2 );
      // return   0.5 * ( myEmbedder( lo ) + myEmbedder( hi ) );
    }
    
    
    /// Computes the d-dimensional Hausdorff measure of a d-dimensional cell.
    Scalar Hmeasure( Dimension d ) const
    {
      Scalar H = 1.0;
      for ( Dimension k = d; k != 0; --k ) H *= myH;
      return H;
    }

    Scalar sRelativeIntersection( RealPoint p, Scalar r, SCell c )
    {
      return relativeIntersection( p, r, space().unsigns( c ) );
    }

    /// Computes an approximation of the relative intersection of the
    /// ball of radius \a r and center \a p with the given cell. Cells
    /// are embedded naturally in the grid of step h.
    ///
    /// @return the relative intersection as a scalar between 0 (no
    /// intersection) and 1 (inclusion).
    Scalar relativeIntersection( RealPoint p, Scalar r, Cell c )
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
      if      ( d_max <= r     ) return 1.0;
      else if ( r     <= d_min ) return 0.0;
      return ( r - d_min ) / ( d_max - d_min );
    }

    /// \f$ \mu_0 \f$ Lipschitz-Killing measure. It corresponds to a
    /// corrected area measure, and is non null only on 2-cells (or Vertex).
    ///
    /// @param[in] c any 2-dimensional cell (or a Vertex in 3D digital
    /// surfaces).
    ///
    /// @return the corrected area measure \f$ \mu_0 := \cos \alpha
    /// d\mathcal{H}^2 \f$.
    Scalar mu0( Vertex c )
    {
      typename Measure0Map::iterator it = myMu0.find( c );
      if ( it != myMu0.end() ) return it->second;
      it->second = computeMu0( c );
      return it->second;
    }
    
    /// \f$ \mu_0 \f$ Lipschitz-Killing measure. It corresponds to a
    /// corrected area measure, and is non null only on 2-cells (or Vertex).
    ///
    /// @param[in] c any 2-dimensional cell (or a Vertex in 3D digital
    /// surfaces).
    ///
    /// @return the corrected area measure \f$ \mu_0 := \cos \alpha
    /// d\mathcal{H}^2 \f$.
    Scalar computeMu0( Vertex c )
    {
      if ( space().sDim( c ) != 2 ) return 0.0;
      return Hmeasure( 2 ) * myTrivialNormals[ c ].dot( myCorrectedNormals[ c ] );
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
    Scalar mu1( Arc arc )
    {
      typename Measure1Map::iterator it = myMu1.find( arc );
      if ( it != myMu1.end() ) return it->second;
      it = myMu1.find( theSurface->opposite( arc ) );
      if ( it != myMu1.end() ) return it->second;
      it->second = computeMu1( arc );
      return it->second;
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
    Scalar computeMu1( Arc arc )
    {
      const KSpace & K = space();
      Surfel    s_plus = theSurface->tail( arc );
      Surfel   s_minus = theSurface->head( arc );
      Cell       linel = K.unsigns( theSurface->separator( arc ) ); // oriented 1-cell
      Dimension      l = *( K.uDirs( linel ) );
      Cell         pta = K.uIncident( linel, l, true );
      Cell         ptb = K.uIncident( linel, l, false );
      auto       faces = theSurface->facesAroundArc( arc );
      if ( faces.size() != 1 )
	faces = theSurface->facesAroundArc( theSurface->opposite( arc ) );	
      Cell       pivot = K.unsigns( theSurface->pivot( faces[ 0 ] ) );
      if ( pivot != pta ) std::swap( pta, ptb );
      RealPoint      a = uCentroid( pta );
      RealVector     e = ( uCentroid( ptb ) - a ).getNormalized();
      RealPoint     s0 = sCentroid( s_plus );
      RealPoint     s1 = sCentroid( s_minus );
      // s_plus must be to the left of e.
      // if ( e.crossProduct( s0 - a ).dot( myTrivialNormals[ s_plus ] ) < 0.0 )
      // 	   std::swap( s_plus, s_minus );
      // Computes u_+ and u_-, then their cross product.
      RealVector    u_p = myCorrectedNormals[ s_plus ];
      RealVector    u_m = myCorrectedNormals[ s_minus ];
      RealVector psi_e1 = u_p.crossProduct( u_m );
      Scalar        ne1 = psi_e1.norm();
      Scalar       npsi = (ne1 == 0.0) ? 0.0 : ( asin( ne1 ) / ne1 );
      return  Hmeasure( 1 ) * npsi * e.dot( psi_e1 );
      // RealVector u_cross = u_p.crossProduct( u_m );
      // Scalar         psi = asin( fabs( u_cross ) );
      // RealVector      e1 = u_cross.getNormalized();
      // return  psi * e.dot( e1 );
    }

    Scalar mu0( RealPoint p, Scalar r, Vertex c )
    {
      Scalar ri = sRelativeIntersection( p, r, c );
      return ri != 0.0 ? ri * mu0( c ) : 0.0;
    }

    Scalar mu1( RealPoint p, Scalar r, Arc a )
    {
      SCell linel = theSurface->separator( a ); // oriented 1-cell
      Scalar   ri = sRelativeIntersection( p, r, linel );
      return ri != 0.0 ? ri * mu1( a ) : 0.0;
    }

    struct SquaredDistance2Point {
      typedef Scalar Value;
      CorrectedNormalCurrent& current;
      RealPoint center;
      SquaredDistance2Point( CorrectedNormalCurrent& aCurrent,
			     const RealPoint&              aPoint )
	: current( aCurrent ), center( aPoint ) {}
      Value operator() ( Vertex v ) const
      {
	RealPoint  x = current.sCentroid( v );
	RealVector w = x - center;
	return w.dot( w );
      }
    };
    
    VertexRange getSurfelsInBall( SCell c, Scalar r )
    {
      typedef DistanceBreadthFirstVisitor
	< Surface, SquaredDistance2Point >     DistanceVisitor;
      typedef typename DistanceVisitor::Node   MyNode;
      typedef typename DistanceVisitor::Scalar MySize;
      
      VertexRange output;
      RealPoint   center = sCentroid( c );
      Scalar       limit = r*r;
      SquaredDistance2Point d2pfct( *this, center );
      DistanceVisitor       visitor( *theSurface, d2pfct, c );
      while ( ! visitor.finished() )
	{
	  MyNode node = visitor.current();
	  Vertex    v = node.first;
	  Scalar    d = node.second;
	  output.push_back( v );
	  if ( d <= limit ) visitor.expand();
	  else              visitor.ignore();
	}
      return output;
    }

    ArcRange getArcsInBall( SCell c, Scalar r )
    {
      const KSpace&    K = space();
      VertexRange scells = getSurfelsInBall( c, r );
      std::set<Arc> arcs;
      for ( auto s : scells )
	{
	  auto l_arcs = theSurface->outArcs( s );
	  for ( auto a : l_arcs )
	    {
	      auto b = theSurface->opposite( a );
	      if ( arcs.find( a ) == arcs.end() && arcs.find( b ) == arcs.end() )
		arcs.insert( a );
	    }
	}
      return ArcRange( arcs.begin(), arcs.end() );
    }

    Scalar mu0Ball( SCell c, Scalar r )
    {
      VertexRange vtcs = getSurfelsInBall( c, r );
      Scalar        m0 = 0.0;
      RealPoint      x = sCentroid( c );
      for ( auto v : vtcs ) {
	m0      += mu0( x, r, v ); 
      }
      return m0;
    }

    Scalar mu1Ball( SCell c, Scalar r )
    {
      ArcRange    arcs = getArcsInBall( c, r );
      Scalar        m1 = 0.0;
      RealPoint      x = sCentroid( c );
      for ( auto a : arcs ) {
	m1      += mu1( x, r, a ); 
      }
      return m1;
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
    /// The map cell -> centroid (to limit computations).
    CentroidMap              myCellCentroids;
    /// The map Vertex -> mu0
    Measure0Map              myMu0;
    /// The map Arc -> mu1
    Measure1Map              myMu1;
    
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
