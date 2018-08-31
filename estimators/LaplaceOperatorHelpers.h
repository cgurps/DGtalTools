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
 * @file LaplaceOperatorHelpers.h
 * @author Thomas Caissard (\c thomas.caissard@gmail.com )
 *
 * @date 2018/06/17
 *
 * Header file for module LaplaceOperatorHelpers.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(LaplaceOperatorHelpers_RECURSES)
#error Recursive header files inclusion detected in LaplaceOperatorHelpers.h
#else // defined(LaplaceOperatorHelpers_RECURSES)
/** Prevents recursive inclusion of headers. */
#define LaplaceOperatorHelpers_RECURSES

#if !defined LaplaceOperatorHelpers_h
/** Prevents repeated inclusion of headers. */
#define LaplaceOperatorHelpers_h

#include "EstimatorHelpers.h"

#include "DGtal/base/BasicFunctors.h"

#include <DGtal/topology/IndexedDigitalSurface.h>

#include <DGtal/shapes/PolygonalSurface.h>
#include <DGtal/shapes/MeshHelpers.h>

#include <DGtal/io/boards/Board3D.h>
#include <DGtal/io/DrawWithDisplay3DModifier.h>
#include "DGtal/io/colormaps/GradientColorMap.h"

#if defined(WITH_EIGEN)

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace DGtal
{
  namespace po  = boost::program_options;

  template <typename TScalar, typename TKSpace>
  struct LaplaceOperatorHelpers
  {
    // typedefs
    typedef TScalar Scalar;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> DenseMatrix;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> DenseVector;
    typedef Eigen::SparseMatrix<Scalar, Eigen::RowMajor, int> SparseMatrix;
    typedef Eigen::Triplet<Scalar> Triplet;

    typedef TKSpace KSpace;
    typedef typename EstimatorHelpers<KSpace>::RealPoint RealPoint;
    typedef typename EstimatorHelpers<KSpace>::Surface Surface;
    typedef typename EstimatorHelpers<KSpace>::SurfaceContainer SurfaceContainer;
    typedef typename EstimatorHelpers<KSpace>::SCell SCell;
    typedef typename EstimatorHelpers<KSpace>::ImplicitShape ImplicitShape;

    typedef typename IndexedDigitalSurface<SurfaceContainer>::Vertex Vertex;

    IndexedDigitalSurface<SurfaceContainer> myIndexedDigitalSurface;
    TriangulatedSurface<RealPoint> myTriangulatedSurface;

    std::vector<RealPoint> myNormals;

    po::variables_map myVM;

    std::map< typename Surface::Vertex, typename TriangulatedSurface<RealPoint>::Index > triangulationMap;

    void
    generateIndexedDigitalSurface( CountedPtr<Surface> surface )
    {
      trace.beginBlock("Constructing Indexed Digital Surface");
      bool fds_ok = myIndexedDigitalSurface.build( surface.get()->container() );
      trace.info() << "- Building IndexedDigitalSurface is "
        << (fds_ok ? "Success" : "Error" ) << std::endl;
      trace.endBlock();
    }

    void generateTriangulatedSurface( CountedPtr<Surface> surface )
    {
      trace.beginBlock("Converting Surface to Triangulated Surface");
      DGtal::CanonicCellEmbedder<KSpace> canonicCellEmbedder( myIndexedDigitalSurface.space() );
      MeshHelpers::digitalSurface2DualTriangulatedSurface( *surface.get(),
          canonicCellEmbedder,
          myTriangulatedSurface,
          triangulationMap );
      trace.endBlock();
    }

    template <typename SurfelIterator, typename RealPointIterator>
    void
    setNormals( SurfelIterator itS, SurfelIterator itSEnd, RealPointIterator itRV )
    {
      trace.beginBlock("Loading Normals");
      myNormals.resize( myIndexedDigitalSurface.nbVertices() );
      for ( ; itS != itSEnd; ++itS )
      {
        const Vertex v = myIndexedDigitalSurface.getVertex( *itS );
        if ( v != myIndexedDigitalSurface.INVALID_FACE )
          myNormals[ v ] = *itRV++;
        else
        trace.warning() << "[LaplaceOperatorHelpers::setNormals]"
          << " Surfel " << *itS << " is not in the surface." << std::endl;
      }
      trace.endBlock();
    }

    void setVariableMap( po::variables_map& vm ) { myVM = vm; }

    DenseVector
    computeHeatCurvatureFromFunctor() const
    {
      DenseVector curvature( myIndexedDigitalSurface.nbVertices() );
      #pragma omp parallel for
      for( int v = 0; v < myIndexedDigitalSurface.nbVertices(); ++v )
        curvature( v ) = 0.5 * heatLaplaceCurvatureFunctor( v ).norm();

      return curvature;
    }

    DenseVector
    computeHeatCurvatureFromTriangulationFunctor() const
    {
      const Scalar h  = myVM["gridstep"].as<Scalar>();
      const Scalar h_mesh  = computeShapeRegularity( h );

      DenseVector curvature( myTriangulatedSurface.nbVertices() );
      DenseVector areas = computeAreaWeights();
      #pragma omp parallel for
      for( int v = 0; v < myTriangulatedSurface.nbVertices(); ++v )
        curvature( v ) = 0.5 * heatLaplaceCurvatureFunctorTriangulation( v, areas, h_mesh ).norm();

      return curvature;
    }

    RealPoint
    heatLaplaceCurvatureFunctor( const typename IndexedDigitalSurface<SurfaceContainer>::Index& vertexIndex ) const
    {
      typedef typename IndexedDigitalSurface<SurfaceContainer>::Index Index;

      const Scalar h  = myVM["gridstep"].as<Scalar>();
      const Scalar mcoef = myVM["m-coef"].as<Scalar>();
      const Scalar mpow  = myVM["m-pow"].as<Scalar>();
      const Scalar t = mcoef * std::pow( h, mpow );
      const Scalar cut_locus = ( log( - log( t ) + 2. ) + 2. ) * sqrt( 2. * t );

      std::vector<unsigned char> visited( myIndexedDigitalSurface.nbVertices(), false );
      visited[ vertexIndex ] = true;

      std::queue<Index> toVisit;
      toVisit.push( vertexIndex );

      const RealPoint p = h * myIndexedDigitalSurface.position( vertexIndex );

      RealPoint accum( 0., 0. ,0. );
      while( ! toVisit.empty() )
      {
        const Index currentIndex = toVisit.front();
        toVisit.pop();

        std::vector<Index> neigh;
        std::back_insert_iterator<std::vector<Index>> backIter( neigh );
        myIndexedDigitalSurface.writeNeighbors( backIter, currentIndex );

        for( const auto& v : neigh )
        {
          if( ! visited[v] && ( p - h * myIndexedDigitalSurface.position( v ) ).norm() < cut_locus )
          {
            visited[v] = true;
            toVisit.push(v);
          }
        }

        const Scalar l2_distance = ( h * myIndexedDigitalSurface.position( currentIndex ) - p ).norm();
        const SCell sCell = myIndexedDigitalSurface.surfel( currentIndex );
        const Scalar measure = std::pow( h, 2 ) * std::abs( myNormals[ currentIndex ][ myIndexedDigitalSurface.space().sOrthDir( sCell ) ] );

        accum += ( measure / ( 4. * M_PI * t * t ) )
          * std::exp( - l2_distance * l2_distance / ( 4. * t ) )
          * ( h * myIndexedDigitalSurface.position( currentIndex ) - p );
      }

      return accum;
    }

    RealPoint
    heatLaplaceCurvatureFunctorTriangulation( const typename TriangulatedSurface<RealPoint>::Index& vertexIndex, const DenseVector& areas, const Scalar& h_mesh ) const
    {
      typedef typename TriangulatedSurface<RealPoint>::Index Index;

      const Scalar h  = myVM["gridstep"].as<Scalar>();
      const Scalar mcoef = myVM["m-coef"].as<Scalar>();
      const Scalar mpow  = myVM["m-pow"].as<Scalar>();
      const Scalar t = mcoef * std::pow( h_mesh, mpow );
      const Scalar cut_locus = ( log( - log( t ) + 2. ) + 2. ) * sqrt( 2. * t );

      std::vector<unsigned char> visited( myTriangulatedSurface.nbVertices(), false );
      visited[ vertexIndex ] = true;

      std::queue<Index> toVisit;
      toVisit.push( vertexIndex );

      const RealPoint p = h * myTriangulatedSurface.position( vertexIndex );

      RealPoint accum( 0., 0. ,0. );
      while( ! toVisit.empty() )
      {
        const Index currentIndex = toVisit.front();
        toVisit.pop();

        std::vector<Index> neigh;
        std::back_insert_iterator<std::vector<Index>> backIter( neigh );
        myTriangulatedSurface.writeNeighbors( backIter, currentIndex );

        for( const auto& v : neigh )
        {
          if( ! visited[v] && ( p - h * myTriangulatedSurface.position( v ) ).norm() < cut_locus )
          {
            visited[v] = true;
            toVisit.push(v);
          }
        }

        const Scalar l2_distance = ( h * myTriangulatedSurface.position( currentIndex ) - p ).norm();
        const Scalar measure = 1. / areas( currentIndex );

        accum += ( measure / ( 4. * M_PI * t * t ) )
          * std::exp( - l2_distance * l2_distance / ( 4. * t ) )
          * ( h * myTriangulatedSurface.position( currentIndex ) - p );
      }

      return accum;
    }

    SparseMatrix
    computeHeatLaplace() const
    {
      trace.beginBlock("Computing Heat Laplace Operator Matrix");

      typedef typename IndexedDigitalSurface<SurfaceContainer>::Index Index;

      const Scalar h  = myVM["gridstep"].as<Scalar>();
      const Scalar mcoef = myVM["m-coef"].as<Scalar>();
      const Scalar mpow  = myVM["m-pow"].as<Scalar>();
      const Scalar t = mcoef * std::pow( h, mpow );
      const Scalar cut_locus = 3.5 * sqrt( 2. * t );

      trace.info() << "t=" << t << std::endl;
      trace.info() << "cut_locus=" << cut_locus << std::endl;

      std::vector<Triplet> laplaceTriplets;

      #pragma omp parallel for
      for( int v_i = 0; v_i < myIndexedDigitalSurface.nbVertices(); ++v_i )
      {
        std::vector<Triplet> localTriplets;
        const RealPoint p_i = h * myIndexedDigitalSurface.position( v_i );
        Scalar accum = 0.;

        for( int v_j = 0; v_j < myIndexedDigitalSurface.nbVertices(); ++v_j )
        {
          if( v_i == v_j ) continue;

          const RealPoint p_j = h * myIndexedDigitalSurface.position( v_j );
          const Scalar l2_distance = ( p_i - p_j ).norm();

          if( l2_distance < cut_locus )
          {
            const SCell surfel = myIndexedDigitalSurface.surfel( v_j );
            const Scalar measure = std::pow( h, 2. ) * std::abs( myNormals[ v_j ][ myIndexedDigitalSurface.space().sOrthDir( surfel ) ] );

            const Scalar laplace_value = ( measure / ( 4. * M_PI * t * t ) ) * std::exp( - l2_distance * l2_distance / ( 4. * t ) );

            localTriplets.push_back( Triplet( v_i, v_j, laplace_value ) );
            accum -= laplace_value;
          }
        }

        localTriplets.push_back( Triplet( v_i, v_i, accum ) );

        {
          #pragma omp critical
          laplaceTriplets.insert( laplaceTriplets.end(), localTriplets.begin(), localTriplets.end() );
        }
      }

      SparseMatrix op( myIndexedDigitalSurface.nbVertices(), myIndexedDigitalSurface.nbVertices() );
      op.setFromTriplets( laplaceTriplets.begin(), laplaceTriplets.end() );

      trace.endBlock();

      return op;
    }

    SparseMatrix
    computeHeatLaplaceTriangulation() const
    {
      trace.beginBlock("Computing Heat Laplace Operator Matrix");

      typedef typename TriangulatedSurface<RealPoint>::Index Index;

      const Scalar h  = myVM["gridstep"].as<Scalar>();
      const Scalar h_mesh  = computeShapeRegularity( h );
      const Scalar mcoef = myVM["m-coef"].as<Scalar>();
      const Scalar mpow  = myVM["m-pow"].as<Scalar>();
      const Scalar t = mcoef * std::pow( h_mesh, mpow );
      const Scalar cut_locus = 2.5 * std::sqrt( 2. * t );

      trace.info() << "t=" << t << std::endl;
      trace.info() << "cut_locus=" << cut_locus << std::endl;

      std::vector<Triplet> laplaceTriplets;

      DenseVector areaWeights = computeAreaWeights();
      #pragma omp parallel for
      for( int v_i = 0; v_i < myTriangulatedSurface.nbVertices(); ++v_i )
      {
        std::vector<Triplet> localTriplets;
        const RealPoint p_i = h * myTriangulatedSurface.position( v_i );
        Scalar accum = 0.;

        for( int v_j = 0; v_j < myTriangulatedSurface.nbVertices(); ++v_j )
        {
          if( v_i == v_j ) continue;

          const RealPoint p_j = h * myTriangulatedSurface.position( v_j );
          const Scalar l2_distance = ( p_i - p_j ).norm();

          if( l2_distance < cut_locus )
          {
            const Scalar measure = 1. / areaWeights( v_j );
            const Scalar laplace_value = ( measure / ( 4. * M_PI * t * t ) ) * std::exp( - l2_distance * l2_distance / ( 4. * t ) );

            localTriplets.push_back( Triplet( v_i, v_j, laplace_value ) );
            accum -= laplace_value;
          }
        }

        localTriplets.push_back( Triplet( v_i, v_i, accum ) );

        {
          #pragma omp critical
          laplaceTriplets.insert( laplaceTriplets.end(), localTriplets.begin(), localTriplets.end() );
        }
      }

      SparseMatrix op( myTriangulatedSurface.nbVertices(), myTriangulatedSurface.nbVertices() );
      op.setFromTriplets( laplaceTriplets.begin(), laplaceTriplets.end() );

      trace.endBlock();

      return op;
    }

    SparseMatrix
    computeCombinatorialLaplace() const
    {
      trace.beginBlock("Computing Combinatorial Laplace Matrix");

      std::vector<Triplet> laplaceTriplets;
      for( auto const& v : myIndexedDigitalSurface.allVertices() )
      {
        laplaceTriplets.push_back( Triplet( v, v, - 1. ) );

        const auto outArcs = myIndexedDigitalSurface.outArcs( v );
        for( auto const& a : outArcs )
          laplaceTriplets.push_back( Triplet( v, myIndexedDigitalSurface.head( a ), 1. / outArcs.size() ) );
      }

      SparseMatrix op( myIndexedDigitalSurface.nbVertices(), myIndexedDigitalSurface.nbVertices() );
      op.setFromTriplets( laplaceTriplets.begin(), laplaceTriplets.end() );

      trace.endBlock();

      return op;
    }

    SparseMatrix
    computeCotanLaplace () const
    {
      trace.beginBlock("Computing Cotan Laplace Matrix");

      const Scalar h = myVM["gridstep"].as<Scalar>();
      std::vector<Triplet> cotanTriplets;
      const auto cotan
        = []( const Scalar& angle ) { return cos( angle ) / sin( angle ); };
      const auto vertex_weight
        = [&]( const typename TriangulatedSurface<RealPoint>::Arc& a )
        {
          const RealPoint p_i = h * myTriangulatedSurface.position( myTriangulatedSurface.head( a ) );
          const RealPoint p_j = h * myTriangulatedSurface.position( myTriangulatedSurface.tail( a ) );
          const RealPoint p_k = h * myTriangulatedSurface.position( myTriangulatedSurface.head( myTriangulatedSurface.next( a ) ) );
          const RealPoint p_l = h * myTriangulatedSurface.position( myTriangulatedSurface.head( myTriangulatedSurface.next( myTriangulatedSurface.opposite( a ) ) ) );

          const Scalar alpha = acos( ( ( p_i - p_k ) / ( p_i - p_k ).norm() ).dot( ( p_j - p_k ) / ( p_j - p_k ).norm() ) );
          const Scalar beta  = acos( ( ( p_i - p_l ) / ( p_i - p_l ).norm() ).dot( ( p_j - p_l ) / ( p_j - p_l ).norm() ) );

          return 0.5 * ( cotan( alpha ) + cotan( beta ) );
        };

      for( auto const& v : myTriangulatedSurface.allVertices() )
      {
        for( auto const& a : myTriangulatedSurface.outArcs( v ) )
        {
          const Scalar laplaceValue = vertex_weight( a );

          if( isnan( laplaceValue ) ) trace.error() << "[computeCotanLaplace] NAN HERE " << v << " -> " << a << std::endl;

          cotanTriplets.push_back( Triplet( v, myTriangulatedSurface.head( a ), laplaceValue ) );
          cotanTriplets.push_back( Triplet( v, v, - laplaceValue ) );
        }
      }

      SparseMatrix op( myTriangulatedSurface.nbVertices(), myTriangulatedSurface.nbVertices() );
      op.setFromTriplets( cotanTriplets.begin(), cotanTriplets.end() );
      DenseVector areaWeights = computeAreaWeights();

      trace.endBlock();

      return areaWeights.asDiagonal() * op;
    }

    SparseMatrix
    computeRLocalLaplace() const
    {
      trace.beginBlock("Computing RLocal Laplace Matrix");

      typedef typename TriangulatedSurface<RealPoint>::Index Index;

      const Scalar h = myVM["gridstep"].as<Scalar>();
      const Scalar h_mesh  = computeShapeRegularity( h );
      const Scalar mcoef = myVM["m-coef"].as<Scalar>();
      const Scalar mpow  = myVM["m-pow"].as<Scalar>();
      const Scalar r = mcoef * std::pow( h_mesh, mpow );

      trace.info() << "h_mesh=" << h_mesh << std::endl;
      trace.info() << "h=" << h << std::endl;
      trace.info() << "r=" << r << std::endl;

      std::vector<Triplet> rTriplets;

      #pragma omp parallel for
      for( int v = 0; v < myTriangulatedSurface.nbVertices(); ++v )
      {
        DenseVector values( myTriangulatedSurface.nbVertices() );
        const RealPoint p_v = h * myTriangulatedSurface.position( v );

        std::vector<unsigned char> visited( myTriangulatedSurface.nbVertices(), false );
        std::queue<Index> pointToVisit;
        pointToVisit.push( v );
        visited[ v ] = true;

        while( ! pointToVisit.empty() )
        {
          const auto v_current = pointToVisit.front();
          pointToVisit.pop();
          std::vector<Index> neigh;
          neigh.reserve( 8 );
          std::back_insert_iterator<std::vector<Index>> back_iter( neigh );
          myTriangulatedSurface.writeNeighbors( back_iter, v_current );

          for( auto const& vv : neigh )
          {
            if( ( p_v - h * myTriangulatedSurface.position( vv ) ).norm() < r && ! visited[ vv ] )
            {
              pointToVisit.push( vv );
              visited[vv] = true;
            }
          }

          const Scalar l2_distance = ( p_v - h * myTriangulatedSurface.position( v_current ) ).norm();
          values( v_current ) = std::max( 0., 1. - ( l2_distance / r ) );
        }

        values /= values.template lpNorm<1>();

        std::fill( visited.begin(), visited.end(), false );
        pointToVisit.push( v );
        visited[ v ] = true;

        std::vector<Triplet> localTriplets;
        while( ! pointToVisit.empty() )
        {
          const auto v_current = pointToVisit.front();
          pointToVisit.pop();
          std::vector<Index> neigh;
          neigh.reserve( 8 );
          std::back_insert_iterator<std::vector<Index>> back_iter( neigh );
          myTriangulatedSurface.writeNeighbors( back_iter, v_current );

          for( auto const& vv : neigh )
          {
            if( ( p_v - h * myTriangulatedSurface.position( vv ) ).norm() < r && ! visited[ vv ] )
            {
              pointToVisit.push( vv );
              visited[vv] = true;
            }
          }

          localTriplets.push_back( Triplet( v, v_current, values( v_current ) ) );
        }

        {
          #pragma omp critical
          rTriplets.insert( rTriplets.end(), localTriplets.begin(), localTriplets.end() );
        }
      }

      SparseMatrix op( myTriangulatedSurface.nbVertices(), myTriangulatedSurface.nbVertices() );
      op.setFromTriplets( rTriplets.begin(), rTriplets.end() );

      op = op * computeCotanLaplace();
      trace.endBlock();

      return op;
    }

    DenseMatrix
    computePrimalEmbedding() const
    {
      trace.beginBlock("Computing Primal Embedding");

      const Scalar h = myVM["gridstep"].as<Scalar>();
      DenseMatrix primalEmbedding( myIndexedDigitalSurface.nbVertices(), 3 );

      //TODO

      trace.endBlock();

      return primalEmbedding;
    }

    DenseMatrix
    computeDualDigitalEmbedding() const
    {
      trace.beginBlock("Computing Dual Embedding");

      const Scalar h = myVM["gridstep"].as<Scalar>();
      DenseMatrix dualEmbedding( myIndexedDigitalSurface.nbVertices(), 3 );

      for( int v = 0; v < myIndexedDigitalSurface.nbVertices(); ++v )
        for( int dim = 0; dim < 3; ++dim )
          dualEmbedding( v, dim ) = h * myIndexedDigitalSurface.position( v )[ dim ];

      trace.endBlock();

      return dualEmbedding;
    }

    DenseMatrix
    computeDualTriangulationEmbedding() const
    {
      trace.beginBlock("Computing Dual Embedding");

      const Scalar h = myVM["gridstep"].as<Scalar>();
      DenseMatrix dualEmbedding( myTriangulatedSurface.nbVertices(), 3 );

      for( auto const& v : myTriangulatedSurface.allVertices() )
        for( int dim = 0; dim < 3; ++dim )
          dualEmbedding( v, dim ) = h * myTriangulatedSurface.position( v )[ dim ];

      trace.endBlock();

      return dualEmbedding;
    }

    void
    exportDualQuantityToOBJ( DenseVector& quantity, std::string filename )
    {
      trace.beginBlock("Exporting Dual Quantity to OBJ");

      trace.info() << "Filename is " << filename << std::endl;
      trace.info() << "Min / Max : " << quantity.minCoeff() << " / " << quantity.maxCoeff() << std::endl;

      DGtal::Board3D<> board;
      DGtal::GradientColorMap<Scalar> gradient( quantity.minCoeff(), quantity.maxCoeff() );
      gradient.addColor( DGtal::Color( 50 , 50 , 255 ) );
      gradient.addColor( DGtal::Color( 255, 0  , 0   ) );
      gradient.addColor( DGtal::Color( 255, 255, 10  ) );

      for( auto const& v : myIndexedDigitalSurface.allVertices() )
      {
        const typename Surface::SCell sCell = myIndexedDigitalSurface.surfel( v );
        board << SetMode3D( sCell.className(), "Basic" );
        board << CustomColors3D( gradient( quantity( v ) ), gradient( quantity( v ) ) );
        board << sCell;
      }

      board.saveOBJ( filename, true );

      trace.endBlock();
    }

    void
    exportTriangulationQuantityToOBJ( DenseVector& quantity, std::string filename )
    {
      trace.beginBlock("Exporting Dual Quantity to OBJ");

      trace.info() << "Filename is " << filename << std::endl;
      trace.info() << "Min / Max : " << quantity.minCoeff() << " / " << quantity.maxCoeff() << std::endl;

      DGtal::Board3D<> board;

      DGtal::GradientColorMap<Scalar> gradient( quantity.minCoeff(), quantity.maxCoeff() );
      gradient.addColor( DGtal::Color( 50 , 50 , 255 ) );
      gradient.addColor( DGtal::Color( 255, 0  , 0   ) );
      gradient.addColor( DGtal::Color( 255, 255, 10  ) );

      for( const auto& v : myIndexedDigitalSurface.allVertices() )
      {
        const typename Surface::SCell sCell = myIndexedDigitalSurface.surfel( v );
        const auto tIndex = triangulationMap.find( sCell )->second;
        board << SetMode3D( sCell.className(), "Basic" );
        board << CustomColors3D( gradient( quantity( tIndex ) ), gradient( quantity( tIndex ) ) );
        board << sCell;
      }

      board.saveOBJ( filename, true );

      trace.endBlock();
    }

    template <typename TIterator>
    void
    exportDualQuantityToOBJ( DenseVector& quantity, std::string filename, TIterator iter )
    {
      trace.beginBlock("Exporting Dual Quantity to OBJ");

      trace.info() << "Filename is " << filename << std::endl;
      trace.info() << "Min / Max : " << quantity.minCoeff() << " / " << quantity.maxCoeff() << std::endl;

      DGtal::Board3D<> board;
      DGtal::GradientColorMap<Scalar> gradient( quantity.minCoeff(), quantity.maxCoeff() );
      gradient.addColor( DGtal::Color( 50 , 50 , 255 ) );
      gradient.addColor( DGtal::Color( 255, 0  , 0   ) );
      gradient.addColor( DGtal::Color( 255, 255, 10  ) );

      for( int i = 0; i < quantity.size(); ++i )
      {
        board << SetMode3D( iter->className(), "Basic" );
        board << CustomColors3D( gradient( quantity( i ) ), gradient( quantity( i ) ) );
        board << *iter++;
      }

      board.saveOBJ( filename, true );

      trace.endBlock();
    }

    void
    exportErrors( DenseVector& groundTruth, DenseVector& estimation, std::string filename )
    {
      const Scalar h = myVM["gridstep"].as<Scalar>();
      trace.info() << filename << " " << h << " " << ( groundTruth - estimation ).cwiseAbs().mean() << " " << ( groundTruth - estimation ).cwiseAbs().maxCoeff() << std::endl;

      std::ofstream errors( filename, std::ofstream::app );
      errors << h << " " << ( groundTruth - estimation ).cwiseAbs().mean() << " " << ( groundTruth - estimation ).cwiseAbs().maxCoeff() << std::endl;
    }

    template <typename TIterator>
    DenseVector
    reorder( DenseVector& quantity, TIterator begin, bool triangulation = false )
    {
      DenseVector reorderedVector( quantity.size() );
      for( int i = 0; i < quantity.size(); ++i )
      {
        const auto cellIndex = ( triangulation ) ?
          triangulationMap.find( *begin++ )->second :
          myIndexedDigitalSurface.getVertex( *begin++ );

        reorderedVector( cellIndex ) = quantity( i );
      }

      return reorderedVector;
    }

    void
    signCurvature( DenseVector& input, DenseVector& output )
    {
      for( int i = 0; i < input.size(); ++i )
        output( i ) *= ( output(i) * input( i ) < 0 ) ? -1. : 1.;
    }

    Scalar
    computeShapeRegularity( const Scalar& h ) const
    {
      Scalar sr = 0.;

      for( auto const& f : myTriangulatedSurface.allFaces() )
      {
        const auto vr = myTriangulatedSurface.verticesAroundFace( f );
        ASSERT( vr.size() == 3 );

        const RealPoint p = h * myTriangulatedSurface.position( vr[0] );
        const RealPoint q = h * myTriangulatedSurface.position( vr[1] );
        const RealPoint r = h * myTriangulatedSurface.position( vr[2] );

        const Scalar a = ( p - q ).norm();
        const Scalar b = ( p - r ).norm();
        const Scalar c = ( r - q ).norm();

        const Scalar s = 0.5 * ( a + b + c );
        const Scalar circumradius = a * b * c / ( 4. * std::sqrt( s * ( a + b - s ) * ( a + c - s ) * ( b + c - s ) ) );

        if( circumradius > sr ) sr = circumradius;
      }

      return sr;
    }

    DenseVector
    computeRealCurvatureTriangulation( CountedPtr<ImplicitShape> shape )
    {
      DenseVector realCurv( myTriangulatedSurface.nbVertices() );
      const Scalar h = myVM["gridstep"].as<Scalar>();

      for( int v = 0; v < myTriangulatedSurface.nbVertices(); ++v )
        realCurv( v ) = shape->meanCurvature( shape->nearestPoint( h * myTriangulatedSurface.position( v ) ) );

      return realCurv;
    }

    DenseVector
    computeAreaWeights() const
    {
      DenseVector areaWeights = DenseVector::Zero( myTriangulatedSurface.nbVertices() );
      for( auto const& v : myTriangulatedSurface.allVertices() )
      {
        const auto facesAround = myTriangulatedSurface.facesAroundVertex( v );
        for( auto const& f : facesAround )
          areaWeights( v ) += triangleArea( f ) / 3.;

        if( isnan( areaWeights( v ) ) ) trace.error() << "[computeAreaWeights] NAN HERE " << v << std::endl;
      }

      return areaWeights.cwiseInverse();
    }

    Scalar
    triangleArea( const typename TriangulatedSurface<RealPoint>::Face& f ) const
    {
      const Scalar h = myVM["gridstep"].as<Scalar>();

      const auto verticesAround = myTriangulatedSurface.verticesAroundFace( f );
      const RealPoint p = h * myTriangulatedSurface.position( verticesAround[0] );
      const RealPoint q = h * myTriangulatedSurface.position( verticesAround[1] );
      const RealPoint r = h * myTriangulatedSurface.position( verticesAround[2] );

      const Scalar a = ( p - q ).norm();
      const Scalar b = ( q - r ).norm();
      const Scalar c = ( r - q ).norm();

      const Scalar s = Scalar( 0.5 ) * ( a + b + c );

      return sqrt( s * ( s -a ) * ( s - b ) * ( s - c ) );
    }

  }; // End Struct
} // namespace DGtal

#endif // EIGEN

#endif // !defined LaplaceOperatorHelpers_h

#undef LaplaceOperatorHelpers_RECURSES
#endif // else defined(LaplaceOperatorHelpers_RECURSES)
