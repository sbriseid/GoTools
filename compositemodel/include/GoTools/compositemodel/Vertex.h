//===========================================================================
//                                                                           
// File: Vertex.h
//                                                                           
// Created: August 25, 2008
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _VERTEX_H
#define _VERTEX_H

#include "GoTools/utils/Point.h"
#include <vector>
#include <memory>

namespace Go
{

    class ftEdge;
    class ftSurface;
    class Body;

/// \brief The vertex class represents the vertex entity in
/// a boundary represented solid or face set

class Vertex
{
 public:
  /// Constructor. Give the geometric position of the vertex
    Vertex(Point vertex_point);

    /// Constructor. Give the geometric position of the vertex and
    /// associated edges
    Vertex(Point vertex_point, std::vector<ftEdge*> edges);

    /// Constructor. Give the geometric position of the vertex and
    /// one associated edge
    Vertex(Point vertex_point, ftEdge* edges);

    /// Constructor. Give one associated edge and an indication on
    /// which of the two vertices belonging to this edge should be
    /// constructed
     Vertex(ftEdge* edge, bool at_start);

     /// Destructor
    ~Vertex();

    /// Vertices belonging to two adjacent edges are represented
    /// as one entity. Used in topology build
    void joinVertex(shared_ptr<Vertex> other);

    /// Add a new edge to this vertex. Used in topology build
    void addEdge(ftEdge* edge);

    /// Remove an edge from this vertex. Used in connetion with topology 
    /// changes in the associated model
    void removeEdge(ftEdge* edge);

    /// Given an edge connected to this vertex, remove twin information
    /// about the edge in the vertex. Used in connetion with topology 
    /// changes in the associated model
    void disconnectTwin(ftEdge* edge);

    /// Returns all edges meeting in this vertex including twins
    std::vector<ftEdge*> allEdges() const;

    /// Returns all geometrically unique edges meeting in this
    /// vertex. One edge is return for a pair of twins.
    std::vector<ftEdge*> uniqueEdges();

    std::vector<ftEdge*> uniqueEdges(Body *bd);

    /// Number of unique edges meeting in this vertex, twin edges
    /// are counted only once
    int nmbUniqueEdges()
    {
      return (int)edges_.size();
    }

    /// Get edges which are not associated a face
    std::vector<ftEdge*> freeEdges();

    /// Get the specified edge, twin edges are counted once and it is
    /// arbitrary which twin edge is returned
    ftEdge* getEdge(int idx)
    {
      return edges_[idx].first;
    }

    /// Get the geometrical position associated to this vertex
    Point getVertexPoint()
	{
	    return vertex_point_;
	}

    /// Set the geometrical position associated to this vertex
    void setVertexPoint(Point vertex_point)
    {
      vertex_point_ = vertex_point;
    }
      
    /// Get all faces meeting in this vertex
    std::vector<ftSurface*> faces() const;

    std::vector<ftSurface*> faces(Body *bd) const;

    /// Get all faces meeting in this vertex and the parameter value in
    /// the face corresponding to the vertex
    std::vector<std::pair<ftSurface*, Point> > getFaces();

    std::vector<std::pair<ftSurface*, Point> > getFaces(Body *bd);

    /// Average corners of spline surfaces corresponding to this vertex
    void averageVertexPos();

    /// Get parameter of associated face corresponding to vertex
    Point getFacePar(ftSurface* face);

    /// Get all boides meeting in this vertex
    std::vector<Body*> getBodies();

    /// The distance between this vertex an another vertex
    double getDist(shared_ptr<Vertex> other_point)
	{
	    return vertex_point_.dist(other_point->getVertexPoint());
	}

    /// Check if the vertex is connected to the given edge
    bool hasEdge(ftEdge *edge) const;

    /// Check if the vertex is connected to the given edge, and this edge
    /// is represented in the vertex with no twin
    bool hasEdgeSingle(ftEdge *edge) const;

    /// Check if two given edges meet in this vertex
    bool meetInVertex(ftEdge *e1, ftEdge *e2) const;

    /// Check if this vertex lies at a model boundary, i.e. is connected to
    /// edges with no twin
    bool isBoundaryVertex() const;

    /// Check if this vertex and the other vertex belongs to the same edge
    bool sameEdge(Vertex* other) const;


    /// Check if this vertex and the other vertex belongs to the same face
    bool sameFace(Vertex* other) const;

    /// Get the edge associated with two vertices, if any
    ftEdge* getCommonEdge(Vertex* other) const;

    ftEdge* getCommonEdgeInFace(Vertex* other,
				ftSurface* face) const;

    /// Get the faces associated with two vertices, if any
    std::vector<ftSurface*> getCommonFaces(Vertex* other) const;

    /// Get the edges meeting in this vertex associated with a given face
    std::vector<ftEdge*> getFaceEdges(ftSurface *face) const;

    /// Collect attached edges where the distance between the endpoints
    /// are larger than the specified tolerance or where the curves meet 
    /// with an angle that are more than the kink tolerance, but less than
    /// the corner tolerance
    void getEdgeDiscontinuities(std::vector<std::pair<ftEdge*, ftEdge*> >& gaps, double tol, 
				std::vector<std::pair<ftEdge*, ftEdge*> >& kinks, double angtol,
				double angmax) const;

     /// Reorganize edges according to twin information
    void reOrganize();

   /// Check the consistency of the edge information in this vertex
    // Debug functionality
    bool checkVertexTopology();

 private:
    /// The spacial position of the vertex
    Point vertex_point_;  

    /// Edges meeting in this vertex. Twin edges are collected in pairs.
    std::vector<std::pair<ftEdge*, ftEdge*> > edges_;
};

} // namespace Go


#endif // _VERTEX_H_