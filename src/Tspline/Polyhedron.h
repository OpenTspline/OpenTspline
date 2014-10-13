/*
  Copyright (c) <2014> <Thomas Mörwald, Vienna University of Technology>

  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted (subject to the limitations in the disclaimer
  below) provided that the following conditions are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.

   * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.

   * Neither the name of <Owner Organization> nor the names of its
     contributors may be used to endorse or promote products derived from this
     software without specific prior written permission.

  NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY THIS
  LICENSE.  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef _TSPLINE_POLYHEDRON_H_
#define _TSPLINE_POLYHEDRON_H_

// This file defines vertex, halfedge and face for CGAL::Polyhedron.
// These classes may not be confused with the Tspline vertex, halfedge and face.
// It implements functions related to domain creation embedded on the polyhedron manifold.
// It implements shortest path computation using boost::dijkstra_shortest_paths,
// gaussian curvature.

#include "Math.hpp"
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/Parameterization_polyhedron_adaptor_3.h>
#include <CGAL/Parameterization_mesh_patch_3.h>
#include <boost/graph/dijkstra_shortest_paths.hpp>

namespace tspline {

typedef CGAL::Tag_true Tag_true;
typedef CGAL::Tag_false Tag_false;

class HalfedgeVector;

// **************************************************
// Polyhedron VERTEX
template < class Refs, class P, class ID>
class HalfedgeDS_vertex : public CGAL::HalfedgeDS_vertex_base< Refs, Tag_true, P>
{
public:
  typedef CGAL::HalfedgeDS_vertex_base< Refs, Tag_true, P> Base ;
  typedef ID size_type ;
  typedef P Point ;

private:
  size_type mID;
  double gauss_curvature;
  std::vector<HalfedgeVector*> incident_halfedges;  /// customized incident halfedges (domain specific halfedges)

public:
  HalfedgeDS_vertex() : mID ( size_type(-1) )  {}
  HalfedgeDS_vertex( Point const& p) : Base(p), mID ( size_type(-1) ) {}
  HalfedgeDS_vertex( Point const& p, size_type i ) : Base(p), mID(i) {}

  size_type&       id()       { return mID; }
  size_type const& id() const { return mID; }

  void set_gauss_curvature(const double& k)
  {
    gauss_curvature = k;
  }
  double get_gauss_curvature() const
  {
    return gauss_curvature;
  }

  void add_incident_vertex(HalfedgeVector* hv)
  {
    incident_halfedges.push_back(hv);
  }
  const std::vector<HalfedgeVector*>& get_incident_halfedges()
  {
    return incident_halfedges;
  }
};

// **************************************************
// Polyhedron HALFEDGE
template < class Refs, class ID>
class HalfedgeDS_halfedge : public CGAL::HalfedgeDS_halfedge_base< Refs, Tag_true, Tag_true, Tag_true >
{
public:
  typedef CGAL::HalfedgeDS_halfedge_base< Refs, Tag_true, Tag_true, Tag_true> Base ;
  typedef typename Base::Base_base Base_base ;
  typedef ID size_type ;


private:
  size_type mID ;
  double m_weight;  /// weights used for dijkstra shortest path computation
  int m_domID;      /// ID of domain, if this halfedge is part of the border of that domain

public:
  HalfedgeDS_halfedge( size_type i = size_type(-1) ) : mID(i), m_domID(-1) {}

  size_type&       id()       { return mID; }
  size_type const& id() const { return mID; }

  void weight(const double& w) { m_weight=w; }
  double weight() const { return m_weight; }

  void domID(int v) { m_domID=v; }
  int domID() const { return m_domID; }
};

// **************************************************
// Polyhedron FACE
template < class Refs, class Pln, class ID>
class HalfedgeDS_face : public CGAL::HalfedgeDS_face_base< Refs, Tag_true, Pln>
{
public:
  typedef CGAL::HalfedgeDS_face_base< Refs, Tag_true, Pln> Base ;
  typedef ID size_type ;

private:
  size_type mID ;
  int m_domID;

public:
  HalfedgeDS_face() : mID ( size_type(-1) ), m_domID(-1) {}
  HalfedgeDS_face( Pln const& p) : Base(p), mID ( size_type(-1) ), m_domID(-1)  {}
  HalfedgeDS_face( Pln const& p, size_type i ) : Base(p), mID (i), m_domID(-1) {}

  size_type&       id()       { return mID; }
  size_type const& id() const { return mID; }

  void domID(int v){ m_domID=v; }
  int domID() const { return m_domID; }
};

// **************************************************
// POLYHEDRON
/** @brief  customized CGAL::Polyhedron_3 class
 *          for operations on polygonal meshes (like quad-domain creation) */
class Polyhedron_with_user_data {
public:
  template < class Refs, class Traits>
  struct Vertex_wrapper {
    typedef typename Traits::Point_3 Point;
    typedef HalfedgeDS_vertex< Refs, Point, std::size_t> Vertex;
  };
  template < class Refs, class Traits>
  struct Halfedge_wrapper {
    typedef HalfedgeDS_halfedge<Refs, std::size_t> Halfedge;
  };
  template < class Refs, class Traits>
  struct Face_wrapper {
    typedef HalfedgeDS_face< Refs, Tag_false, std::size_t>  Face;
  };


};
typedef CGAL::Polyhedron_3<Kernel, Polyhedron_with_user_data> Polyhedron;

// **************************************************
// Polyhedron types for boost::dijkstra
typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator vertex_iterator;
typedef boost::graph_traits<Polyhedron>::edge_descriptor edge_descriptor;
typedef boost::graph_traits<Polyhedron>::edge_iterator edge_iterator;
typedef boost::property_map<Polyhedron, boost::vertex_index_t>::type VertexIdPMap;
typedef boost::property_map<Polyhedron, boost::edge_weight_t>::type WeightPMap;
typedef boost::associative_property_map< std::map<edge_descriptor, double> > WeightAPMap;
typedef boost::iterator_property_map<std::vector<vertex_descriptor>::iterator, VertexIdPMap> PredecessorPMap;
typedef boost::iterator_property_map<std::vector<double>::iterator, VertexIdPMap> DistancePMap;

// **************************************************
// Polyhedron types for parametrization
typedef CGAL::Parameterization_polyhedron_adaptor_3<Polyhedron> ParamPolyAdaptor;
typedef CGAL::Parameterization_mesh_patch_3<ParamPolyAdaptor> ParamMeshPatch;
typedef std::list<ParamPolyAdaptor::Vertex_handle> Seam;
typedef std::list<Point2d> SeamUV;

// **************************************************
// Domain border
/** @brief  Halfedge structure for Domain, defining a full domain border (one side)
 *          derived from a vector of Polyhedron::Halfedge
 *          traversal of domain boundaries
 *          insertion and removal of halfeges from domain boundary*/
class HalfedgeVector : public std::vector<Polyhedron::Halfedge_iterator>
{
private:
  int m_domID;
  HalfedgeVector* m_opposite;
  HalfedgeVector* m_next;
  HalfedgeVector* m_prev;
  std::vector<vertex_descriptor> m_intermediate;

protected:
  void opposite(HalfedgeVector* v) { m_opposite = v; }
  void next(HalfedgeVector* v) { m_next = v; }
  void prev(HalfedgeVector* v) { m_prev = v; }
  void intermediate(vertex_descriptor v) { m_intermediate.push_back(v); }

  bool has_opposite() const { return m_opposite != NULL; }

  bool is_opposite(const HalfedgeVector& v) const
  {
    Polyhedron::Vertex_const_iterator v1a = front()->opposite()->vertex();
    Polyhedron::Vertex_const_iterator v1b = back()->vertex();
    Polyhedron::Vertex_const_iterator v2a = v.front()->opposite()->vertex();
    Polyhedron::Vertex_const_iterator v2b = v.back()->vertex();

    if(v1a==v2b && v1b==v2a)
      return true;

    return false;
  }

  void append(HalfedgeVector& v)
  {
    v.set_halfedges_domID();
    v.set_halfedges_weights(DBL_MAX);
    v.set_halfedges_opposite_weights(DBL_MAX);
    insert(this->end(),v.begin(),v.end());
  }

  void remove()
  {
    for(size_t i=0; i<size(); i++)
    {
      Polyhedron::Halfedge_iterator h = at(i);
      h->domID(-1);
      if(h->opposite()->domID()==-1)
      {
        double w = sqrt( (h->vertex()->point() - h->opposite()->vertex()->point()).squared_length() );
        h->weight(w);
        h->opposite()->weight(w);
      }
    }
  }

  void set_halfedges_domID() { for(size_t i=0; i<size(); i++) at(i)->domID(m_domID); }
  void set_halfedges_weights(double v) { for(size_t i=0; i<size(); i++) at(i)->weight(v); }
  void set_halfedges_opposite_weights(double v) { for(size_t i=0; i<size(); i++) at(i)->opposite()->weight(v); }

  HalfedgeVector(int domID) : std::vector<Polyhedron::Halfedge_iterator>(),
    m_domID(domID), m_opposite(NULL), m_next(NULL), m_prev(NULL)
  {}

  friend class Domain;

public:
  HalfedgeVector* opposite() { return m_opposite; }
  HalfedgeVector* next() { return m_next; }
  HalfedgeVector* prev() { return m_prev; }
  const std::vector<vertex_descriptor>& intermediate() const { return m_intermediate; }

};

// *************************
// Polyhedron functions
Polyhedron::Halfedge_iterator GetHalfedge(Polyhedron::Vertex_iterator a, Polyhedron::Vertex_iterator b);

/** @brief compute gaussian curvature of a vertex */
void ComputeGaussCurvature(Polyhedron::Vertex_iterator v1);

/** @brief compute gaussian curvature of all vertices of the polyhedron */
void ComputeGaussCurvature(Polyhedron &P);

/** @brief implements dijkstra shortest path from one vertex to another
 *  @param p in: polyhedron
 *  @param a in: start vertex
 *  @param b in: end vertex
 *  @param hv out: path given as halfedge sequence (HalfedgeVector)
*/
void PolyhedronShortestPath(Polyhedron& P,
                            vertex_descriptor &a, vertex_descriptor &b,
                            HalfedgeVector &hv);

/** @param set weight of halfedges according to euclidean distance of adjascent vertices */
void SetEdgeWeightsEuclidean(Polyhedron& P);

} // namespace tspline

#endif
