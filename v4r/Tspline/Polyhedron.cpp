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

#include "Polyhedron.h"

namespace tspline{

Polyhedron::Halfedge_iterator GetHalfedge(Polyhedron::Vertex_iterator a, Polyhedron::Vertex_iterator b)
{
  Polyhedron::Halfedge_around_vertex_circulator first, curr;
  first = curr = b->vertex_begin();
  do
  {
    if(a == curr->opposite()->vertex())
      return curr;
  }
  while(++curr!=first);

  throw std::runtime_error("[Polyhedron::GetHalfedge] Error, a is not connected to b.");
}

void ComputeGaussCurvature(Polyhedron::Vertex_iterator v1)
{
  double theta(0.0);
  Point3d p1 = v1->point();
  Point3d p2,p3;
  Vector3d v12, v13;
  Polyhedron::Halfedge_around_vertex_circulator first, curr, next;
  first = curr = v1->vertex_begin();
  do
  {
    next = curr;
    next++;

    p2 = curr->opposite()->vertex()->point();
    p3 = next->opposite()->vertex()->point();

    v12 = p2-p1;
    v13 = p3-p1;

    theta += std::acos( dot(v12,v13) / (norm(v12)*norm(v13)) );
  }
  while(++curr!=first);

  v1->set_gauss_curvature(2.0*M_PI - theta);
}

void ComputeGaussCurvature(Polyhedron &P)
{
  Polyhedron::Vertex_iterator vit;
  for(vit=P.vertices_begin(); vit!=P.vertices_end(); vit++)
    ComputeGaussCurvature(vit);
}

void PolyhedronShortestPath(Polyhedron& P,
                            vertex_descriptor &a, vertex_descriptor &b,
                            HalfedgeVector &hv)
{
  // create indices, predecessor and distance map
  VertexIdPMap vertex_index_pmap = get(boost::vertex_index, P);

  std::vector<vertex_descriptor> predecessor(boost::num_vertices(P));
  PredecessorPMap predecessor_pmap = PredecessorPMap(predecessor.begin(), vertex_index_pmap);

  std::vector<double> distance(boost::num_vertices(P));
  DistancePMap distance_pmap = DistancePMap(distance.begin(), vertex_index_pmap);

  // boost: set weights as edge length (default weight is squared length between vertices)
  std::map<edge_descriptor, double> edge2weight;
  WeightAPMap weight_apmap = WeightAPMap(edge2weight);
  edge_iterator ei, ei_end;
  for(boost::tie(ei,ei_end)=boost::edges(P); ei!=ei_end; ++ei)
    edge2weight.insert(make_pair(*ei,(*ei)->opposite()->weight())); // opposite: we compute from b to a

  // run dijkstra
  boost::detail::dijkstra_dispatch2(P, b, distance_pmap, weight_apmap, vertex_index_pmap,
                                    distance_map(distance_pmap).predecessor_map(predecessor_pmap));

  // extract edges
  vertex_descriptor v = a;
  for(vertex_descriptor u = predecessor_pmap[v]; u!=v; v=u, u=predecessor_pmap[v])
    hv.push_back(GetHalfedge(v,u));
}

void SetEdgeWeightsEuclidean(Polyhedron& P)
{
  edge_descriptor e;
  for(e=P.edges_begin(); e!=P.edges_end(); e++)
  {
    e->weight (sqrt( (e->vertex()->point() - e->opposite()->vertex()->point()).squared_length() ));
    e->opposite()->weight(e->weight());
  }
}

} // namespace tspline
