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

#ifndef _TSPLINE_DOMAIN_H_
#define _TSPLINE_DOMAIN_H_

// This file defines a rectangular domain embedded on a CGAL::Polyhedron manifold
// It implements insertion and removal of edges

// This file defines an arrangement of domains called MultiDomain
// It implements adding and removal of boundaries and domains and provides functions
// for domain editing (4-point clicking of domain corners, see DomainEditor)

#include "Types.h"
#include "Polyhedron.h"

namespace tspline {


/** @brief  CGAL::Polyhedron Domain (not T-spline related)
 *          for creation of a rectangular domain embedded on a Polyhedron structure
 */
class Domain
{
public:
  struct Transition
  {
    int domain1;
    int border1;
    int domain2;
    int border2;
    unsigned char multiplicity;
    Transition(int d1, int b1, int d2, int b2, unsigned char mult) :
      domain1(d1), border1(b1), domain2(d2), border2(b2), multiplicity(mult) {}
    bool operator==(const Transition& t)
    {
      return (domain1==t.domain1 &&
              border1==t.border1 &&
              domain2==t.domain2 &&
              border1==t.border1);
    }
  };

private:
  int mID;
  Polyhedron* polyhedron;
  std::vector<HalfedgeVector> borders;
  std::vector<vertex_descriptor> vertices;
  bool complete;
  bool filled;
  bool parameterized;

  Polyhedron polydomain;

  ParamPolyAdaptor* adaptor;
  ParamMeshPatch* patch;

  typedef std::pair<Polyhedron::Face_iterator,Polyhedron::Face_iterator> FacePair;
  typedef std::list< FacePair > Queue;

  std::vector<Transition> transitions;

public:
  Domain(Polyhedron* p, int id);
  ~Domain();

  /** @brief Adds a vertex and the resulting HalfedgeVector to the Domain
   *  @param vd  the vertex
   *  @param is_corner true if vertex should be treated as corner of the domain border */
  void add(vertex_descriptor vd, bool is_corner);

  bool add_border(Domain *domain, const unsigned char &border_idx, bool flip, unsigned char multiplicity);

  void remove_domain_from_transitions(int dom_id);
  //  void connect(Domain* d);

  bool is_complete() const{ return complete; }
  bool is_filled() const{ return filled; }
  bool is_parameterized() const { return parameterized; }
  std::vector<vertex_descriptor> get_selectables(size_t vertex_id) const;
  const std::vector<vertex_descriptor>& get_vertices() const { return vertices; }
  const HalfedgeVector& get_border(unsigned char i) const { return borders[i]; }

  void fill(Queue& queue);
  void fill();

  void create_seam(Seam& seam, SeamUV& seam_uv);
  bool parameterize();

  void get_points_and_params(vector_vec3d& points, vector_vec2d& params) const;

  Polyhedron& get_polydomain() { return polydomain; }

  const std::vector<Transition>& get_transitions() const { return transitions; }

  const ParamPolyAdaptor& get_adaptor() const
  { if(adaptor!=NULL) return *adaptor;
    throw std::runtime_error("[Domain::get_adaptor] Error, no adaptor available."); }
  const ParamMeshPatch& get_patch() const
  { if(patch!=NULL) return *patch;
    throw std::runtime_error("[Domain::get_patch] Error, no patch available."); }

  const int& id() const { return mID; }
};

/** @brief  CGAL object (not T-spline related)
 *          for creation of multiple connected quadratic domains embedded on a Polyhedron structure */
class MultiDomain
{
private:
  Polyhedron* polyhedron;
  std::vector<Domain*> domains;

  bool prev_vertex_selectable;
  vertex_descriptor prev_selectable;

  unsigned char border_multiplicity;

public:
  MultiDomain(Polyhedron* p) : polyhedron(p), prev_vertex_selectable(false), border_multiplicity(0) {}
  ~MultiDomain();

  /** @brief Adds a vertex and the resulting HalfedgeVector to the Domain
   *  @param vd  the vertex
   *  @param is_corner true if vertex should be treated as corner of the domain border
   *  @return true if domain is complete, false if domain is in-complete (i.e. borders/vertices missing) */
  void add(vertex_descriptor vd, bool is_corner, bool selectable);

  bool add_border(Domain* dom, vertex_descriptor v1, vertex_descriptor v2);

  //  bool exists(vertex_descriptor vd, size_t& domain_id, size_t& vertex_id);

  std::vector<vertex_descriptor> get_selectables() const;
  std::vector<vertex_descriptor> get_vertices() const;
  std::vector<HalfedgeVector> get_borders() const;

  void set_multiplicity(unsigned char m) { border_multiplicity = m; }
  const unsigned char& get_multiplicity() { return border_multiplicity; }

  //  void connect();
  void fill();
  void parametrize();
  void pop_back();

  const Domain& get_domain(size_t i) const
  { if(i>=0 && i<domains.size()) return *domains[i];
    throw std::runtime_error("[Domain::get_domain] Error, index out ouf bounds."); }


  size_t size(){ return domains.size(); }

  bool has_prev_selectable() { return prev_vertex_selectable; }
  const vertex_descriptor& get_prev_selectable() const { return prev_selectable; }
};



} // namespace tspline

#endif
