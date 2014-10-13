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

#include "Domain.h"
#include "Given_border_parameterizer_3.h"
#include <CGAL/parameterize.h>
#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/Discrete_authalic_parameterizer_3.h>

using namespace tspline;

// **********************************************************
// Domain

Domain::Domain(Polyhedron* p, int id) :
  mID(id), polyhedron(p), complete(false), filled(false), parameterized(false),
  adaptor(NULL), patch(NULL)
{
  borders.assign(4, HalfedgeVector(mID));

  borders[0].prev(&borders[3]);
  borders[1].prev(&borders[0]);
  borders[2].prev(&borders[1]);
  borders[3].prev(&borders[2]);

  borders[0].next(&borders[1]);
  borders[1].next(&borders[2]);
  borders[2].next(&borders[3]);
  borders[3].next(&borders[0]);
}

Domain::~Domain()
{
  for(size_t i=0; i<borders.size(); i++)
    borders[i].remove();

  if(adaptor!=NULL)
    delete adaptor;
  if(patch!=NULL)
    delete patch;
}

void Domain::add(vertex_descriptor vd, bool is_corner)
{
  // first vertex inserted
  if(vertices.empty())
  {
    vertices.push_back(vd);
    return;
  }

  // intermediate vertices / borders
  size_t vs = vertices.size();
  if(vs>0 && vs<4)
  {
    vertex_descriptor vd_last;
    if(borders[vs-1].empty())
      vd_last = vertices.back();
    else
      vd_last = borders[vs-1].back()->vertex();

    HalfedgeVector hv(mID);
    PolyhedronShortestPath(*polyhedron, vd_last, vd, hv);
    borders[vs-1].append(hv);

    if(is_corner)
      vertices.push_back(vd);
    else
      borders[vs-1].intermediate(vd);
  }

  // last vertex and last border
  if(vs==4)
  {
    vertex_descriptor vd_last;
    if(borders[vs-1].empty())
      vd_last = vertices.back();
    else
      vd_last = borders[vs-1].back()->vertex();

    if(is_corner)
    {
      HalfedgeVector hv(mID);
      PolyhedronShortestPath(*polyhedron, vd_last, vertices.front(), hv);
      borders[vs-1].append(hv);
      complete = true;
    }
    else
    {
      HalfedgeVector hv(mID);
      PolyhedronShortestPath(*polyhedron, vd_last, vd, hv);
      borders[vs-1].append(hv);
      borders[vs-1].intermediate(vd);
    }
  }

  return;
}

bool Domain::add_border(Domain* domain, const unsigned char& border_idx, bool flip, unsigned char multiplicity)
{
  HalfedgeVector& border = domain->borders.at(border_idx);

  if(vertices.empty())
  {
    printf("[Domain::add_border] [%d] Warning, vertices empty.\n", mID);
    return false;
  }

  size_t vs = vertices.size();
  if(!borders[vs-1].empty())
  {
    printf("[Domain::add_border] [%d] Warning, border not empty.\n", mID);
    return false;
  }

  HalfedgeVector hv(mID);
  if(flip)
  {
    HalfedgeVector::const_reverse_iterator rev;
    for(rev = border.rbegin(); rev!=border.rend(); ++rev)
      hv.push_back( (*rev)->opposite() );
  }
  else
  {
    hv.insert(hv.end(), border.begin(), border.end());
  }

  if(hv.front()->opposite()->vertex()!=vertices.back())
  {
    printf("[Domain::add_border] Warning, starting vertex of border does not coincide with last vertex of domain.\n");
    return false;
  }

  borders[vs-1].append(hv);

  borders[vs-1].opposite(&border);
  border.opposite(&borders[vs-1]);
  transitions.push_back(Transition(mID, vs-1, domain->id(), border_idx, multiplicity));

  if(vs<4)
    vertices.push_back(hv.back()->vertex());
  else
    complete = true;
  return true;
}

void Domain::remove_domain_from_transitions(int dom_id)
{
  std::vector<Transition>::iterator it;
  for(it=transitions.begin(); it!=transitions.end(); it++)
  {
    if((*it).domain1 == dom_id || (*it).domain2 == dom_id)
    {
      it=transitions.erase(it);
      if(it!=transitions.begin())
        it--;
    }
  }
}

//void Domain::connect(Domain *dom)
//{
//  if(!is_complete() || !dom->is_complete())
//    throw std::runtime_error("[Domain::connect] Error, domain not complete.");

//  for(size_t i=0; i<4; i++)
//  {
//    HalfedgeVector& v1 = borders[i];
//    if(v1.has_opposite())
//      continue;

//    for(size_t j=0; j<4; j++)
//    {
//      HalfedgeVector& v2 = dom->borders[j];
//      if(v2.has_opposite())
//        continue;

//      if(v1.is_opposite(v2))
//      {
//        printf("[Domain::connect] is opposite domain[%d]->borders[%lu] domain[%d]->borders[%lu]\n",
//               mID, i, dom->id(), j);
//        v1.opposite(&v2);
//        v2.opposite(&v1);
//        transitions.push_back(Transition(mID, i, dom->id(), j));
//      }

//    }
//  }

//}

std::vector<vertex_descriptor> Domain::get_selectables(size_t vertex_id) const
{
  std::vector<vertex_descriptor> v;
  int v1 = static_cast<int>(vertex_id) - 1;
  int v2 = static_cast<int>(vertex_id) + 1;

  if(v1 < 0)    v1 = 3;
  if(v1 > 3)    v1 = 0;
  if(v2 < 0)    v2 = 3;
  if(v2 > 3)    v2 = 0;

  v.push_back(vertices[v1]);
  v.push_back(vertices[v2]);

  return v;
}

void Domain::fill(Queue &queue)
{
  Polyhedron::Face_iterator f0 = queue.front().first;
  Polyhedron::Face_iterator f1 = queue.front().second;
  queue.pop_front();

  // fill bordering faces
  Polyhedron::Halfedge_around_facet_circulator first0, curr0, first1, curr1;
  first0 = curr0 = f0->facet_begin();
  first1 = curr1 = f1->facet_begin();

  size_t i(0);
  for(curr1=f1->facet_begin();
      !equal(curr1->vertex()->point(),curr0->vertex()->point()) && i<3;
      curr1++) i++;

  if(i>=3)
    throw std::runtime_error("[Domain::fill] Error, start point not found.");

  do
  {
    // if halfedge 'curr' is not a border of domain
    if(curr0->domID()!=mID)
    {
      Polyhedron::Halfedge_iterator hb0 = curr0->opposite();
      Polyhedron::Halfedge_iterator hb1 = curr1->opposite();
      Polyhedron::Face_iterator fb0 = hb0->facet();

      // bordering facet is not labled yet
      if(fb0->domID()==-1)
      {
        if(!hb1->is_border())
          throw std::runtime_error("[Domain::fill] Error, not a border.");

        // add face to polydomain and label it
        fb0->domID(mID);
        Polyhedron::Halfedge_iterator hn1 = polydomain.add_vertex_and_facet_to_border(hb1->prev(), hb1);

        // set new vertex
        Point3d pn0 = hb0->next()->vertex()->point();
        hn1->vertex()->point() = pn0;
        hn1->vertex()->id() = polydomain.size_of_vertices() - 1;

        queue.push_back(FacePair(fb0, hb1->face()));
      }
    }

    ++curr1;
  }
  while(++curr0!=first0);
}

void Domain::fill()
{
  if(filled)
    return;

  // make sure all halfedges have the correct domain id
  for(size_t i=0; i<borders.size(); i++)
    borders.at(i).set_halfedges_domID();

  if(borders.empty())
    throw std::runtime_error("[Domain::fill] Error, borders empty.");

  // find a halfedge inside the domain
  HalfedgeVector::iterator he;
  HalfedgeVector& hv0 = borders[0];
  for( he=hv0.begin();
       he!=hv0.end() && (*he)->opposite()->domID()==(*he)->domID();
       he++);

  if(he==hv0.end())
    throw std::runtime_error("[Domain::fill] Error, no starting halfedge found.");

  // clear polydomain
  polydomain = Polyhedron();

  // create first triangle
  Polyhedron::Face_iterator f0 = (*he)->facet();
  size_t i(0);
  Point3d p[3];
  Polyhedron::Halfedge_around_facet_circulator first, curr;
  first = curr = f0->facet_begin();
  do
  {
    p[i] = curr->vertex()->point();
  }
  while(++curr!=first && ++i<3);
  f0->domID(mID);
  Polyhedron::Face_iterator f1 = polydomain.make_triangle(p[0],p[1],p[2])->facet();

  // assign indices
  Polyhedron::Vertex_iterator vit;
  size_t index(0);
  for(vit=polydomain.vertices_begin(); vit!=polydomain.vertices_end(); vit++)
    vit->id() = index++;

  // add triangles
  Queue queue;
  queue.push_back(FacePair(f0,f1));
  while(!queue.empty())
    fill(queue);

  filled = true;
}

void Domain::create_seam(Seam& seam, SeamUV &seam_uv)
{
  if(!complete)
    throw std::runtime_error("[Domain::create_seam] Error, patch not complete.");

  seam.clear();
  seam_uv.clear();

  // create seam and assign uv coordinates
  double u,v,d;
  for(size_t i=0; i<borders.size(); i++)
  {
    HalfedgeVector hv = borders[i];

    for(size_t j=0; j<hv.size(); j++)
    {
      d = double(j+1) / (hv.size());

      if(i==0) { u = d; v = 0.0; } // bottom / south
      if(i==1) { u = 1.0; v = d; } // right / east
      if(i==2) { u = 1.0-d; v = 1.0; } // top / north
      if(i==3) { u = 0.0, v = 1.0-d; } // left / west

      seam.push_front(hv[j]->vertex());
      seam_uv.push_front(Point2d(u,v));
    }
  }
}

bool Domain::parameterize()
{
  if(parameterized)
    return true;

  if(!complete)
  {
    printf("[Domain::parameterize()] Warning, domain %d not complete.\n", mID);
    return false;
  }

  // create seam from borders
  Seam seam;
  SeamUV seam_uv;
  create_seam(seam, seam_uv);
  if (seam.empty())
    throw std::runtime_error("[Domain::parameterize()] Error, seam not created.");

  // Create Polyhedron adaptor
  adaptor = new ParamPolyAdaptor(*polyhedron);

  // Create a second adaptor that virtually "cuts" the mesh following the 'seam' path
  patch = new ParamMeshPatch(*adaptor, seam.begin(), seam.end());
  if (!patch->is_valid())
  {
    printf("[Domain::parameterize()] Error, non manifold shape or invalid cutting");
    return false;
  }

  // assign uv coordinates of seam
  Seam::iterator sit = seam.begin();
  SeamUV::iterator suv = seam_uv.begin();
  for(sit  =  seam.begin(); sit !=  seam.end(); sit++)
  {
    patch->get_decorated_mesh().set_vertex_uv((*sit), (*suv));
    suv++;
  }

  // Border parametrization
  typedef CGAL::Given_border_parameterizer_3<ParamMeshPatch> BorderParam;

  // Eigen solver
  typedef CGAL::Eigen_solver_traits<> Solver;

  // Discrete Authalic Parameterization (square border) with Eigen solver
  typedef CGAL::Discrete_authalic_parameterizer_3<ParamMeshPatch,BorderParam,Solver> Parameterizer;

  Parameterizer::Error_code err = CGAL::parameterize(*patch, Parameterizer());

  switch(err)
  {
  case Parameterizer::OK: // Success
    break;

  case Parameterizer::ERROR_EMPTY_MESH:
  case Parameterizer::ERROR_NON_TRIANGULAR_MESH:
  case Parameterizer::ERROR_NO_TOPOLOGICAL_DISC:
  case Parameterizer::ERROR_BORDER_TOO_SHORT:
  default:
    std::cerr << "[Domain::parameterize] Error: " << Parameterizer::get_error_message(err) << std::endl;
    printf("[Domain::parameterize] Error, parameterization failed.\n");
    return false;
  };

  //  ParamMeshPatch::Vertex_iterator pmp_it;
  //  for(pmp_it  = patch->mesh_vertices_begin();
  //      pmp_it != patch->mesh_vertices_end();
  //      pmp_it++)
  //  {
  //    Point2d uv = patch->get_vertex_uv(pmp_it);
  //    std::cout << "[Domain::parameterize] " << uv.x() << ", " << uv.y();
  //    uv = adjust(uv);
  //    std::cout << " | " << uv.x() << ", " << uv.y() << std::endl;
  //  }

  parameterized = true;
  return true;
}

void Domain::get_points_and_params(vector_vec3d& points, vector_vec2d& params) const
{
  if(!parameterized)
  {
    printf("[Domain::get_points_and_params] Warning, domain %d not parameterized\n", mID);
    return;
  }

  ParamMeshPatch::Vertex_iterator pmp_it;
  for(pmp_it  = patch->mesh_vertices_begin();
      pmp_it != patch->mesh_vertices_end();
      pmp_it++)
  {
    Point3d p = patch->get_vertex_position(pmp_it);
    Point2d uv = adjust(patch->get_vertex_uv(pmp_it));
    //    std::cout << " | " << uv.x() << ", " << uv.y() << std::endl;

    points.push_back(Eigen::Vector3d(p.x(),p.y(),p.z()));
    params.push_back(Eigen::Vector2d(uv.x(),uv.y()));
  }
}

// **********************************************************
// Multi Domain

MultiDomain::~MultiDomain()
{
  for(size_t i=0; i<domains.size(); i++)
    delete domains[i];
}

void MultiDomain::add(vertex_descriptor vd, bool is_corner, bool selectable)
{
  if(domains.empty())
    domains.push_back(new Domain(polyhedron,0));

  Domain* dom = domains.back();

  if(domains.back()->is_complete())
  {
    int dom_id = static_cast<int>(domains.size());
    printf("[MultiDomain::add] adding new domain: %d\n", dom_id);
    domains.push_back(new Domain(polyhedron, dom_id));
    dom = domains.back();
    prev_vertex_selectable = false;
  }

  //  if( selectable && prev_vertex_selectable )
  //    add_border(dom, prev_selectable, vd);
  //  else if( dom->get_vertices().size()==4 && is_corner && prev_vertex_selectable )
  //    add_border(dom, prev_selectable, dom->get_vertices().front());
  //  else
  //    dom->add(vd, is_corner);

  if(!selectable ||
     !prev_vertex_selectable ||
     !add_border(dom, prev_selectable, vd))
    dom->add(vd, is_corner);

  prev_vertex_selectable = selectable;
  if(prev_vertex_selectable)
    prev_selectable = vd;

  return;
}

bool MultiDomain::add_border(Domain* dom, vertex_descriptor v1, vertex_descriptor v2)
{
  for(size_t i=0; i<domains.size(); i++)
  {
    Domain* dom_i = domains[i];
    if(dom_i->is_complete())
    {
      const std::vector<vertex_descriptor>& vertices = dom_i->get_vertices();
      for(size_t j=0; j<vertices.size(); j++)
      {
        if( v1==vertices.at(j%4) && v2==vertices.at((j+1)%4) )
        {
          printf("[MultiDomain::add_border] Error, wrong direction.\n");
          return false;
        }
        //          return dom->add_border(dom_i->get_border(j), false);

        if( v1==vertices.at((j+1)%4) && v2==vertices.at(j%4) )
          return dom->add_border(dom_i, j, true, border_multiplicity);
      }
    }
  }

  printf("[MultiDomain::add_border] Warning, no border added.\n");
  return false;
}


//bool MultiDomain::exists(vertex_descriptor vd, size_t &domain_id, size_t &vertex_id)
//{
//  size_t& i = domain_id;
//  size_t& j = vertex_id;

//  for(i=0; i<domains.size(); i++)
//  {
//    const std::vector<vertex_descriptor>& vertices = domains[i]->get_vertices();
//    for(j=0; j<vertices.size(); j++)
//      if(vd==vertices.at(j))
//        return true;
//  }
//  return false;
//}

std::vector<vertex_descriptor> MultiDomain::get_selectables() const
{
  std::vector<vertex_descriptor> vertices;
  for(size_t i=0; i<domains.size(); i++)
  {
    if(domains[i]->is_complete())
    {
      const std::vector<vertex_descriptor>& dv = domains[i]->get_vertices();
      for(size_t j=0; j<dv.size(); j++)
        if(find(vertices.begin(), vertices.end(), dv.at(j))==vertices.end())
          vertices.push_back(dv.at(j));
    }
  }
  return vertices;
}

std::vector<vertex_descriptor> MultiDomain::get_vertices() const
{
  std::vector<vertex_descriptor> vertices;
  for(size_t i=0; i<domains.size(); i++)
  {
    const std::vector<vertex_descriptor>& dv = domains[i]->get_vertices();
    for(size_t j=0; j<dv.size(); j++)
      if(find(vertices.begin(), vertices.end(), dv.at(j))==vertices.end())
        vertices.push_back(dv.at(j));
  }
  return vertices;
}

std::vector<HalfedgeVector> MultiDomain::get_borders() const
{
  std::vector<HalfedgeVector> borders;
  for(size_t i=0; i<domains.size(); i++)
  {
    for(size_t j=0; j<4; j++)
    {
      const HalfedgeVector& b = domains[i]->get_border(j);
      if(!b.empty())
        borders.push_back(b);
    }
  }
  return borders;
}

//void MultiDomain::connect()
//{
//  for(size_t i=0; i<domains.size()-1; i++)
//  {
//    for(size_t j=(i+1); j<domains.size(); j++)
//    {
//      if(domains[j]->is_complete())
//        domains[i]->connect(domains[j]);
//    }
//  }
//}

void MultiDomain::fill()
{
  for(size_t i=0; i<domains.size(); i++)
    if(domains[i]->is_complete() && !domains[i]->is_filled())
      domains[i]->fill();
}

void MultiDomain::parametrize()
{
  for(size_t i=0; i<domains.size(); i++)
    if(domains[i]->is_complete() && !domains[i]->is_parameterized())
    {
      printf("[MultiDomain::parametrize] parameterizing domain %d ...\n", domains[i]->id());
      bool success = domains[i]->parameterize();
      if(success)
        printf("[MultiDomain::parametrize] parameterizing domain %d done\n", domains[i]->id());
      else
        printf("[MultiDomain::parametrize] parameterizing domain %d failed\n", domains[i]->id());
    }
}

void MultiDomain::pop_back()
{
  if(!domains.empty())
  {
    int id = domains.back()->id();
    delete(domains.back());
    domains.pop_back();
    for(size_t i=0; i<domains.size(); i++)
      domains[i]->remove_domain_from_transitions(id);
  }
}
