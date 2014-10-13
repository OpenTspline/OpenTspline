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

#include "Tspline.h"

using namespace std;
using namespace tspline;

void Tspline::update_params(Vertex_iterator &v1)
{
  TVertex vext1 = v1->data();
  if (vext1.updated)
    return;

  vext1.updated = true;
  v1->set_data(vext1);

  Halfedge_around_vertex_circulator first, curr;
  first = curr = v1->incident_halfedges();

  do {
    Vertex_iterator v0 = curr->source();
    TVertex vext0 = v0->data();
    THalfedge hext = curr->data();

    Vector2d v = v1->point() - v0->point();
    double n = v.squared_length();
    v = v / sqrt(n);
    //if (hext.d == 0.0)
    if(equal(hext.d, 0.0))
      vext0.param = adjust(vext1.param);
    else
      vext0.param = adjust(vext1.param - v * hext.d); // minus, because of direction of v
    v0->set_data(vext0);


    if(smaller(vext0.param.x(), param_min.x()))
      param_min = Point2d(vext0.param.x(), param_min.y());
    if (greater(vext0.param.x(), param_max.x()))
      param_max = Point2d(vext0.param.x(), param_max.y());
    if (smaller(vext0.param.y(), param_min.y()))
      param_min = Point2d(param_min.x(), vext0.param.y());
    if (greater(vext0.param.y(), param_max.y()))
      param_max = Point2d(param_max.x(), vext0.param.y());

    if(smaller(v0->point().x(), point_min.x()))
      point_min = Point2d(v0->point().x(), point_min.y());
    if (greater(v0->point().x(), point_max.x()))
      point_max = Point2d(v0->point().x(), point_max.y());
    if (smaller(v0->point().y(), point_min.y()))
      point_min = Point2d(point_min.x(), v0->point().y());
    if (greater(v0->point().y(), point_max.y()))
      point_max = Point2d(point_max.x(), v0->point().y());

    update_params(v0);
  } while (++curr != first);
}

void Tspline::update_params()
{
  for (Vertex_iterator vit = vertices_begin(); vit != vertices_end(); vit++)
  {
    TVertex& vext = vit->data();
    vext.updated = false;
  }

  Vertex_iterator v = vertices_begin();
  TVertex vext = v->data();
  vext.param = v->point();
  v->set_data(vext);
  param_min = Point2d(DBL_MAX, DBL_MAX);
  param_max = Point2d(-DBL_MAX, -DBL_MAX);
  point_min = Point2d(DBL_MAX, DBL_MAX);
  point_max = Point2d(-DBL_MAX, -DBL_MAX);

  update_params(v);
}

void Tspline::update_knot_vectors()
{
  Vertex_iterator vit;
  for (vit = vertices_begin(); vit != vertices_end(); vit++)
  {
    TVertex& vext = vit->data();
    compute_knot_vectors(vit, vext.s, vext.t);
  }
}

void Tspline::update_vertex_ids()
{
  int id(0);
  Vertex_iterator vit;
  for (vit = vertices_begin(); vit != vertices_end(); vit++)
    vit->data().id = id++;
}

unsigned Tspline::refine_average_error(double avg_error, unsigned min_points)
{

  double avg_error_sqr = avg_error * avg_error;

  std::vector<Face_handle> faces;
  Face_iterator fit;
  for(fit=faces_begin(); fit!=faces_end(); fit++)
    if(!fit->is_unbounded() && fit->has_outer_ccb())
      faces.push_back(fit);

  unsigned refined(0);
  std::vector<Face_handle>::reverse_iterator rit;
  for(rit=faces.rbegin(); rit!=faces.rend(); rit++)
  {
    Face_handle fh = (*rit);
    double err = fh->data().error_sqr;
    unsigned pts = fh->data().point_count;
    double avg = err / pts;

    if(avg > avg_error_sqr)
      if(refine_point_count(fh, min_points))
        refined++;
  }

  insert_missing_edges();
  update_knot_vectors();

  return refined;
}

void Tspline::refine(Face_handle f, bool insert_edges, bool update_knots) // assumes that only rectangular faces exist in grid
{
  if (f->is_unbounded() || !f->has_outer_ccb())
    throw std::runtime_error("[Tspline::refine] Error, face is unbounded or has no outer contour\n");

  Face_const_iterator fc = f;
  Eigen::Vector4d bb;
  face_bounding_box(fc, bb);

  double s = (bb(1) + bb(0)) * 0.5;
  double t = (bb(3) + bb(2)) * 0.5;

  Point2d param(s, t);
  insert_vertex(f, param);

  if(insert_edges)
    insert_missing_edges();
  if(update_knots)
    update_knot_vectors();
}

void Tspline::refine_congruent(Face_handle f) // assumes that only rectangular faces exist in grid
{
  if (f->is_unbounded() || !f->has_outer_ccb())
    throw std::runtime_error("[Tspline::refine_congruent] Error, face is unbounded or has no outer contour\n");

  Face_const_iterator fc = f;
  Eigen::Vector4d bb;
  face_bounding_box(fc, bb);

  double w = bb(1) - bb(0);
  double h = bb(3) - bb(2);

  if(w>h)
    split_vertical_congruent(f, bb(0)+0.5*w);
  else
    split_horizontal_congruent(f, bb(2)+0.5*h);
}

unsigned Tspline::refine_point_count(unsigned max_num_points)
{
  std::vector<Face_handle> faces;
  Face_iterator fit;
  for(fit=faces_begin(); fit!=faces_end(); fit++)
    if(!fit->is_unbounded() && fit->has_outer_ccb())
      faces.push_back(fit);

  unsigned refined(0);
  std::vector<Face_handle>::reverse_iterator rit;
  for(rit=faces.rbegin(); rit!=faces.rend(); rit++)
    if(refine_point_count((*rit), max_num_points))
      refined++;

  insert_missing_edges();
  update_knot_vectors();

  return refined;
}

bool Tspline::refine_point_count(Face_handle f, unsigned max_num_points)
{
  const TFace& fext = f->data();
  if(fext.point_count < max_num_points)
    return false;

  Eigen::Vector4d bb;
  face_bounding_box(f, bb);
  double ds, dt;
//  face_dimensions_squared(f, ds, dt, bb);
    ds = bb(1)-bb(0);
    dt = bb(3)-bb(2);

  if(ds > dt)
  {
    // split vertical
    split_vertical( f, 0.5*(bb(0)+bb(1)) );
    return true;
  }
  else
  {
    // split horizontal
    split_horizontal(f, 0.5*(bb(2)+bb(3)) );
    return true;
  }

  return false;
}

unsigned Tspline::refine_chord(double chord_length_threshold)
{
  std::vector<Face_handle> faces;
  Face_iterator fit;
  for(fit=faces_begin(); fit!=faces_end(); fit++)
    if(!fit->is_unbounded() && fit->has_outer_ccb())
      faces.push_back(fit);

  unsigned refined(0);
  std::vector<Face_handle>::reverse_iterator rit;
  for(rit=faces.rbegin(); rit!=faces.rend(); rit++)
    if(refine_chord((*rit), chord_length_threshold))
      refined ++;

  insert_missing_edges();
  update_knot_vectors();

  return refined;
}

bool Tspline::refine_chord(Face_handle f, double chord_length_threshold)
{
  Eigen::Vector4d bb;
  face_bounding_box(f, bb);
  double ds, dt;
  face_dimensions_squared(f, ds, dt, bb);

  double clt_sqr = chord_length_threshold * chord_length_threshold;

  if(ds > clt_sqr || dt > clt_sqr)
  {
    if(ds > dt)
    {
      // split vertical
      split_vertical( f, 0.5*(bb(0)+bb(1)) );
      return true;
    }
    else
    {
      // split horizontal
      split_horizontal(f, 0.5*(bb(2)+bb(3)) );
      return true;
    }
  }

  return false;
}

Tspline::Vertex_iterator Tspline::insert_vertex(const double &s, const double &t,
                                                bool insert_edges, bool update_knots)
{
  Point2d param(s, t);

  CGAL::Object obj = locate_param(s, t);

  Face_iterator f;
  Halfedge_iterator e;
  Vertex_iterator v;

  if (!obj.empty()) {
    if (CGAL::assign(f, obj))
      v=insert_vertex(f, param);  // re-implement using split-edge and split-face operations
    else if (CGAL::assign(e, obj))
      v=insert_vertex(e, param);
    else if (CGAL::assign(v, obj))
      ;//printf("[Tspline::insert_vertex] Warning, vertex already existing at param position (%e %e)\n", s, t);
    else
      throw std::runtime_error("[Tspline::insert_vertex] Error, no object found");

    if(insert_edges)
      insert_missing_edges();
    if(update_knots)
      update_knot_vectors();

    return v;
  }

  return vertices_end();
}

Tspline::Halfedge_iterator Tspline::split_horizontal(Face_iterator f)
{
  if (f->is_unbounded () || f==faces_end())
    throw std::runtime_error("[Tspline::split_horizontal] Error, face not valid.");

  Eigen::Vector4d bb;
  face_bounding_box(f, bb);

  double y = (bb(3) + bb(2)) * 0.5;

  Halfedge_iterator h = split_horizontal(f,y);

  if(h==halfedges_end())
    throw std::runtime_error("[Tspline::split_horizontal] Error, split not valid.");

  insert_missing_edges();
  update_knot_vectors();

  return h;
}

Tspline::Halfedge_iterator Tspline::split_vertical(Face_iterator f)
{
  if (f->is_unbounded () || f==faces_end())
    throw std::runtime_error("[Tspline::split_vertical] Error, face not valid.");

  Eigen::Vector4d bb;
  face_bounding_box(f, bb);

  double x = (bb(1) + bb(0)) * 0.5;

  Halfedge_iterator h = split_vertical(f,x);

  if(h==halfedges_end())
    throw std::runtime_error("[Tspline::split_vertical] Error, split not valid.");

  insert_missing_edges();
  update_knot_vectors();

  return h;
}

Tspline::Halfedge_iterator Tspline::split_horizontal_congruent(Face_iterator &f, double y)
{
  if (f->is_unbounded () || f==faces_end())
    throw std::runtime_error("[Tspline::split_horizontal_congruent] Error, face not valid.");

  int num_vertices_prev = number_of_vertices();

  Halfedge_iterator h = split_horizontal(f,y);

  if(h==halfedges_end())
    throw std::runtime_error("[Tspline::split_horizontal_congruent] Error, split not valid.");

  insert_missing_edges();
  make_congruent(num_vertices_prev);
  insert_missing_edges();

  return h;
}

Tspline::Halfedge_iterator Tspline::split_vertical_congruent(Face_iterator &f, double x)
{
  if (f->is_unbounded () || f==faces_end())
    throw std::runtime_error("[Tspline::split_vertical_congruent] face not valid.");

  int num_vertices_prev = number_of_vertices();

  Halfedge_iterator h = split_vertical(f,x);

  if(h==halfedges_end())
    throw std::runtime_error("[Tspline::split_vertical_congruent] Error, split not valid.");

  insert_missing_edges();
  make_congruent(num_vertices_prev);
  insert_missing_edges();

  return h;
}

Tspline::Vertex_iterator Tspline::split_edge_congruent(Halfedge_iterator &h, Point2d param)
{
  if(h==halfedges_end())
    return vertices_end();

  int num_vertices_prev = number_of_vertices();

  Vertex_iterator v = insert_vertex(h, param);

  if(v==vertices_end())
    return v;

  printf("[Tspline::split_edge_congruent] insert vertex: %d\n", v->data().id);
  v->data().PrintKnotVectors();

  insert_missing_edges();
  make_congruent(num_vertices_prev);
  insert_missing_edges();

  return v;
}

namespace CGAL{
template <class GeomTraits, class TopTraits>
bool remove_vertex
(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
 typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Vertex_handle v,
 typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Halfedge_handle& h)
{
  typedef Arrangement_on_surface_2<GeomTraits, TopTraits>   Arr;
  typedef Arr_traits_adaptor_2<GeomTraits>                  Traits_adaptor_2;

  // Notify the arrangement observers that a global operation is about to
  // take place.
  Arr_accessor<Arr>    arr_access (arr);

  arr_access.notify_before_global_change();

  // Act according to the number of edges incident to v.
  bool      removed = false;

  if (v->degree() == 2)
  {
    // If the vertex has now two incident edges, and the curves associated
    // with these edges can be merged, merge the two edges and remove the
    // vertex.
    const Traits_adaptor_2 * traits =
        static_cast<const Traits_adaptor_2*>(arr.geometry_traits());
    typename Arr::Halfedge_around_vertex_circulator  circ;
    typename Arr::Halfedge_handle                    e1, e2;

    circ = v->incident_halfedges();
    e1 = circ;
    ++circ;
    e2 = circ;

    if (traits->are_mergeable_2_object() (e1->curve(), e2->curve()))
    {
      // Merge the two curves.
      typename GeomTraits::X_monotone_curve_2   cv;
      traits->merge_2_object() (e1->curve(), e2->curve(),
                                cv);

      // Merge the two edges in the arrangement.
      h = arr.merge_edge (e1, e2, cv);
      removed = true;
    }
  }

  // Notify the arrangement observers that the global operation has been
  // completed.
  arr_access.notify_after_global_change();

  // Return an indication whether the vertex has been removed or not.
  return (removed);
}
}

void Tspline::remove_vertex(Vertex_handle vh,
                            bool remove_isolated,
                            bool insert_edges,
                            bool update_knots,
                            bool remove_I_junctions)
{
  if(vh->is_isolated())
  {
    CGAL::remove_vertex(*this, vh);
    return;
  }

  if(is_boundary(vh))
  {
    printf("[Tspline::remove_vertex] Error, vertex %d is boundary\n", vh->data().id);
    throw std::runtime_error("[Tspline::remove_vertex] Error, vertex is boundary");
  }

  Halfedge_around_vertex_circulator vFirst, vCurr;
  vFirst=vCurr=vh->incident_halfedges();
  std::vector<Halfedge_handle> halfedges;
  std::vector<Vertex_handle> neighbor_vertices;
  do
  {
    halfedges.push_back(vCurr);
    if(remove_I_junctions)
      neighbor_vertices.push_back(vCurr->source());
    //        printf("[Tspline::remove_vertex] he: %d %d\n", vCurr->source()->data().id, vCurr->target()->data().id);
  }while(++vCurr!=vFirst);

  Face_iterator fit;
  std::vector<Halfedge_handle>::reverse_iterator rhit;
  for(rhit=halfedges.rbegin(); rhit!=halfedges.rend(); rhit++)
    fit = this->remove_edge((*rhit), remove_isolated, true); // removes vh when it becomes isolated

  if(insert_edges)
    insert_missing_edges(fit);

  for(size_t i=0; i<neighbor_vertices.size(); i++)
  {
    Halfedge_handle h;
    Vertex_handle v = neighbor_vertices[i];
    if(remove_I_junctions && is_Ijunction(v))
      CGAL::remove_vertex(*this, v, h);
  }

  if(update_knots)
    update_knot_vectors();
}

size_t Tspline::remove_redundant_vertices(bool Ijunctions, bool loose, bool isolated)
{
  size_t removed(0);
  std::vector<Vertex_handle> vertices;
  Vertex_iterator vit;
  for(vit=vertices_begin(); vit!=vertices_end(); vit++)
  {
    if(!is_corner(vit)) // don't remove corner vertices
    {
      if(vit->is_isolated())
        vertices.push_back(vit); // isolated vertices

      else if(vit->degree()==1)
        vertices.push_back(vit); // loose edges

      else if(is_Ijunction(vit))
        vertices.push_back(vit); // I junction
    }
  }

  if(vertices.empty())
    return 0;

  std::vector<Vertex_handle>::reverse_iterator rit;
  Vertex_handle v;
  for(rit=vertices.rbegin(); rit!=vertices.rend(); rit++)
  {
    v = (*rit);
    if(v->is_isolated()) // remove isolated vertices
    {
      if(CGAL::remove_vertex(*this, v))
        removed++;
    }
    else if(v->degree()==1) // remove loose edges
    {
      Halfedge_handle h = v->incident_halfedges();
      if(this->remove_edge(h, false, true)!=faces_end())
        removed++;
    }
    else if(v->degree()==2) // remove I junction
    {
      Halfedge_around_vertex_circulator circ = v->incident_halfedges();
      double d = circ->data().d;
      circ++;
      d += circ->data().d;
      Halfedge_handle h;
      if(CGAL::remove_vertex(*this, v, h))
      {
        h->data().d = d;
        h->twin()->data().d = d;
        removed++;
      }
    }
  }

  return removed;
}

typedef std::pair<int, double> ID_Curvature;
typedef std::vector< ID_Curvature > ID_Curvature_Vector;

static bool ascending_curvature(const ID_Curvature& a, const ID_Curvature& b)
{
  return a.second < b.second;
}

void Tspline::reduce_vertices(double min_cage_curvature)
{
  ID_Curvature_Vector curvature;
  Tspline::Vertex_iterator vit;
  for(vit=vertices_begin(); vit!=vertices_end(); vit++)
  {
    if(!is_boundary(vit))
    {
      Eigen::Vector4d n = grid_normal_curvature(vit);
      curvature.push_back(ID_Curvature(vit->data().id, std::abs<double>(n[3])));
    }
  }

  std::sort(curvature.begin(), curvature.end(), ascending_curvature);

  for(size_t i=0; i<curvature.size(); i++)
  {
    ID_Curvature& c = curvature[i];
    Tspline::Vertex_handle vh = get_vertex(c.first);

    if(vh==vertices_end())
      throw std::runtime_error("[DomainCreator::SimplifyTsplineMultiPatch] Error, vertex not valid.");

    if(c.second<min_cage_curvature)
    {
      //      printf("[DomainCreator::SimplifyTsplineMultiPatch] removing: %d );\n", c.id);
      remove_vertex(vh, false, false, false, false);
    }
  }

  while(remove_redundant_vertices()>0) {};
  update_vertex_ids();
  insert_edges_at_L_junctions();
  insert_missing_edges();
  update_knot_vectors();
}

bool Tspline::reduce_vertices(double fitting_threshold, double sampling_resolution)
{
  double s0 = param_min.x();
  double s1 = param_max.x();
  double t0 = param_min.y();
  double t1 = param_max.y();
  double w = s1 - s0;
  double h = t1 - t0;
  double stepS = w * 0.25;
  double stepT = h * 0.25;
  double epsS = stepS * 0.1;
  double epsT = stepT * 0.1;

  Vertex_iterator v;
  for(double j=t0; j<t1+epsT; j+=stepT)
  {
    for(double i=s0; i<s1+epsS; i+=stepS)
    {
      CGAL::Object obj = locate_param(i,j);
      if(!CGAL::assign(v, obj))
      {
        printf("[Tspline::reduce_vertices] Error, T-spline not reduceable.\n");
        return false;
      }
    }
  }



  return true;
}

void Tspline::evaluate_basis(const double &s, const double &t, const TVertex &vext, double &B) const
{
  std::vector<double> Ns, Nt, vs, vt;
  tspline::cox(s, degree, vext.s, Ns);
  tspline::cox(t, degree, vext.t, Nt);

  size_t idx = vext.s.size() * degree; // idx of basis function that has peak at middle knot[2]

  B = Ns[idx] * Nt[idx];

  if (clamped) {
    std::vector<double> Nsc, Ntc;
    bool cl_left(false), cl_right(false), cl_bottom(false), cl_top(false);

    if (equal(vext.param.x(), param_min.x()) && smaller(s, vext.s[3])) // left
    {
      vs = vext.s;
      vs[4] = vs[3];
      vs[3] = vs[2];
      tspline::cox(s, degree, vs, Nsc);
      B += Nsc[idx] * Nt[idx];
      cl_left = true;
    } else if (equal(vext.param.x(), param_max.x()) && greater(s, vext.s[1])) // right
    {
      vs = vext.s;
      vs[0] = vs[1];
      vs[1] = vs[2];
      tspline::cox(s, degree, vs, Nsc);
      B += Nsc[idx] * Nt[idx];
      cl_right = true;
    }

    if (equal(vext.param.y(), param_min.y()) && smaller(t, vext.t[3])) // bottom
    {
      vt = vext.t;
      vt[4] = vt[3];
      vt[3] = vt[2];
      tspline::cox(t, degree, vt, Ntc);
      B += Ns[idx] * Ntc[idx];
      cl_bottom = true;
    } else if (equal(vext.param.y(), param_max.y()) && greater(t, vext.t[1])) // top
    {
      vt = vext.t;
      vt[0] = vt[1];
      vt[1] = vt[2];
      tspline::cox(t, degree, vt, Ntc);
      B += Ns[idx] * Ntc[idx];
      cl_top = true;
    }

    if ((cl_left || cl_right) && (cl_bottom || cl_top)) // corners
      B += Nsc[idx] * Ntc[idx];

  } // if(clamped)
}

void Tspline::evaluate_basis(const double &s, const double &t, const TVertex &vext, double &B, double &Bs, double &Bt) const
{
  std::vector<double> Ns, Nt, Nds, Ndt, vs, vt;
  tspline::cox(s, degree, vext.s, Ns);
  tspline::cox(t, degree, vext.t, Nt);
  tspline::coxder(degree, vext.s, Ns, Nds);
  tspline::coxder(degree, vext.t, Nt, Ndt);

  size_t idx = vext.s.size() * degree; // idx of basis function that has peak at middle knot[2]

  B = Ns[idx] * Nt[idx];
  Bs = Nds[idx] * Nt[idx];
  Bt = Ns[idx] * Ndt[idx];

  if (clamped) {
    std::vector<double> Nsc, Ntc;
    bool cl_left(false), cl_right(false), cl_bottom(false), cl_top(false);

    //if (vext.param.x() == param_min.x() && s < vext.s[3]) // left
    if (equal(vext.param.x(), param_min.x()) && smaller(s, vext.s[3])) // left
    {
      vs = vext.s;
      vs[4] = vs[3];
      vs[3] = vs[2];
      tspline::cox(s, degree, vs, Nsc);
      B += Nsc[idx] * Nt[idx];
      Bt += Nsc[idx] * Ndt[idx];
      cl_left = true;
    } //else if (vext.param.x() == param_max.x() && s > vext.s[1]) // right
    else if (equal(vext.param.x(), param_max.x()) && greater(s, vext.s[1])) // right
    {
      vs = vext.s;
      vs[0] = vs[1];
      vs[1] = vs[2];
      tspline::cox(s, degree, vs, Nsc);
      B += Nsc[idx] * Nt[idx];
      Bt += Nsc[idx] * Ndt[idx];
      cl_right = true;
    }

    //if (vext.param.y() == param_min.y() && t < vext.t[3]) // bottom
    if (equal(vext.param.y(), param_min.y()) && smaller(t, vext.t[3])) // bottom
    {
      vt = vext.t;
      vt[4] = vt[3];
      vt[3] = vt[2];
      tspline::cox(t, degree, vt, Ntc);
      B += Ns[idx] * Ntc[idx];
      Bs += Nds[idx] * Ntc[idx];
      cl_bottom = true;
    } //else if (vext.param.y() == param_max.y() && t > vext.t[1]) // top
    else if (equal(vext.param.y(), param_max.y()) && greater(t, vext.t[1])) // top
    {
      vt = vext.t;
      vt[0] = vt[1];
      vt[1] = vt[2];
      tspline::cox(t, degree, vt, Ntc);
      B += Ns[idx] * Ntc[idx];
      Bs += Nds[idx] * Ntc[idx];
      cl_top = true;
    }

    if ((cl_left || cl_right) && (cl_bottom || cl_top)) // corners
      B += Nsc[idx] * Ntc[idx];

    //if (s == param_min.x() || s == param_max.x())
    if (equal(s, param_min.x()) || equal(s, param_max.x()))
      Bs = 0.0;
    //if (t == param_min.y() || t == param_max.y())
    if (equal(t, param_min.y()) || equal(t, param_max.y()))
      Bt = 0.0;
  }
}

Eigen::Vector3d Tspline::evaluate(const double &s, const double &t) const
{
  Eigen::Vector3d p(0.0, 0.0, 0.0);
  double B(0.0);
  double b(0.0);

  // get supporting vertices
  for (Vertex_const_iterator v = vertices_begin(); v != vertices_end(); v++)
  {
    const TVertex &vext = v->data();

    if ((greater(s, vext.s[0]) && smaller(s, vext.s[4])) ||
        (equal(s, vext.s[1]) && equal(s, vext.s[0])) ||
        (equal(s, vext.s[3]) && equal(s, vext.s[4])))
    {
      if ((greater(t, vext.t[0]) && smaller(t, vext.t[4])) ||
          (equal(t, vext.t[1]) && equal(t, vext.t[0])) ||
          (equal(t, vext.t[3]) && equal(t, vext.t[4])))
      {
        const Point4d& cpw = vext.GetCP();
        Eigen::Vector3d cp(cpw.x(),cpw.y(),cpw.z());
        const double& w = cpw.w();

        evaluate_basis(s, t, vext, b);
        B += b * w;
        p += cp * b;
      }
    }

  } // for

  p /= B;

  return p;
}

Eigen::Vector3d Tspline::evaluate(const double &s, const double &t, Eigen::Vector3d &ts, Eigen::Vector3d &tt) const
{
  Eigen::Vector3d p(0.0, 0.0, 0.0);
  ts = Eigen::Vector3d(0.0, 0.0, 0.0);
  tt = Eigen::Vector3d(0.0, 0.0, 0.0);
  double B(0.0), Bs(0.0), Bt(0.0);
  double b(0.0), bs(0.0), bt(0.0);

  // get supporting vertices
  for (Vertex_const_iterator v = vertices_begin(); v != vertices_end(); v++)
  {
    const TVertex &vext = v->data();

    if ((greater(s, vext.s[0]) && smaller(s, vext.s[4])) ||
        (equal(s, vext.s[1]) && equal(s, vext.s[0])) ||
        (equal(s, vext.s[3]) && equal(s, vext.s[4])))
    {
      if ((greater(t, vext.t[0]) && smaller(t, vext.t[4])) ||
          (equal(t, vext.t[1]) && equal(t, vext.t[0])) ||
          (equal(t, vext.t[3]) && equal(t, vext.t[4])))
      {
        const Point4d& cpw = vext.GetCP();
        Eigen::Vector3d cp(cpw.x(),cpw.y(),cpw.z());
        const double& w = cpw.w();

        evaluate_basis(s, t, vext, b, bs, bt);
        B += b * w;
        p += cp * b;
        Bs += bs * w;
        Bt += bt * w;
        ts += cp * bs;
        tt += cp * bt;
      }
    }
  } // for

  ts = (ts * B - p * Bs) / (B * B); // quotient rule
  tt = (tt * B - p * Bt) / (B * B); // quotient rule
  p /= B;

  return p;
}

Eigen::Vector3d Tspline::evaluate(const double &s, const double &t, std::vector<tspline::CPbasis> &cpbasis) const
{
  Eigen::Vector3d p(0.0, 0.0, 0.0);
  double B(0.0);
  double b(0.0);

  //    printf ("[Tspline::evaluate] param_min: %f %f\n", param_min.x (), param_min.y ());
  //    printf ("[Tspline::evaluate] param_max: %f %f\n", param_max.x (), param_max.y ());
  //    printf ("[Tspline::evaluate] (%f %f)\n", s, t);

  // get supporting vertices
  for (Vertex_const_iterator v = vertices_begin(); v != vertices_end(); v++)
  {
    const TVertex &vext = v->data();

    if(!quiet && vext.id==8)
    {
      printf("> (%.1f,%.1f) id: %d\n", s, t, vext.id);
      printf(">    s: %f %f %f %f %f\n", vext.s[0], vext.s[1], vext.s[2], vext.s[3], vext.s[4]);
      printf(">    t: %f %f %f %f %f\n", vext.t[0], vext.t[1], vext.t[2], vext.t[3], vext.t[4]);
    }

    if ((greater(s, vext.s[0]) && smaller(s, vext.s[4])) ||
        (equal(s, vext.s[1]) && equal(s, vext.s[0])) ||
        (equal(s, vext.s[3]) && equal(s, vext.s[4])))
    {
      if ((greater(t, vext.t[0]) && smaller(t, vext.t[4])) ||
          (equal(t, vext.t[1]) && equal(t, vext.t[0])) ||
          (equal(t, vext.t[3]) && equal(t, vext.t[4])))
      {
        const Point4d& cpw = vext.GetCP();
        Eigen::Vector3d cp(cpw.x(),cpw.y(),cpw.z());
        const double& w = cpw.w();

        evaluate_basis(s, t, vext, b);
        B += b * w;
        p += cp * b;

        CPbasis cpb;
        cpb.vit = v;
        cpb.b = b;
        cpbasis.push_back(cpb);

      }
    }

  } // for

  p /= B;

  return p;
}

Eigen::Vector3d Tspline::evaluate(const double &s, const double &t,
                                  Eigen::Vector3d &ts, Eigen::Vector3d &tt,
                                  std::vector<tspline::CPbasis> &cpbasis) const
{
  Eigen::Vector3d p(0.0, 0.0, 0.0);
  ts = Eigen::Vector3d(0.0, 0.0, 0.0);
  tt = Eigen::Vector3d(0.0, 0.0, 0.0);
  double B(0.0), Bs(0.0), Bt(0.0);
  double b(0.0), bs(0.0), bt(0.0);

  // get supporting vertices
  for (Vertex_const_iterator v = vertices_begin(); v != vertices_end(); v++)
  {
    const TVertex &vext = v->data();

    if ((greater(s, vext.s[0]) && smaller(s, vext.s[4])) ||
        (equal(s, vext.s[1]) && equal(s, vext.s[0])) ||
        (equal(s, vext.s[3]) && equal(s, vext.s[4])))
    {
      if ((greater(t, vext.t[0]) && smaller(t, vext.t[4])) ||
          (equal(t, vext.t[1]) && equal(t, vext.t[0])) ||
          (equal(t, vext.t[3]) && equal(t, vext.t[4])))
      {
        const Point4d& cpw = vext.GetCP();
        Eigen::Vector3d cp(cpw.x(),cpw.y(),cpw.z());
        const double& w = cpw.w();

        evaluate_basis(s, t, vext, b, bs, bt);
        B += b * w;
        p += cp * b;
        Bs += bs * w;
        Bt += bt * w;
        ts += cp * bs;
        tt += cp * bt;

        CPbasis cpb;
        cpb.vit = v;
        cpb.b = b;
        cpb.bs = bs;
        cpb.bt = bt;
        cpbasis.push_back(cpb);

      }
    }
  } // for

  ts = (ts * B - p * Bs) / (B * B); // quotient rule
  tt = (tt * B - p * Bt) / (B * B); // quotient rule

  p /= B;

  return p;
}

//Eigen::Vector3d Tspline::evaluate(const double &s, const double &t, std::vector<double> &basis) const
//{
//  Eigen::Vector3d p(0.0, 0.0, 0.0);
//  double B(0.0);

//  basis.resize(number_of_vertices(), 0.0);

//  // get supporting vertices
//  size_t i(0);
//  for (Vertex_const_iterator v = vertices_begin(); v != vertices_end(); v++)
//  {
//    const TVertex &vext = v->data();

//    if ((greater(s, vext.s[0]) && smaller(s, vext.s[4])) ||
//        (equal(s, param_min.x()) && equal(s, vext.s[0])) ||
//        (equal(s, param_max.x()) && equal(s, vext.s[4])))
//    {
//      if ((greater(t, vext.t[0]) && smaller(t, vext.t[4])) ||
//          (equal(t, param_min.y()) && equal(t, vext.t[0])) ||
//          (equal(t, param_max.y()) && equal(t, vext.t[4])))
//      {
//        Eigen::Vector3d cp;
//        vext.GetCP(cp);

//        double &b = basis[i];

//        evaluate_basis(s, t, vext, b);
//        B += b;
//        p += (cp * b);

//      }
//    }

//    i++;
//  } // for

//  p /= B;

//  return p;
//}

