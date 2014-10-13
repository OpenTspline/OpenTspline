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

#include "Types.h"
#include "Tspline.h"

using namespace tspline;

TVertex::TVertex() :
  cp(0.0, 0.0, 0.0),
  id(-1), patch_id(-1),
  param(0.0, 0.0), s(5, 0.0), t(5, 0.0),
  updated(false), is_primary(true)
{

}

void TVertex::SetCP(const Point4d& p)
{
  cp = p;
  for(size_t i=0; i<cp_linked.size(); i++)
    cp_linked[i]->cp = cp;
}

void TVertex::SetCP(const Point3d& p)
{
  cp = p;
  for(size_t i=0; i<cp_linked.size(); i++)
    cp_linked[i]->cp = cp;
}

void TVertex::SetCP(const Eigen::Vector4d& p)
{
  cp = Point4d(p(0),p(1),p(2),p(3));
  for(size_t i=0; i<cp_linked.size(); i++)
    cp_linked[i]->cp = cp;
}

void TVertex::SetCP(const Eigen::Vector3d& p)
{
  cp = Point3d(p(0),p(1),p(2));
  for(size_t i=0; i<cp_linked.size(); i++)
    cp_linked[i]->cp = cp;
}

const Point4d& TVertex::GetCP() const
{
  return cp;
}

void TVertex::GetCP(Eigen::Vector4d &p) const
{
  p(0) = cp.x();
  p(1) = cp.y();
  p(2) = cp.z();
  p(3) = cp.w();
}

void TVertex::GetCP(Eigen::Vector3d &p) const
{
  p(0) = cp.x();
  p(1) = cp.y();
  p(2) = cp.z();
}

void TVertex::LinkCP(TVertex *vext)
{
  if(vext->patch_id==patch_id && vext->id==id)
    return;

  if(IsLinked(vext->patch_id, vext->id))
    return;

  cp_linked.push_back(vext);
  vext->cp_linked.push_back(this);
  cp = vext->GetCP();

  // link new vertex with vertices in cp_linked
  for(size_t i=0; i<cp_linked.size()-1; i++)
    cp_linked[i]->LinkCP(vext);

  // link this vertex with vertices in vext->cp_linked
  for(size_t i=0; i<vext->cp_linked.size()-1; i++)
    vext->cp_linked[i]->LinkCP(this);

  // determine dominant control point
  ComputeDominant(patch_id, id);
}

// computes the vertex with the lowest pid,vid and sets it dominant
TVertex* TVertex::ComputeDominant(int pid, int vid)
{
  is_primary = false;
  TVertex* vdom = this;
  for(size_t i=0; i<cp_linked.size(); i++)
  {
    TVertex* v2 = cp_linked[i];
    if(v2->patch_id < pid ||
       (v2->patch_id == pid && v2->id < vid))
    {
      pid = v2->patch_id;
      vid = v2->id;
      vdom = v2;
    }
    v2->is_primary = false;
  }
  vdom->is_primary = true;
  return vdom;
}

void TVertex::UnLink()
{
  if(cp_linked.empty())
    return;

  // remove this vertex from all link-lists of vertices from this->cp_linked
  for(size_t i=0; i<cp_linked.size(); i++)
  {
    TVertex* vext = cp_linked[i];
    std::vector<TVertex*>::iterator vit;
    for(vit=vext->cp_linked.begin(); vit!=vext->cp_linked.end(); vit++)
    {
      if((*vit)->patch_id == patch_id && (*vit)->id == id)
      {
        vext->cp_linked.erase(vit);
        break;
      }
    }
  }

  // if this one was dominant, recompute it from cp_linked
  if(is_primary)
    ComputeDominant(INT_MAX, INT_MAX);
  is_primary = true; // since this vertex is not linked any more it is dominant

  cp_linked.clear();  
}

bool TVertex::IsLinked() const
{
  return !cp_linked.empty();
}

bool TVertex::IsLinked(int pid, int vid) const
{
  for(size_t i=0; i<cp_linked.size(); i++)
    if(cp_linked[i]->patch_id==pid && cp_linked[i]->id==vid)
      return true;
  return false;
}

bool TVertex::IsIn(std::vector<const TVertex*>& links) const
{
  for(size_t i=0; i<links.size(); i++)
    if(patch_id==links[i]->patch_id && id==links[i]->id)
      return true;

  return false;
}

bool TVertex::IsLinked(std::vector<const TVertex*>& links) const
{
  if(cp_linked.empty())
    return false;

  for(size_t i=0; i<links.size(); i++)
    if(links[i]->IsLinked(patch_id, id))
      return true;

  return false;
}

size_t TVertex::GetNumLinks() const
{
  return cp_linked.size();
}

void TVertex::GetLinks(std::vector<int> &_patch_ids) const
{
  for(size_t i=0; i<cp_linked.size(); i++)
    _patch_ids.push_back(cp_linked[i]->patch_id);
}
void TVertex::GetLinks(std::vector<int> &_patch_ids, std::vector<int> &_vertex_ids) const
{
  for(size_t i=0; i<cp_linked.size(); i++)
  {
    _patch_ids.push_back(cp_linked[i]->patch_id);
    _vertex_ids.push_back(cp_linked[i]->id);
  }
}

int TVertex::GetPrimaryID() const
{
  int pid, vid;
  GetPrimaryID(pid, vid);
  return vid;
}

void TVertex::GetPrimaryID(int& pid, int& vid) const
{
  pid = patch_id;
  vid = id;

  for(size_t i=0; i<cp_linked.size(); i++)
  {
    TVertex* v2 = cp_linked[i];
    if(v2->is_primary)
    {
      pid = patch_id;
      vid = id;
      return;
    }
  }
}

void TVertex::PrintKnotVectors() const
{
  printf("  s: %f %f %f %f %f\n", s[0], s[1], s[2], s[3], s[4]);
  printf("  t: %f %f %f %f %f\n", t[0], t[1], t[2], t[3], t[4]);
}

BasisFunction::~BasisFunction()
{
  if(B[0]!=NULL)
    delete B[0];
  if(B[1]!=NULL)
    delete B[1];
  B[0] = NULL;
  B[1] = NULL;
}

int BasisFunction::id()
{
  return vit->data().id;
}

void BasisFunction::split(const std::vector<double> &N, const double &k,
                          std::vector<double> &N0, std::vector<double> &N1,
                          double& c0, double& c1)
{
  if(sequal(k,N[0]) || gequal(k,N[4]))
    throw std::runtime_error("[BasisFunction::split] Error, k out of bounds.");

  if(equal(k,N[1]) || equal(k,N[2]) || equal(k,N[3]))
  {
    printf("[BasisFunction::split] k: %f  N[1]: %f  N[2]: %f  N[3]: %f\n", k, N[1], N[2], N[3]);
    throw std::runtime_error("[BasisFunction::split] Error, k lies on existing knot.");
  }

  N0.assign(5,0.0);
  N1.assign(5,0.0);

  if(greater(k,N[0]) && smaller(k,N[1])) // Eq. (6)
  {
    N0[0] = N[0];
    N0[1] = k;
    N0[2] = N[1];
    N0[3] = N[2];
    N0[4] = N[3];

    N1[0] = k;
    N1[1] = N[1];
    N1[2] = N[2];
    N1[3] = N[3];
    N1[4] = N[4];

    c0 = (k-N[0]) / (N[3]-N[0]);
    c1 = 1.0;
  }

  else if(greater(k,N[1]) && smaller(k,N[2])) // Eq. (7)
  {
    N0[0] = N[0];
    N0[1] = N[1];
    N0[2] = k;
    N0[3] = N[2];
    N0[4] = N[3];

    N1[0] = N[1];
    N1[1] = k;
    N1[2] = N[2];
    N1[3] = N[3];
    N1[4] = N[4];

    c0 = (k-N[0]) / (N[3]-N[0]);

    if(N[2]==N[3])
      c1 = 1.0; // double knot (e.g. at the boundaries of the T-spline)
    else
      c1 = (N[4]-k) / (N[4]-N[1]);
  }

  else if(greater(k,N[2]) && smaller(k,N[3])) // Eq. (8)
  {
    N0[0] = N[0];
    N0[1] = N[1];
    N0[2] = N[2];
    N0[3] = k;
    N0[4] = N[3];

    N1[0] = N[1];
    N1[1] = N[2];
    N1[2] = k;
    N1[3] = N[3];
    N1[4] = N[4];

    if(N[1]==N[2])
      c0 = 1.0; // double knot (e.g. at the boundaries of the T-spline)
    else
      c0 = (k-N[0]) / (N[3]-N[0]);

    c1 = (N[4]-k) / (N[4]-N[1]);
  }

  else if(greater(k,N[3]) && smaller(k,N[4])) // Eq. (8)
  {
    N0[0] = N[0];
    N0[1] = N[1];
    N0[2] = N[2];
    N0[3] = N[3];
    N0[4] = k;

    N1[0] = N[1];
    N1[1] = N[2];
    N1[2] = N[3];
    N1[3] = k;
    N1[4] = N[4];

    c0 = 1.0;
    c1 = (N[4]-k) / (N[4]-N[1]);
  }

  else throw std::runtime_error("[BasisFunction::split] Error, unhandled case.");
}

bool BasisFunction::split_s(const double& k)
{
//  printf("[BasisFunction::split_s(%f)] id: %d  is_basis: %d\n", k, id(), is_basis_function());
//  printf("  s:  %f %f %f %f %f\n", s[0], s[1], s[2], s[3], s[4]);
//  printf("  t:  %f %f %f %f %f\n", t[0], t[1], t[2], t[3], t[4]);

  B[0] = new BlendingFunction(tsp, this);
  B[1] = new BlendingFunction(tsp, this);

  B[0]->t = t;
  B[1]->t = t;

  double c0, c1;
  split(s, k, B[0]->s, B[1]->s, c0, c1);

  if(B[0]->s[2] == s[2])
    B[0]->vit = vit;
  else
  {
    Tspline::Halfedge_iterator he = tsp->get_left_halfedge(vit);
    if(he==tsp->halfedges_end())
    {
//      printf("[BasisFunction::split_s] left\n");
      Tspline::Face_iterator f = tsp->get_bottom_halfedge(vit)->face();
      tsp->split_horizontal(f, t[2]);
      he = tsp->get_left_halfedge(vit);
    }
    B[0]->vit = he->source();
  }
  B[0]->scale.c = c0;
  B[0]->scale.i = vit->data().id;
  B[0]->scale.j = B[0]->vit->data().id;

  if(B[1]->s[2] == s[2])
    B[1]->vit = vit;
  else
  {
    Tspline::Halfedge_iterator he = tsp->get_right_halfedge(vit);
    if(he==tsp->halfedges_end())
    {
//      printf("[BasisFunction::split_s] right\n");
      Tspline::Face_iterator f = tsp->get_top_halfedge(vit)->face();
      tsp->split_horizontal(f, t[2]);
      he = tsp->get_right_halfedge(vit);
    }
    B[1]->vit = he->source();
  }
  B[1]->scale.c = c1;
  B[1]->scale.i = vit->data().id;
  B[1]->scale.j = B[1]->vit->data().id;

  is_split = true;

//  printf("  s0: %f %f %f %f %f\n", B[0]->s[0], B[0]->s[1], B[0]->s[2], B[0]->s[3], B[0]->s[4]);
//  printf("  s1: %f %f %f %f %f\n", B[1]->s[0], B[1]->s[1], B[1]->s[2], B[1]->s[3], B[1]->s[4]);
//  printf("[BasisFunction::split_s(%f)] id: %d(%f) %d(%f)\n",
//         k, B[0]->vit->data().id, B[0]->scale.c, B[1]->vit->data().id, B[1]->scale.c);
//  std::cout << std::endl;

  return true;
}

bool BasisFunction::split_t(const double& k)
{
//  printf("[BasisFunction::split_t(%f)] id: %d  is_basis: %d\n", k, id(), is_basis_function());
//  printf("  s:  %f %f %f %f %f\n", s[0], s[1], s[2], s[3], s[4]);
//  printf("  t:  %f %f %f %f %f\n", t[0], t[1], t[2], t[3], t[4]);

  B[0] = new BlendingFunction(tsp, this);
  B[1] = new BlendingFunction(tsp, this);

  B[0]->s = s;
  B[1]->s = s;

  double c0, c1;
  split(t, k, B[0]->t, B[1]->t, c0, c1);

  if(B[0]->t[2] == t[2])
    B[0]->vit = vit;
  else
  {
    Tspline::Halfedge_iterator he = tsp->get_bottom_halfedge(vit);
    if(he==tsp->halfedges_end())
    {
//      printf("[BasisFunction::split_t] bottom\n");
      Tspline::Face_iterator f = tsp->get_right_halfedge(vit)->face();
      tsp->split_vertical(f, s[2]);
      he = tsp->get_bottom_halfedge(vit);
    }
    B[0]->vit = he->source();
  }
  B[0]->scale.c = c0;
  B[0]->scale.i = vit->data().id;
  B[0]->scale.j = B[0]->vit->data().id;

  if(B[1]->t[2] == t[2])
    B[1]->vit = vit;
  else
  {
    Tspline::Halfedge_iterator he = tsp->get_top_halfedge(vit);
    if(he==tsp->halfedges_end())
    {
//      printf("[BasisFunction::split_t] top\n");
      Tspline::Face_iterator f = tsp->get_left_halfedge(vit)->face();
      tsp->split_vertical(f, s[2]);
      he = tsp->get_top_halfedge(vit);
    }
    B[1]->vit = he->source();
  }
  B[1]->scale.c = c1;
  B[1]->scale.i = vit->data().id;
  B[1]->scale.j = B[1]->vit->data().id;

  is_split = true;

//  printf("  t0: %f %f %f %f %f\n", B[0]->t[0], B[0]->t[1], B[0]->t[2], B[0]->t[3], B[0]->t[4]);
//  printf("  t1: %f %f %f %f %f\n", B[1]->t[0], B[1]->t[1], B[1]->t[2], B[1]->t[3], B[1]->t[4]);
//  printf("[BasisFunction::split_t(%f)] id: %d(%f) %d(%f)\n",
//         k, B[0]->vit->data().id, B[0]->scale.c, B[1]->vit->data().id, B[1]->scale.c);
//  std::cout << std::endl;

  return true;
}

bool BasisFunction::compute_violation_1()
{
  if(is_split)
  {
    bool b0 = B[0]->compute_violation_1();
    bool b1 = B[1]->compute_violation_1();
    return (b0 || b1);
  }
  return false;
}

bool BasisFunction::compute_violation_2()
{
  if(is_split)
  {
    bool b0 = B[0]->compute_violation_2();
    bool b1 = B[1]->compute_violation_2();
    return (b0 || b1);
  }
  return false;
}

void BasisFunction::get_scale(std::vector<Scale>& g_scale)
{
  B[0]->get_scale(g_scale);
  B[1]->get_scale(g_scale);
}

namespace tspline{
std::vector<BasisFunction*>::iterator find(std::vector<BasisFunction*>& f, int id)
{
  std::vector<BasisFunction*>::iterator it;
  for(it=f.begin(); it!=f.end(); it++)
  {
    if((*it)->id() == id)
      return it;
  }
  return f.end();
}

Scale operator*(const Scale& s1, const Scale& s2)
{
  Scale result;

  result.c = s1.c * s2.c;
  result.i = s1.i;
  result.j = s2.j;

  return result;
}
} // namespace tspline

bool BlendingFunction::compute_violation_1()
{
  if(is_split)
  {
    bool b0 = B[0]->compute_violation_1();
    bool b1 = B[1]->compute_violation_1();
    return (b0 || b1);
  }
  else
  {
    std::vector<double> knots;

    // left
    tsp->shoot_left(vit, 2, knots);
    if(smaller(s[1],knots[0]))
    {
//      printf("[BlendingFunction::compute_violation_1] A id: %d parent: %d\n", id(), parent->id());
      return split_s(knots[0]);
    }
    else if(smaller(s[0],knots[1]))
    {
//      printf("[BlendingFunction::compute_violation_1] B id: %d parent: %d\n", id(), parent->id());
      return split_s(knots[1]);
    }

    // right
    tsp->shoot_right(vit, 2, knots);
    if(greater(s[3],knots[0]))
    {
//      printf("[BlendingFunction::compute_violation_1] C id: %d parent: %d\n", id(), parent->id());
      return split_s(knots[0]);
    }
    else if(greater(s[4], knots[1]))
    {
//      printf("[BlendingFunction::compute_violation_1] D id: %d parent: %d\n", id(), parent->id());
      return split_s(knots[1]);
    }

    // down
    tsp->shoot_down(vit, 2, knots);
    if(smaller(t[1],knots[0]))
    {
//      printf("[BlendingFunction::compute_violation_1] E id: %d parent: %d\n", id(), parent->id());
      return split_t(knots[0]);
    }
    else if(smaller(t[0],knots[1]))
    {
//      printf("[BlendingFunction::compute_violation_1] F id: %d parent: %d\n", id(), parent->id());
      return split_t(knots[1]);
    }

    // up
    tsp->shoot_up(vit, 2, knots);
    if(greater(t[3],knots[0]))
    {
//      printf("[BlendingFunction::compute_violation_1] G id: %d parent: %d\n", id(), parent->id());
      return split_t(knots[0]);
    }
    else if(greater(t[4], knots[1]))
    {
//      printf("[BlendingFunction::compute_violation_1] H id: %d parent: %d\n", id(), parent->id());
      return split_t(knots[1]);
    }

    return false;
  }

}

bool BlendingFunction::compute_violation_2()
{
  if(is_split)
  {
    bool b0 = B[0]->compute_violation_2();
    bool b1 = B[1]->compute_violation_2();
    return (b0 || b1);
  }
  else
  {
    double x,y;
    bool violation(false);
    std::vector<double> knots;
    Tspline::Face_iterator f;

    // left
    tsp->shoot_left(vit, 2, knots);
    if(greater(s[1],knots[0]))
    {
//      printf("[BlendingFunction::compute_violation_2] Violation_2 (%f). A id: %d parent: %d\n",
//             s[1], id(), parent->id());

      if(parent->t[2] > t[2])
        f = tsp->get_top_halfedge(vit)->twin()->face();
      else
        f = tsp->get_bottom_halfedge(vit)->face();
      tsp->split_vertical(f, s[1]);
      violation = true;
    }
    else if(greater(s[0],knots[1]))
    {
//      printf("[BlendingFunction::compute_violation_2] Violation_2 (%f). B id: %d parent: %d\n",
//             s[0], id(), parent->id());
      Tspline::Face_iterator f;
      if(parent->t[2] > t[2])
        f = tsp->get_top_halfedge(vit)->twin()->face();
      else
        f = tsp->get_bottom_halfedge(vit)->face();
      Tspline::Halfedge_iterator h = tsp->get_left_halfedge_at_param(f, t[2], x)->twin();
      f = h->face();
      tsp->split_vertical(f, s[0]);
      violation = true;
    }

    // right
    tsp->shoot_right(vit, 2, knots);
    if(smaller(s[3],knots[0]))
    {
//      printf("[BlendingFunction::compute_violation_2] Violation_2 (%f). C id: %d parent: %d\n",
//             s[3], id(), parent->id());
      Tspline::Face_iterator f;
      if(parent->t[2] > t[2])
        f = tsp->get_top_halfedge(vit)->face();
      else
        f = tsp->get_bottom_halfedge(vit)->twin()->face();
      tsp->split_vertical(f, s[3]);
      violation = true;
    }
    else if(smaller(s[4], knots[1]))
    {
//      printf("[BlendingFunction::compute_violation_2] Violation_2 (%f). D id: %d parent: %d\n",
//             s[4], id(), parent->id());
      Tspline::Face_iterator f;
      if(parent->t[2] > t[2])
        f = tsp->get_top_halfedge(vit)->face();
      else
        f = tsp->get_bottom_halfedge(vit)->twin()->face();
      Tspline::Halfedge_iterator h = tsp->get_right_halfedge_at_param(f, t[2], x)->twin();
      f = h->face();
      tsp->split_vertical(f, s[4]);
      violation = true;
    }

    // down
    tsp->shoot_down(vit, 2, knots);
    if(greater(t[1],knots[0]))
    {
//      printf("[BlendingFunction::compute_violation_2] Violation_2 (%f). E id: %d parent: %d\n",
//             t[1], id(), parent->id());
      Tspline::Face_iterator f;
      if(parent->s[2] > s[2])
        f = tsp->get_right_halfedge(vit)->face();
      else
        f = tsp->get_left_halfedge(vit)->twin()->face();
      tsp->split_horizontal(f, t[1]);
      violation = true;
    }
    else if(greater(t[0],knots[1]))
    {
//      printf("[BlendingFunction::compute_violation_2] Violation_2 (%f). F id: %d parent: %d\n",
//             t[0], id(), parent->id());
      Tspline::Face_iterator f;
      if(parent->s[2] > s[2])
        f = tsp->get_right_halfedge(vit)->face();
      else
        f = tsp->get_left_halfedge(vit)->twin()->face();
      Tspline::Halfedge_iterator h = tsp->get_bottom_halfedge_at_param(f, s[2], y)->twin();
      f = h->face();
      tsp->split_horizontal(f, t[0]);
      violation = true;
    }

    // up
    tsp->shoot_up(vit, 2, knots);
    if(smaller(t[3],knots[0]))
    {
//      printf("[BlendingFunction::compute_violation_2] Violation_2 (%f). G id: %d parent: %d\n",
//             t[3], id(), parent->id());
      Tspline::Face_iterator f;
      if(parent->s[2] > s[2])
        f = tsp->get_right_halfedge(vit)->twin()->face();
      else
        f = tsp->get_left_halfedge(vit)->face();
      tsp->split_horizontal(f, t[3]);
      violation = true;
    }
    else if(smaller(t[4], knots[1]))
    {
//      printf("[BlendingFunction::compute_violation_2] Violation_2 (%f). H id: %d parent: %d\n",
//             t[4], id(), parent->id());
      Tspline::Face_iterator f;
      if(parent->s[2] > s[2])
        f = tsp->get_right_halfedge(vit)->twin()->face();
      else
        f = tsp->get_left_halfedge(vit)->face();
      Tspline::Halfedge_iterator h = tsp->get_top_halfedge_at_param(f, s[2], y)->twin();
      f = h->face();
      tsp->split_horizontal(f, t[4]);
      violation = true;
    }

    return violation;
  }
}

void BlendingFunction::get_scale(std::vector<Scale>& g_scale)
{
  if(is_split)
  {
//    printf("  c_%d^%d(%f) * ", scale.i, scale.j, scale.c);
    B[0]->get_scale(g_scale, scale);
//    printf("  c_%d^%d(%f) * ", scale.i, scale.j, scale.c);
    B[1]->get_scale(g_scale, scale);
  }
  else
  {
//    printf("  c_%d^%d(%f)\n", scale.i, scale.j, scale.c);
    g_scale.push_back(scale);
  }
}

void BlendingFunction::get_scale(std::vector<Scale>& g_scale, Scale parent_scale)
{
  Scale s = parent_scale * scale;
  if(is_split)
  {
//    printf("c_%d^%d(%f) * ", scale.i, scale.j, scale.c);
    B[0]->get_scale(g_scale, s);
//    printf("c_%d^%d(%f) * ", scale.i, scale.j, scale.c);
    B[1]->get_scale(g_scale, s);
  }
  else
  {
//    printf("c_%d^%d(%f) = c_%d^%d(%f)\n", scale.i, scale.j, scale.c, s.i, s.j, s.c);
    g_scale.push_back(s);
  }
}
