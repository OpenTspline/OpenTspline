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

#include "TsplineMultiPatch.h"
#include "TsplineCreator.h"

using namespace tspline;

TsplinePatch::~TsplinePatch()
{
  for(size_t i=0; i<m_north.size(); i++)
    delete m_north[i];
  for(size_t i=0; i<m_east.size(); i++)
    delete m_east[i];
  for(size_t i=0; i<m_south.size(); i++)
    delete m_south[i];
  for(size_t i=0; i<m_west.size(); i++)
    delete m_west[i];
}

Tspline::Vertex_iterator TsplinePatch::insert_vertex_no_transition(const double &s, const double &t)
{
  Point2d param(s, t);

  CGAL::Object obj = locate_param(s, t);

  Face_iterator f;
  Halfedge_iterator e;
  Vertex_iterator v;

  if (!obj.empty()) {
    if (CGAL::assign(f, obj))
      // this cannot be done because it would subsequently introduce new vertices which would not be labelled (as below)
      throw std::runtime_error("[TsplinePatch::insert_vertex_no_transition] Error, cannot insert_vertex on face in TsplinePatch. "
                               "Use TsplineMultiPatch::InsertVertex instead.");
    else if (CGAL::assign(e, obj))
      // calls vertex insertion on an edge, without taking care about transition (necessary when creating transitions)
      v=Tspline::insert_vertex(e, param);
    else if (CGAL::assign(v, obj))
      ;//printf("[TsplinePatch::insert_vertex] Warning, vertex already existing at param position (%e %e)\n", s, t);
    else
      throw std::runtime_error("[TsplinePatch::insert_vertex] Error, no object found");

    TVertex& vext = v->data();
    vext.patch_id = id;

    return v;
  }

  return vertices_end();
}

Tspline::Vertex_iterator TsplinePatch::insert_vertex(Halfedge_iterator &e, Point2d &param)
{
  Vertex_iterator vn = Tspline::insert_vertex(e, param);

  // insert vertex at adjacent patch
  Vertex_const_iterator vnc = vn;
  Halfedge_const_iterator hec;
  if(this->is_boundary(vnc, hec))
  {
    Transition* trans;
    const TVertex& vext0 = hec->source()->data();
    const TVertex& vext1 = hec->target()->data();
    if(equal(vext0.param.x(), vext1.param.x()))
    {
      if(smaller(vext0.param.y(), vext1.param.y())){ // WEST
        //        if(m_west.Valid()){
        trans = GetWest(param.y());
        if(trans){
          Point2d p1 = trans->ConvertParam(param);
          trans->target->insert_vertex_no_transition(p1.x(), p1.y());
        }
      }else{ // EAST
        //        if(m_east.Valid()){
        trans = GetEast(param.y());
        if(trans){
          Point2d p1 = trans->ConvertParam(param);
          trans->target->insert_vertex_no_transition(p1.x(), p1.y());
        }
      }
    }
    if(equal(vext0.param.y(), vext1.param.y()))
    {
      if(smaller(vext0.param.x(), vext1.param.x())){ // NORTH
        //        if(m_north.Valid()){
        trans = GetNorth(param.x());
        if(trans){
          Point2d p1 = trans->ConvertParam(param);
          trans->target->insert_vertex_no_transition(p1.x(), p1.y());
        }
      }else{ // SOUTH
        //        if(m_south.Valid()){
        trans = GetSouth(param.x());
        if(trans){
          Point2d p1 = trans->ConvertParam(param);
          trans->target->insert_vertex_no_transition(p1.x(), p1.y());
        }
      }
    }

  } // if(this->is_boundary(vnc, hec))

  return vn;
}

bool  TsplinePatch::IsNorth(double s) const
{
  for(size_t i=0; i<m_north.size(); i++)
    if(m_north[i]->InSourceRange(s))
      return true;
  return false;
}

bool TsplinePatch::IsEast(double t) const
{
  for(size_t i=0; i<m_east.size(); i++)
    if(m_east[i]->InSourceRange(t))
      return true;
  return false;
}

bool TsplinePatch::IsSouth(double s) const
{
  for(size_t i=0; i<m_south.size(); i++)
    if(m_south[i]->InSourceRange(s))
      return true;
  return false;
}

bool TsplinePatch::IsWest(double t) const
{
  for(size_t i=0; i<m_west.size(); i++)
    if(m_west[i]->InSourceRange(t))
      return true;
  return false;
}

Transition* TsplinePatch::GetNorth(double s)
{
  for(size_t i=0; i<m_north.size(); i++)
    if(m_north[i]->InSourceRange(s))
      return m_north[i];
  return 0;
}

Transition *TsplinePatch::GetEast(double t)
{
  for(size_t i=0; i<m_east.size(); i++)
    if(m_east[i]->InSourceRange(t))
      return m_east[i];
  return 0;
}

Transition *TsplinePatch::GetSouth(double s)
{
  for(size_t i=0; i<m_south.size(); i++)
    if(m_south[i]->InSourceRange(s))
      return m_south[i];
  return 0;
}

Transition *TsplinePatch::GetWest(double t)
{
  for(size_t i=0; i<m_west.size(); i++)
    if(m_west[i]->InSourceRange(t))
      return m_west[i];
  return 0;
}

const Transition* TsplinePatch::GetNorthConst(double s) const
{
  for(size_t i=0; i<m_north.size(); i++)
    if(m_north[i]->InSourceRange(s))
      return m_north[i];
  return 0;
}

const Transition *TsplinePatch::GetEastConst(double t) const
{
  for(size_t i=0; i<m_east.size(); i++)
    if(m_east[i]->InSourceRange(t))
      return m_east[i];
  return 0;
}

const Transition *TsplinePatch::GetSouthConst(double s) const
{
  for(size_t i=0; i<m_south.size(); i++)
    if(m_south[i]->InSourceRange(s))
      return m_south[i];
  return 0;
}

const Transition *TsplinePatch::GetWestConst(double t) const
{
  for(size_t i=0; i<m_west.size(); i++)
    if(m_west[i]->InSourceRange(t))
      return m_west[i];
  return 0;
}

bool TsplinePatch::HasTransition(Direction source_dir, int target_id)
{
  if(source_dir==NORTH)
  {
    for(size_t i=0; i<m_north.size(); i++)
      if(m_north[i]->target->id==target_id)
        return true;
    return false;
  }

  if(source_dir==EAST)
  {
    for(size_t i=0; i<m_east.size(); i++)
      if(m_east[i]->target->id==target_id)
        return true;
    return false;
  }

  if(source_dir==SOUTH)
  {
    for(size_t i=0; i<m_south.size(); i++)
      if(m_south[i]->target->id==target_id)
        return true;
    return false;
  }

  if(source_dir==WEST)
  {
    for(size_t i=0; i<m_west.size(); i++)
      if(m_west[i]->target->id==target_id)
        return true;
    return false;
  }

  return false;
}

bool TsplinePatch::HasTransition(int target_id)
{
  for(size_t i=0; i<m_north.size(); i++)
    if(m_north[i]->target->id==target_id)
      return true;

  for(size_t i=0; i<m_east.size(); i++)
    if(m_east[i]->target->id==target_id)
      return true;

  for(size_t i=0; i<m_south.size(); i++)
    if(m_south[i]->target->id==target_id)
      return true;

  for(size_t i=0; i<m_west.size(); i++)
    if(m_west[i]->target->id==target_id)
      return true;

  return false;
}


Transition*  TsplinePatch::AddTransition(Direction source_dir, TsplinePatch* target, Direction target_dir,
                                         unsigned char mult)
{
  if(target==NULL)
    throw std::runtime_error("[TsplinePatch::AddTransition] Error, null pointer exception");

  if(source_dir==NORTH && m_north.empty()){
    m_north.push_back(new NorthTransition(this, source_dir, target, target_dir, mult));
    return m_north.back();
  }
  if(source_dir==EAST && m_east.empty()){
    this->m_east.push_back(new EastTransition(this, source_dir, target, target_dir, mult));
    return m_east.back();
  }
  if(source_dir==SOUTH && m_south.empty()){
    this->m_south.push_back(new SouthTransition(this, source_dir, target, target_dir, mult));
    return m_south.back();
  }
  if(source_dir==WEST && m_west.empty()){
    this->m_west.push_back(new WestTransition(this, source_dir, target, target_dir, mult));
    return m_west.back();
  }

  printf("[TsplinePatch::AddTransition] Error, transition could not be added id0: %d dir0: %d  id1: %d dir1: %d\n",
         id, source_dir, target->id, target_dir);
  throw std::runtime_error("[TsplinePatch::AddTransition] Error, transition could not be added");
}

Transition* TsplinePatch::AddTransition(Direction d0, double r0a, double r0b,
                                        TsplinePatch* p1, Direction d1, double r1a, double r1b,
                                        unsigned char mult)
{
  if(p1==NULL)
    throw std::runtime_error("[TsplinePatch::AddTransition] Error, null pointer exception");

  if(d0==NORTH){
    m_north.push_back(new NorthTransition(this, d0, r0a, r0b, p1, d1, r1a, r1b, mult));
    return m_north.back();
  }
  if(d0==EAST){
    m_east.push_back(new EastTransition(this, d0, r0a, r0b, p1, d1, r1a, r1b, mult));
    return m_east.back();
  }
  if(d0==SOUTH){
    m_south.push_back(new SouthTransition(this, d0, r0a, r0b, p1, d1, r1a, r1b, mult));
    return m_south.back();
  }
  if(d0==WEST){
    m_west.push_back(new WestTransition(this, d0, r0a, r0b, p1, d1, r1a, r1b, mult));
    return m_west.back();
  }

  printf("[TsplinePatch::AddTransition] Error, transition could not be added "
         "id0: %d dir0: %d r: (%f, %f) "
         "id1: %d dir1: %d r: (%f, %f)\n",
         id, d0, r0a, r0b, p1->id, d1, r1a, r1b);
  throw std::runtime_error("[TsplinePatch::AddTransition] Error, transition could not be added");
}

//void TsplinePatch::evaluate_basis(const double &s, const double &t, const TVertex &vext, double &B) const
//{
//  std::vector<double> Ns, Nt, vs, vt;
//  tspline::cox(s, degree, vext.s, Ns);
//  tspline::cox(t, degree, vext.t, Nt);

//  size_t idx = vext.s.size() * degree;

//  B = Ns[idx] * Nt[idx];

//  //  if (clamped) {
//  std::vector<double> Nsc, Ntc;
//  bool cl_left(false), cl_right(false), cl_bottom(false), cl_top(false);
//  //  if (!m_west.Valid() && equal(vext.param.x(), param_min.x()) && smaller(s, vext.s[3])) // left
//  if (equal(vext.param.x(), param_min.x()) && !IsWest(t) && smaller(s, vext.s[3])) // left
//  {
//    vs = vext.s;
//    vs[4] = vs[3];
//    vs[3] = vs[2];
//    tspline::cox(s, degree, vs, Nsc);
//    B += Nsc[idx] * Nt[idx];
//    cl_left = true;
//    //  } else if (!m_east.Valid() && equal(vext.param.x(), param_max.x()) && greater(s, vext.s[1])) // right
//  } else if (equal(vext.param.x(), param_max.x()) && !IsEast(t) && greater(s, vext.s[1])) // right
//  {
//    vs = vext.s;
//    vs[0] = vs[1];
//    vs[1] = vs[2];
//    tspline::cox(s, degree, vs, Nsc);
//    B += Nsc[idx] * Nt[idx];
//    cl_right = true;
//  }

//  //  if (!m_south.Valid() && equal(vext.param.y(), param_min.y()) && smaller(t, vext.t[3])) // bottom
//  if (equal(vext.param.y(), param_min.y()) && !IsSouth(s) && smaller(t, vext.t[3])) // bottom
//  {
//    vt = vext.t;
//    vt[4] = vt[3];
//    vt[3] = vt[2];
//    tspline::cox(t, degree, vt, Ntc);
//    B += Ns[idx] * Ntc[idx];
//    cl_bottom = true;
//    //  } else if (!m_north.Valid(s) && equal(vext.param.y(), param_max.y()) && greater(t, vext.t[1])) // top
//  } else if (equal(vext.param.y(), param_max.y()) && !IsNorth(s) && greater(t, vext.t[1])) // top
//  {
//    vt = vext.t;
//    vt[0] = vt[1];
//    vt[1] = vt[2];
//    tspline::cox(t, degree, vt, Ntc);
//    B += Ns[idx] * Ntc[idx];
//    cl_top = true;
//  }

//  if ((cl_left || cl_right) && (cl_bottom || cl_top)) // corners
//    B += Nsc[idx] * Ntc[idx];

//  //  } // if(clamped)
//}

void TsplinePatch::evaluate_basis(const double &s, const double &t, const TVertex &vext,
                                  CPbasis &basis) const
//                                  double &B, double &Bs, double &Bt) const
{
  std::vector<double> Ns, Nt, Nds, Ndt, vs, vt;
  tspline::cox(s, degree, vext.s, Ns);
  tspline::cox(t, degree, vext.t, Nt);
  tspline::coxder(degree, vext.s, Ns, Nds);
  tspline::coxder(degree, vext.t, Nt, Ndt);

  size_t idx = vext.s.size() * degree;

  basis.b = Ns[idx] * Nt[idx];
  basis.bs = Nds[idx] * Nt[idx];
  basis.bt = Ns[idx] * Ndt[idx];

  //  if (clamped) {
  std::vector<double> Nsc, Ntc;
  bool cl_left(false), cl_right(false), cl_bottom(false), cl_top(false);
  if (m_west.empty() && equal(vext.param.x(), param_min.x()) && smaller(s, vext.s[3])) // left
  {
    vs = vext.s;
    vs[4] = vs[3];
    vs[3] = vs[2];
    tspline::cox(s, degree, vs, Nsc);
    basis.b += Nsc[idx] * Nt[idx];
    basis.bt += Nsc[idx] * Ndt[idx];
    cl_left = true;
  } else if (m_east.empty() && equal(vext.param.x(), param_max.x()) && greater(s, vext.s[1])) // right
  {
    vs = vext.s;
    vs[0] = vs[1];
    vs[1] = vs[2];
    tspline::cox(s, degree, vs, Nsc);
    basis.b += Nsc[idx] * Nt[idx];
    basis.bt += Nsc[idx] * Ndt[idx];
    cl_right = true;
  }

  if (m_south.empty() && equal(vext.param.y(), param_min.y()) && smaller(t, vext.t[3])) // bottom
  {
    vt = vext.t;
    vt[4] = vt[3];
    vt[3] = vt[2];
    tspline::cox(t, degree, vt, Ntc);
    basis.b += Ns[idx] * Ntc[idx];
    basis.bs += Nds[idx] * Ntc[idx];
    cl_bottom = true;
  } else if (m_north.empty() && equal(vext.param.y(), param_max.y()) && greater(t, vext.t[1])) // top
  {
    vt = vext.t;
    vt[0] = vt[1];
    vt[1] = vt[2];
    tspline::cox(t, degree, vt, Ntc);
    basis.b += Ns[idx] * Ntc[idx];
    basis.bs += Nds[idx] * Ntc[idx];
    cl_top = true;
  }

  if ((cl_left || cl_right) && (cl_bottom || cl_top)) // corners
    basis.b += Nsc[idx] * Ntc[idx];

  if ((m_west.empty() && equal(s, param_min.x())) || (m_east.empty() && equal(s, param_max.x())))
    basis.bs = 0.0;
  if ((m_south.empty() && equal(t, param_min.y())) || (m_north.empty() && equal(t, param_max.y())))
    basis.bt = 0.0;
  //  } // if(clamped)
}

void TsplinePatch::evaluate_basis_sum(const double &s, const double &t, BasisSum &basis,
                                      std::vector<const TVertex*>& links,
                                      std::vector<tspline::CPbasis> &cpbasis) const
{
  for (Vertex_const_iterator v = vertices_begin(); v != vertices_end(); v++)
  {
    const TVertex &vext = v->data();

    if ((greater(s, vext.s[0]) && smaller(s, vext.s[4])) ||
        (equal(s, param_min.x()) && equal(vext.s[0], vext.s[2])) ||
        (equal(s, param_max.x()) && equal(vext.s[2], vext.s[4])))
    {
      if ((greater(t, vext.t[0]) && smaller(t, vext.t[4])) ||
          (equal(t, param_min.y()) && equal(vext.t[0], vext.t[2])) ||
          (equal(t, param_max.y()) && equal(vext.t[2], vext.t[4])))
      {
        Eigen::Vector3d cp;
        vext.GetCP(cp);
        CPbasis b;
        evaluate_basis(s, t, vext, b);

        basis.B += b.b;
        basis.cp += (cp * b.b);
        basis.Bs += b.bs;
        basis.Bt += b.bt;
        basis.ts += cp * b.bs;
        basis.tt += cp * b.bt;

        b.vit = v;
        cpbasis.push_back(b);

        links.push_back(&vext);

      }
    }

  } // for
}

void TsplinePatch::evaluate_basis_sum(const double &s, const double &t, BasisSum &basis, const Transition *trans,
                                      std::vector<const TVertex*>& links,
                                      std::vector<tspline::CPbasis> &cpbasis) const
{
  Point2d p(s,t);
  p = trans->ConvertParam(Point2d(s,t));

  Point2d pb;
  for (Vertex_const_iterator v = vertices_begin(); v != vertices_end(); v++)
  {
    const TVertex &vext = v->data();

    if(!vext.IsIn(links) && !vext.IsLinked(links))
    {
      if(!quiet && id==1)
      {
        printf("  [%d,%d] (%.2f,%.2f)\n", id, vext.id, p.x(),p.y());
        printf("    s: %f  |  %f %f %f %f %f\n", p.x(), vext.s[0], vext.s[1], vext.s[2], vext.s[3], vext.s[4]);
        printf("    t: %f  |  %f %f %f %f %f\n", p.y(), vext.t[0], vext.t[1], vext.t[2], vext.t[3], vext.t[4]);

      }


      if ((greater(p.x(), vext.s[0]) && smaller(p.x(), vext.s[4])) ||
          (equal(p.x(), param_min.x()) && equal(vext.s[0], vext.s[2])) ||
          (equal(p.x(), param_max.x()) && equal(vext.s[2], vext.s[4]))) // if(s in range)
      {
        if ((greater(p.y(), vext.t[0]) && smaller(p.y(), vext.t[4])) ||
            (equal(p.y(), param_min.y()) && equal(vext.t[0], vext.t[2])) ||
            (equal(p.y(), param_max.y()) && equal(vext.t[2], vext.t[4])))  // if(t in range)
        {
          Eigen::Vector3d cp;
          vext.GetCP(cp);
          CPbasis b;
          evaluate_basis(p.x(), p.y(), vext, b);
          pb = trans->RotateParamInverse(Point2d(b.bs, b.bt));
          b.bs = pb.x();
          b.bt = pb.y();

          basis.B += b.b;
          basis.cp += (cp * b.b);
          basis.Bs += b.bs;
          basis.Bt += b.bt;
          basis.ts += cp * b.bs;
          basis.tt += cp * b.bt;

          b.vit = v;
          cpbasis.push_back(b);

          links.push_back(&vext);

        } // if(t in range)
      } // if(s in range)
    }// if(!vext.IsLinked())
  } // for
}

Eigen::Vector3d TsplinePatch::evaluate(const double &s, const double &t) const
{
  std::vector<tspline::CPbasis> cpbasis;
  Eigen::Vector3d ts, tt;
  return evaluate(s,t,ts,tt,cpbasis);
}

Eigen::Vector3d TsplinePatch::evaluate(const double &s, const double &t,
                                       Eigen::Vector3d &ts, Eigen::Vector3d &tt) const
{
  std::vector<tspline::CPbasis> cpbasis;
  return evaluate(s,t,ts,tt,cpbasis);
}

Eigen::Vector3d TsplinePatch::evaluate(const double &s, const double &t,
                                       std::vector<tspline::CPbasis> &cpbasis) const
{\
  Eigen::Vector3d ts, tt;
  return evaluate(s,t,ts,tt,cpbasis);
}

Eigen::Vector3d TsplinePatch::evaluate(const double &s, const double &t,
                                       Eigen::Vector3d &ts, Eigen::Vector3d &tt,
                                       std::vector<tspline::CPbasis> &cpbasis) const
{
  BasisSum basis;
  std::vector<const TVertex*> links;

  evaluate_basis_sum(s, t, basis, links, cpbasis);

  if(!m_north.empty() || !m_east.empty() || !m_south.empty() || !m_west.empty())
  {
    // sort linkage according to shortest distance to patch border
    std::vector<LinkSort> linksort;
    linksort.push_back(LinkSort(param_max.y()-t, NORTH));
    linksort.push_back(LinkSort(param_max.x()-s, EAST));
    linksort.push_back(LinkSort(t-param_min.y(), SOUTH));
    linksort.push_back(LinkSort(s-param_min.x(), WEST));

    std::sort(linksort.begin(), linksort.end());

    for(size_t j=0; j<linksort.size(); j++)
    {
      if(linksort[j].dir==NORTH)
        for(size_t i=0; i<m_north.size(); i++)
          m_north[i]->target->evaluate_basis_sum(s, t, basis, m_north[i], links, cpbasis);
      if(linksort[j].dir==EAST)
        for(size_t i=0; i<m_east.size(); i++)
          m_east[i]->target->evaluate_basis_sum(s, t, basis, m_east[i], links, cpbasis);
      if(linksort[j].dir==SOUTH)
        for(size_t i=0; i<m_south.size(); i++)
          m_south[i]->target->evaluate_basis_sum(s, t, basis, m_south[i], links, cpbasis);
      if(linksort[j].dir==WEST)
        for(size_t i=0; i<m_west.size(); i++)
          m_west[i]->target->evaluate_basis_sum(s, t, basis, m_west[i], links, cpbasis);
    }
  }

  ts = (basis.ts * basis.B - basis.cp * basis.Bs) / (basis.B * basis.B); // quotient rule
  tt = (basis.tt * basis.B - basis.cp * basis.Bt) / (basis.B * basis.B); // quotient rule
  basis.cp /= basis.B;

  return basis.cp;
}

// ------------------------------------------------------------------------------
// operators
// ------------------------------------------------------------------------------

//void TsplinePatch::update_knot_vectors()
//{
//  int id(0);
//  Vertex_iterator vit;
//  for (vit = vertices_begin(); vit != vertices_end(); vit++) {
//    Halfedge_iterator he;
//    Vertex_iterator v_curr;
//    TVertex vext = vit->data();

//    vext.id = id++;

//    vext.s[2] = vext.param.x();
//    vext.t[2] = vext.param.y();

//    std::vector<double> knots;

//    // left
//    this->shoot_left(vit, 2, knots);
//    vext.s[0] = knots[1];
//    vext.s[1] = knots[0];

//    // right
//    this->shoot_right(vit, 2, knots);
//    vext.s[3] = knots[0];
//    vext.s[4] = knots[1];

//    // up
//    this->shoot_up(vit, 2, knots);
//    vext.t[3] = knots[0];
//    vext.t[4] = knots[1];

//    // down
//    this->shoot_down(vit, 2, knots);
//    vext.t[1] = knots[0];
//    vext.t[0] = knots[1];

//    vit->set_data(vext);
//  }
//}

void TsplinePatch::shoot_left (Vertex_iterator &v, const unsigned &n, std::vector<double> &knots)
{
  knots.clear ();
  TVertex vext = v->data ();
  double x_curr = vext.param.x ();
  double y_curr = vext.param.y ();
  Vertex_iterator v_curr = v;
  Halfedge_iterator he_v = get_left_halfedge (v_curr);
  Face_iterator f = faces_end ();
  Transition* west;

  while (knots.size () < n)
  {
    if (he_v == halfedges_end ())
    { // find edge on the left of the left face of the vertex
      if (f == faces_end ())
      { // get face left of vertex
        Halfedge_iterator he_top = get_top_halfedge (v_curr);
        if (he_top == halfedges_end ())
        {
          Halfedge_iterator he_bot = get_bottom_halfedge (v_curr);
          if (he_bot == halfedges_end ())
            throw std::runtime_error ("[TsplinePatch::shoot_left] Error, invalid vertex");
          f = he_bot->face ();
        }
        else
          f = he_top->twin ()->face ();
      }

      if (f->is_unbounded ())
      {
        //        if(m_west.Valid())
        west = GetWest(y_curr);
        if(west)
        {
          for(unsigned char m=0; m<west->multiplicity && knots.size()<n; m++)
            knots.push_back(x_curr);
          if(knots.size()<n)
            west->shoot(Point2d(x_curr, y_curr), n-knots.size(), knots);
        }else
          knots.push_back (x_curr); // multiply border knots
      }else
      {
        Halfedge_iterator he_left = get_left_halfedge_at_param (f, y_curr, x_curr);
        knots.push_back (x_curr);
        f = he_left->twin ()->face ();
      }
    }
    else
    { // parse along edges
      v_curr = he_v->source ();
      TVertex vext_curr = v_curr->data ();
      x_curr = vext_curr.param.x ();
      knots.push_back (x_curr);
      he_v = get_left_halfedge (v_curr);
    }
  }
}

void TsplinePatch::shoot_right (Vertex_iterator &v, const unsigned &n, std::vector<double> &knots)
{
  knots.clear ();
  TVertex vext = v->data ();
  double x_curr = vext.param.x ();
  double y_curr = vext.param.y ();
  Vertex_iterator v_curr = v;
  Halfedge_iterator he_v = get_right_halfedge (v_curr);
  Face_iterator f = faces_end ();
  Transition* east;

  while (knots.size () < n)
  {
    if (he_v == halfedges_end ())
    { // find edge on the right of the right face of the vertex
      if (f == faces_end ())
      { // get face right of vertex
        Halfedge_iterator he_top = get_top_halfedge (v_curr);
        if (he_top == halfedges_end ())
        {
          Halfedge_iterator he_bot = get_bottom_halfedge (v_curr);
          if (he_bot == halfedges_end ())
          {
            TVertex vext = v_curr->data();
            printf("[TsplinePatch::shoot_right] param: %f %f\n", vext.param.x(), vext.param.y());
            throw std::runtime_error ("[TsplinePatch::shoot_right] Error, invalid vertex");
          }
          f = he_bot->twin ()->face ();
        }
        else
          f = he_top->face ();
      }

      if (f->is_unbounded ())
      {
        //        if(m_east.Valid())
        east = GetEast(y_curr);
        if(east)
        {
          for(unsigned char m=0; m<east->multiplicity && knots.size()<n; m++)
            knots.push_back(x_curr);
          if(knots.size()<n)
            east->shoot(Point2d(x_curr, y_curr), n-knots.size(), knots);
        }else
          knots.push_back (x_curr); // multiply border knots
      }else
      {
        Halfedge_iterator he_right = get_right_halfedge_at_param (f, y_curr, x_curr);
        knots.push_back (x_curr);
        f = he_right->twin ()->face ();
      }
    }
    else
    { // parse along edges
      v_curr = he_v->source ();
      TVertex vext_curr = v_curr->data ();
      x_curr = vext_curr.param.x ();
      knots.push_back (x_curr);
      he_v = get_right_halfedge (v_curr);
    }
  }
}

void TsplinePatch::shoot_up (Vertex_iterator &v, const unsigned &n, std::vector<double> &knots)
{
  knots.clear ();
  TVertex vext = v->data ();
  double x_curr = vext.param.x ();
  double y_curr = vext.param.y ();
  Vertex_iterator v_curr = v;
  Halfedge_iterator he_v = get_top_halfedge (v_curr);
  Face_iterator f = faces_end ();
  Transition* north;

  while (knots.size () < n)
  {
    if (he_v == halfedges_end ())
    { // find edge on the top of the top face of the vertex
      if (f == faces_end ())
      { // get face top of vertex
        Halfedge_iterator he_left = get_left_halfedge (v_curr);
        if (he_left == halfedges_end ())
        {
          Halfedge_iterator he_right = get_right_halfedge (v_curr);
          if (he_right == halfedges_end ()){
            printf("[TsplinePatch::shoot_up] param: %e %e\n", v_curr->data().param.x(), v_curr->data().param.y());
            throw std::runtime_error ("[TsplinePatch::shoot_up] Error, invalid vertex");
          }
          f = he_right->twin ()->face ();
        }
        else
          f = he_left->face ();
      }

      if (f->is_unbounded ())
      {
        //        if(GetNorth(x_curr, north)){
        north=GetNorth(x_curr);
        if(north)
        {
          for(unsigned char m=0; m<north->multiplicity && knots.size()<n; m++)
            knots.push_back(y_curr);
          if(knots.size()<n)
            north->shoot(Point2d(x_curr, y_curr), n-knots.size(), knots);
        }else
          knots.push_back (y_curr); // multiply border knots
      }else
      {
        Halfedge_iterator he_top = get_top_halfedge_at_param (f, x_curr, y_curr);
        knots.push_back (y_curr);
        f = he_top->twin ()->face ();
      }
    }
    else
    { // parse along edges
      v_curr = he_v->source ();
      TVertex vext_curr = v_curr->data ();
      y_curr = vext_curr.param.y ();
      knots.push_back (y_curr);
      he_v = get_top_halfedge (v_curr);
    }
  }
}

void TsplinePatch::shoot_down (Vertex_iterator &v, const unsigned &n, std::vector<double> &knots)
{
  knots.clear ();
  TVertex vext = v->data ();
  double x_curr = vext.param.x ();
  double y_curr = vext.param.y ();
  Vertex_iterator v_curr = v;
  Halfedge_iterator he_v = get_bottom_halfedge (v_curr);
  Face_iterator f = faces_end ();
  Transition* south;

  while (knots.size () < n)
  {
    if (he_v == halfedges_end ())
    { // find edge on the bottom of the bottom face of the vertex
      if (f == faces_end ())
      { // get face bottom of vertex
        Halfedge_iterator he_left = get_left_halfedge (v_curr);
        if (he_left == halfedges_end ())
        {
          Halfedge_iterator he_right = get_right_halfedge (v_curr);
          if (he_right == halfedges_end ())
            throw std::runtime_error ("[TsplinePatch::shoot_down] Error, invalid vertex");
          f = he_right->face ();
        }
        else
          f = he_left->twin ()->face ();
      }

      if (f->is_unbounded ())
      {
        //        if(m_south.Valid())
        south = GetSouth(x_curr);
        if(south)
        {
          for(unsigned char m=0; m<south->multiplicity && knots.size()<n; m++)
            knots.push_back(y_curr);
          if(knots.size()<n)
            south->shoot(Point2d(x_curr, y_curr), n-knots.size(), knots);
        }else
          knots.push_back (y_curr); // multiply border knots
      }else
      {
        Halfedge_iterator he_bottom = get_bottom_halfedge_at_param (f, x_curr, y_curr);
        knots.push_back (y_curr);
        f = he_bottom->twin ()->face ();
      }
    }
    else
    { // parse along edges
      v_curr = he_v->source ();
      TVertex vext_curr = v_curr->data ();
      y_curr = vext_curr.param.y ();
      knots.push_back (y_curr);
      he_v = get_bottom_halfedge (v_curr);
    }
  }
}
