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

TsplineMultiPatch::TsplineMultiPatch()
  : quiet(true)
{

}

//TsplineMultiPatch::TsplineMultiPatch(const TsplineMultiPatch& tsmp)
//{
//  std::list<TsplinePatch*>::const_iterator pit0;
//  for(pit0=tsmp.patchlist.begin(); pit0!=tsmp.patchlist.end(); pit0++)
//  {
//    TsplinePatch* np = new TsplinePatch((*pit0)->id);
//    *np = *(*pit0);
//    patchlist.push_back(np);

//    tspline::Tspline::Vertex_iterator vit;
//    for(vit=np->vertices_begin(); vit!=np->vertices_end(); vit++)
//    {
//      TVertex& vext = vit->data();
//      vext.UnLink();
//    }
//  }

//  CopyTransitions(tsmp);

//  quiet = tsmp.quiet;
//}

TsplineMultiPatch::~TsplineMultiPatch()
{
  std::list<TsplinePatch*>::iterator pit;
  for(pit=patchlist.begin(); pit!=patchlist.end(); pit++)
    delete(*pit);
}

int TsplineMultiPatch::AddPatch(unsigned segX, unsigned segY, double s0, double s1, double t0, double t1)
{
  patchlist.push_back(new TsplinePatch(patchlist.size()));
  TsplinePatch* tsp = patchlist.back();

  TsplineCreator::CreatePlaneXY(*tsp, s0, t0, 0.0, s1-s0, t1-t0, segX, segY);

  Tspline::Vertex_iterator vit;
  for(vit=tsp->vertices_begin(); vit!=tsp->vertices_end(); vit++)
  {
    TVertex& vext = vit->data();
    vext.patch_id = tsp->GetID();
  }

  return tsp->GetID();
}

int TsplineMultiPatch::AddPatch(const Tspline& _tsp)
{
  patchlist.push_back(new TsplinePatch(patchlist.size(), _tsp));
  TsplinePatch* tsp = patchlist.back();

  Tspline::Vertex_iterator vit;
  for(vit=tsp->vertices_begin(); vit!=tsp->vertices_end(); vit++)
  {
    TVertex& vext = vit->data();
    vext.patch_id = tsp->GetID();
  }

  return tsp->GetID();
}

void TsplineMultiPatch::AddTransition(const size_t& id0, const Direction& ed0,
                                      const size_t& id1, const Direction& ed1,
                                      unsigned char mult)
{
  AddTransition(id0,ed0,0.0,1.0, id1,ed1,0.0,1.0,mult);

  //  TsplinePatch* patch0 = NULL;
  //  TsplinePatch* patch1 = NULL;
  //  std::list<TsplinePatch*>::iterator it;
  //  for(it=patchlist.begin(); it!=patchlist.end(); it++)
  //  {
  //    if((*it)->GetID()==id0)
  //      patch0 = (*it);
  //    if((*it)->GetID()==id1)
  //      patch1 = (*it);
  //  }
  //  if(patch0==NULL || patch1==NULL)
  //    throw std::runtime_error("[TsplineMultiPatch::AddTransition] Error, patches not found");

  //  Transition* trans0 = patch0->AddTransition(ed0, patch1, ed1);
  //  Transition* trans1 = patch1->AddTransition(ed1, patch0, ed0);

  //  // link control points
  //  Tspline::Vertex_iterator vit0, vit1;
  //  for(vit0=patch0->vertices_begin(); vit0!=patch0->vertices_end(); vit0++)
  //  {
  //    TVertex& vext0 = vit0->data();
  //    if(ed0==NORTH && equal(vext0.param.y(),1.0))
  //    {
  //      Point2d param1 = trans0->ConvertParam(vext0.param);
  //      vit1 = patch1->insert_vertex_no_transition(param1.x(), param1.y());
  //      TVertex &vext1 = vit1->data();
  //      vext0.LinkCP(&vext1);
  //    }
  //    if(ed0==EAST && equal(vext0.param.x(),1.0))
  //    {
  //      Point2d param1 = trans0->ConvertParam(vext0.param);
  //      vit1 = patch1->insert_vertex_no_transition(param1.x(), param1.y());
  //      TVertex &vext1 = vit1->data();
  //      vext0.LinkCP(&vext1);
  //    }
  //    if(ed0==SOUTH && equal(vext0.param.y(),0.0))
  //    {
  //      Point2d param1 = trans0->ConvertParam(vext0.param);
  //      vit1 = patch1->insert_vertex_no_transition(param1.x(), param1.y());
  //      TVertex &vext1 = vit1->data();
  //      vext0.LinkCP(&vext1);
  //    }
  //    if(ed0==WEST && equal(vext0.param.x(),0.0))
  //    {
  //      Point2d param1 = trans0->ConvertParam(vext0.param);
  //      vit1 = patch1->insert_vertex_no_transition(param1.x(), param1.y());
  //      TVertex &vext1 = vit1->data();
  //      vext0.LinkCP(&vext1);
  //    }
  //  }

  //  for(vit1=patch1->vertices_begin(); vit1!=patch1->vertices_end(); vit1++)
  //  {
  //    TVertex& vext1 = vit1->data();
  //    if(ed1==NORTH && equal(vext1.param.y(),1.0))
  //    {
  //      Point2d param0 = trans1->ConvertParam(vext1.param);
  //      vit0 = patch0->insert_vertex_no_transition(param0.x(), param0.y());
  //      TVertex &vext0 = vit0->data();
  //      vext1.LinkCP(&vext0);
  //    }
  //    if(ed1==EAST && equal(vext1.param.x(),1.0))
  //    {
  //      Point2d param0 = trans1->ConvertParam(vext1.param);
  //      vit0 = patch0->insert_vertex_no_transition(param0.x(), param0.y());
  //      TVertex &vext0 = vit0->data();
  //      vext1.LinkCP(&vext0);
  //    }
  //    if(ed1==SOUTH && equal(vext1.param.y(),0.0))
  //    {
  //      Point2d param0 = trans1->ConvertParam(vext1.param);
  //      vit0 = patch0->insert_vertex_no_transition(param0.x(), param0.y());
  //      TVertex &vext0 = vit0->data();
  //      vext1.LinkCP(&vext0);
  //    }
  //    if(ed1==WEST && equal(vext1.param.x(),0.0))
  //    {
  //      Point2d param0 = trans1->ConvertParam(vext1.param);
  //      vit0 = patch0->insert_vertex_no_transition(param0.x(), param0.y());
  //      TVertex &vext0 = vit0->data();
  //      vext1.LinkCP(&vext0);
  //    }
  //  }

  //  for(it=patchlist.begin(); it!=patchlist.end(); it++)
  //    (*it)->update_knot_vectors();

  //  UpdateVertices();
}

void TsplineMultiPatch::AddTransition(const int& id0, const Direction& ed0, double r0a, double r0b,
                                      const int& id1, const Direction& ed1, double r1a, double r1b,
                                      unsigned char mult)
{
  TsplinePatch* patch0 = NULL;
  TsplinePatch* patch1 = NULL;
  std::list<TsplinePatch*>::iterator it;
  for(it=patchlist.begin(); it!=patchlist.end(); it++)
  {
    if((*it)->GetID()==id0)
      patch0 = (*it);
    if((*it)->GetID()==id1)
      patch1 = (*it);
  }
  if(patch0==NULL || patch1==NULL)
    throw std::runtime_error("[TsplineMultiPatch::AddTransition] Error, patches not found");

  Transition* trans0 = patch0->AddTransition(ed0, r0a, r0b, patch1, ed1, r1a, r1b, mult);
  Transition* trans1 = patch1->AddTransition(ed1, r1a, r1b, patch0, ed0, r0a, r0b, mult);
  trans0->twin = trans1;
  trans1->twin = trans0;

  bool quiet(true);

  // link control points
  Tspline::Vertex_iterator vit0, vit1;
  for(vit0=patch0->vertices_begin(); vit0!=patch0->vertices_end(); vit0++)
  {
    TVertex& vext0 = vit0->data();
    if(ed0==NORTH && equal(vext0.param.y(),patch0->param_max.y()) && trans0->InSourceRange(vext0.param.x()))
    {
      if(!quiet)
        printf("[TsplineMultiPatch::AddTransition] ed0==NORTH\n");
      Point2d param1 = trans0->ConvertParam(vext0.param);
      if(!quiet)
        printf("  param0: %f %f -> param1: %f %f\n", vext0.param.x(), vext0.param.y(), param1.x(), param1.y());
      vit1 = patch1->insert_vertex_no_transition(param1.x(), param1.y());
      TVertex &vext1 = vit1->data();
      vext1.LinkCP(&vext0);
    }
    if(ed0==EAST && equal(vext0.param.x(),patch0->param_max.x()) && trans0->InSourceRange(vext0.param.y()))
    {
      if(!quiet)
        printf("[TsplineMultiPatch::AddTransition] ed0==EAST\n");
      Point2d param1 = trans0->ConvertParam(vext0.param);
      if(!quiet)
        printf("  param0: %f %f -> param1: %f %f\n", vext0.param.x(), vext0.param.y(), param1.x(), param1.y());
      vit1 = patch1->insert_vertex_no_transition(param1.x(), param1.y());
      TVertex &vext1 = vit1->data();
      vext1.LinkCP(&vext0);
    }
    if(ed0==SOUTH && equal(vext0.param.y(),patch0->param_min.y()) && trans0->InSourceRange(vext0.param.x()))
    {
      if(!quiet)
        printf("[TsplineMultiPatch::AddTransition] ed0==SOUTH  dir1: %d\n", trans0->GetTargetDir());
      Point2d param1 = trans0->ConvertParam(vext0.param);
      if(!quiet)
        printf("  param0: %f %f -> param1: %f %f\n", vext0.param.x(), vext0.param.y(), param1.x(), param1.y());
      vit1 = patch1->insert_vertex_no_transition(param1.x(), param1.y());
      TVertex &vext1 = vit1->data();
      vext1.LinkCP(&vext0);
    }
    if(ed0==WEST && equal(vext0.param.x(),patch0->param_min.y()) && trans0->InSourceRange(vext0.param.y()))
    {
      if(!quiet)
        printf("[TsplineMultiPatch::AddTransition] ed0==WEST\n");
      Point2d param1 = trans0->ConvertParam(vext0.param);
      if(!quiet)
        printf("  param0: %f %f -> param1: %f %f\n", vext0.param.x(), vext0.param.y(), param1.x(), param1.y());
      vit1 = patch1->insert_vertex_no_transition(param1.x(), param1.y());
      TVertex &vext1 = vit1->data();
      vext1.LinkCP(&vext0);
    }
  }

  for(vit1=patch1->vertices_begin(); vit1!=patch1->vertices_end(); vit1++)
  {
    TVertex& vext1 = vit1->data();
    if(ed1==NORTH && equal(vext1.param.y(),patch1->param_max.y()) && trans1->InSourceRange(vext1.param.x()))
    {
      if(!quiet)
        printf("[TsplineMultiPatch::AddTransition] ed1==NORTH  dir0: %d\n", trans1->GetTargetDir());
      Point2d param0 = trans1->ConvertParam(vext1.param);
      if(!quiet)
        printf("  param1: %f %f -> param0: %f %f\n", vext1.param.x(), vext1.param.y(), param0.x(), param0.y());
      vit0 = patch0->insert_vertex_no_transition(param0.x(), param0.y());
      TVertex &vext0 = vit0->data();
      vext0.LinkCP(&vext1);
    }
    if(ed1==EAST && equal(vext1.param.x(),patch1->param_max.x()) && trans1->InSourceRange(vext1.param.y()))
    {
      if(!quiet)
        printf("[TsplineMultiPatch::AddTransition] ed1==EAST\n");
      Point2d param0 = trans1->ConvertParam(vext1.param);
      if(!quiet)
        printf("  param1: %f %f -> param0: %f %f\n", vext1.param.x(), vext1.param.y(), param0.x(), param0.y());
      vit0 = patch0->insert_vertex_no_transition(param0.x(), param0.y());
      TVertex &vext0 = vit0->data();
      vext0.LinkCP(&vext1);
    }
    if(ed1==SOUTH && equal(vext1.param.y(),patch1->param_min.y()) && trans1->InSourceRange(vext1.param.x()))
    {
      if(!quiet)
        printf("[TsplineMultiPatch::AddTransition] ed1==SOUTH\n");
      Point2d param0 = trans1->ConvertParam(vext1.param);
      if(!quiet)
        printf("  param1: %f %f -> param0: %f %f\n", vext1.param.x(), vext1.param.y(), param0.x(), param0.y());
      vit0 = patch0->insert_vertex_no_transition(param0.x(), param0.y());
      TVertex &vext0 = vit0->data();
      vext0.LinkCP(&vext1);
    }
    if(ed1==WEST && equal(vext1.param.x(),patch1->param_min.x()) && trans1->InSourceRange(vext1.param.y()))
    {
      if(!quiet)
        printf("[TsplineMultiPatch::AddTransition] ed1==WEST\n");
      Point2d param0 = trans1->ConvertParam(vext1.param);
      if(!quiet)
        printf("  param1: %f %f -> param0: %f %f\n", vext1.param.x(), vext1.param.y(), param0.x(), param0.y());
      vit0 = patch0->insert_vertex_no_transition(param0.x(), param0.y());
      TVertex &vext0 = vit0->data();
      vext0.LinkCP(&vext1);
    }
  }

  // todo: this is very exhaustive; only update neighboring patches of p0, and p1 without p0, p1
  for(it=patchlist.begin(); it!=patchlist.end(); it++)
    (*it)->update_knot_vectors();
}

void TsplineMultiPatch::InsertVertex(const int &id, const Point2d& param)
{
  TsplinePatch* patch = NULL;
  std::list<TsplinePatch*>::iterator it;
  for(it=patchlist.begin(); it!=patchlist.end(); it++)
  {
    if((*it)->GetID()==id)
      patch = (*it);
  }
  if(patch==NULL)
    throw std::runtime_error("[TsplineMultiPatch::InsertVertex] Error, patch not found");

  patch->insert_vertex(param.x(), param.y());
  patch->update_knot_vectors();

  for(size_t i=0; i<patch->m_north.size(); i++)
    patch->m_north[i]->target->update_knot_vectors();
  for(size_t i=0; i<patch->m_east.size(); i++)
    patch->m_east[i]->target->update_knot_vectors();
  for(size_t i=0; i<patch->m_south.size(); i++)
    patch->m_south[i]->target->update_knot_vectors();
  for(size_t i=0; i<patch->m_west.size(); i++)
    patch->m_west[i]->target->update_knot_vectors();
}

vector_vec3d TsplineMultiPatch::compute_cp_normals_by_footpoints() const
{
  vector_vec3d normals;
  std::vector<tspline::Tspline::Vertex_const_iterator> controlpoints = get_controlpoints();

  for(size_t i=0; i<controlpoints.size(); i++)
  {
    const TVertex& vext = controlpoints[i]->data();
    Eigen::Vector3d ts, tt, n;

    std::list<TsplinePatch*>::const_iterator pit;
    for(pit=patchlist.begin(); pit!=patchlist.end(); pit++)
      if((*pit)->id==vext.patch_id)
        break;

    double s=vext.param.x();
    double t=vext.param.y();
    //    if(equal(s,(*pit)->param_min.x()))
    //      s+=0.01;//(10.0*epsilon);
    //    if(equal(s,(*pit)->param_max.x()))
    //      s-=0.01;//(10.0*epsilon);
    //    if(equal(t,(*pit)->param_min.y()))
    //      t+=0.01;//(10.0*epsilon);
    //    if(equal(t,(*pit)->param_max.y()))
    //      t-=0.01;//(10.0*epsilon);

    (*pit)->evaluate(s, t, ts, tt);
    n = ts.cross(tt);
    n.normalize();
    normals.push_back(n);

    //    if(viewer!=NULL)
    //    {
    //      Eigen::Vector3d a;
    //      vext.GetCP(a);
    //      float s(0.01);
    //      viewer->AddPoint3D(fp(0),fp(1),fp(2),255, 0, 255, 5.0);
    //      viewer->AddLine3D(a(0),a(1),a(2),a(0)+n(0)*s,a(1)+n(1)*s,a(2)+n(2)*s, 255, 0, 255, 4.0f);
    //    }
  }
  return normals;
}

void TsplineMultiPatch::CopyTransitions(const TsplineMultiPatch &tsmp)
{
  std::list<TsplinePatch*>::const_iterator pit0, pit1;

  // relink vertices
  // simultaneous loop over new and old patchlist
  pit1 = patchlist.begin();
  for(pit0=tsmp.patchlist.begin(); pit0!=tsmp.patchlist.end(); pit0++)
  {
    TsplinePatch* p0 = (*pit0);
    TsplinePatch* p1 = (*pit1);
    for(size_t i=0; i<p0->m_north.size(); i++) //p0->m_north.Valid())
      AddTransition(p1->GetID(), NORTH, p0->m_north[i]->target->GetID(), p0->m_north[i]->target_dir, p0->m_north[i]->multiplicity);
    for(size_t i=0; i<p0->m_east.size(); i++) //if(p0->m_east.Valid())
      AddTransition(p1->GetID(), EAST, p0->m_east[i]->target->GetID(), p0->m_east[i]->target_dir, p0->m_east[i]->multiplicity);
    for(size_t i=0; i<p0->m_south.size(); i++) //if(p0->m_south.Valid())
      AddTransition(p1->GetID(), SOUTH, p0->m_south[i]->target->GetID(), p0->m_south[i]->target_dir, p0->m_south[i]->multiplicity);
    for(size_t i=0; i<p0->m_west.size(); i++) //if(p0->m_west.Valid())
      AddTransition(p1->GetID(), WEST, p0->m_west[i]->target->GetID(), p0->m_west[i]->target_dir, p0->m_west[i]->multiplicity);

    pit1++;
  }
}

std::vector<tspline::Tspline::Vertex_iterator> TsplineMultiPatch::get_controlpoints()
{
  std::vector<tspline::Tspline::Vertex_iterator> controlpoints;
  std::list<TsplinePatch*>::iterator pit;
  for(pit=patchlist.begin(); pit!=patchlist.end(); pit++)
    for(Tspline::Vertex_iterator vit=(*pit)->vertices_begin(); vit!=(*pit)->vertices_end(); vit++)
      if(vit->data().is_primary)
        controlpoints.push_back(vit);

  return controlpoints;
}

std::vector<tspline::Tspline::Vertex_const_iterator> TsplineMultiPatch::get_controlpoints() const
{
  std::vector<tspline::Tspline::Vertex_const_iterator> controlpoints;
  std::list<TsplinePatch*>::const_iterator pit;
  for(pit=patchlist.begin(); pit!=patchlist.end(); pit++)
    for(Tspline::Vertex_const_iterator vit=(*pit)->vertices_begin(); vit!=(*pit)->vertices_end(); vit++)
      if(vit->data().is_primary)
        controlpoints.push_back(vit);

  return controlpoints;
}

std::vector<tspline::Tspline::Vertex_iterator> TsplineMultiPatch::get_vertices()
{
  std::vector<tspline::Tspline::Vertex_iterator> vertices;
  std::list<TsplinePatch*>::iterator pit;
  for(pit=patchlist.begin(); pit!=patchlist.end(); pit++)
    for(Tspline::Vertex_iterator vit=(*pit)->vertices_begin(); vit!=(*pit)->vertices_end(); vit++)
      vertices.push_back(vit);

  return vertices;
}

std::vector<tspline::Tspline::Vertex_const_iterator> TsplineMultiPatch::get_vertices() const
{
  std::vector<tspline::Tspline::Vertex_const_iterator> vertices;
  std::list<TsplinePatch*>::const_iterator pit;
  for(pit=patchlist.begin(); pit!=patchlist.end(); pit++)
    for(Tspline::Vertex_const_iterator vit=(*pit)->vertices_begin(); vit!=(*pit)->vertices_end(); vit++)
      vertices.push_back(vit);

  return vertices;
}

size_t TsplineMultiPatch::number_of_controlpoints() const
{
  return get_controlpoints().size();
}
size_t TsplineMultiPatch::number_of_vertices() const
{
  return get_vertices().size();
}

Eigen::Vector3d TsplineMultiPatch::center() const
{
  Eigen::Vector3d center(0.0,0.0,0.0), cp;
  std::vector<tspline::Tspline::Vertex_const_iterator> cps = get_controlpoints();
  for(size_t i=0; i<cps.size(); i++)
  {
    cps[i]->data().GetCP(cp);
    center += cp;
  }
  center /= cps.size();
  return center;
}
