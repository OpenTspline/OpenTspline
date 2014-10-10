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

Arrangement_2::Halfedge_iterator Tspline::get_left_halfedge_at_param (Face_iterator &f, const double &y, double &x)
{
  if (f->is_unbounded())
    return halfedges_end();

  Ccb_halfedge_circulator first, curr;
  first = curr = f->outer_ccb();

  do
  {
    Vertex_iterator p0 = curr->source();
    Vertex_iterator p1 = curr->target();

    if (greater(p0->point().y(), p1->point().y())) // is left
    {
      TVertex v0 = p0->data();
      TVertex v1 = p1->data();
      if(gequal(v0.param.y(), y) && sequal(v1.param.y(), y))
      {
        x = v0.param.x();
        return curr;
      }
    }
  }
  while (++curr != first);

  return halfedges_end();
}

Arrangement_2::Halfedge_iterator Tspline::get_right_halfedge_at_param (Face_iterator &f, const double &y, double &x)
{
  if (f->is_unbounded())
    return halfedges_end();

  Ccb_halfedge_circulator first, curr;
  first = curr = f->outer_ccb();

  do
  {
    Vertex_iterator p0 = curr->source();
    Vertex_iterator p1 = curr->target();

    if (smaller(p0->point().y(), p1->point().y()))
    {
      TVertex v0 = p0->data();
      TVertex v1 = p1->data();
      if (sequal(v0.param.y(), y) && gequal(v1.param.y(), y))
      {
        x = v0.param.x();
        return curr;
      }
    }
  }
  while (++curr != first);

  return halfedges_end();
}

Arrangement_2::Halfedge_iterator Tspline::get_top_halfedge_at_param (Face_iterator &f, const double &x, double &y)
{
  if (f->is_unbounded())
    return halfedges_end();

  Ccb_halfedge_circulator first, curr;
  first = curr = f->outer_ccb();

  do
  {
    Vertex_iterator p0 = curr->source();
    Vertex_iterator p1 = curr->target();

    if (greater(p0->point().x(), p1->point().x()))
    {
      TVertex v0 = p0->data();
      TVertex v1 = p1->data();
      if (gequal(v0.param.x(), x) && sequal(v1.param.x(), x))
      {
        y = v0.param.y();
        return curr;
      }
    }
  }
  while (++curr != first);

  return halfedges_end();
}

Arrangement_2::Halfedge_iterator Tspline::get_bottom_halfedge_at_param (Face_iterator &f, const double &x, double &y)
{
  if (f->is_unbounded())
    return halfedges_end();

  Ccb_halfedge_circulator first, curr;
  first = curr = f->outer_ccb();

  do
  {
    Vertex_iterator p0 = curr->source();
    Vertex_iterator p1 = curr->target();

    if(smaller(p0->point().x(), p1->point().x()))
    {
      TVertex v0 = p0->data();
      TVertex v1 = p1->data();

      if(sequal(v0.param.x(), x) && gequal(v1.param.x(), x))
      {
        y = v0.param.y();
        return curr;
      }
    }

  }
  while (++curr != first);

  return halfedges_end();
}

Arrangement_2::Halfedge_iterator Tspline::get_left_halfedge_at_point (Face_iterator f, const Point2d &p)
{
  if (f->is_unbounded())
    return halfedges_end();

  Ccb_halfedge_circulator first, curr;
  Halfedge_iterator edge;
  bool edge_found(false);
  double edge_dist(DBL_MAX);

  for(Hole_iterator hole=f->holes_begin(); hole!=f->holes_end(); hole++)
  {
    first = curr = (*hole);
    do
    {
      const Point2d& p0 = curr->source()->point();
      const Point2d& p1 = curr->target()->point();
      if (greater(p0.y(), p1.y()) && // edge is vertical and pointing downwards
          smaller(p0.x(), p.x())  && // edge is left of p
          gequal(p0.y(), p.y()) && sequal(p1.y(), p.y())) // edge covers range of p
      {
        double d = p.x() - p0.x();
        if(d<edge_dist)
        {
          edge = curr;
          edge_found = true;
          edge_dist = d;
        }
      }
    }while (++curr != first);
  }

  first = curr = f->outer_ccb();
  do
  {
    const Point2d& p0 = curr->source()->point();
    const Point2d& p1 = curr->target()->point();
    if (greater(p0.y(), p1.y()) && // edge is vertical and pointing downwards
        smaller(p0.x(), p.x())  && // edge is left of p
        gequal(p0.y(), p.y()) && sequal(p1.y(), p.y())) // edge covers range of p
    {
      double d = p.x() - p0.x();
      if(d<edge_dist)
      {
        edge = curr;
        edge_found = true;
        edge_dist = d;
      }
    }
  }while (++curr != first);

  if(edge_found)
    return edge;

  return halfedges_end();
}


Arrangement_2::Halfedge_iterator Tspline::get_right_halfedge_at_point (Face_iterator f, const Point2d &p)
{
  if (f->is_unbounded())
    return halfedges_end();

  Ccb_halfedge_circulator first, curr;
  Halfedge_iterator edge;
  bool edge_found(false);
  double edge_dist(DBL_MAX);

  for(Hole_iterator hole=f->holes_begin(); hole!=f->holes_end(); hole++)
  {
    first = curr = (*hole);
    do
    {
      const Point2d& p0 = curr->source()->point();
      const Point2d& p1 = curr->target()->point();
      if (smaller(p0.y(), p1.y()) && // edge is vertical and pointing upwards
          greater(p0.x(), p.x())  && // edge is right of point
          sequal(p0.y(), p.y()) && gequal(p1.y(), p.y())) // edge covers range of p
      {
        double d = p0.x() - p.x();
        if(d<edge_dist)
        {
          edge = curr;
          edge_found = true;
          edge_dist = d;
        }
      }
    }while (++curr != first);
  }

  first = curr = f->outer_ccb();
  do
  {
    const Point2d& p0 = curr->source()->point();
    const Point2d& p1 = curr->target()->point();
    if (smaller(p0.y(), p1.y()) && // edge is vertical and pointing upwards
        greater(p0.x(), p.x())  && // edge is right of point
        sequal(p0.y(), p.y()) && gequal(p1.y(), p.y())) // edge covers range of p
    {
      double d = p0.x() - p.x();
      if(d<edge_dist)
      {
        edge = curr;
        edge_found = true;
        edge_dist = d;
      }
    }
  }while (++curr != first);

  if(edge_found)
    return edge;

  return halfedges_end();
}

Arrangement_2::Halfedge_iterator Tspline::get_top_halfedge_at_point (Face_iterator f, const Point2d &p)
{
  if (f->is_unbounded())
    return halfedges_end();

  Ccb_halfedge_circulator first, curr;
  Halfedge_iterator edge;
  bool edge_found(false);
  double edge_dist(DBL_MAX);

  for(Hole_iterator hole=f->holes_begin(); hole!=f->holes_end(); hole++)
  {
    first = curr = (*hole);
    do
    {
      const Point2d& p0 = curr->source()->point();
      const Point2d& p1 = curr->target()->point();
      if (greater(p0.x(), p1.x()) && // edge is horizontal and a top edge
          greater(p0.y(),p.y()) &&  // edge is above p
          gequal(p0.x(), p.x()) && sequal(p1.x(), p.x())) // edge covers range of p
      {
        double d = p0.y() - p.y();
        if(d<edge_dist)
        {
          edge = curr;
          edge_found = true;
          edge_dist = d;
        }
      }
    }while (++curr != first);
  }

  first = curr = f->outer_ccb();
  do
  {
    const Point2d& p0 = curr->source()->point();
    const Point2d& p1 = curr->target()->point();
    if (greater(p0.x(), p1.x()) && // edge is horizontal and a top edge
        greater(p0.y(),p.y()) &&  // edge is above p
        gequal(p0.x(), p.x()) && sequal(p1.x(), p.x())) // edge covers range of p
    {
      double d = p0.y() - p.y();
      if(d<edge_dist)
      {
        edge = curr;
        edge_found = true;
        edge_dist = d;
      }
    }
  }while (++curr != first);

  if(edge_found)
    return edge;

  return halfedges_end();
}

Arrangement_2::Halfedge_iterator Tspline::get_bottom_halfedge_at_point (Face_iterator f, const Point2d& p)
{
  if (f->is_unbounded())
    return halfedges_end();

  Ccb_halfedge_circulator first, curr;
  Halfedge_iterator edge;
  bool edge_found(false);
  double edge_dist(DBL_MAX);

  for(Hole_iterator hole=f->holes_begin(); hole!=f->holes_end(); hole++)
  {
    first = curr = (*hole);
    do
    {
      const Point2d& p0 = curr->source()->point();
      const Point2d& p1 = curr->target()->point();
      if( smaller(p0.x(), p1.x()) && // edge is horizontal and a bottom edge
          smaller(p0.y(),p.y()) &&  // edge is below point
          sequal(p0.x(), p.x()) && gequal(p1.x(), p.x()) ) // edge covers range of point
      {
        double d = p.y() - p0.y();
        if(d<edge_dist)
        {
          edge = curr;
          edge_found = true;
          edge_dist = d;

        }
      }
    }while (++curr != first);
  }

  first = curr = f->outer_ccb();
  do
  {
    const Point2d& p0 = curr->source()->point();
    const Point2d& p1 = curr->target()->point();
    if( smaller(p0.x(), p1.x()) && // edge is horizontal and a bottom edge
        smaller(p0.y(),p.y()) &&  // edge is below point
        sequal(p0.x(), p.x()) && gequal(p1.x(), p.x()) ) // edge covers range of point
    {
      double d = p.y() - p0.y();
      if(d<edge_dist)
      {
        edge = curr;
        edge_found = true;
        edge_dist = d;

      }
    }
  }while (++curr != first);

  if(edge_found)
    return edge;

  return halfedges_end();
}

void Tspline::shoot_left (Vertex_iterator &v, const unsigned &n, std::vector<double> &knots)
{
  knots.clear();
  TVertex vext = v->data();
  double x_curr = vext.param.x();
  double y_curr = vext.param.y();
  Vertex_iterator v_curr = v;
  Halfedge_iterator he_v = get_left_halfedge (v_curr);
  Face_iterator f = faces_end();

  while (knots.size() < n)
  {
    if (he_v == halfedges_end())
    { // find edge on the left of the left face of the vertex
      if (f == faces_end())
      { // get face left of vertex
        Halfedge_iterator he_top = get_top_halfedge (v_curr);
        if (he_top == halfedges_end())
        {
          Halfedge_iterator he_bot = get_bottom_halfedge (v_curr);
          if (he_bot == halfedges_end())
            throw std::runtime_error ("[Tspline::shoot_left] Error, invalid vertex");
          f = he_bot->face();
        }
        else
          f = he_top->twin()->face();
      }

      if (f->is_unbounded())
        knots.push_back (x_curr); // multiply border knots
      else
      {
        Halfedge_iterator he_left = get_left_halfedge_at_param (f, y_curr, x_curr);
        knots.push_back (x_curr);
        f = he_left->twin()->face();
      }
    }
    else
    { // parse along edges
      v_curr = he_v->source();
      TVertex vext_curr = v_curr->data();
      x_curr = vext_curr.param.x();
      knots.push_back (x_curr);
      he_v = get_left_halfedge (v_curr);
    }
  }
}

void Tspline::shoot_right (Vertex_iterator &v, const unsigned &n, std::vector<double> &knots)
{
  knots.clear();
  TVertex vext = v->data();
  double x_curr = vext.param.x();
  double y_curr = vext.param.y();
  Vertex_iterator v_curr = v;
  Halfedge_iterator he_v = get_right_halfedge (v_curr);
  Face_iterator f = faces_end();

  while (knots.size() < n)
  {
    if (he_v == halfedges_end())
    { // find edge on the right of the right face of the vertex
      if (f == faces_end())
      { // get face right of vertex
        Halfedge_iterator he_top = get_top_halfedge (v_curr);
        if (he_top == halfedges_end())
        {
          Halfedge_iterator he_bot = get_bottom_halfedge (v_curr);
          if (he_bot == halfedges_end())
            throw std::runtime_error ("[Tspline::shoot_right] Error, invalid vertex");
          f = he_bot->twin()->face();
        }
        else
          f = he_top->face();
      }

      if (f->is_unbounded())
        knots.push_back (x_curr); // multiply border knots
      else
      {
        Halfedge_iterator he_right = get_right_halfedge_at_param (f, y_curr, x_curr);
        knots.push_back (x_curr);
        f = he_right->twin()->face();
      }
    }
    else
    { // parse along edges
      v_curr = he_v->source();
      TVertex vext_curr = v_curr->data();
      x_curr = vext_curr.param.x();
      knots.push_back (x_curr);
      he_v = get_right_halfedge (v_curr);
    }
  }
}

void Tspline::shoot_up (Vertex_iterator &v, const unsigned &n, std::vector<double> &knots)
{
  knots.clear();
  TVertex vext = v->data();
  double x_curr = vext.param.x();
  double y_curr = vext.param.y();
  Vertex_iterator v_curr = v;
  Halfedge_iterator he_v = get_top_halfedge (v_curr);
  Face_iterator f = faces_end();

  while (knots.size() < n)
  {
    if (he_v == halfedges_end())
    { // find edge on the top of the top face of the vertex
      if (f == faces_end())
      { // get face top of vertex
        Halfedge_iterator he_left = get_left_halfedge (v_curr);
        if (he_left == halfedges_end())
        {
          Halfedge_iterator he_right = get_right_halfedge (v_curr);
          if (he_right == halfedges_end())
            throw std::runtime_error ("[Tspline::shoot_up] Error, invalid vertex");
          f = he_right->twin()->face();
        }
        else
          f = he_left->face();
      }

      if (f->is_unbounded())
        knots.push_back (y_curr); // multiply border knots
      else
      {
        Halfedge_iterator he_top = get_top_halfedge_at_param (f, x_curr, y_curr);
        knots.push_back (y_curr);
        f = he_top->twin()->face();
      }
    }
    else
    { // parse along edges
      v_curr = he_v->source();
      TVertex vext_curr = v_curr->data();
      y_curr = vext_curr.param.y();
      knots.push_back (y_curr);
      he_v = get_top_halfedge (v_curr);
    }
  }
}

void Tspline::shoot_down (Vertex_iterator &v, const unsigned &n, std::vector<double> &knots)
{
  knots.clear();
  TVertex vext = v->data();
  double x_curr = vext.param.x();
  double y_curr = vext.param.y();
  Vertex_iterator v_curr = v;
  Halfedge_iterator he_v = get_bottom_halfedge (v_curr);
  Face_iterator f = faces_end();

  while (knots.size() < n)
  {
    if (he_v == halfedges_end())
    { // find edge on the bottom of the bottom face of the vertex
      if (f == faces_end())
      { // get face bottom of vertex
        Halfedge_iterator he_left = get_left_halfedge (v_curr);
        if (he_left == halfedges_end())
        {
          Halfedge_iterator he_right = get_right_halfedge (v_curr);
          if (he_right == halfedges_end())
            throw std::runtime_error ("[Tspline::shoot_down] Error, invalid vertex");
          f = he_right->face();
        }
        else
          f = he_left->twin()->face();
      }

      if (f->is_unbounded())
        knots.push_back (y_curr); // multiply border knots
      else
      {
        Halfedge_iterator he_bottom = get_bottom_halfedge_at_param (f, x_curr, y_curr);
        knots.push_back (y_curr);
        f = he_bottom->twin()->face();
      }
    }
    else
    { // parse along edges
      v_curr = he_v->source();
      TVertex vext_curr = v_curr->data();
      y_curr = vext_curr.param.y();
      knots.push_back (y_curr);
      he_v = get_bottom_halfedge(v_curr);
    }
  }
}

void Tspline::shoot_left(Halfedge_iterator& he, const double& y_curr, const unsigned &n, std::vector<double> &knots)
{
  double x_curr = he->source()->data().param.x();
  if(!equal(x_curr, he->target()->data().param.x()))
    throw std::runtime_error("[Tspline::shoot_left_he] Error, halfedge is not vertical");

  knots.clear();
  if(smaller(he->target()->data().param.y(), he->source()->data().param.y()))
    he = he->twin();

  Face_iterator f = he->face();
  if(f->is_unbounded())
    throw std::runtime_error("[Tspline::shoot_left_he] Error, invalid halfedge");

  while (knots.size() < n)
  {
    if (f->is_unbounded())
      knots.push_back (x_curr); // multiply border knots
    else
    {
      Halfedge_iterator he_left = get_left_halfedge_at_param (f, y_curr, x_curr);
      knots.push_back (x_curr);
      f = he_left->twin()->face();
    }
  }
}

void Tspline::shoot_right(Halfedge_iterator& he, const double& y_curr, const unsigned &n, std::vector<double> &knots)
{
  double x_curr = he->source()->data().param.x();
  if(!equal(x_curr, he->target()->data().param.x()))
    throw std::runtime_error("[Tspline::shoot_right_he] Error, halfedge is not vertical");

  knots.clear();
  if(greater(he->target()->data().param.y(), he->source()->data().param.y()))
    he = he->twin();

  Face_iterator f = he->face();
  if(f->is_unbounded())
    throw std::runtime_error("[Tspline::shoot_right_he] Error, invalid halfedge");

  while (knots.size() < n)
  {
    if (f->is_unbounded())
      knots.push_back (x_curr); // multiply border knots
    else
    {
      Halfedge_iterator he_right = get_right_halfedge_at_param (f, y_curr, x_curr);
      knots.push_back (x_curr);
      f = he_right->twin()->face();
    }
  }
}

void Tspline::shoot_up(Halfedge_iterator& he, const double& x_curr, const unsigned &n, std::vector<double> &knots)
{
  double y_curr = he->source()->data().param.y();
  if(!equal(y_curr, he->target()->data().param.y()))
    throw std::runtime_error("[Tspline::shoot_up_he] Error, halfedge is not horizontal");

  knots.clear();
  if(smaller(he->target()->data().param.x(), he->source()->data().param.x()))
    he = he->twin();

  Face_iterator f = he->face();
  if(f->is_unbounded())
    throw std::runtime_error("[Tspline::shoot_up_he] Error, invalid halfedge");

  while (knots.size() < n)
  {
    if (f->is_unbounded())
      knots.push_back (y_curr); // multiply border knots
    else
    {
      Halfedge_iterator he_top = get_top_halfedge_at_param (f, x_curr, y_curr);
      knots.push_back (y_curr);
      f = he_top->twin()->face();
    }
  }
}

void Tspline::shoot_down(Halfedge_iterator& he, const double& x_curr, const unsigned &n, std::vector<double> &knots)
{
  double y_curr = he->source()->data().param.y();
  if(!equal(y_curr, he->target()->data().param.y()))
    throw std::runtime_error("[Tspline::shoot_down_he] Error, halfedge is not horizontal");

  knots.clear();
  if(greater(he->target()->data().param.x(), he->source()->data().param.x()))
    he = he->twin();

  Face_iterator f = he->face();
  if(f->is_unbounded())
    throw std::runtime_error("[Tspline::shoot_down_he] Error, invalid halfedge");

  while (knots.size() < n)
  {
    if (f->is_unbounded())
      knots.push_back (y_curr); // multiply border knots
    else
    {
      Halfedge_iterator he_bottom = get_bottom_halfedge_at_param (f, x_curr, y_curr);
      knots.push_back (y_curr);
      f = he_bottom->twin()->face();
    }
  }
}

Arrangement_2::Halfedge_iterator Tspline::get_right_halfedge (Vertex_iterator &v)
{
  Halfedge_around_vertex_circulator first, curr;
  first = curr = v->incident_halfedges();

  do
  {
    Vertex_iterator u = curr->source();
    Vector2d vu = u->point() - v->point();
    //if (vu.x() > 0)
    if(greater(vu.x(), 0.0))
      return curr;
  }
  while (++curr != first);

  return halfedges_end();
}

Arrangement_2::Halfedge_iterator Tspline::get_left_halfedge (Vertex_iterator &v)
{
  Halfedge_around_vertex_circulator first, curr;
  first = curr = v->incident_halfedges();

  do
  {
    Vertex_iterator u = curr->source();
    Vector2d vu = u->point() - v->point();
    //if (vu.x() < 0)
    if (smaller(vu.x(), 0.0))
      return curr;
  }
  while (++curr != first);

  return halfedges_end();
}

Arrangement_2::Halfedge_iterator Tspline::get_top_halfedge (Vertex_iterator &v)
{
  Halfedge_around_vertex_circulator first, curr;
  first = curr = v->incident_halfedges();

  do
  {
    Vertex_iterator u = curr->source();
    Vector2d vu = u->point() - v->point();
    //if (vu.y() > 0)
    if (greater(vu.y(), 0.0))
      return curr;
  }
  while (++curr != first);

  return halfedges_end();
}

Arrangement_2::Halfedge_iterator Tspline::get_bottom_halfedge (Vertex_iterator &v)
{
  Halfedge_around_vertex_circulator first, curr;
  first = curr = v->incident_halfedges();

  do
  {
    Vertex_iterator u = curr->source();
    Vector2d vu = u->point() - v->point();
    if (smaller(vu.y(), 0.0))
      return curr;
  }
  while (++curr != first);

  return halfedges_end();
}

void Tspline::get_halfedges (Face_iterator &f, Point2d &param, Halfedge_iterator &right, Halfedge_iterator &left,
                             Halfedge_iterator &top, Halfedge_iterator &bottom)
{
  if (!f->has_outer_ccb() || f->is_unbounded())
    throw std::runtime_error ("[Tspline::get_halfedges] Error, face is unbound");

  Ccb_halfedge_circulator first, curr;
  first = curr = f->outer_ccb();

  right = left = top = bottom = halfedges_end();

  do
  {
    Point2d p0 = curr->source()->data().param;
    Point2d p1 = curr->target()->data().param;

    Vector2d v01 = p1 - p0;
    //if (v01.y() > 0 && param.y() >= p0.y() && param.y() <= p1.y()) // right
    if(greater(v01.y(), 0.0) && gequal(param.y(), p0.y()) && sequal(param.y(), p1.y()))
      right = curr;
    //if (v01.y() < 0 && param.y() <= p0.y() && param.y() >= p1.y()) // left
    if(smaller(v01.y(), 0.0) && sequal(param.y(), p0.y()) && gequal(param.y(), p1.y()))
      left = curr;
    //if (v01.x() < 0 && param.x() <= p0.x() && param.x() >= p1.x()) // top
    if(smaller(v01.x(), 0.0) && sequal(param.x(), p0.x()) && gequal(param.x(), p1.x()))
      top = curr;
    //if (v01.x() > 0 && param.x() >= p0.x() && param.x() <= p1.x()) // bottom
    if(greater(v01.x(), 0.0) && gequal(param.x(), p0.x()) && sequal(param.x(), p1.x()))
      bottom = curr;
  }
  while (++curr != first);

}


Eigen::Vector2d Tspline::closest_point(const Eigen::Vector3d &point, const Eigen::Vector2d &hint, unsigned iterations,
                                       double accuracy) const
{
  Eigen::Vector3d p, ts, tt;
  return closest_point(point, hint, p, ts, tt, iterations, accuracy);
}

// todo replace relations with tolerant relation funtions in Math.hpp
Eigen::Vector2d Tspline::closest_point(const Eigen::Vector3d &point, const Eigen::Vector2d &hint, Eigen::Vector3d &p,
                                       Eigen::Vector3d &ts, Eigen::Vector3d &tt, unsigned iterations, double accuracy) const
{
  Eigen::Vector2d curr = hint;
  Eigen::Vector3d r;
  Eigen::Vector2d b, d;
  Eigen::Matrix2d A;

  for (unsigned i = 0; i < iterations; i++) {
    p = evaluate(curr(0), curr(1), ts, tt);

    if (isnan(p(0)) || isnan(ts(0)) || isnan(tt(0))) {
      printf("[Tspline::closest_point] Error, method did not converge (nan) (%e %d)\n", accuracy, iterations);
      printf("  point: %f %f %f\n", point(0), point(1), point(2));
      printf("  %f %f ... %f %f  p: %f %f %f\n", hint(0), hint(1), curr(0), curr(1), p(0), p(1), p(2));
      throw std::runtime_error("[Tspline::closest_point] Error, method did not converge");
    }

    r = p - point;

    Eigen::Vector3d null(0.0, 0.0, 0.0);
    if (equal(ts, null) && equal(tt, null)) {
      return curr;
    } else if (equal(ts, null)) {
      d(0) = 0.0;
      d(1) = -r.dot(tt) / tt.dot(tt);
    } else if (equal(tt, null)) {
      d(0) = -r.dot(ts) / ts.dot(ts);
      d(1) = 0.0;
    } else {
      b(0) = -r.dot(ts);
      b(1) = -r.dot(tt);

      A(0, 0) = ts.dot(ts);
      A(0, 1) = ts.dot(tt);
      A(1, 0) = A(0, 1);
      A(1, 1) = tt.dot(tt);

      d = A.ldlt().solve(b);
    }

    if (d.norm() < accuracy) {
      return curr;
    } else {
      curr = curr + d * 0.25; // todo step-width control

      if (curr(0) < param_min.x())
        curr(0) = param_min.x();
      else if (curr(0) > param_max.x())
        curr(0) = param_max.x();

      if (curr(1) < param_min.y())
        curr(1) = param_min.y();
      else if (curr(1) > param_max.y())
        curr(1) = param_max.y();
    }

  }

  if (!quiet) {
    printf("[Tspline::closest_point] Warning: Method did not converge (%e %d)\n", accuracy, iterations);
    printf("  %f %f ... %f %f  p: %f %f %f\n", hint(0), hint(1), curr(0), curr(1), p(0), p(1), p(2));
  }

  return hint;
}

vector_vec3d Tspline::compute_cp_normals_by_cage()
{
  vector_vec3d normals;

  for(Tspline::Vertex_iterator vit=vertices_begin(); vit!=vertices_end(); vit++)
  {
    const TVertex &vext = vit->data();
    Eigen::Vector3d cp;
    vext.GetCP(cp);
    if(isnan(cp(0)))
      printf("[TsplineReprojection::compute_cp_normals_by_cage A] %f %f %f\n", cp(0),cp(1),cp(2));
  }

  for(Tspline::Vertex_iterator vitA=vertices_begin(); vitA!=vertices_end(); vitA++)
  {
    Eigen::Vector3d a;
    vitA->data().GetCP(a);
    Eigen::Vector3d normal(0,0,0);
    unsigned num_normals(0);

    Tspline::Halfedge_around_vertex_circulator first, curr, next;
    first = curr = vitA->incident_halfedges();
    do
    {
      Eigen::Vector3d b, c, e0, e1, n;

      next = curr;
      next++;

      if((!curr->face()->is_unbounded() || !next->twin()->face()->is_unbounded()))
      {
        curr->source()->data().GetCP(b);
        next->source()->data().GetCP(c);

        e0 = b - a;
        e1 = c - a;

        n = e1.cross(e0);

        double norm = n.norm();
        if(!equal(norm,0.0))
        {
          n /= norm;
          n.normalize();
          normal += n;
          num_normals++;
        }
      }

    }
    while (++curr != first);

    normal /= num_normals;
    normal.normalize();
    normals.push_back(normal);

    //    if(viewer!=NULL)
    //    {
    //      Eigen::Vector3d &n = normal;
    //      float s(0.01);
    //      viewer->AddLine3D(a(0),a(1),a(2),a(0)+n(0)*s,a(1)+n(1)*s,a(2)+n(2)*s, 255, 255, 255, 4.0f);
    //    }
  }

  return normals;
}

vector_vec3d Tspline::compute_cp_normals_by_footpoints()
{
  vector_vec3d normals;

  for(Tspline::Vertex_iterator vitA=vertices_begin(); vitA!=vertices_end(); vitA++)
  {
    const TVertex& vext = vitA->data();
    Eigen::Vector3d ts, tt, n;

    double s=vext.param.x();
    double t=vext.param.y();
    if(equal(s,param_min.x()))
      s+=0.01;//(10.0*epsilon);
    if(equal(s,param_max.x()))
      s-=0.01;//(10.0*epsilon);
    if(equal(t,param_min.y()))
      t+=0.01;//(10.0*epsilon);
    if(equal(t,param_max.y()))
      t-=0.01;//(10.0*epsilon);

    evaluate(s, t, ts, tt);
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

CGAL::Object Tspline::locate_param(const double &s, const double &t)
{
  return locate_param(Point2d(s,t));
}

CGAL::Object Tspline::locate_param (Point2d param)
{
  if (smaller(param.x(), param_min.x()) || greater(param.x(), param_max.x()))
    printf("[Tspline::locate_param] Warning, parameter.x() out of bounds (%e; %e %e)\n", param.x(), param_min.x(), param_max.x());
  if (smaller(param.y(), param_min.y()) || greater(param.y(), param_max.y()))
    printf("[Tspline::locate_param] Warning, parameter.y() out of bounds (%e; %e %e)\n", param.y(), param_min.y(), param_max.y());

  // find closest vertex
  Vertex_iterator v_close = get_closest_vertex (param);

  const TVertex& v_close_ext = v_close->data();

  //printf("[Tspline::locate_param] v_close: %e %e v: %e %e\n", v_close_ext.param.x(), v_close_ext.param.y(), s, t);
  //printf("[Tspline::locate_param] equal: %d %d\n", equal(v_close_ext.param.x(), param.x()), equal(v_close_ext.param.y(), param.y()));
  if(equal(v_close_ext.param.x(), param.x()) && equal(v_close_ext.param.y(), param.y()))
  {
    //printf("[Tspline::locate_param] vertex found\n");
    return CGAL::make_object (v_close);
  }



  // iterate around vertex to find face, where param point is lying in
  Halfedge_around_vertex_circulator first, curr;
  first = curr = v_close->incident_halfedges();
  Halfedge_iterator he_right, he_left, he_top, he_bottom;
  Face_iterator f;
  do
  {
    f = curr->face();
    if (!f->is_unbounded())
    {
      //printf("[Tspline::locate_param] f->is_bounded\n");

      get_halfedges (f, param, he_right, he_left, he_top, he_bottom);
      if (he_right != halfedges_end() && he_left != halfedges_end() && he_top != halfedges_end() && he_bottom
          != halfedges_end())
      {
        //printf("[Tspline::locate_param] halfedges found\n");

        // Check if param point eventually lies on edge or vertex
        // right
        Point2d p0 = he_right->source()->data().param;
        Point2d p1 = he_right->target()->data().param;

        if (equal(param.x(), p0.x())) //if (param.x() == p0.x())
        {
          if (equal(param.y(), p0.y())) //if (param.y() == p0.y())
            return CGAL::make_object (he_right->source()); // vertex
          else if (equal(param.y(), p1.y())) //else if (param.y() == p1.y())
            return CGAL::make_object (he_right->target()); // vertex
          else
            return CGAL::make_object (he_right); // edge
        }
        else
        {
          // left
          p0 = he_left->source()->data().param;
          p1 = he_left->target()->data().param;
          if (equal(param.x(), p0.x())) //if (param.x() == p0.x())
          {
            if (equal(param.y(), p0.y())) //if (param.y() == p0.y())
              return CGAL::make_object (he_left->source()); // vertex
            else if (equal(param.y(), p1.y())) //else if (param.y() == p1.y())
              return CGAL::make_object (he_left->target()); // vertex
            else
              return CGAL::make_object (he_left); // edge
          }
          else
          {
            // top
            p0 = he_top->source()->data().param;
            p1 = he_top->target()->data().param;
            if (equal(param.y(), p0.y())) //if (param.y() == p0.y())
            {
              if (equal(param.x(), p0.x())) //if (param.x() == p0.x())
                return CGAL::make_object (he_top->source()); // vertex
              else if (equal(param.x(), p1.x())) //if (param.x() == p1.x())
                return CGAL::make_object (he_top->target()); // vertex
              else
                return CGAL::make_object (he_top); // edge
            }
            else
            {
              // bottom
              p0 = he_bottom->source()->data().param;
              p1 = he_bottom->target()->data().param;
              if (equal(param.y(), p0.y())) //if (param.y() == p0.y())
              {
                if (equal(param.x(), p0.x())) //if (param.x() == p0.x())
                  return CGAL::make_object (he_bottom->source()); // vertex
                else if (equal(param.x(), p1.x())) // if (param.x() == p1.x())
                  return CGAL::make_object (he_bottom->target()); // vertex
                else
                  return CGAL::make_object (he_bottom); // edge
              }
            }
          }
        }

        return CGAL::make_object (f); // face
      }
    }
  }
  while (++curr != first);

  return CGAL::Object();
}

CGAL::Object Tspline::locate_point(const Point2d& p)
{
  NaivePointLocation npl(*this);

  return npl.locate(p);
}

Eigen::Vector3d Tspline::center() const
{
  Eigen::Vector3d meanCP(0,0,0);
  Tspline::Vertex_const_iterator vit;
  for(vit=vertices_begin(); vit!=vertices_end(); vit++)
  {
    Eigen::Vector3d cp;
    vit->data().GetCP(cp);
    meanCP += cp;
  }
  meanCP /= number_of_vertices();
  return meanCP;
}

bool Tspline::is_rational() const
{
  for(Vertex_const_iterator v=vertices_begin(); v!=vertices_end(); v++)
    if(!equal(1.0, v->data().GetCP().w()))
      return true;

  return false;
}

bool Tspline::is_NURBS() const
{
  for(Vertex_const_iterator v=vertices_begin(); v!=vertices_end(); v++)
  {
    if(is_corner(v))
      continue;

    bool bnd = is_boundary(v);

    if(bnd && v->degree()!=3)
      return false;

    if(!bnd && v->degree()!=4)
      return false;
  }

  return true;
}

size_t Tspline::number_of_controlpoints() const
{
  size_t ncps(0);

  Vertex_const_iterator vit;
  for(vit=vertices_begin(); vit!=vertices_end(); vit++)
    if(vit->data().is_primary)
      ncps++;

  return ncps;
}

void Tspline::get_bounding_box(Eigen::Vector3d& bb_min, Eigen::Vector3d& bb_max) const
{
  bb_min = Eigen::Vector3d(DBL_MAX, DBL_MAX, DBL_MAX);
  bb_max = Eigen::Vector3d(-DBL_MAX, -DBL_MAX, -DBL_MAX);

  Vertex_const_iterator vit;
  for(vit=vertices_begin(); vit!= vertices_end(); vit++)
  {
    Eigen::Vector3d cp;
    vit->data().GetCP(cp);
    for(int i=0; i<3; i++)
    {
      if(cp(i)>bb_max(i))
        bb_max(i)=cp(i);
      if(cp(i)<bb_min(i))
        bb_min(i)=cp(i);
    }
  }
}

void Tspline::get_NURBS_dimension(unsigned &width, unsigned &height)
{
  if(!is_NURBS())
    throw std::runtime_error("[Tspline::get_NURBS_dimension] Error, T-spline is not a NURBS.");

  Face_iterator fit;
  for(fit=faces_begin(); fit!=faces_end(); fit++)
    if(face_is_bottom_left_corner(fit))
      break;

  if(fit==faces_end())
    throw std::runtime_error("[Tspline::get_NURBS_dimension] Error, bottom-left face not found.");

  Eigen::Vector2d ec = face_center(fit);
  Point2d c(ec(0),ec(1));

  Face_iterator f_right=fit;
  width = 0;
  while(!f_right->is_unbounded())
  {
    f_right = get_right_halfedge_at_point(f_right, c)->twin()->face();
    width++;
  }

  Face_iterator f_top=fit;
  height = 0;
  while(!f_top->is_unbounded())
  {
    f_top = get_top_halfedge_at_point(f_top, c)->twin()->face();
    height++;
  }
}

Tspline::Vertex_iterator Tspline::get_NURBS_vertex(unsigned col, unsigned row)
{
  Vertex_iterator vit;
  for(vit=vertices_begin(); vit!=vertices_end(); vit++)
    if(is_bottom_left_vertex(vit))
      break;

  if(vit==vertices_end())
    throw std::runtime_error("[Tspline::get_NURBS_vertex] Error, no bottom left vertex found.");

  Vertex_iterator v;
  unsigned c(0);
  for(v=vit; v!=vertices_end() && c<col; c++)
  {
    Halfedge_iterator h = get_right_halfedge(v);
    if(h==halfedges_end())
      throw std::runtime_error("[Tspline::get_NURBS_vertex] Error, index out of bounds (right).");
    v = h->source();
  }

  unsigned r(0);
  for(v=v; v!=vertices_end() && r<row; r++)
  {
    Halfedge_iterator h = get_top_halfedge(v);
    if(h==halfedges_end())
      throw std::runtime_error("[Tspline::get_NURBS_vertex] Error, index out of bounds (top).");
    v = h->source();
  }

  return v;
}

Tspline::Face_iterator Tspline::locate_face(const double &s, const double &t)
{
  CGAL::Object obj = locate_param(s, t);

  Face_iterator f;

  if (!obj.empty())
    if (CGAL::assign(f, obj))
      return f;

  return faces_end();
}

Eigen::Vector2d Tspline::face_center(Face_const_handle f) const
{
  if(f->is_unbounded())
    throw std::runtime_error("[Tspline::face_center] Error, face is unbounded\n");
  if(!f->has_outer_ccb())
    throw std::runtime_error("[Tspline::face_center] Error, has no outer contour\n");

  Tspline::Ccb_halfedge_const_circulator curr, first;
  curr = first = f->outer_ccb();
  Eigen::Vector2d param = Eigen::Vector2d(0,0);

  size_t n(0);
  do
  {
    const Point2d &p = curr->source()->data().param;
    param += Eigen::Vector2d(p.x(), p.y());
    n++;
  }while(++curr!=first);

  param /= n;
  return param;
}

Eigen::Vector2d Tspline::face_center(Face_const_handle f, Eigen::Vector3d &point) const
{
  Eigen::Vector2d param = face_center(f);
  point = evaluate(param[0], param[1]);
  return param;
}

void Tspline::face_bounding_box(Face_const_handle f,
                                Point2d &bbmin, Point2d &bbmax) const
{
  Eigen::Vector4d bb;
  face_bounding_box(f, bb);
  bbmin = Point2d(bb(0),bb(2));
  bbmax = Point2d(bb(1),bb(3));
}

void Tspline::face_bounding_box(Face_const_handle f,
                                Eigen::Vector4d &bb) const
{
  if (f->is_unbounded() || !f->has_outer_ccb())
    throw std::runtime_error("[Tspline::face_bounding_box] Error, face is unbounded or has no outer contour\n");

  Ccb_halfedge_const_circulator first, curr;
  first = curr = f->outer_ccb();

  bb(0) = DBL_MAX;
  bb(1) = -DBL_MAX;
  bb(2) = DBL_MAX;
  bb(3) = -DBL_MAX;

  do
  {
    Point2d p0 = curr->source()->data().param;

    if (p0.x() < bb(0))
      bb(0) = p0.x();
    if (p0.x() > bb(1))
      bb(1) = p0.x();

    if (p0.y() < bb(2))
      bb(2) = p0.y();
    if (p0.y() > bb(3))
      bb(3) = p0.y();
  }
  while (++curr != first);
}

void Tspline::face_bounding_box_by_point(Face_const_handle f,
                                         Point2d &bbmin, Point2d &bbmax) const
{
  Eigen::Vector4d bb;
  face_bounding_box_by_point(f, bb);
  bbmin = Point2d(bb(0),bb(2));
  bbmax = Point2d(bb(1),bb(3));
}

void Tspline::face_bounding_box_by_point(Face_const_handle f,
                                         Eigen::Vector4d &bb) const
{
  if (f->is_unbounded() || !f->has_outer_ccb())
    throw std::runtime_error("[Tspline::face_bounding_box] Error, face is unbounded or has no outer contour\n");

  Ccb_halfedge_const_circulator first, curr;
  first = curr = f->outer_ccb();

  bb(0) = DBL_MAX;
  bb(1) = -DBL_MAX;
  bb(2) = DBL_MAX;
  bb(3) = -DBL_MAX;

  do
  {
    Point2d p0 = curr->source()->point();

    if (p0.x() < bb(0))
      bb(0) = p0.x();
    if (p0.x() > bb(1))
      bb(1) = p0.x();

    if (p0.y() < bb(2))
      bb(2) = p0.y();
    if (p0.y() > bb(3))
      bb(3) = p0.y();
  }
  while (++curr != first);
}

bool Tspline::face_is_collapsed(Face_const_handle f) const
{
  Eigen::Vector4d bb;
  face_bounding_box(f,bb);

  if(equal(bb(1)-bb(0),0.0) || equal(bb(3)-bb(2),0.0))
    return true;

  return false;
}

bool Tspline::face_is_corner(Face_const_handle f) const
{
  if(f->is_unbounded())
    return false;

  Tspline::Ccb_halfedge_const_circulator first, curr;
  first = curr = f->outer_ccb();

  bool has_vertical_boundary(false);
  bool has_horizontal_boundary(false);

  do
  {
    if(curr->twin()->face()->is_unbounded())
    {
      const Point2d& a = curr->source()->point();
      const Point2d& b = curr->target()->point();

      if( equal(a.x(),b.x()) )
        has_horizontal_boundary = true;

      if( equal(a.y(),b.y()) )
        has_vertical_boundary = true;

      if(has_vertical_boundary && has_horizontal_boundary)
        return true;
    }

  }while(++first!=curr);

  return false;
}


bool Tspline::face_is_bottom_left_corner(Face_const_handle f) const
{
  if(f->is_unbounded())
    return false;

  Tspline::Ccb_halfedge_const_circulator first, curr;
  first = curr = f->outer_ccb();

  bool left(false);
  bool bottom(false);

  do
  {
    if(curr->twin()->face()->is_unbounded())
    {
      const Point2d& a = curr->source()->point();
      const Point2d& b = curr->target()->point();

      if( greater(a.y(),b.y()) )
        left = true;

      if( smaller(a.x(),b.x()) )
        bottom = true;

      if(left && bottom)
        return true;
    }

  }while(++curr!=first);

  return false;
}


void Tspline::face_dimensions_squared(Face_handle f, double& ds, double& dt, const Eigen::Vector4d &bb)
{
  // midpoints of edges of bounding-box
  Point2d a( bb(0), 0.5*(bb(2)+bb(3)) );
  Point2d b( bb(1), 0.5*(bb(2)+bb(3)) );
  Point2d c( 0.5*(bb(0)+bb(1)), bb(2) );
  Point2d d( 0.5*(bb(0)+bb(1)), bb(3) );

  Eigen::Vector3d pa = evaluate(a.x(),a.y());
  Eigen::Vector3d pb = evaluate(b.x(),b.y());
  Eigen::Vector3d pc = evaluate(c.x(),c.y());
  Eigen::Vector3d pd = evaluate(d.x(),d.y());

  ds = (pa-pb).squaredNorm();
  dt = (pc-pd).squaredNorm();
}

void Tspline::face_dimensions_squared(Face_handle f, double& ds, double& dt)
{
  Eigen::Vector4d bb;
  face_bounding_box(f, bb);

  // midpoints of edges of bounding-box
  Point2d a( bb(0), 0.5*(bb(2)+bb(3)) );
  Point2d b( bb(1), 0.5*(bb(2)+bb(3)) );
  Point2d c( 0.5*(bb(0)+bb(1)), bb(2) );
  Point2d d( 0.5*(bb(0)+bb(1)), bb(3) );

  Eigen::Vector3d pa = evaluate(a.x(),a.y());
  Eigen::Vector3d pb = evaluate(b.x(),b.y());
  Eigen::Vector3d pc = evaluate(c.x(),c.y());
  Eigen::Vector3d pd = evaluate(d.x(),d.y());

  ds = (pa-pb).squaredNorm();
  dt= (pc-pd).squaredNorm();
}


Tspline::Face_iterator Tspline::face(Face_const_iterator fc)
{
  if (fc->is_unbounded() || !fc->has_outer_ccb())
    throw std::runtime_error("[Tspline::face] Error, face is unbounded or has no outer contour\n");

  Face_iterator fit;
  for(fit=faces_begin(); fit!=faces_end(); fit++)
  {
    if(fc==Face_const_iterator(fit))
      return fit;
  }

  throw std::runtime_error("[Tspline::face] Error, face not found\n");
}

Tspline::Vertex_iterator Tspline::get_vertex(const int& id)
{
  for(Vertex_iterator vit=vertices_begin(); vit!=vertices_end(); vit++)
  {
    if(vit->data().id==id)
      return vit;
  }
  return vertices_end();
}

bool Tspline::is_boundary (Vertex_const_iterator v) const
{
  Halfedge_around_vertex_const_circulator first, curr;
  first = curr = v->incident_halfedges();

  do
  {
    if (curr->face()->is_unbounded())
      return true;
  }
  while (++curr != first);

  return false;
}

bool Tspline::is_boundary (Vertex_const_iterator v, Halfedge_const_iterator &he) const
{
  Halfedge_around_vertex_const_circulator first, curr;
  first = curr = v->incident_halfedges();

  do
  {
    if (curr->face()->is_unbounded())
    {
      he = curr;
      return true;
    }
  }
  while (++curr != first);

  return false;
}

bool Tspline::is_corner(Vertex_const_iterator v) const
{
  if(v->degree()!=2)
    return false;

  if(!is_Ljunction(v))
    return false;

  Halfedge_around_vertex_const_circulator h1, h2;
  h1 = h2 = v->incident_halfedges();
  h2++;

  if(h1->face()->is_unbounded() && h2->twin()->face()->is_unbounded())
    return true;

  if(h1->twin()->face()->is_unbounded() && h2->face()->is_unbounded())
    return true;

  return false;
}

bool Tspline::is_Xjunction(Vertex_const_iterator v) const
{
  if(v->degree()==4)
    return true;

  return false;
}

bool Tspline::is_Tjunction(Vertex_const_iterator v) const
{
  if(v->degree()!=3)
    return false;

  if(is_boundary(v))
    return false;

  return true;
}

bool Tspline::is_Ljunction(Vertex_const_iterator v) const
{
  if(v->degree()!=2)
    return false;

  Halfedge_around_vertex_const_circulator circ = v->incident_halfedges();
  const Point2d& p1a = circ->source()->point();
  const Point2d& p1b = circ->target()->point();
  circ++;
  const Point2d& p2a = circ->source()->point();
  const Point2d& p2b = circ->target()->point();

  // if edge 1 is vertical and edge 2 is horizontal
  if( equal(p1a.x(),p1b.x()) && equal(p2a.y(),p2b.y()) )
    return true;

  // if edge 1 is horizontal and edge 2 is vertical
  if( equal(p1a.y(),p1b.y()) && equal(p2a.x(),p2b.x()) )
    return true;

  return false;
}

bool Tspline::is_Ijunction(Vertex_const_iterator v) const
{
  if(v->degree()!=2)
    return false;

  Halfedge_around_vertex_const_circulator circ = v->incident_halfedges();
  const Point2d& p1a = circ->source()->point();
  const Point2d& p1b = circ->target()->point();
  circ++;
  const Point2d& p2a = circ->source()->point();
  const Point2d& p2b = circ->target()->point();

  // if edge 1 is vertical and edge 2 is vertical
  if( equal(p1a.x(),p1b.x()) && equal(p2a.x(),p2b.x()) )
    return true;

  // if edge 1 is horizontal and edge 2 is horizontal
  if( equal(p1a.y(),p1b.y()) && equal(p2a.y(),p2b.y()) )
    return true;

  return false;
}

bool Tspline::is_bottom_left_vertex(Vertex_const_iterator v) const
{
  if(v->degree()!=2)
    return false;

  if(!is_Ljunction(v))
    return false;

  Halfedge_around_vertex_const_circulator h1, h2;
  h1 = h2 = v->incident_halfedges();
  h2++;

  const Point2d& h1s = h1->source()->point();
  const Point2d& h1t = h1->target()->point();
  const Point2d& h2s = h2->source()->point();
  const Point2d& h2t = h2->target()->point();

  if(h1->face()->is_unbounded() && h2->twin()->face()->is_unbounded())
    if( greater(h1s.x(),h1t.x()) && greater(h2s.y(),h2t.y()) )
      return true;

  if(h1->twin()->face()->is_unbounded() && h2->face()->is_unbounded())
    if( greater(h2s.x(),h2t.x()) && greater(h1s.y(),h1t.y()) )
      return true;

  return false;
}

Tspline::Vertex_iterator Tspline::get_closest_vertex_by_point(const Point2d &point, double& d_sqr)
{
  Vertex_iterator v_close;
  d_sqr = DBL_MAX;
  for (Vertex_iterator v = vertices_begin(); v != vertices_end(); v++)
  {
    double d = (v->point() - point).squared_length();
    if (smaller(d, d_sqr))
    {
      d_sqr = d;
      v_close = v;
    }
  }
  return v_close;
}

Tspline::Halfedge_iterator Tspline::get_closest_edge_by_point(const Point2d &point, double& d_sqr)
{
  Halfedge_iterator he_close = halfedges_end();
  d_sqr = DBL_MAX;
  for(Halfedge_iterator h=edges_begin(); h!=edges_end(); h++)
  {
    Point2d p0 = h->source()->point();
    Point2d p1 = h->target()->point();

    double d(DBL_MAX);
    if(equal(p0.x(),p1.x()))
      if( (smaller(point.y(),p0.y()) && greater(point.y(),p1.y())) ||
          (smaller(point.y(),p1.y()) && greater(point.y(),p0.y())))
        d = std::abs<double>(p0.x()-point.x());

    if(equal(p0.y(),p1.y()))
      if( (smaller(point.x(),p0.x()) && greater(point.x(),p1.x())) ||
          (smaller(point.x(),p1.x()) && greater(point.x(),p0.x())))
        d = std::abs<double>(p0.y()-point.y());

    if(d<d_sqr)
    {
      d_sqr = d;
      he_close = h;
    }
  }

  d_sqr *= d_sqr;
  return he_close;
}

Tspline::Face_const_iterator Tspline::get_face_by_point(const Point2d &point)
{
  CGAL::Object obj = locate_point(point);

  Face_const_iterator f;

  if (!obj.empty())
    if (CGAL::assign(f, obj))
      return f;

  return faces_end();
}

Tspline::Vertex_iterator Tspline::get_closest_vertex (const Point2d &param)
{
  Vertex_iterator v_close;
  double d_sqr (DBL_MAX);
  for (Vertex_iterator v = vertices_begin(); v != vertices_end(); v++)
  {
    double d = (v->data().param - param).squared_length();
    if (smaller(d, d_sqr))
    {
      d_sqr = d;
      v_close = v;
    }
  }
  return v_close;
}

Eigen::Vector4d Tspline::grid_normal_curvature(Vertex_iterator v) const
{
  if(v->degree()<2)
    return Eigen::Vector4d(0.0,0.0,0.0,0.0);

  Halfedge_around_vertex_const_circulator curr, first;
  curr = first = v->incident_halfedges();

  Point3d p0 = v->data().GetCP();
  vector_vec3d edges;

  do
  {
    Point3d p1 = curr->source()->data().GetCP();
    Eigen::Vector3d e( p1.x()-p0.x(), p1.y()-p0.y(), p1.z()-p0.z() );
    e.normalize();
    edges.push_back(e);
  }
  while(++curr!=first);

  double sum_theta(0.0);
  Eigen::Vector3d normal(0.0,0.0,0.0);
  size_t es = edges.size();
  for(size_t i=0; i<es; i++)
  {
    const Eigen::Vector3d& e1 = edges[i];
    const Eigen::Vector3d& e2 = edges[(i+1)%es];

    normal += e2.cross(e1);
    sum_theta += acos( e2.dot(e1) );
  }

  normal /= es;

  double curvature = 2.0 * M_PI - sum_theta;
  if(v->degree()<3)
    curvature = 0.0;

  return Eigen::Vector4d(normal(0),normal(1),normal(2),curvature);
}

double Tspline::edge_length_by_point(Halfedge_const_iterator h)
{
  const Point2d& a = h->source()->point();
  const Point2d& b = h->target()->point();
  return sqrt((a-b).squared_length());
}

Tspline::Halfedge_handle Tspline::get_shortest_edge()
{
  double sqr_length(DBL_MAX);

  Halfedge_iterator hit;
  Halfedge_handle hh;
  for(hit=edges_begin(); hit!=edges_end(); hit++)
  {
    const Point2d& a = hit->source()->data().param;
    const Point2d& b = hit->source()->data().param;

    double d = (b-a).squared_length();
    if(d<sqr_length)
    {
      sqr_length=d;
      hh = hit;
    }
  }

  return hh;
}

bool Tspline::is_connectable(Vertex_iterator v0, Vertex_iterator v1)
{
  const Point2d& p0 = v0->point();
  const Point2d& p1 = v1->point();

  if(p0.x()==p1.x() && p0.y()==p1.y())
  {
    printf("[Tspline::is_connectable] Error, point identical (%d %d)\n.", v0->data().id, v1->data().id);
    throw std::runtime_error("[Tspline::is_connectable] Error, point identical.");
  }

  bool connectable(false);
  if(p0.y()==p1.y())
  {
    if(p1.x()>p0.x() && get_right_halfedge(v0)==halfedges_end())
      connectable = true;
    else if(p0.x()>p1.x() && get_left_halfedge(v0)==halfedges_end())
      connectable = true;
  }
  else if(p0.x()==p1.x())
  {
    if(p1.y()>p0.y() && get_top_halfedge(v0)==halfedges_end())
      connectable = true;
    else if(p0.y()>p1.y() && get_bottom_halfedge(v0)==halfedges_end())
      connectable = true;
  }

  return connectable;
}

Arrangement_2::Halfedge_iterator Tspline::insert_missing_edges(Face_iterator fit)
{
  if(fit->is_unbounded() || !fit->has_outer_ccb())
    return halfedges_end();

  Hole_iterator hole; size_t num_holes(0);
  for(hole=fit->holes_begin(); hole!=fit->holes_end(); hole++)
    num_holes++;

  if(num_holes>0)
    printf("[Tspline::insert_missing_edges] Warning, T-spline face has holes.\n");

  //  {
  //    printf("[Tspline::insert_missing_edges(fit)]");
  //    Ccb_halfedge_circulator first, i, j;
  //    first = i = fit->outer_ccb();
  //    do
  //    {
  //      printf(" %d", i->source()->data().id);

  //    }while(++i!=first);
  //    printf("   .\n");
  //  }

  Ccb_halfedge_circulator first, i;
  first = i = fit->outer_ccb();

  VertexPair doubleLedge, Ledge, edge;
  double doubleLdist(DBL_MAX), Ldist(DBL_MAX), dist(DBL_MAX);
  bool doubleLedge_found(false), Ledge_found(false), edge_found(false);

  do
  {
    Vertex_iterator v0 = i->source();
    Ccb_halfedge_circulator j = i;
    while(++j!=first)
    {
      Vertex_iterator v1 = j->source();

      if(is_connectable(v0,v1))
      {
        double d = (v1->point()-v0->point()).squared_length();
        bool v0_is_L = is_Ljunction(v0);
        bool v1_is_L = is_Ljunction(v1);
        if(v0_is_L && v1_is_L && d<doubleLdist)
        {
          doubleLedge = VertexPair(v0,v1);
          doubleLedge_found = true;
          doubleLdist = d;
        }
        else if( (v0_is_L || v1_is_L) && d<Ldist )
        {
          Ledge = VertexPair(v0,v1);
          Ledge_found = true;
          Ldist = d;
        }
        else if( d<dist )
        {
          edge = VertexPair(v0,v1);
          edge_found = true;
          dist = d;
        }
      }

    } // while(++j!=first)
  }while (++i!=first);

  // connect edges with two Ljunctions (priority 1)
  if(doubleLedge_found)
  {
    Vertex_iterator& v0 = doubleLedge.first;
    Vertex_iterator& v1 = doubleLedge.second;
    Halfedge_iterator he = insert_at_vertices(Segment2(v0->point(),v1->point()), v0, v1);
    he->data().d = sqrt((v1->data().param-v0->data().param).squared_length());
    he->twin()->data().d = he->data().d;
    return he;
  }

  // connect edges with one Ljunction (priority 2)
  if(Ledge_found)
  {
    Vertex_iterator& v0 = Ledge.first;
    Vertex_iterator& v1 = Ledge.second;
    Halfedge_iterator he = insert_at_vertices(Segment2(v0->point(),v1->point()), v0, v1);
    he->data().d = sqrt((v1->data().param-v0->data().param).squared_length());
    he->twin()->data().d = he->data().d;
    return he;
  }

  // connect edges (priority 3)
  if(edge_found)
  {
    Vertex_iterator& v0 = edge.first;
    Vertex_iterator& v1 = edge.second;
    Halfedge_iterator he = insert_at_vertices(Segment2(v0->point(),v1->point()), v0, v1);
    he->data().d = sqrt((v1->data().param-v0->data().param).squared_length());
    he->twin()->data().d = he->data().d;
    return he;
  }

  return halfedges_end();
}

void Tspline::insert_missing_edges()
{
  bool fit_prev_valid(false);
  Face_iterator fit_prev; // last face without edge insertion
  for(Face_iterator fit=faces_begin(); fit!=faces_end(); fit++)
  {
    if(fit->is_unbounded())
      continue;

    if(insert_missing_edges(fit)==halfedges_end()) // no edge inserted
    {
      fit_prev = fit;
      fit_prev_valid = true;
    }
    else // edge inserted
    {
      if(fit_prev_valid)
        fit = fit_prev;
      else
        fit = faces_begin();
    }
  }
}

void Tspline::insert_edges_at_L_junctions()
{
  Vertex_iterator vit;
  for(vit=vertices_begin(); vit!=vertices_end(); vit++)
  {
    if(is_corner(vit) || !is_Ljunction(vit))
      continue;

    const Point2d& p = vit->point();

    //    printf("[Tspline::insert_edges_at_L_junctions] insert vertex %d\n", vit->data().id);

    // get incident halfedges
    Halfedge_iterator h_bottom = get_bottom_halfedge(vit);
    Halfedge_iterator h_top = get_top_halfedge(vit);
    Halfedge_iterator h_left = get_left_halfedge(vit);
    Halfedge_iterator h_right = get_right_halfedge(vit);

    if( (h_bottom==halfedges_end() && h_top==halfedges_end()) ||
        (h_bottom!=halfedges_end() && h_top!=halfedges_end()) ||
        (h_left==halfedges_end() && h_right==halfedges_end()) ||
        (h_left!=halfedges_end() && h_right!=halfedges_end()) )
      throw std::runtime_error("[Tspline::insert_edges_at_L_junctions] Error, no L junction.");

    Halfedge_iterator edge;
    bool edge_found(false);
    double edge_dist(DBL_MAX);
    Point2d param;
    Vertex_iterator vn = vertices_end();
    Face_iterator fit;

    // bottom
    if(h_bottom==halfedges_end())
    {
      if(h_left==halfedges_end())
        fit = h_right->face();
      else
        fit = h_left->twin()->face();

      Halfedge_iterator e = get_bottom_halfedge_at_point(fit, p);
      double d = vit->point().y() - e->source()->point().y();
      if(d<edge_dist)
      {
        edge = e;
        edge_dist = d;
        param = Point2d(vit->data().param.x(), e->source()->data().param.y());
        edge_found = true;
        if(equal(e->source()->point().x(), p.x()))
          vn = e->source();
        else if(equal(e->target()->point().x(),p.x()))
          vn = e->target();
      }
    }

    // top
    if(h_top==halfedges_end())
    {
      if(h_left==halfedges_end())
        fit = h_right->twin()->face();
      else
        fit = h_left->face();

      Halfedge_iterator e = get_top_halfedge_at_point(fit, p);
      double d = e->source()->point().y() - vit->point().y();
      if(d<edge_dist)
      {
        edge = e;
        edge_dist = d;
        param = Point2d(vit->data().param.x(), e->source()->data().param.y());
        edge_found = true;
        if(equal(e->source()->point().x(), p.x()))
          vn = e->source();
        else if(equal(e->target()->point().x(),p.x()))
          vn = e->target();
      }
    }

    // left
    if(h_left==halfedges_end())
    {
      if(h_top==halfedges_end())
        fit = h_bottom->face();
      else
        fit = h_top->twin()->face();

      Halfedge_iterator e = get_left_halfedge_at_point(fit, p);
      double d = vit->point().x() - e->source()->point().x();
      if(d<edge_dist)
      {
        edge = e;
        edge_dist = d;
        param = Point2d(e->source()->data().param.x(), vit->data().param.y());
        edge_found = true;
        if(equal(e->source()->point().y(), p.y()))
          vn = e->source();
        else if(equal(e->target()->point().y(),p.y()))
          vn = e->target();
      }
    }

    // right
    if(h_right==halfedges_end())
    {
      if(h_top==halfedges_end())
        fit = h_bottom->twin()->face();
      else
        fit = h_top->face();

      Halfedge_iterator e = get_right_halfedge_at_point(fit, p);
      double d = e->source()->point().x() - vit->point().x();
      if(d<edge_dist)
      {
        edge = e;
        edge_dist = d;
        param = Point2d(e->source()->data().param.x(), vit->data().param.y());
        edge_found = true;
        if(equal(e->source()->point().y(), p.y()))
          vn = e->source();
        else if(equal(e->target()->point().y(),p.y()))
          vn = e->target();
      }
    }

    if(edge_found)
    {
      //      printf("[Tspline::insert_edges_at_L_junctions] %d %d\n",
      //             edge->source()->data().id, edge->target()->data().id);

      if(vn==vertices_end())
        vn = insert_vertex( edge, param );

      Halfedge_iterator hn = insert_at_vertices(Segment2(vit->point(), vn->point()), vit, vn);
      hn->data().d = sqrt( (vit->data().param - vn->data().param).squared_length() );
      hn->twin()->data().d = hn->data().d;
    }

  } // for(vit=vertices_begin(); vit!=vertices_end(); vit++)
}

void Tspline::compute_knot_vectors(Vertex_iterator vit)
{
  compute_knot_vectors(vit, vit->data().s, vit->data().t);
}

void Tspline::compute_knot_vectors(Tspline::Vertex_iterator vit,
                                   std::vector<double>& s, std::vector<double>& t)
{
  s.assign(5,0.0);
  t.assign(5,0.0);

  s[2] = vit->data().param.x();
  t[2] = vit->data().param.y();

  std::vector<double> knots;

  // left
  this->shoot_left(vit, 2, knots);
  s[0] = knots[1];
  s[1] = knots[0];

  // right
  this->shoot_right(vit, 2, knots);
  s[3] = knots[0];
  s[4] = knots[1];

  // up
  this->shoot_up(vit, 2, knots);
  t[3] = knots[0];
  t[4] = knots[1];

  // down
  this->shoot_down(vit, 2, knots);
  t[1] = knots[0];
  t[0] = knots[1];
}

bool Tspline::compute_violation_1(std::vector<BasisFunction*>& blend, int num_vertices_prev)
{
  bool violation(false);

  for (Vertex_iterator vit = vertices_begin(); vit != vertices_end(); vit++)
  {
    const TVertex& vext = vit->data();

    if(vext.id >= num_vertices_prev) // do not treat newly inserted control points
      continue;

    std::vector<double> knots;
    std::vector<BasisFunction*>::iterator it = find(blend, vext.id);
    if(it!=blend.end())
    {
      bool v1 = (*it)->compute_violation_1();
      violation = violation || v1;
      continue;
    }

    // left
    this->shoot_left(vit, 2, knots);
    if(smaller(vext.s[1],knots[0]))
    {
      if(!quiet) printf("[Tspline::compute_violation_1] Violation 1 occured. A id: %d\n", vext.id);
      blend.push_back(new BasisFunction(this,vext.s,vext.t,vit));
      blend.back()->split_s(knots[0]);
      blend.back()->compute_violation_1();
      violation = true;
      continue;
    }else if(smaller(vext.s[0],knots[1]))
    {
      if(!quiet) printf("[Tspline::compute_violation_1] Violation 1 occured. B id: %d\n", vext.id);
      blend.push_back(new BasisFunction(this,vext.s,vext.t,vit));
      blend.back()->split_s(knots[1]);
      blend.back()->compute_violation_1();
      violation = true;
      continue;
    }

    // right
    this->shoot_right(vit, 2, knots);
    if(greater(vext.s[3],knots[0]))
    {
      if(!quiet) printf("[Tspline::compute_violation_1] Violation 1 occured. C id: %d\n", vext.id);
      blend.push_back(new BasisFunction(this,vext.s,vext.t,vit));
      blend.back()->split_s(knots[0]);
      blend.back()->compute_violation_1();
      violation = true;
      continue;
    }else if(greater(vext.s[4], knots[1]))
    {
      if(!quiet) printf("[Tspline::compute_violation_1] Violation 1 occured. D id: %d\n", vext.id);
      blend.push_back(new BasisFunction(this,vext.s,vext.t,vit));
      blend.back()->split_s(knots[1]);
      blend.back()->compute_violation_1();
      violation = true;
      continue;
    }

    // down
    this->shoot_down(vit, 2, knots);
    if(smaller(vext.t[1],knots[0]))
    {
      if(!quiet) printf("[Tspline::compute_violation_1] Violation 1 occured. E id: %d\n", vext.id);
      blend.push_back(new BasisFunction(this,vext.s,vext.t,vit));
      blend.back()->split_t(knots[0]);
      blend.back()->compute_violation_1();
      violation = true;
      continue;
    }else if(smaller(vext.t[0],knots[1]))
    {
      if(!quiet) printf("[Tspline::compute_violation_1] Violation 1 occured. F id: %d\n", vext.id);
      blend.push_back(new BasisFunction(this,vext.s,vext.t,vit));
      blend.back()->split_t(knots[1]);
      blend.back()->compute_violation_1();
      violation = true;
      continue;
    }

    // up
    this->shoot_up(vit, 2, knots);
    if(greater(vext.t[3],knots[0]))
    {
      if(!quiet) printf("[Tspline::compute_violation_1] Violation 1 occured. G\n");
      blend.push_back(new BasisFunction(this,vext.s,vext.t,vit));
      blend.back()->split_t(knots[0]);
      blend.back()->compute_violation_1();
      violation = true;
      continue;
    }else if(greater(vext.t[4], knots[1]))
    {
      if(!quiet) printf("[Tspline::compute_violation_1] Violation 1 occured. H\n");
      blend.push_back(new BasisFunction(this,vext.s,vext.t,vit));
      blend.back()->split_t(knots[1]);
      blend.back()->compute_violation_1();
      violation = true;
      continue;
    }
  }
  return violation;
}

bool Tspline::compute_violation_2(std::vector<BasisFunction*> &blend)
{
  // DEBUG ON
  for (Vertex_iterator vit = vertices_begin(); vit != vertices_end(); vit++)
  {
    const TVertex& vext = vit->data();
    std::vector<double> knots;

    // left
    this->shoot_left(vit, 2, knots);
    if(greater(vext.s[1],knots[0]))
      throw std::runtime_error("[Tspline::compute_violation_2] Violation_2 cannot occur for top-level basis functions. A");
    else if(greater(vext.s[0],knots[1]))
      throw std::runtime_error("[Tspline::compute_violation_2] Violation_2 cannot occur for top-level basis functions. B");

    // right
    this->shoot_right(vit, 2, knots);
    if(smaller(vext.s[3],knots[0]))
      throw std::runtime_error("[Tspline::compute_violation_2] Violation_2 cannot occur for top-level basis functions. C");
    else if(smaller(vext.s[4], knots[1]))
      throw std::runtime_error("[Tspline::compute_violation_2] Violation_2 cannot occur for top-level basis functions. D");

    // down
    this->shoot_down(vit, 2, knots);
    if(greater(vext.t[1],knots[0]))
      throw std::runtime_error("[Tspline::compute_violation_2] Violation_2 cannot occur for top-level basis functions. E");
    else if(greater(vext.t[0],knots[1]))
      throw std::runtime_error("[Tspline::compute_violation_2] Violation_2 cannot occur for top-level basis functions. F");

    // up
    this->shoot_up(vit, 2, knots);
    if(smaller(vext.t[3],knots[0]))
      throw std::runtime_error("[Tspline::compute_violation_2] Violation_2 cannot occur for top-level basis functions. G");
    else if(smaller(vext.t[4], knots[1]))
      throw std::runtime_error("[Tspline::compute_violation_2] Violation_2 cannot occur for top-level basis functions. H");
  }
  // DEBUG OFF

  bool violation(false);
  for(size_t i=0; i<blend.size(); i++)
  {
    bool v2 = blend[i]->compute_violation_2();
    violation = violation || v2;
  }
  return violation;
}

void Tspline::scale_2_matrix(const std::vector<Scale>& scales, Eigen::MatrixXd& M, int num_vertices_prev)
{
  // Pn = M * P

  // determine size of matrix M
  int cols(num_vertices_prev);
  int rows(number_of_vertices());

  if(cols>rows)
    throw std::runtime_error("[Tspline::scale_2_matrix] Error, more columns than rows.");

  // init matrix
  M = Eigen::MatrixXd::Zero(rows,cols);

  // assemble matrix
  for(size_t n=0; n<scales.size(); n++)
  {
    const Scale& s = scales[n];
    M(s.j,s.i) += s.c;
  }

  for(int i=0; i<cols; i++)
  {
    double& c = M(i,i);
    if(c==0.0)
      c = 1.0;
  }

//  if(!quiet) printf("[Tspline::scale_2_matrix] check if the following code works (enforcing standardT-splines)\n");
//  for(int r=0; r<M.rows(); r++)
//  {
//    double w(0.0);
//    for(int c=0; c<M.cols(); c++)
//      w += M(r,c);

//    if(equal(w,0.0))
//      throw std::runtime_error("[Tspline::scale_2_matrix] Error, zero row detected.");

//    if(!equal(w,1.0))
//      for(int c=0; c<M.cols(); c++)
//        M(r,c) /= w;
//  }
}

void Tspline::update_control_points(const std::vector<Scale>& scales, int num_vertices_prev, Eigen::MatrixXd& M)
{
  if(!quiet) printf("\n[Tspline::update_control_points]\n");

  scale_2_matrix(scales, M, num_vertices_prev);

  vector_vec4d vertices(number_of_vertices());
  for(size_t j=0; j<number_of_vertices(); j++)
  {
    Eigen::Vector4d& cpn = vertices[j];
    cpn = Eigen::Vector4d(0,0,0,0);
    for(int i=0; i<num_vertices_prev; i++)
    {
      Eigen::Vector4d cp;
      get_vertex(i)->data().GetCP(cp);

      double& c = M(j,i);
      cpn += cp * c;
      if(!quiet && c!=0.0) printf("  Pn%lu = P%d * %f\n", j, i, M(j,i));
    }
  }

  for(size_t i=0; i<vertices.size(); i++)
  {
    Eigen::Vector4d& cpn = vertices[i];
    get_vertex(i)->data().SetCP(cpn);
    if(!quiet)
    {
      if(static_cast<int>(i)==num_vertices_prev)
        printf("-----------------------------------------------\n");
      printf("  Pn%lu  (%f %f %f %f)\n", i, cpn(0), cpn(1), cpn(2), cpn(3));
    }
  }

  update_knot_vectors();

  // ################################################################
  // OLD
  //    std::map<size_t, Eigen::Vector4d> new_cps;
  //    std::map<size_t, bool> complete; // i.e. true if new CP contains its old version, or is an inserted CP

  //    if(!quiet) printf("[Tspline::update_control_points] Scales:\n");
  //    for(size_t n=0; n<scales.size(); n++)
  //    {
  //      const Scale& s = scales[n];

  //      Eigen::Vector4d cp;
  //      get_vertex(s.i)->data().GetCP(cp);

  //      cp *= s.c;

  //      // check if sum for CP is complete
  //      if(complete.count(s.j)==0)
  //        complete[s.j] = false;
  //      if(s.i==s.j || s.j >= num_vertices_prev)
  //        complete[s.j] = true;

  //      if(new_cps.count(s.j)==0)
  //      {
  //        new_cps[s.j] = cp;
  //        if(!quiet) printf("  Pn%d = P%d * %f\n", s.j, s.i, s.c);

  //      }else
  //      {
  //        new_cps[s.j] += cp;
  //        if(!quiet) printf("  Pn%d += P%d * %f\n", s.j, s.i, s.c);
  //      }

  //    }

  //    if(!quiet) printf("[Tspline::update_control_points] control-point weights:\n");
  //    std::map<size_t, Eigen::Vector4d >::iterator it;
  //    std::map<size_t, bool>::iterator itb=complete.begin();
  //    for(it=new_cps.begin(); it!=new_cps.end(); it++)
  //    {
  //      Tspline::Vertex_iterator vit = get_vertex(it->first);
  //      TVertex& vext = vit->data();

  //      Eigen::Vector4d& cp = it->second;

  //      if(!itb->second) // if CP does not contain itself, insert itself
  //      {
  //        Eigen::Vector4d cp_old;
  //        vext.GetCP(cp_old);
  //        cp += cp_old;
  //      }

  //      vext.SetCP(cp);
  //      compute_knot_vectors(vit, vext.s, vext.t); // note that the originial knot vectors have not been touched yet

  //      //    printf("  Pn%d  (%f %f %f %f) %d\n", vext.id, cp(0), cp(1), cp(2), cp(3), itb->second);
  //      itb++;
  //    }

  //    for(size_t i=0; i<number_of_vertices() && !quiet; i++)
  //    {
  //      Eigen::Vector4d cp;
  //      get_vertex(i)->data().GetCP(cp);
  //      if(static_cast<int>(i)==num_vertices_prev)
  //        printf("-----------------------------------------------\n");
  //      printf("  Pn%lu  (%f %f %f %f)\n", i, cp(0), cp(1), cp(2), cp(3));
  //    }

}

void Tspline::make_congruent(int num_vertices_prev)
{
  Eigen::MatrixXd M;
  make_congruent(num_vertices_prev, M);
}

void Tspline::make_congruent(int num_vertices_prev, Eigen::MatrixXd& M)
{
  std::vector<BasisFunction*> blend;

  bool v1(true);
  bool v2(true);

  while(v1 || v2)
  {
    if(!quiet) printf("\n[Tspline::make_congruent] compute violation 1\n");
    v1 = compute_violation_1(blend, num_vertices_prev);

    if(!quiet) printf("[Tspline::make_congruent] compute violation 2\n");
    v2 = compute_violation_2(blend);
    if(!quiet) printf("\n[Tspline::make_congruent] v1: %d  v2: %d  blend.size: %lu\n", v1, v2, blend.size());
  }

  if(!quiet) printf("\n[Tspline::make_congruent] Scales:\n");
  std::vector<Scale> scales;
  for(size_t i=0; i<blend.size(); i++)
    blend[i]->get_scale(scales);

  if(!scales.empty())
    update_control_points(scales, num_vertices_prev, M);
}

Arrangement_2::Vertex_iterator Tspline::insert_vertex(Halfedge_iterator e, Point2d param)
{
  Vertex_iterator v0 = e->source();
  Vertex_iterator v1 = e->target();
  const Point2d &p0 = v0->data().param;
  const Point2d &p1 = v1->data().param;

  // if vertex already exists, return
  if (equal(p0, param))
    return v0;
  if (equal(p1, param))
    return v1;

  // find point location in grid
  double d0 = sqrt ( (p0 - param).squared_length());
  double d1 = sqrt ( (p1 - param).squared_length());
  double r0 = d0 / (d0 + d1);
  const Point2d &pt0 = v0->point();
  const Point2d &pt1 = v1->point();
  Vector2d v01 = pt1 - pt0;
  Point2d pn = pt0 + v01 * r0; //r0;

  // insert vertex
  Halfedge_iterator hen = this->split_edge(e, Segment2(pt0, pn), Segment2(pt1, pn));
  Vertex_iterator vn = hen->target(); //CGAL::insert_point (*this, pn);

  if(vn->degree()!=2)
    throw std::runtime_error("[Tspline::insert_vertex] Error, vertex insertion failed.");

  Point2d &ptn = vn->point();
  ptn = Point2d(adjust(ptn.x()), adjust(ptn.y()));

  TVertex& vext = vn->data();
  vext.id = (number_of_vertices() - 1);
  vext.param = param;
  vext.SetCP(v0->data().GetCP() + ( (v1->data().GetCP() - v0->data().GetCP()) * r0)); // todo create vertex correct (Sederberg 1998)
  if(!quiet)
    printf("*[Tspline::insert_vertex(Halfedge)] vertex inserted: %d\n", vext.id);


  // update edge distances d
  Halfedge_around_vertex_circulator first, curr;
  first = curr = vn->incident_halfedges();
  do
  {
    THalfedge hext = curr->data();
    if(equal(curr->source()->data().param, p0) || equal(curr->target()->data().param, p0))
      hext.d = d0;
    if(equal(curr->source()->data().param, p1) || equal(curr->target()->data().param, p1))
      hext.d = d1;
    curr->set_data (hext);
    curr->twin()->set_data (hext);
  }
  while(++curr != first);

  compute_knot_vectors(vn, vext.s, vext.t);

  return vn;
}

Arrangement_2::Vertex_iterator Tspline::insert_vertex(Face_iterator f, Point2d param)
{
  Halfedge_iterator he_right, he_left, he_top, he_bottom;
  get_halfedges (f, param, he_right, he_left, he_top, he_bottom);

  if (he_right == halfedges_end() || he_left == halfedges_end() || he_top == halfedges_end() || he_bottom
      == halfedges_end())
  {
    printf ("[Tspline::insert_vertex] Error, param point does not lie on the face f\n");
    return vertices_end();
  }

  // insert vertices at four edges
  Point2d param_left = Point2d(he_left->source()->data().param.x(), param.y());
  Vertex_iterator vertex_left = insert_vertex(he_left, param_left);
  Point2d param_right = Point2d(he_right->source()->data().param.x(), param.y());
  Vertex_iterator vertex_right = insert_vertex(he_right, param_right);
  Point2d param_top = Point2d(param.x(), he_top->source()->data().param.y());
  Vertex_iterator vertex_top = insert_vertex (he_top, param_top);
  Point2d param_bottom = Point2d(param.x(), he_bottom->source()->data().param.y());
  Vertex_iterator vertex_bottom = insert_vertex(he_bottom, param_bottom);

  // insert edges
  Point2d pn(vertex_bottom->point().x(), vertex_left->point().y());
  pn = adjust(pn);
  THalfedge hext;

  // insert vertex in face
  Vertex_iterator vn = this->insert_in_face_interior(pn, f);
  vn->data().id = (number_of_vertices() - 1);

  Halfedge_iterator hn_left = this->insert_at_vertices(Segment2(vertex_left->point(), pn), vertex_left, vn);
  hext.d = param.x() - param_left.x();
  hn_left->set_data (hext);
  hn_left->twin()->set_data (hext);

  Halfedge_iterator hn_right = this->insert_at_vertices(Segment2(vertex_right->point(), pn), vertex_right, vn);
  hext.d = param_right.x() - param.x();
  hn_right->set_data (hext);
  hn_right->twin()->set_data (hext);

  Halfedge_iterator hn_top = this->insert_at_vertices(Segment2(vertex_top->point(), pn), vertex_top, vn);
  hext.d = param_top.y() - param.y();
  hn_top->set_data (hext);
  hn_top->twin()->set_data (hext);

  Halfedge_iterator hn_bottom = this->insert_at_vertices(Segment2(vertex_bottom->point(), pn), vertex_bottom, vn);
  hext.d = param.y() - param_bottom.y();
  hn_bottom->set_data (hext);
  hn_bottom->twin()->set_data (hext);

  // create vertex data
  TVertex& vext = vn->data();
  vext.param = param;
  const Point3d &g0 = vertex_left->data().GetCP();
  const Point3d &g1 = vertex_right->data().GetCP();
  const Point3d &g2 = vertex_top->data().GetCP();
  const Point3d &g3 = vertex_bottom->data().GetCP();
  vext.SetCP(Point3d ( (g0.x() + g1.x() + g2.x() + g3.x()) * 0.25, (g0.y() + g1.y() + g2.y() + g3.y()) * 0.25,
                       (g0.z() + g1.z() + g2.z() + g3.z()) * 0.25)); // todo create vertex correct (Sederberg 1998)

  compute_knot_vectors(vn, vext.s, vext.t);

  return vn;
}

Arrangement_2::Halfedge_iterator Tspline::split_horizontal(Face_iterator f, double y)
{
  if (f->is_unbounded())
    throw std::runtime_error("[Tspline::split_horizontal] Error, face not valid.");

  double x_left, x_right;
  Halfedge_iterator h_left  = get_left_halfedge_at_param(f, y, x_left);
  Halfedge_iterator h_right = get_right_halfedge_at_param(f, y, x_right);

  if( h_left==halfedges_end() || h_right==halfedges_end() )
    throw std::runtime_error("[Tspline::split_horizontal] Error, halfedges not found.");

  Vertex_iterator v_left = vertices_end();
  Vertex_iterator v_right = vertices_end();

  if( equal(h_left->source()->data().param.y(), y) )
    v_left = h_left->source();
  if( equal(h_left->target()->data().param.y(), y) )
    v_left = h_left->target();

  if( equal(h_right->source()->data().param.y(), y) )
    v_right = h_right->source();
  if( equal(h_right->target()->data().param.y(), y) )
    v_right = h_right->target();

  if( v_left!=vertices_end() && v_right!=vertices_end() )
  {
    if(!quiet) printf("[Tspline::split_horizontal] Warning, opposing vertices not connected.\n");
    Halfedge_iterator he = insert_at_vertices(Segment2(v_left->point(),v_right->point()), v_left, v_right);
    he->data().d = sqrt((v_right->data().param-v_left->data().param).squared_length());
    he->twin()->data().d = he->data().d;
    return he;
  }

  if( v_left==vertices_end() )
    v_left = insert_vertex(h_left, Point2d(x_left, y));

  if( v_right==vertices_end() )
    v_right = insert_vertex(h_right, Point2d(x_right, y));

  Halfedge_iterator hn = insert_at_vertices(Segment2(v_left->point(),v_right->point()), v_left, v_right);
  hn->data().d = x_right - x_left;
  hn->twin()->data().d = x_right - x_left;

  //  printf("[Tspline::split_horizontal] id: %d\n", v_left->data().id);
  //  std::vector<double> s = v_left->data().s;
  //  std::vector<double> t = v_left->data().t;
  //  printf("  s:  %f %f %f %f %f\n", s[0], s[1], s[2], s[3], s[4]);
  //  printf("  t:  %f %f %f %f %f\n", t[0], t[1], t[2], t[3], t[4]);

  //  printf("[Tspline::split_horizontal] id: %d\n", v_right->data().id);
  //  s = v_right->data().s;
  //  t = v_right->data().t;
  //  printf("  s:  %f %f %f %f %f\n", s[0], s[1], s[2], s[3], s[4]);
  //  printf("  t:  %f %f %f %f %f\n", t[0], t[1], t[2], t[3], t[4]);

  return hn;
}

Arrangement_2::Halfedge_iterator Tspline::split_vertical(Face_iterator f, double x)
{
  if (f->is_unbounded())
    throw std::runtime_error("[Tspline::split_vertical] Error, face not valid.");

  double y_bottom, y_top;
  Halfedge_iterator h_bottom = get_bottom_halfedge_at_param(f, x, y_bottom);
  Halfedge_iterator h_top = get_top_halfedge_at_param(f, x, y_top);

  if( h_bottom==halfedges_end() || h_top==halfedges_end() )
    throw std::runtime_error("[Tspline::split_vertical] Error, halfedges not found.");

  Vertex_iterator v_bottom = vertices_end();
  Vertex_iterator v_top = vertices_end();

  if( equal(h_bottom->source()->data().param.x(), x) )
    v_bottom = h_bottom->source();
  if( equal(h_bottom->target()->data().param.x(), x) )
    v_bottom = h_bottom->target();

  if( equal(h_top->source()->data().param.x(), x) )
    v_top = h_top->source();
  if( equal(h_top->target()->data().param.x(), x) )
    v_top = h_top->target();

  if( v_bottom!=vertices_end() && v_top!=vertices_end() )
  {
    if(!quiet) printf("[Tspline::split_vertical] Warning, opposing vertices not connected.\n");
    Halfedge_iterator he = insert_at_vertices(Segment2(v_bottom->point(),v_top->point()), v_bottom, v_top);
    he->data().d = sqrt((v_top->data().param-v_bottom->data().param).squared_length());
    he->twin()->data().d = he->data().d;
    return he;
  }

  if( v_bottom==vertices_end() )
    v_bottom = insert_vertex(h_bottom, Point2d(x, y_bottom));

  if( v_top==vertices_end() )
    v_top = insert_vertex(h_top, Point2d(x, y_top));

  Halfedge_iterator hn = insert_at_vertices(Segment2(v_bottom->point(),v_top->point()), v_bottom, v_top);
  hn->data().d = y_top - y_bottom;
  hn->twin()->data().d = y_top - y_bottom;

  //  printf("[Tspline::split_vertical] id: %d\n", v_bottom->data().id);
  //  std::vector<double> s = v_bottom->data().s;
  //  std::vector<double> t = v_bottom->data().t;
  //  printf("  s:  %f %f %f %f %f\n", s[0], s[1], s[2], s[3], s[4]);
  //  printf("  t:  %f %f %f %f %f\n", t[0], t[1], t[2], t[3], t[4]);

  //  printf("[Tspline::split_vertical] id: %d\n", v_top->data().id);
  //  s = v_top->data().s;
  //  t = v_top->data().t;
  //  printf("  s:  %f %f %f %f %f\n", s[0], s[1], s[2], s[3], s[4]);
  //  printf("  t:  %f %f %f %f %f\n", t[0], t[1], t[2], t[3], t[4]);

  return hn;
}

