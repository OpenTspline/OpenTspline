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

bool Transition::InSourceRange(const double &p) const
{
  if(gequal(p, source_range[0]))
  {
    if(smaller(p, source_range[1]))
      return true;
    else
    {
      if(source_dir==NORTH || source_dir==SOUTH)
        if(equal(source_range[1], source->param_max.x()))
          return true;
      if(source_dir==EAST || source_dir==WEST)
        if(equal(source_range[1], source->param_max.y()))
          return true;
    }
  }

  return false;
}

bool Transition::InTargetRange(const double &p) const
{
  if(gequal(p, target_range[0]))
  {
    if(smaller(p, target_range[1]))
      return true;
    else
    {
      if(target_dir==NORTH || target_dir==SOUTH)
        if(equal(target_range[1], target->param_max.x()))
          return true;
      if(target_dir==EAST || target_dir==WEST)
        if(equal(target_range[1], target->param_max.y()))
          return true;
    }
  }

  return false;
}

void Transition::shoot(const Point2d& p0, unsigned n, std::vector<double> &knots)
{
  Point2d p1 = ConvertParam(p0);
  CGAL::Object obj = target->locate_param(p1.x(), p1.y());
  Tspline::Vertex_iterator v;
  Tspline::Halfedge_iterator he;
  bool edge = false;
  if (obj.empty())
    throw std::runtime_error("[Transition::shoot] Error, no object found");
  if (!CGAL::assign(v, obj))
  {
    if(!CGAL::assign(he, obj))
      throw std::runtime_error("[NorthTransition::shoot] Error, object not a vertex or edge");
    else
      edge = true;
  }

  if(target_dir==NORTH)
  {
    std::vector<double> tmp;
    if(edge)
      target->shoot_down(he, p1.x(), n, tmp);
    else
      target->shoot_down(v, n, tmp);

    for(size_t i=0; i<n; i++)
    {
      Point2d p0n = twin->ConvertParam(Point2d(p1.x(), tmp[i]));
      if(source_dir==NORTH || source_dir==SOUTH)
        knots.push_back(p0n.y());
      if(source_dir==EAST || source_dir==WEST)
        knots.push_back(p0n.x());
    }

  }
  if(target_dir==EAST)
  {
    std::vector<double> tmp;
    if(edge)
      target->shoot_left(he, p1.y(), n, tmp);
    else
      target->shoot_left(v, n, tmp);
    for(size_t i=0; i<n; i++)
    {
      Point2d p0n = twin->ConvertParam(Point2d(tmp[i],p1.y()));
      if(source_dir==NORTH || source_dir==SOUTH)
        knots.push_back(p0n.y());
      if(source_dir==EAST || source_dir==WEST)
        knots.push_back(p0n.x());
    }
  }
  if(target_dir==SOUTH)
  {
    std::vector<double> tmp;
    if(edge)
      target->shoot_up(he, p1.x(), n, tmp);
    else
      target->shoot_up(v, n, tmp);

    for(size_t i=0; i<n; i++)
    {
      Point2d p0n = twin->ConvertParam(Point2d(p1.x(),tmp[i]));
      if(source_dir==NORTH || source_dir==SOUTH)
        knots.push_back(p0n.y());
      if(source_dir==EAST || source_dir==WEST)
        knots.push_back(p0n.x());
    }
  }
  if(target_dir==WEST)
  {
    std::vector<double> tmp;
    if(edge)
      target->shoot_right(he, p1.y(), n, tmp);
    else
      target->shoot_right(v, n, tmp);
    for(size_t i=0; i<n; i++)
    {
      Point2d p0n = twin->ConvertParam(Point2d(tmp[i],p1.y()));
      if(source_dir==NORTH || source_dir==SOUTH)
        knots.push_back(p0n.y());
      if(source_dir==EAST || source_dir==WEST)
        knots.push_back(p0n.x());
    }
  }

}

Point2d NorthTransition::ConvertParam(const Point2d& p) const
{
  const double& a0 = source_range[0];
  const double& b0 = source_range[1];
  const double& a1 = target_range[0];

  if(target_dir==NORTH)
    return Point2d( a1+(b0-p.x()), target->param_max.y()+(source->param_max.y()-p.y()) );
  if(target_dir==EAST)
    return Point2d( target->param_max.x()+(source->param_max.y()-p.y()), a1+(p.x()-a0) );
  if(target_dir==SOUTH)
    return Point2d( a1+(p.x()-a0), target->param_min.y()-(source->param_max.y()-p.y()) );
  if(target_dir==WEST)
    return Point2d( target->param_min.x()-(source->param_max.y()-p.y()), a1+(b0-p.x()) );

  throw std::runtime_error("[NorthTransition::ConvertParam] Error, unknown transition");
}

Point2d NorthTransition::RotateParamInverse(const Point2d& p) const
{
  double s = (target_range[1] - target_range[0]) / (source_range[1] - source_range[0]);

  if(s!=1.0)
    throw std::runtime_error("[NorthTransition::RotateParamInverse] Error, parameter scaling");

  if(target_dir==NORTH)
    return Point2d( -p.x()*s, -p.y() );
  if(target_dir==EAST)
    return Point2d( p.y()*s, -p.x() );
  if(target_dir==SOUTH)
    return Point2d( p.x()*s, p.y() );
  if(target_dir==WEST)
    return Point2d( -p.y()*s, p.x() );

  throw std::runtime_error("[NorthTransition::RotateParamInverse] Error, unknown transition");
}

Point2d EastTransition::ConvertParam(const Point2d& p) const
{
  const double& a0 = source_range[0];
  const double& b0 = source_range[1];
  const double& a1 = target_range[0];

  if(target_dir==NORTH)
    return Point2d( a1+(p.y()-a0), target->param_max.y()+(source->param_max.x()-p.x()) );
  if(target_dir==EAST)
    return Point2d( target->param_max.x()+(source->param_max.x()-p.x()), a1+(b0-p.y()) );
  if(target_dir==SOUTH)
    return Point2d( a1+(b0-p.y()), target->param_min.y()-(source->param_max.x()-p.x()) );
  if(target_dir==WEST)
    return Point2d( target->param_min.x()-(source->param_max.x()-p.x()), a1+(p.y()-a0) );

  throw std::runtime_error("[EastTransition::ConvertParam] Error, unknown transition");
}

Point2d EastTransition::RotateParamInverse(const Point2d& p) const
{
  double s = (target_range[1] - target_range[0]) / (source_range[1] - source_range[0]);

  if(s!=1.0)
    throw std::runtime_error("[NorthTransition::RotateParamInverse] Error, parameter scaling");

  if(target_dir==NORTH)
    return Point2d( -p.y(), p.x()*s );
  if(target_dir==EAST)
    return Point2d( -p.x(), -p.y()*s );
  if(target_dir==SOUTH)
    return Point2d( p.y(), -p.x()*s );
  if(target_dir==WEST)
    return Point2d( p.x(), p.y()*s );

  throw std::runtime_error("[EastTransition::RotateParamInverse] Error, unknown transition");
}

Point2d SouthTransition::ConvertParam(const Point2d& p) const
{
  const double& a0 = source_range[0];
  const double& b0 = source_range[1];
  const double& a1 = target_range[0];

  if(target_dir==NORTH)
    return Point2d( a1+(p.x()-a0), target->param_max.y()+(p.y()-source->param_min.y()) );
  if(target_dir==EAST)
    return Point2d( target->param_max.x()+(p.y()-source->param_min.y()), a1+(b0-p.x()) );
  if(target_dir==SOUTH)
    return Point2d( a1+(b0-p.x()), target->param_min.y()-(p.y()-source->param_min.y()) );
  if(target_dir==WEST)
    return Point2d( target->param_min.x()-(p.y()-source->param_min.y()), a1+(p.x()-a0) );

  throw std::runtime_error("[SouthTransition::ConvertParam] Error, unknown transition");
}

Point2d SouthTransition::RotateParamInverse(const Point2d& p) const
{
  double s = (target_range[1] - target_range[0]) / (source_range[1] - source_range[0]);

  if(s!=1.0)
    throw std::runtime_error("[NorthTransition::RotateParamInverse] Error, parameter scaling");

  if(target_dir==NORTH)
    return Point2d( p.x()*s, p.y() );
  if(target_dir==EAST)
    return Point2d( -p.y()*s, p.x() );
  if(target_dir==SOUTH)
    return Point2d( -p.x()*s, -p.y() );
  if(target_dir==WEST)
    return Point2d( p.y()*s, -p.x() );

  throw std::runtime_error("[SouthTransition::RotateParamInverse] Error, unknown transition");
}

Point2d WestTransition::ConvertParam(const Point2d& p) const
{
  const double& a0 = source_range[0];
  const double& b0 = source_range[1];
  const double& a1 = target_range[0];

  if(target_dir==NORTH)
    return Point2d( a1+(b0-p.y()), target->param_max.y()+(p.x()-source->param_min.x()) );
  if(target_dir==EAST)
    return Point2d( target->param_max.x()+(p.x()-source->param_min.x()), a1+(p.y()-a0) );
  if(target_dir==SOUTH)
    return Point2d( a1+(p.y()-a0), target->param_min.y()-(p.x()-source->param_min.x()) );
  if(target_dir==WEST)
    return Point2d( target->param_min.x()-(p.x()-source->param_min.x()), a1+(b0-p.y()) );

  throw std::runtime_error("[WestTransition::ConvertParam] Error, unknown transition");
}

Point2d WestTransition::RotateParamInverse(const Point2d& p) const
{
  double s = (target_range[1] - target_range[0]) / (source_range[1] - source_range[0]);

  if(s!=1.0)
    throw std::runtime_error("[NorthTransition::RotateParamInverse] Error, parameter scaling");

  if(target_dir==NORTH)
    return Point2d( p.y(), -p.x()*s );
  if(target_dir==EAST)
    return Point2d( p.x(), p.y()*s );
  if(target_dir==SOUTH)
    return Point2d( -p.y(), p.x()*s );
  if(target_dir==WEST)
    return Point2d( -p.x(), -p.y()*s );

  throw std::runtime_error("[WestTransition::RotateParamInverse] Error, unknown transition");
}
