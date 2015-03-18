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

#ifndef _TSPLINE_UTILS_HPP_
#define _TSPLINE_UTILS_HPP_

#include "Tspline/Tspline.h"
#include "Tspline/TsplineMultiPatch.h"
#include "v4r/TomGine/tgTomGineThread.h"
#include "v4r/TomGine/tgShapeCreator.h"

#include <opencv2/highgui/highgui.hpp>

static float DIV_255 = 1.0f / 255.0f;
static float DIV_3 = 1.0f / 3.0f;

static TomGine::vec3 getColor(const cv::Mat &image, const int &x, const int &y)
{
  cv::Vec3b c = image.at<cv::Vec3b> (y, x);
  return TomGine::vec3(float(c[0]), float(c[1]), float(c[2]));
}

static TomGine::vec3 getColor(const cv::Mat &image, const float &x, const float &y)
{
  cv::Vec3b c = image.at<cv::Vec3b> ((int) round(y), (int) round(x));
  return TomGine::vec3(float(c[0]), float(c[1]), float(c[2]));
}

static TomGine::vec3 getColorBilinear(const cv::Mat &image, const float &x, const float &y, bool quiet=true)
{
  int x1 = floor(x);
  int y1 = floor(y);

  float dx = x - x1;
  float dy = y - y1;

  TomGine::vec3 a = getColor(image, x1+0, y1);
  TomGine::vec3 b = getColor(image, x1+1, y1);
  TomGine::vec3 c = getColor(image, x1+0, y1+1);
  TomGine::vec3 d = getColor(image, x1+1, y1+1);


  TomGine::vec3 e = a * (1.0 - dx) + b * dx;
  TomGine::vec3 f = c * (1.0 - dx) + d * dx;

  if(!quiet)
  {
    printf("getColorBilinear (%d %d)\n", x1, y1);
    printf("a: %f %f %f  b: %f %f %f dx: %f\n", a.x, a.y, a.z, b.x, b.y, b.z, dx);
    printf("c: %f %f %f  d: %f %f %f dx: %f\n", c.x, c.y, c.z, d.x, d.y, d.z, dx);
    printf("e: %f %f %f  f: %f %f %f dy: %f\n", e.x, e.y, e.z, f.x, f.y, f.z, dy);
  }

  return e * (1.0 - dy) + f * dy;
}


static TomGine::vec3 getGradient(const cv::Mat &gradient,
                                 const int &x,
                                 const int &y)
{
  cv::Vec3b gb = gradient.at<cv::Vec3b> (y, x);
  TomGine::vec3 g;
  g.x = ((float(gb[0]) + 0.5f) * DIV_255 - 0.5f); // fixme: rounding is not correct (+0.5f)
  g.y = ((float(gb[1]) + 0.5f) * DIV_255 - 0.5f);
  g.z = ((float(gb[2]) + 0.5f) * DIV_255 - 0.5f);
  return g;
}

static TomGine::vec3 getGradient(const cv::Mat &image,
                                 const float &x,
                                 const float &y)
{
  cv::Vec3b gb = image.at<cv::Vec3b> ((int) round(y), (int) round(x));
  TomGine::vec3 g;
  g.x = ((float(gb[0]) + 0.5f) * DIV_255 - 0.5f);
  g.y = ((float(gb[1]) + 0.5f) * DIV_255 - 0.5f);
  g.z = ((float(gb[2]) + 0.5f) * DIV_255 - 0.5f);
  return g;
}

static TomGine::vec3 getGradientBilinear(const cv::Mat &image, const float &x, const float &y)
{
  int x1 = floor(x);
  int y1 = floor(y);

  float dx = x - x1;
  float dy = y - y1;

  TomGine::vec3 a = getGradient(image, x1+0, y1);
  TomGine::vec3 b = getGradient(image, x1+1, y1);
  TomGine::vec3 c = getGradient(image, x1+0, y1+1);
  TomGine::vec3 d = getGradient(image, x1+1, y1+1);

  TomGine::vec3 e = a * (1.0 - dx) + b * dx;
  TomGine::vec3 f = c * (1.0 - dx) + d * dx;

  return e * (1.0 - dy) + f * dy;
}

static void convertTspline2tgModel(const tspline::Tspline &tsp,
                                   TomGine::tgModel &model,
                                   unsigned resU, unsigned resV,
                                   bool per_face=false)
{
  model.Clear();
  double param_w = tsp.param_max.x() - tsp.param_min.x();
  double param_h = tsp.param_max.y() - tsp.param_min.y();
  double div_w = 1.0 / param_w;
  double div_h = 1.0 / param_h;

  if(per_face)
  {
    tspline::Tspline::Face_const_iterator fit;
    for(fit=tsp.faces_begin(); fit!=tsp.faces_end(); fit++)
    {
      if(fit->is_unbounded() || !fit->has_outer_ccb())
        continue;

      Eigen::Vector4d bb;
      tsp.face_bounding_box(fit, bb);

      TomGine::tgShapeCreator::CreatePlaneXY(model,
                                             bb(0), bb(2), 0.0f,
                                             bb(1)-bb(0), bb(3)-bb(2),
                                             resU, resV);
    }
  }
  else
  {
    TomGine::tgShapeCreator::CreatePlaneXY(model, tsp.param_min.x(), tsp.param_min.y(), 0.0f, param_w, param_h,
                                           resU, resV);
  }

#pragma omp parallel for
  for (size_t i = 0; i < model.m_vertices.size(); i++)
  {
    TomGine::tgVertex &v = model.m_vertices[i];

    v.texCoord.x = (v.pos.x - tsp.param_min.x()) * div_w;
    v.texCoord.y = (v.pos.y - tsp.param_min.y()) * div_h;

    Eigen::Vector3d ts, tt;
    Eigen::Vector3d s = tsp.evaluate(v.pos.x, v.pos.y, ts, tt);
    v.pos = TomGine::vec3(s(0), s(1), s(2));

    Eigen::Vector3d n = ts.cross(tt);
    n.normalize();
    v.normal = TomGine::vec3(n(0), n(1), n(2));
  }
}

static void convertTsplineFace2tgModel(tspline::Tspline::Face_iterator &fit, tspline::Tspline &tsp,
                                       TomGine::tgModel &model, unsigned resU, unsigned resV)
{
  if(fit->is_unbounded())
    return;

  size_t idx = model.m_vertices.size();

  tspline::Tspline::Ccb_halfedge_circulator first, curr;
  first = curr = fit->outer_ccb ();

  tspline::Point2d p_max(-DBL_MAX, -DBL_MAX);
  tspline::Point2d p_min(DBL_MAX, DBL_MAX);

  do
  {
    tspline::Point2d p = curr->source()->data().param;

    if(p.x()<p_min.x())
      p_min = tspline::Point2d(p.x(),p_min.y());
    if(p.x()>p_max.x())
      p_max = tspline::Point2d(p.x(),p_max.y());
    if(p.y()<p_min.y())
      p_min = tspline::Point2d(p_min.x(),p.y());
    if(p.y()>p_max.y())
      p_max = tspline::Point2d(p_max.x(),p.y());
  }
  while (++curr != first);

  tspline::Point2d p_size( p_max.x()-p_min.x(), p_max.y()-p_min.y() );

  //  printf("p_min: %f %f   p_max: %f %f\n", p_min.x(), p_min.y(), p_max.x(), p_max.y());

  TomGine::tgShapeCreator::CreatePlaneXY(model, p_min.x(),p_min.y(), 0.0f,p_size.x(),p_size.y(), resU,resV);

  //  double div_w = 1.0 / p_size.x();
  //  double div_h = 1.0 / p_size.y();

#pragma omp parallel for
  for (size_t i = idx; i < model.m_vertices.size(); i++)
  {
    TomGine::tgVertex &v = model.m_vertices[i];

    //v.texCoord.x = (v.pos.x - p_min.x()) * div_w;
    //v.texCoord.y(v.pos.y - p_min.y()) * div_h;
    v.texCoord.x = v.pos.x;
    v.texCoord.y = v.pos.y;

    Eigen::Vector3d ts, tt;
    Eigen::Vector3d s = tsp.evaluate(v.pos.x, v.pos.y, ts, tt);
    v.pos = TomGine::vec3(s(0), s(1), s(2));

    Eigen::Vector3d n = ts.cross(tt);
    n.normalize();
    v.normal = TomGine::vec3(n(0), n(1), n(2));
  }

}

static void convertTsplineControl2tgModel(const tspline::Tspline &tsp, TomGine::tgModel &model)
{
  // control points
  tspline::Tspline::Vertex_const_iterator vit;
  for (vit = tsp.vertices_begin(); vit != tsp.vertices_end(); ++vit)
  {
    const tspline::Point3d &g = vit->data().GetCP();
    model.m_points.push_back(TomGine::vec3(g.x(), g.y(), g.z()));
  }
  model.m_point_size = 10.0f;
  model.m_point_color = TomGine::vec3(0.0f, 1.0f, 0.0f);

  // edges
  tspline::Tspline::Edge_const_iterator eit;
  for (eit = tsp.edges_begin(); eit != tsp.edges_end(); ++eit)
  {
    const tspline::Point3d &geom1 = eit->source()->data().GetCP();
    const tspline::Point3d &geom2 = eit->target()->data().GetCP();

    Eigen::Vector3d p1(geom1.x(), geom1.y(), geom1.z());
    Eigen::Vector3d p2(geom2.x(), geom2.y(), geom2.z());

    TomGine::tgLine line;
    line.start = TomGine::vec3(p1(0), p1(1), p1(2));
    line.end = TomGine::vec3(p2(0), p2(1), p2(2));
    model.m_lines.push_back(line);
  }
  model.m_line_color = TomGine::vec3(0.0, 0.8, 0.0);
  model.m_line_width = 2.0f;
}

static void convertTsplineEdges2tgModel(const tspline::Tspline &tsp, TomGine::tgModel &model,
                                        unsigned subres=8, double offset=0.0)
{
  if(subres<=0)
    throw std::runtime_error("[Utils::convertTsplineEdges2tgModel] resolution too low\n");

  Eigen::Vector3d p1, p2, ts1, tt1, ts2, tt2, n1, n2;


  // edges
  tspline::Tspline::Edge_const_iterator eit;
  for (eit = tsp.edges_begin(); eit != tsp.edges_end(); ++eit)
  {
    const tspline::Point2d &param1 = eit->source()->data().param;
    const tspline::Point2d &param2 = eit->target()->data().param;

    Eigen::Vector2d a(param1.x(),param1.y());
    Eigen::Vector2d b(param2.x(),param2.y());
    Eigen::Vector2d d = b-a;

    for(unsigned i=0; i<subres; i++)
    {
      Eigen::Vector2d sp = a + d * double(i) / subres;
      Eigen::Vector2d ep = a + d * double(i+1) / subres;


      p1 = tsp.evaluate(sp[0],sp[1], ts1, tt1);
      p2 = tsp.evaluate(ep[0],ep[1], ts2, tt2);

      n1 = ts1.cross(tt1);
      n2 = ts2.cross(tt2);
      n1 /= n1.norm();
      n2 /= n2.norm();

      TomGine::tgLine line;
      line.start = TomGine::vec3(p1[0], p1[1], p1[2]) + TomGine::vec3(n1[0],n1[1],n1[2]) * offset;
      line.end = TomGine::vec3(p2[0], p2[1], p2[2]) + TomGine::vec3(n2[0],n2[1],n2[2]) * offset;
      model.m_lines.push_back(line);
    }
  }
  model.m_line_color = TomGine::vec3(0.0, 0.8, 0.0);
  model.m_line_width = 2.0f;
}

static void convertTspline2tgModel(const tspline::TsplineMultiPatch &tsmp,
                                   TomGine::tgModel &model,
                                   unsigned resU, unsigned resV,
                                   bool per_face=false)
{
  model.Clear();
  size_t vidx(0);
  std::list<tspline::TsplinePatch*>::const_iterator it;
  for(it=tsmp.patchlist.begin(); it!=tsmp.patchlist.end(); it++)
  {
    const tspline::TsplinePatch& tsp = *(*it);
    double param_w = tsp.param_max.x() - tsp.param_min.x();
    double param_h = tsp.param_max.y() - tsp.param_min.y();
    double div_w = 1.0 / param_w;
    double div_h = 1.0 / param_h;
    vidx = model.m_vertices.size();

    if(per_face)
    {
      tspline::TsplinePatch::Face_const_iterator fit;
      for(fit=tsp.faces_begin(); fit!=tsp.faces_end(); fit++)
      {
        if(fit->is_unbounded() || !fit->has_outer_ccb())
          continue;

        Eigen::Vector4d bb;
        tsp.face_bounding_box(fit, bb);

        TomGine::tgShapeCreator::CreatePlaneXY(model,
                                               bb(0), bb(2), 0.0f,
                                               bb(1)-bb(0), bb(3)-bb(2),
                                               resU, resV);
      }
    }
    else
    {
      TomGine::tgShapeCreator::CreatePlaneXY(model,
                                             tsp.param_min.x(), tsp.param_min.y(), 0.0f,
                                             param_w, param_h,
                                             resU, resV);
    }

#pragma omp parallel for
    for (size_t i = vidx; i < model.m_vertices.size(); i++)
    {
      TomGine::tgVertex &v = model.m_vertices[i];

      v.texCoord.x = (v.pos.x - tsp.param_min.x()) * div_w;
      v.texCoord.y = (v.pos.y - tsp.param_min.y()) * div_h;

      v.color[0] = v.texCoord.x * 255;
      v.color[1] = v.texCoord.y * 255;
      v.color[2] = 0.0;

      Eigen::Vector3d ts, tt;
      Eigen::Vector3d p = tsp.evaluate(v.pos.x, v.pos.y, ts, tt);
      v.pos = TomGine::vec3(p(0), p(1), p(2));
      Eigen::Vector3d n = ts.cross(tt);
      n.normalize();
      v.normal = TomGine::vec3(n(0), n(1), n(2));

      //      Eigen::Vector3d p = tsp.evaluate(v.pos.x, v.pos.y);
      //      v.pos = TomGine::vec3(p(0), p(1), p(2));
    }

  }
}

static void convertTspline2tgTextureModel(const tspline::TsplineMultiPatch &tsmp,
                                          std::vector<cv::Mat>& textures,
                                          TomGine::tgTextureModel &model,
                                          unsigned resU, unsigned resV,
                                          bool per_face=false)
{
  if(textures.size() != tsmp.patchlist.size())
    throw std::runtime_error("[Utils::convertTspline2tgTextureModel] Error, size of patches does not match size of textures.");

  model.Clear();
  model.m_tex_cv.resize(tsmp.patchlist.size());

  size_t vidx(0), fidx(0), i(0);
  std::list<tspline::TsplinePatch*>::const_iterator it;
  for(it=tsmp.patchlist.begin(); it!=tsmp.patchlist.end(); it++)
  {
    const tspline::TsplinePatch& tsp = *(*it);
    double param_w = tsp.param_max.x() - tsp.param_min.x();
    double param_h = tsp.param_max.y() - tsp.param_min.y();
    double div_w = 1.0 / param_w;
    double div_h = 1.0 / param_h;
    vidx = model.m_vertices.size();
    fidx = model.m_faces.size();

    if(per_face)
    {
      tspline::TsplinePatch::Face_const_iterator fit;
      for(fit=tsp.faces_begin(); fit!=tsp.faces_end(); fit++)
      {
        if(fit->is_unbounded() || !fit->has_outer_ccb())
          continue;

        Eigen::Vector4d bb;
        tsp.face_bounding_box(fit, bb);

        TomGine::tgShapeCreator::CreatePlaneXY(model,
                                               bb(0), bb(2), 0.0f,
                                               bb(1)-bb(0), bb(3)-bb(2),
                                               resU, resV);
      }
    }
    else
    {
      TomGine::tgShapeCreator::CreatePlaneXY(model,
                                             tsp.param_min.x(), tsp.param_min.y(), 0.0f,
                                             param_w, param_h,
                                             resU, resV);
    }

    model.m_face_tex_id.insert(model.m_face_tex_id.end(),model.m_faces.size()-fidx, i);

    if(!textures[i].empty())
      textures[i].copyTo(model.m_tex_cv[i]);

#pragma omp parallel for
    for (size_t j = vidx; j < model.m_vertices.size(); j++)
    {
      TomGine::tgVertex &v = model.m_vertices[j];

      v.texCoord.x = (v.pos.x - tsp.param_min.x()) * div_w;
      v.texCoord.y = (v.pos.y - tsp.param_min.y()) * div_h;

      v.color[0] = v.texCoord.x * 255;
      v.color[1] = v.texCoord.y * 255;
      v.color[2] = 0.0;

      Eigen::Vector3d ts, tt;
      Eigen::Vector3d p = tsp.evaluate(v.pos.x, v.pos.y, ts, tt);
      v.pos = TomGine::vec3(p(0), p(1), p(2));
      Eigen::Vector3d n = ts.cross(tt);
      n.normalize();
      v.normal = TomGine::vec3(n(0), n(1), n(2));

      //      Eigen::Vector3d p = tsp.evaluate(v.pos.x, v.pos.y);
      //      v.pos = TomGine::vec3(p(0), p(1), p(2));
    }

    i++;
  }

  model.SetColor(1.0f,1.0f,1.0f);
}

static void convertTsplineControl2tgModel(const tspline::TsplineMultiPatch &tsmp, TomGine::tgModel &model)
{
  std::list<tspline::TsplinePatch*>::const_iterator it;
  for(it=tsmp.patchlist.begin(); it!=tsmp.patchlist.end(); it++)
  {
    const tspline::TsplinePatch& tsp = *(*it);
    //    // control points
    //    tspline::Tspline::Vertex_const_iterator vit;
    //    for (vit = tsp.vertices_begin(); vit != tsp.vertices_end(); ++vit)
    //    {
    //      const tspline::Point3d &g = vit->data().GetCP();
    //      model.m_points.push_back(TomGine::vec3(g.x(), g.y(), g.z()));
    //    }
    //    model.m_point_size = 10.0f;
    //    model.m_point_color = TomGine::vec3(0.0f, 0.0f, 1.0f);

    // edges
    tspline::Tspline::Edge_const_iterator eit;
    for (eit = tsp.edges_begin(); eit != tsp.edges_end(); ++eit)
    {
      const tspline::Point3d &geom1 = eit->source()->data().GetCP();
      const tspline::Point3d &geom2 = eit->target()->data().GetCP();

      Eigen::Vector3d p1(geom1.x(), geom1.y(), geom1.z());
      Eigen::Vector3d p2(geom2.x(), geom2.y(), geom2.z());

      TomGine::tgLine line;
      line.start = TomGine::vec3(p1(0), p1(1), p1(2));
      line.end = TomGine::vec3(p2(0), p2(1), p2(2));
      model.m_lines.push_back(line);
    }
    model.m_line_color = TomGine::vec3(0.0, 0.8, 0.0);
    model.m_line_width = 2.0f;
  }
}

static void convertTsplineEdges2tgModel(const tspline::TsplineMultiPatch &tsmp, TomGine::tgModel &model,
                                        unsigned subres=8, TomGine::vec3 offset = TomGine::vec3(0,0,0))
{
  if(subres<=0)
    throw std::runtime_error("[Utils::convertTsplineEdges2tgModel] resolution too low\n");

  std::list<tspline::TsplinePatch*>::const_iterator it;
  for(it=tsmp.patchlist.begin(); it!=tsmp.patchlist.end(); it++)
  {
    const tspline::TsplinePatch& tsp = *(*it);

    // edges
    tspline::Tspline::Edge_const_iterator eit;
    for (eit = tsp.edges_begin(); eit != tsp.edges_end(); ++eit)
    {
      const tspline::Point2d &p1 = eit->source()->data().param;
      const tspline::Point2d &p2 = eit->target()->data().param;

      Eigen::Vector2d a(p1.x(),p1.y());
      Eigen::Vector2d b(p2.x(),p2.y());
      Eigen::Vector2d d = b-a;

      for(unsigned i=0; i<subres; i++)
      {
        Eigen::Vector2d sp = a + d * double(i) / subres;
        Eigen::Vector2d ep = a + d * double(i+1) / subres;

        Eigen::Vector3d s = tsp.evaluate(sp[0],sp[1]);
        Eigen::Vector3d e = tsp.evaluate(ep[0],ep[1]);

        TomGine::tgLine line;
        line.start = TomGine::vec3(s(0), s(1), s(2)) + offset;
        line.end = TomGine::vec3(e(0), e(1), e(2)) + offset;
        model.m_lines.push_back(line);
      }
    }
  }
  model.m_line_color = TomGine::vec3(0.0, 0.8, 0.0);
  model.m_line_width = 2.0f;
}

//static void cut_file_name(std::string full_file_name, std::string& file_name, std::string& path)
//{
//  size_t c_slash_idx = full_file_name.find_last_of("/\\");
//  size_t c_dot_idx = full_file_name.find_last_of(".");
//  if(c_slash_idx == std::string::npos)
//    c_slash_idx = 0;
//  else
//    c_slash_idx++;
//  if(c_dot_idx == std::string::npos || c_dot_idx < c_slash_idx)
//    c_dot_idx = full_file_name.size();
//  file_name = full_file_name.substr(c_slash_idx,c_dot_idx-c_slash_idx);
//  path = full_file_name.substr(0, c_slash_idx);
//  if(c_slash_idx==0)
//    path= ".";
//}

#endif

