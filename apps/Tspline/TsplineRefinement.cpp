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

#include "Tspline/Utils.hpp"
#include "Tspline/TsplineCreator.h"
#include "Tspline/File.h"

TomGine::tgTomGineThread viewer(800, 600, "T-spline Refinement");

using namespace tspline;

void render(Tspline &tsp, unsigned res=16);

int main()
{
  viewer.SetClearColor(1.0f);

  double x0(-1.0);
  double y0(-1.0);
  double width(2.0);
  double height(2.0);
  double segX(2);
  double segY(2);

  // create T-spline
  Tspline tsp;
  TsplineCreator::CreatePlaneXY(tsp, x0,y0,0.0, width,height, segX,segY);
  //  TsplineCreator::CreateUniformPlaneXY(tsp, -1.0,0.0,-1.0, 2.0,2.0, 4,4);
  //  tspline::File::Load(tsp, "result.tsp");

  TVertex& vext = tsp.get_vertex(4)->data();
  Point3d p = vext.GetCP();
  vext.SetCP(Point3d(p.x(),p.y(),1.0));

  TomGine::tgRenderModel mesh;
  convertTspline2tgModel(tsp, mesh, 16, 16);
  mesh.m_line_width = 2.0f;
  viewer.AddModel3D(mesh);

  TomGine::mat4 ext;
  TomGine::tgCamera cam = viewer.GetCamera();
  cam.SetExtrinsic(ext);
  cam.TranslateF(-3.0);
  cam.ApplyTransform();
  viewer.SetCamera(cam);

  render(tsp, 16);
  viewer.WaitForEvent(TomGine::TMGL_Press, TomGine::TMGL_Space);

  Tspline::Face_iterator f = tsp.locate_face(0.5, 0.5);
  tsp.split_horizontal_congruent(f, 0.5);
  render(tsp, 16);
  viewer.WaitForEvent(TomGine::TMGL_Press, TomGine::TMGL_Space);

  f = tsp.locate_face(-0.5, -0.5);
  tsp.split_vertical_congruent(f, -0.5);
  render(tsp, 16);
  viewer.WaitForEvent(TomGine::TMGL_Press, TomGine::TMGL_Space);

  f = tsp.locate_face(-0.25, 0.25);
  tsp.split_horizontal_congruent(f, 0.25);
  render(tsp, 16);
  viewer.WaitForEvent(TomGine::TMGL_Press, TomGine::TMGL_Space);

  f = tsp.locate_face(0.5, -0.5);
  tsp.split_vertical_congruent(f, 0.5);
  render(tsp, 16);
  viewer.WaitForEvent(TomGine::TMGL_Press, TomGine::TMGL_Space);


  viewer.Update();
  viewer.WaitForEvent(TomGine::TMGL_Press, TomGine::TMGL_Escape);
  return (0);
}

void render(Tspline &tsp, unsigned res)
{
  static int mesh_id = -1;

  TomGine::tgRenderModel mesh;
  convertTspline2tgModel(tsp, mesh, res, res);
  convertTsplineControl2tgModel(tsp, mesh);
  mesh.m_line_width = 2.0f;

  if(mesh_id==-1)
    mesh_id = viewer.AddModel3D(mesh);
  else
    viewer.SetModel3D(mesh_id, mesh);

  viewer.ClearLabels();
  Tspline::Vertex_iterator vit;
  for(vit=tsp.vertices_begin(); vit!=tsp.vertices_end(); vit++)
  {
    std::ostringstream os;
    os << vit->data().id;
    Eigen::Vector3d cp;
    vit->data().GetCP(cp);
    viewer.AddLabel3D(os.str(), 20, TomGine::vec3(cp[0],cp[1],cp[2]), 0.0f,0.0f,0.0f);
  }

  viewer.Update();
}
