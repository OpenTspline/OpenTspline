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

#include "v4r/Tspline/Utils.hpp"
#include "v4r/Tspline/File.h"

#include "v4r/Tspline/TsplineCreator.h"

TomGine::tgTomGineThread viewer (800, 600, "T-spline", false, 0.1, 100.0);

size_t w (11);
size_t h (11);

using namespace tspline;

int main ()
{
  viewer.SetClearColor (1.0);
  viewer.Update ();

  tspline::Tspline tsp;
  tspline::TsplineCreator::CreatePlaneXY(tsp,-1.0,-1.0,0.0,2.0,2.0,4,4);


  tsp.remove_vertex(tsp.get_vertex(6));
  tsp.remove_vertex(tsp.get_vertex(7));
  tsp.remove_vertex(tsp.get_vertex(11));

  tspline::Tspline::Vertex_iterator vit;
  for(vit=tsp.vertices_begin(); vit!=tsp.vertices_end(); vit++)
  {
    tspline::TVertex& vext = vit->data();
    std::ostringstream os;
    os << vext.id; //vext.param.x() << "|" << vext.param.y();
    Point3d p = vext.GetCP();
    p = Point3d(p.x(),p.y(),sin(4.0*p.x())*sin(4.0*p.y()));
//    if(p.y() > -0.5 && p.x()>-1.0 && p.x()<1.0)
//      p = Point3d(0.0,p.y(),1.0);
    vext.SetCP(p);
    viewer.AddLabel3D(os.str(), 20, p.x(),p.y(),p.z());

//    printf("[%d] s: %f %f %f %f %f\n", vext.id, vext.s[0], vext.s[1], vext.s[2], vext.s[3], vext.s[4]);
//    printf("[%d] t: %f %f %f %f %f\n", vext.id, vext.t[0], vext.t[1], vext.t[2], vext.t[3], vext.t[4]);
//    printf("\n");
  }

  TomGine::tgRenderModel model;
  model.SetColor (0.5f, 0.3f, 0.8f, 1.0f);
  convertTspline2tgModel (tsp, model, 30, 30);
  convertTsplineControl2tgModel(tsp, model);
  viewer.AddModel3D (model);

  //  tsp.set_quiet(false);
  //  Eigen::Vector3d ts, tt;
  //  Eigen::Vector3d p = tsp.evaluate(1.1, 1.0, ts, tt);
  //  viewer.AddPoint3D(p(0),p(1),p(2), 255,0,0, 5.0f);
  //  printf("[main]  %e %e %e\n", p(0),p(1),p(2));
  //  printf("[main]  %e %e %e\n", ts(0),ts(1),ts(2));
  //  printf("[main]  %e %e %e\n", tt(0),tt(1),tt(2));


  TomGine::tgCamera cam = viewer.GetCamera();
  TomGine::mat4 e;
  e.identity();
  cam.SetExtrinsic(e);
  cam.TranslateF(-4.0f);
  viewer.SetCamera(cam);

  viewer.Update ();
  viewer.WaitForEvent (TomGine::TMGL_Press, TomGine::TMGL_Space);
  return (0);
}
