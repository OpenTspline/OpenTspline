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

#include "v4r/Tspline/TsplineMultiPatch.h"
#include "v4r/Tspline/File.h"
#include "v4r/Tspline/Utils.hpp"

using namespace tspline;
using namespace TomGine;

tgTomGineThread viewer (800, 600, "T-spline Box");

int main ()
{
  TsplineMultiPatch tsmp;

  int id0 = tsmp.AddPatch(1, 4);
  int id1 = tsmp.AddPatch(3, 4);

  std::list<TsplinePatch*>::iterator it;
  for(it=tsmp.patchlist.begin(); it!=tsmp.patchlist.end(); it++)
  {
    TsplinePatch& tsp = *(*it);

    Tspline::Vertex_iterator vit;
    for(vit=tsp.vertices_begin(); vit!=tsp.vertices_end(); vit++)
    {
      TVertex vext = vit->data();
      double a = M_PI * vext.param.y();

      if(tsp.GetID() == id0)
        vext.SetCP(Point3d(vext.param.x(), cos(a), sin(a)));
      else if(tsp.GetID() == id1)
      {
        if(!equal(vext.param.x(),0.0) && !equal(vext.param.x(),1.0))
          vext.SetCP(Point3d(1.0-vext.param.x(), cos(a), 0.0));
        else
          vext.SetCP(Point3d(1.0-vext.param.x(), cos(a), -sin(a)));
      }

      vit->set_data(vext);
    }
  }

  tsmp.AddTransition(id1, NORTH, id0, NORTH, 0);
  tsmp.AddTransition(id0, SOUTH, id1, SOUTH, 0);

  //tsmp.AddTransition(id0, NORTH, id0, SOUTH);


//  File::Save(tsmp, "cylinder.tmp");

  TomGine::tgRenderModel mesh;
  mesh.SetColor(0.4f, 0.4f, 0.4f);
  convertTspline2tgModel(tsmp, mesh, 16, 16);
  convertTsplineControl2tgModel(tsmp, mesh);
  viewer.AddModel3D(mesh);
  viewer.Update ();

  viewer.WaitForEvent (TMGL_Press, TMGL_Escape);
  return (0);
}
