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

using namespace tspline;
using namespace TomGine;



class SequenceViewer : public tgEventListener
{
public:
  std::vector<TomGine::tgTextureModel> meshes;
  int index;
  int mesh_id;
  int label_id;
  TomGine::tgTomGineThread* viewer;

  void EventFunction(TomGine::Event event)
  {
    if(event.type==TomGine::TMGL_Press && event.input==TomGine::TMGL_Left)
    {
      if(index > 0)
      {
        index--;
        viewer->SetModel3D(mesh_id, meshes[index]);
        std::ostringstream os;
        os << index;
        viewer->SetLabel2D(label_id, os.str(), 20, 5,5);
      }

    }
    if(event.type==TomGine::TMGL_Press && event.input==TomGine::TMGL_Right)
    {
      if(index < static_cast<int>(meshes.size())-1)
      {
        index++;
        viewer->SetModel3D(mesh_id, meshes[index]);
        std::ostringstream os;
        os << index;
        viewer->SetLabel2D(label_id, os.str(), 20, 5,5);
      }
    }
  }

  SequenceViewer(TomGine::tgTomGineThread& tg) : index(0)
  {
    viewer = &tg;
    viewer->RegisterEventListener(this);
    meshes.resize(1);
    mesh_id = viewer->AddModel3D(meshes.back());
    label_id = viewer->AddLabel2D("0", 20, 5,5);
  }
};


int main (int argc, char * argv[])
{
  std::string tsp_file;
  size_t width(800), height(600);
  int start(0), end(1), res(16);

  if (argc > 3)
  {
    tsp_file = argv[1];
    start = atoi(argv[2]);
    end = atoi(argv[3]);
  }
  else
  {
    printf("Usage: TsplineViewer tsp_file_sequence start_index end_index\n");
    printf("   tsp-file-sequence:  filename_%%d.tsp\n\n");
    return 0;
  }

  if(argc > 4)
    res = atoi(argv[4]);

  std::ostringstream os;
  os << "T-spline Sequence (" << tsp_file << ")";
  tgTomGineThread viewer (width, height, os.str());
  viewer.SetClearColor(1.0f);
  viewer.SetInputSpeeds(0.1);

  // schnipp
  TomGine::tgCamera cam = viewer.GetCamera();
//    dino_sparse_front
  float ext_data[16] = {  -0.382363, -0.097398, -0.918871, 0.000105,
                          -0.181148, 0.983042, -0.028819, -0.027788,
                          0.906093, 0.155430, -0.393519, -0.365332,
                          0.000000, 0.000000, 0.000000, 1.000000 };
//     dino
  float int_data[16] = { 10.344999, 0.000000, 0.010219, 0.000000,
                          0.000000, 13.856250, -0.164375, 0.000000,
                          0.000000, 0.000000, -1.000200, -0.020002,
                          0.000000, 0.000000, -1.000000, 0.000000};
  TomGine::mat4 E(ext_data);
  TomGine::mat4 I(int_data);
  cam.SetIntrinsic(I.transpose());
  cam.SetExtrinsic(E.transpose());
  viewer.SetCamera(cam);
  // schnapp

  SequenceViewer sv(viewer);
  sv.meshes.resize(end-start+1);

  for(int i=start; i<=end && !viewer.Stopped(); i++)
  {
    // Tspline
    char file[256];
    sprintf(file, tsp_file.c_str(), i);
    printf("[main] Loading T-spline %d: '%s'\n", i, file);
    Tspline tsp;
    File::Load(tsp, file);

    TomGine::tgTextureModel& mesh = sv.meshes[i];
    mesh.m_material.Color( 0.1,0.1,0.1,1.0,
                           0.6,1.0,0.6,1.0,
                           0.5,0.5,0.5,1.0, 50);
    mesh.SetColor(0.5f,0.5f,0.5f,0.5f);

    convertTspline2tgModel(tsp, mesh, res, res, true);
    convertTsplineControl2tgModel(tsp, mesh);
//    convertTsplineEdges2tgModel(tsp, mesh, 4, 1e-4);
    mesh.m_line_color = TomGine::vec3(0,0,0);
    mesh.m_line_width = 5.0f;
    mesh.m_point_color = TomGine::vec3(0,0,0);
    mesh.m_point_size = 15.0f;

    if(!tsp.texture.empty())
    {
      // texture from tspline
      printf("[main] tseigen normalizep.texture.size: %lu  tsp.texture: %s\n", tsp.texture.size(), tsp.texture.c_str());
      unsigned found = tsp_file.find_last_of("/\\");
      std::string texture = tsp_file.substr(0,found) + "/" + tsp.texture;
      mesh.m_tex_cv.push_back(cv::imread(texture));
      mesh.m_face_tex_id.assign(mesh.m_faces.size(), 0);
      mesh.SetColor(1.0f,1.0f,1.0f);
    }

    //    convertTsplineControl2tgModel(tsp, mesh);

    //    TomGine::vec3 cor(0,0,0);
    //    for(size_t i=0; i<mesh.m_lines.size(); i++)
    //    {
    //      vec3 a = mesh.m_lines[i].start;
    //      vec3 b = mesh.m_lines[i].end;
    //      viewer.AddLine3D(a.x,a.y,a.z, b.x,b.y,b.z, 0,0,128, 1.0);
    //    }

    //    size_t cpi(0);
    //    Tspline::Vertex_iterator vit;
    //    for(vit=tsp.vertices_begin(); vit!=tsp.vertices_end(); vit++)
    //    {
    //      Point3d cp = vit->data().GetCP();
    //      viewer.AddPoint3D(cp.x(),cp.y(),cp.z(), 0,0,255, 5.0f);
    //      std::ostringstream os;
    //      os << cpi;
    //      viewer.AddLabel3D(os.str(), 12, cp.x(),cp.y(),cp.z());
    //      cor+=vec3(cp.x(),cp.y(),cp.z());
    //      cpi++;
    //    }
    //    cor /= cpi;
    //    viewer.LookAt(cor);
    //    viewer.SetRotationCenter(cor);

//    mesh.m_lines.clear();
//    mesh.m_points.clear();

    //    convertTsplineEdges2tgModel(tsp, mesh, 4);
    //    mesh.m_line_color = TomGine::vec3(1,0,0);
    //    mesh.m_line_width = 4.0;

    if(i==start)
    {
      mesh.ComputeBoundingSphere();
      viewer.SetInputSpeeds(1.0f, mesh.m_bs.radius, mesh.m_bs.radius);
      viewer.SetRotationCenter(mesh.m_bs.center);
    }
  }

  printf("[main] Loading T-splines done\n");



  //  viewer.AddModel3D(meshes);


  viewer.Update();
  viewer.WaitForEvent (TMGL_Press, TMGL_Escape);
  return (0);
}
