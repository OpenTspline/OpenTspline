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
#include "Tspline/File.h"

using namespace tspline;
using namespace TomGine;



int main (int argc, char * argv[])
{
  std::string tsp_file;
  std::string ply_file;
  size_t width(800), height(600);
  bool isolines(false);
  bool use_texture(true);
  unsigned resolution(16);

  if(argc > 1)
    tsp_file = argv[1];

  // Input
  for ( int i = 1; i < argc; i++ )
  {
    if( strcmp ( argv[i], "-o" ) == 0 && i<argc-1 )
      ply_file = argv[i+1];
    if( strcmp ( argv[i], "-r" ) == 0 && i<argc-1 )
      resolution = atoi(argv[i+1]);
    if( strcmp ( argv[i], "-isolines" ) == 0 )
      isolines=true;
    if( strcmp ( argv[i], "-no-texture" ) == 0 )
      use_texture = false;
  }

  if(tsp_file.empty())
  {
    printf("Usage: TsplineViewer tsp_file [-o ply_file] [-r resolution] [-isolines] [-no-texture]\n\n");
    return 0;
  }

  std::ostringstream os;
  os << "T-spline Viewer (" << tsp_file << ")";
  tgTomGineThread viewer (width, height, os.str());
  viewer.SetClearColor(1.0f);

  // Tspline
  Tspline tsp;
  File::Load(tsp, tsp_file);

  TomGine::tgTextureModel meshTsp;
  meshTsp.m_material.Color( 0.1,0.1,0.1,1.0,
                            0.6,1.0,0.6,1.0,
                            0.5,0.5,0.5,1.0, 50);

  tspline::Tspline::Face_iterator fit;
  for(fit=tsp.faces_begin(); fit!=tsp.faces_end(); fit++)
    if(!fit->is_unbounded())
      convertTsplineFace2tgModel(fit, tsp, meshTsp, resolution, resolution);

  meshTsp.ComputeBoundingSphere();
  viewer.SetInputSpeeds(1.0f, meshTsp.m_bs.radius, meshTsp.m_bs.radius);

  if(use_texture && !tsp.texture.empty())
  {
    // texture from tspline
    printf("[main] tsp.texture.size: %lu  tsp.texture: %s\n", tsp.texture.size(), tsp.texture.c_str());
    unsigned found = tsp_file.find_last_of("/\\");
    std::string texture = tsp_file.substr(0,found) + "/" + tsp.texture;
    cv::Mat img = cv::imread(texture);
    if(!img.empty())
    {
      meshTsp.m_tex_cv.push_back(img);
      meshTsp.m_face_tex_id.assign(meshTsp.m_faces.size(), 0);
      meshTsp.SetColor(1.0f,1.0f,1.0f);
    }
    else
      printf("[main] Error, could not load texture '%s'\n", texture.c_str());
  }

  convertTsplineControl2tgModel(tsp, meshTsp);

  TomGine::vec3 cor(0,0,0);
  for(size_t i=0; i<meshTsp.m_lines.size(); i++)
  {
    vec3 a = meshTsp.m_lines[i].start;
    vec3 b = meshTsp.m_lines[i].end;
    viewer.AddLine3D(a.x,a.y,a.z, b.x,b.y,b.z, 0,0,0, 2.0);
  }

  size_t cpi(0);
  Tspline::Vertex_iterator vit;
  for(vit=tsp.vertices_begin(); vit!=tsp.vertices_end(); vit++)
  {
    Point3d cp = vit->data().GetCP();
    viewer.AddPoint3D(cp.x(),cp.y(),cp.z(), 0,0,0, 5.0f);
    std::ostringstream os;
    os << cpi;
    viewer.AddLabel3D(os.str(), 12, cp.x(),cp.y(),cp.z());
    cor+=vec3(cp.x(),cp.y(),cp.z());
    cpi++;
  }
  cor /= cpi;
  viewer.LookAt(cor);
  viewer.SetRotationCenter(cor);

  meshTsp.m_lines.clear();
  meshTsp.m_points.clear();

  if(isolines)
  {
    convertTsplineEdges2tgModel(tsp, meshTsp, resolution, meshTsp.m_bs.radius*4e-3);
    meshTsp.m_line_color = TomGine::vec3(1,0,0);
    meshTsp.m_line_width = 4.0;
  }
  //  meshTsp.SetColor(0.5f,0.5f,0.5f);

//  // schnipp
//  viewer.SetClearColor(1.0f);
//  //  mesh.ComputeNormals();
//  meshTsp.SetColor(0.5f,0.5f,0.5f);
//  TomGine::tgCamera cam = viewer.GetCamera();

//  // liver0
//  float ext_data[16] = {  -0.504010, 0.643584, -0.576026, 48.225060,
//                          0.612255, 0.736638, 0.287327, -315.040070,
//                          0.609232, -0.207850, -0.765292, -142.501160,
//                          0.000000, 0.000000, 0.000000, 1.000000 };

//  // liver0
//  float int_data[16] = {  1.344443, 0.000000, 0.000000, 0.000000,
//                          0.000000, 1.792591, 0.000000, 0.000000,
//                          0.000000, 0.000000, -1.000031, -0.020000,
//                          0.000000, 0.000000, -1.000000, 0.000000};

//  TomGine::mat4 E(ext_data);
//  TomGine::mat4 I(int_data);
//  cam.SetIntrinsic(I.transpose());
//  cam.SetExtrinsic(E.transpose());
//  viewer.SetCamera(cam);
//  // schnapp

  viewer.AddModel3D(meshTsp);

  viewer.WaitForEvent (TMGL_Press, TMGL_Space);
  viewer.Snapshot("liver_tspline.png");

  printf("[mesh] vertices: %lu faces: %lu\n", meshTsp.m_vertices.size(), meshTsp.m_faces.size());
  printf("[tspline] vertices: %lu faces: %lu  clamped: %d\n",
         tsp.number_of_vertices(), tsp.number_of_faces(), tsp.clamped);

  if(!ply_file.empty())
    TomGine::tgModelLoader::SavePly(meshTsp, ply_file);

  viewer.Update();
  viewer.WaitForEvent (TMGL_Press, TMGL_Escape);
  return (0);
}
