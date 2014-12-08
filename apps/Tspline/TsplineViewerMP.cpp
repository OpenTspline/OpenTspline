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

#include "Tspline/TsplineMultiPatch.h"
#include "Tspline/TsplineCreator.h"
#include "Tspline/Utils.hpp"
#include "Tspline/File.h"

using namespace tspline;
using namespace TomGine;

tgTomGineThread viewer (800, 600, "T-spline Viewer", true);

int main (int argc, char * argv[])
{
  //  viewer.SetClearColor(1.0f);
  std::string tmp_file;
  bool isolines(false);
  bool use_texture(true);

  for ( int i = 1; i < argc; i++ )
  {
    if( strcmp ( argv[i], "-tf" ) == 0 && i<argc-1 )
      tmp_file = argv[i+1];
    if( strcmp ( argv[i], "-isolines" ) == 0 )
      isolines=true;
    if( strcmp ( argv[i], "-no-texture" ) == 0 )
      use_texture = false;
  }

  if(tmp_file.empty())
  {
    printf("Usage: TsplineMPViewer -tf tmp_file [-isolines] [-no-texture]\n\n");
    return 0;
  }

  viewer.SetClearColor(0.5f);

  std::string tex_name, path;
  cut_file_name(tmp_file, tex_name, path);

  // Load Tspline
  TsplineMultiPatch tsmp;
  File::Load(tsmp, tmp_file);

  // Load textures
  std::vector<cv::Mat> textures(tsmp.patchlist.size());
  std::list<tspline::TsplinePatch*>::iterator pit;
  size_t i(0);
  bool textured(false);
  for(pit=tsmp.patchlist.begin(); pit!=tsmp.patchlist.end(); pit++)
  {
    tspline::TsplinePatch* tsp = (*pit);
    if(!tsp->texture.empty())
    {
      std::string texture = path + "/" + tsp->texture;
      textures[i] = cv::imread(texture);
      textured = true;
    }
    i++;
  }

  TomGine::tgTextureModel mesh;
  if(use_texture && textured)
    convertTspline2tgTextureModel(tsmp, textures, mesh, 8, 8, true);
  else
    convertTspline2tgModel(tsmp, mesh, 8, 8, true);

  convertTsplineControl2tgModel(tsmp, mesh);
  for(size_t i=0; i<mesh.m_lines.size(); i++)
  {
    vec3 a = mesh.m_lines[i].start;
    vec3 b = mesh.m_lines[i].end;
    viewer.AddLine3D(a.x,a.y,a.z, b.x,b.y,b.z, 0,0,128, 1.0);
  }

  std::vector<tspline::Tspline::Vertex_iterator> cps = tsmp.get_controlpoints();
  for(size_t i=0; i<cps.size(); i++)
  {
    tspline::Tspline::Vertex_const_iterator vit = cps[i];
    const TVertex& vext = vit->data();
    const Point3d& cp = vext.GetCP();
    viewer.AddPoint3D(cp.x(),cp.y(),cp.z(), 0,0,255, 5.0f);
    std::ostringstream os;
    os << vext.GetPrimaryID();
    viewer.AddLabel3D(os.str(), 12, cp.x(),cp.y(),cp.z());
  }
  mesh.m_lines.clear();
  mesh.m_points.clear();

  if(isolines)
    convertTsplineEdges2tgModel(tsmp, mesh);

  //  convertTsplineControl2tgModel(tsmp, mesh);
  mesh.ComputeBoundingSphere();

  viewer.SetRotationCenter(mesh.m_bs.center);
  viewer.LookAt(mesh.m_bs.center);
  viewer.AddModel3D(mesh);


  //  // schnipp
  //  viewer.WaitForEvent (TMGL_Press, TMGL_Space);

  //  TsplinePatch* tsp = tsmp.patchlist.front();

  //  std::vector<int> vertex_ids;

  //  vertex_ids.push_back( 227 );
  //  vertex_ids.push_back( 81 );

  //  for(size_t i=0; i<vertex_ids.size(); i++)
  //  {
  //    printf("[main] removing vertex %d\n", vertex_ids[i]);
  //    tspline::Tspline::Vertex_iterator vit = tsp->get_vertex(vertex_ids[i]);
  //    if(vit==tsp->vertices_end())
  //    {
  //      printf("[main] vertex not valid.\n");
  //      continue;
  //    }
  //    tsp->remove_vertex(vit, false, false, false);

  //    viewer.Clear();
  //    for(vit=tsp->vertices_begin(); vit!=tsp->vertices_end(); vit++)
  //    {
  //      std::ostringstream os;
  //      os << vit->data().id;
  //      const Point3d& p = vit->data().GetCP();
  //      viewer.AddLabel3D(os.str(),8,p.x(),p.y(),p.z());
  //    }

  //    TomGine::tgModel mesh_sparse;
  //    convertTspline2tgModel(tsmp, mesh_sparse, 32, 32);
  //    convertTsplineControl2tgModel(tsmp, mesh_sparse);
  //    viewer.AddModel3D(mesh_sparse);

  //    //    viewer.Update();
  //    //    viewer.WaitForEvent (TMGL_Press, TMGL_Space);
  //  }

  //  printf("[main] remove redundant\n");
  //  while(tsp->remove_redundant_vertices()>0) {};
  //  printf("[main] update vertex ids\n");
  //  tsp->update_vertex_ids();
  //  printf("[main] insert edges at L junctions\n");
  //  tsp->insert_edges_at_L_junctions();
  //  printf("[main] insert missing edges\n");
  //  tsp->insert_missing_edges();
  //  printf("[main] update knot vectors\n");
  //  tsp->update_knot_vectors();

  //  printf("[main] drawing\n");
  //  viewer.Clear();
  //  tspline::Tspline::Vertex_iterator vit;
  //  for(vit=tsp->vertices_begin(); vit!=tsp->vertices_end(); vit++)
  //  {
  //    std::ostringstream os;
  //    os << vit->data().id;
  //    const Point3d& p = vit->data().GetCP();
  //    viewer.AddLabel3D(os.str(),8,p.x(),p.y(),p.z());
  //  }

  //  TomGine::tgModel mesh_sparse;
  //  convertTspline2tgModel(tsmp, mesh_sparse, 32, 32);
  //  convertTsplineControl2tgModel(tsmp, mesh_sparse);
  //  viewer.AddModel3D(mesh_sparse);


  //  printf("[main] done\n");

  //  // schnapp


  viewer.Update();
  viewer.WaitForEvent (TMGL_Press, TMGL_Escape);
  return (0);
}
