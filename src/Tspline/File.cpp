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

#include "File.h"

using namespace tspline;

template<typename T>
size_t File::read_vector (std::vector<T> &vec, FILE* pFile)
{
  size_t result (0);
  size_t s;
  result += fread (&s, sizeof(size_t), 1, pFile) * sizeof(size_t);
  vec.assign (s, T ());
  for (size_t i = 0; i < s; i++)
    result += fread (&vec[i], sizeof(T), 1, pFile) * sizeof(T);

  return result;
}

size_t File::read_string(std::string &str, FILE* pFile)
{
  std::vector<char> cvec;
  size_t result = read_vector(cvec, pFile);

  if(cvec.size() > 1) // if cvec contains more than termination char (0)
    str = std::string(cvec.begin(),cvec.end());

  return result;
}

size_t File::read_point(Point2d &p, FILE* pFile)
{
  size_t result(0);
  double x, y;
  result += fread(&x, sizeof(double), 1, pFile) * sizeof(double);
  result += fread(&y, sizeof(double), 1, pFile) * sizeof(double);
  p = Point2d(x,y);
  return result;
}

size_t File::read_point(Point3d &p, FILE* pFile)
{
  size_t result(0);
  double x, y, z;
  result += fread(&x, sizeof(double), 1, pFile) * sizeof(double);
  result += fread(&y, sizeof(double), 1, pFile) * sizeof(double);
  result += fread(&z, sizeof(double), 1, pFile) * sizeof(double);
  p = Point3d(x,y,z);
  return result;
}

size_t File::read_point(Point4d &p, FILE* pFile)
{
  size_t result(0);
  double x, y, z, w;
  result += fread(&x, sizeof(double), 1, pFile) * sizeof(double);
  result += fread(&y, sizeof(double), 1, pFile) * sizeof(double);
  result += fread(&z, sizeof(double), 1, pFile) * sizeof(double);
  result += fread(&w, sizeof(double), 1, pFile) * sizeof(double);
  p = Point4d(x,y,z,w);
  return result;
}

size_t File::read_vertex(Tspline &tsp, FILE* pFile, Link &link)
{
  size_t result(0);

  TVertex vext;

  Point2d p;
  result += read_point(p, pFile);

  Tspline::Vertex_iterator vit = insert_point(tsp, p);

  Point4d cp;
  result += read_point(cp, pFile);
  vext.SetCP(cp);

  result += read_vector<int>(link.patch_ids, pFile);
  result += read_vector<int>(link.vertex_ids, pFile);

  result += fread(&vext.id, sizeof(int), 1, pFile) * sizeof(int);
  result += fread(&vext.patch_id, sizeof(int), 1, pFile) * sizeof(int);
  result += read_point(vext.param, pFile);

  result += read_vector<double>(vext.s, pFile);
  result += read_vector<double>(vext.t, pFile);

  result += fread(&vext.updated, sizeof(bool), 1, pFile) * sizeof(bool);

  link.vertex_id = vext.id;
  link.patch_id = vext.patch_id;

  vit->set_data(vext);

  return result;
}

size_t File::read_edge(Tspline &tsp, FILE* pFile)
{
  size_t result(0);

  int id_source;
  result += fread(&id_source, sizeof(int), 1, pFile) * sizeof(int);

  int id_target;
  result += fread(&id_target, sizeof(int), 1, pFile) * sizeof(int);

  THalfedge hext;
  result += fread(&hext.d, sizeof(double), 1, pFile) * sizeof(double);

  Tspline::Vertex_iterator vs = tsp.get_vertex(id_source);
  Tspline::Vertex_iterator vt = tsp.get_vertex(id_target);

  Tspline::Halfedge_iterator hit = tsp.insert_at_vertices(Segment2(vs->point(), vt->point()), vs, vt);
  hit->set_data(hext);
  hit->twin()->set_data(hext);

  return result;
}

size_t File::read_tspline(Tspline &tsp, FILE* pFile, std::vector<Link> &links)
{
  size_t result(0);

  size_t vs;
  result += fread(&vs, sizeof(size_t), 1, pFile) * sizeof(size_t);
  for(size_t i=0; i<vs; i++)
  {
    Link link;
    result += read_vertex(tsp, pFile, link);
    links.push_back(link);
  }

  size_t es;
  result += fread(&es, sizeof(size_t), 1, pFile) * sizeof(size_t);
  for(size_t i=0; i<es; i++)
  {
    result += read_edge(tsp, pFile);
  }

  result += fread(&tsp.degree, sizeof(unsigned), 1, pFile) * sizeof(unsigned);
  result += fread(&tsp.clamped, sizeof(bool), 1, pFile) * sizeof(bool);
  result += read_point(tsp.param_min, pFile);
  result += read_point(tsp.param_max, pFile);
  result += read_point(tsp.point_min, pFile);
  result += read_point(tsp.point_max, pFile);
  result += read_string(tsp.texture, pFile);

  return result;
}

size_t File::read_tspline_patch(tspline::TsplinePatch &tsp, FILE* pFile)
{
  size_t result(0);

  std::vector<Link> links;
  result += read_tspline(tsp, pFile, links);

  result += fread(&tsp.id, sizeof(size_t), 1, pFile) * sizeof(size_t);

  return result;
}

size_t File::read_tspline_multipatch(TsplineMultiPatch &tsp, FILE* pFile)
{
  size_t result(0);

  size_t ps;
  result += fread(&ps, sizeof(size_t), 1, pFile) * sizeof(size_t);

  for(size_t i=0; i<ps; i++)
  {
    TsplinePatch* patch = new TsplinePatch(0);
    result += read_tspline_patch(*patch, pFile);
    tsp.patchlist.push_back(patch);
  }

  size_t ts;
  result += fread(&ts, sizeof(size_t), 1, pFile) * sizeof(size_t);
  for(size_t i=0; i<ts; i++)
  {
    Trans t;
    result += fread(&t.source_patch_id, sizeof(size_t), 1, pFile) * sizeof(size_t);
    result += fread(&t.source_dir, sizeof(Direction), 1, pFile) * sizeof(Direction);
    result += fread(&t.source_range, sizeof(double), 2, pFile) * sizeof(double);
    result += fread(&t.target_patch_id, sizeof(size_t), 1, pFile) * sizeof(size_t);
    result += fread(&t.target_dir, sizeof(Direction), 1, pFile) * sizeof(Direction);
    result += fread(&t.target_range, sizeof(double), 2, pFile) * sizeof(double);
    result += fread(&t.multiplicity, sizeof(unsigned char), 1, pFile) * sizeof(unsigned char);

    TsplinePatch* p = tsp.GetPatch(t.source_patch_id);
    if(p->HasTransition(t.source_dir, t.target_patch_id))
      continue;

    tsp.AddTransition(t.source_patch_id, t.source_dir, t.source_range[0], t.source_range[1],
        t.target_patch_id, t.target_dir, t.target_range[0], t.target_range[1], t.multiplicity);
  }

  return result;
}

template<typename T>
void File::write_vector(const std::vector<T> &vec, FILE* pFile)
{
  size_t s = vec.size ();
  fwrite (&s, sizeof(size_t), 1, pFile);
  for (size_t i = 0; i < vec.size (); i++)
    fwrite (&vec[i], sizeof(T), 1, pFile);
}

void File::write_string(const std::string& str, FILE* pFile)
{
  std::vector<char> vchar(str.begin(), str.end());
  write_vector<char>(vchar, pFile);
}

void File::write_point(const Point2d &p, FILE *pFile)
{
  double x = p.x();
  double y = p.y();
  fwrite(&x, sizeof(double), 1, pFile);
  fwrite(&y, sizeof(double), 1, pFile);
}

void File::write_point(const Point3d &p, FILE *pFile)
{
  double x = p.x();
  double y = p.y();
  double z = p.z();
  fwrite(&x, sizeof(double), 1, pFile);
  fwrite(&y, sizeof(double), 1, pFile);
  fwrite(&z, sizeof(double), 1, pFile);
}

void File::write_point(const Point4d &p, FILE *pFile)
{
  double x = p.x();
  double y = p.y();
  double z = p.z();
  double w = p.w();
  fwrite(&x, sizeof(double), 1, pFile);
  fwrite(&y, sizeof(double), 1, pFile);
  fwrite(&z, sizeof(double), 1, pFile);
  fwrite(&w, sizeof(double), 1, pFile);
}

void File::write_vertex(Tspline::Vertex_const_iterator &vit, FILE* pFile)
{
  write_point(vit->point(), pFile);

  const TVertex& vext = vit->data();

  write_point(vext.GetCP(), pFile);

  std::vector<int> patch_ids;
  std::vector<int> vertex_ids;
  vext.GetLinks(patch_ids, vertex_ids);
  write_vector<int>(patch_ids, pFile);
  write_vector<int>(vertex_ids, pFile);

  fwrite(&vext.id, sizeof(int), 1, pFile);
  fwrite(&vext.patch_id, sizeof(int), 1, pFile);
  write_point(vext.param, pFile);

  write_vector<double>(vext.s, pFile);
  write_vector<double>(vext.t, pFile);

  fwrite(&vext.updated, sizeof(bool), 1, pFile);
}

void File::write_edge(Tspline::Edge_const_iterator &eit, FILE* pFile)
{
  const TVertex &vext_source = eit->source()->data();
  fwrite(&vext_source.id, sizeof(int), 1, pFile);

  const TVertex &vext_target = eit->target()->data();
  fwrite(&vext_target.id, sizeof(int), 1, pFile);

  const THalfedge &hext = eit->data();
  fwrite(&hext.d, sizeof(double), 1, pFile);
}

void File::write_tspline(const tspline::Tspline &tsp, FILE* pFile)
{
  size_t vs = tsp.number_of_vertices();
  fwrite (&vs, sizeof(size_t), 1, pFile);
  Tspline::Vertex_const_iterator vit;
  for(vit=tsp.vertices_begin(); vit!=tsp.vertices_end(); vit++)
    write_vertex(vit, pFile);

  size_t es = tsp.number_of_edges();
  fwrite (&es, sizeof(size_t), 1, pFile);
  Tspline::Edge_const_iterator eit;
  for(eit=tsp.edges_begin(); eit!=tsp.edges_end(); eit++)
    write_edge(eit, pFile);

  fwrite(&tsp.degree, sizeof(unsigned), 1, pFile);
  fwrite(&tsp.clamped, sizeof(bool), 1, pFile);
  write_point(tsp.param_min, pFile);
  write_point(tsp.param_max, pFile);
  write_point(tsp.point_min, pFile);
  write_point(tsp.point_max, pFile);
  write_string(tsp.texture, pFile);
}

void File::write_tspline_patch(const tspline::TsplinePatch &tsp, FILE* pFile, std::vector<Trans> &t)
{
  write_tspline(tsp, pFile);

  fwrite(&tsp.GetID(), sizeof(size_t), 1, pFile);

  for(size_t i=0; i<tsp.m_north.size(); i++)
  {
    const Transition* trans = tsp.m_north[i];
    double r0a, r0b, r1a, r1b;
    trans->GetSourceRange(r0a, r0b);
    trans->GetTargetRange(r1a, r1b);
    t.push_back( Trans(tsp.GetID(), trans->GetSourceDir(), r0a, r0b,
                       trans->target->GetID(), trans->GetTargetDir(), r1a, r1b, trans->GetMultiplicity()) );
  }

  for(size_t i=0; i<tsp.m_east.size(); i++)
  {
    const Transition* trans = tsp.m_east[i];
    double r0a, r0b, r1a, r1b;
    trans->GetSourceRange(r0a, r0b);
    trans->GetTargetRange(r1a, r1b);
    t.push_back( Trans(tsp.GetID(), trans->GetSourceDir(), r0a, r0b,
                       trans->target->GetID(), trans->GetTargetDir(), r1a, r1b, trans->GetMultiplicity()) );
  }

  for(size_t i=0; i<tsp.m_south.size(); i++)
  {
    const Transition* trans = tsp.m_south[i];
    double r0a, r0b, r1a, r1b;
    trans->GetSourceRange(r0a, r0b);
    trans->GetTargetRange(r1a, r1b);
    t.push_back( Trans(tsp.GetID(), trans->GetSourceDir(), r0a, r0b,
                       trans->target->GetID(), trans->GetTargetDir(), r1a, r1b, trans->GetMultiplicity()) );
  }

  for(size_t i=0; i<tsp.m_west.size(); i++)
  {
    const Transition* trans = tsp.m_west[i];
    double r0a, r0b, r1a, r1b;
    trans->GetSourceRange(r0a, r0b);
    trans->GetTargetRange(r1a, r1b);
    t.push_back( Trans(tsp.GetID(), trans->GetSourceDir(), r0a, r0b,
                       trans->target->GetID(), trans->GetTargetDir(), r1a, r1b, trans->GetMultiplicity()) );
  }
}

void File::write_tspline_multipatch(tspline::TsplineMultiPatch &tsmp, FILE* pFile)
{
  size_t ps = tsmp.patchlist.size();
  fwrite(&ps, sizeof(size_t), 1, pFile);

  std::vector<Trans> trans;
  std::list<TsplinePatch*>::iterator it;
  for(it=tsmp.patchlist.begin(); it!=tsmp.patchlist.end(); it++)
    write_tspline_patch(*(*it), pFile, trans);

  size_t ts = trans.size();
  fwrite(&ts, sizeof(size_t), 1, pFile);
  for(size_t i=0; i<trans.size(); i++)
  {
    const Trans &t = trans[i];
    fwrite(&t.source_patch_id, sizeof(size_t), 1, pFile);
    fwrite(&t.source_dir, sizeof(Direction), 1, pFile);
    fwrite(&t.source_range, sizeof(double), 2, pFile);
    fwrite(&t.target_patch_id, sizeof(size_t), 1, pFile);
    fwrite(&t.target_dir, sizeof(Direction), 1, pFile);
    fwrite(&t.target_range, sizeof(double), 2, pFile);
    fwrite(&t.multiplicity, sizeof(unsigned char), 1, pFile);
  }
}

bool File::Save(const tspline::Tspline &tsp, const std::string& filename)
{
  FILE * pFile;
  pFile = fopen (filename.c_str (), "wb");
  if (pFile == NULL)
  {
    printf ("[tspline::File::Save(const Tspline&,..)] Error: cannot open file for writing: '%s'\n", filename.c_str ());
    return false;
  }
  if (tsp.number_of_vertices()==0)
  {
    printf ("[tspline::File::Save(const Tspline&,..)] Warning: no data to save: '%s'\n", filename.c_str ());
    return false;
  }

  write_tspline(tsp, pFile);

  fclose (pFile);

  printf("[File::Save(tspline::Tspline,..)] Saved file '%s'\n", filename.c_str());

  return true;
}

bool File::Load(tspline::Tspline &tsp, const std::string& filename)
{
  tsp = Tspline();
  FILE * pFile;
  pFile = fopen (filename.c_str (), "rb");
  if (pFile == NULL)
  {
    printf ("[tspline::File::Load(tspline::Tspline &tsp,..)] Error: cannot open file for reading: '%s'\n",
            filename.c_str ());
    return false;
  }
  long lSize;
  fseek (pFile, 0, SEEK_END);
  lSize = ftell (pFile);
  rewind (pFile);

  std::vector<Link> links;
  long result = read_tspline (tsp, pFile, links);

  if (result != lSize)
  {
    printf ("[tspline::File::Load(tspline::Tspline &tsp,..)] %ld %ld\n", result, lSize);
    throw std::runtime_error ("[tspline::File::Load(tspline::Tspline &tsp,..)] Reading error (memory size does not match) ");
  }

  fclose (pFile);

  return true;
}

bool File::Save(tspline::TsplineMultiPatch &tsmp, const std::string& filename)
{
  FILE * pFile;
  pFile = fopen (filename.c_str (), "wb");
  if (pFile == NULL)
  {
    printf ("[tspline::File::Save(const TsplineMultiPatch&,..)] Error: cannot open file for writing: '%s'\n", filename.c_str ());
    return false;
  }
  if (tsmp.patchlist.size()==0)
  {
    printf ("[tspline::File::Save(const TsplineMultiPatch&,..)] Warning: no data to save: '%s'\n", filename.c_str ());
    return false;
  }

  write_tspline_multipatch(tsmp, pFile);

  fclose (pFile);

  printf("[File::Save(tspline::TsplineMultiPatch,..)] Saved file '%s'\n", filename.c_str());

  return true;
}

bool File::Load(tspline::TsplineMultiPatch &tsmp, const std::string& filename)
{
  FILE * pFile;
  pFile = fopen (filename.c_str (), "rb");
  if (pFile == NULL)
  {
    printf ("[tspline::File::Load(tspline::TsplineMultiPatch &tsp,..)] Error: cannot open file for reading: '%s'\n",
            filename.c_str ());
    return false;
  }
  long lSize;
  fseek (pFile, 0, SEEK_END);
  lSize = ftell (pFile);
  rewind (pFile);

  long result = read_tspline_multipatch (tsmp, pFile);

  if (result != lSize)
  {
    printf ("[tspline::File::Load(tspline::TsplineMultiPatch &tsp,..)] %ld %ld\n", result, lSize);
    throw std::runtime_error ("[tspline::File::Load(tspline::TsplineMultiPatch &tsp,..)] Reading error (memory size does not match) ");
  }

  fclose (pFile);

  return true;
}
