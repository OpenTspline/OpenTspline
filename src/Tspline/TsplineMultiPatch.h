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


#ifndef _TSPLINE_TSPLINE_MULTI_PATCH_H_
#define _TSPLINE_TSPLINE_MULTI_PATCH_H_

#include "TsplinePatch.h"

namespace tspline
{

class TsplineMultiPatch
{
protected:
  bool quiet;

  void CopyTransitions(const TsplineMultiPatch& tsmp);

public:
  std::list<TsplinePatch*> patchlist; // todo: move to private

private:
  TsplineMultiPatch(const TsplineMultiPatch& tsmp);// TODO this causes a seg fault at the moment (disable by private)

public:
  TsplineMultiPatch();
  ~TsplineMultiPatch();

  virtual int AddPatch(unsigned segX, unsigned segY,
                          double s0=0.0, double s1=1.0, double t0=0.0, double t1=1.0);

  virtual int AddPatch(const Tspline& _tsp);

  void AddTransition(const size_t& id0, const Direction& ed0,
                     const size_t& id1, const Direction& ed1,
                     unsigned char mult);

  void AddTransition(const int &id0, const Direction& ed0, double r0a, double r0b,
                     const int &id1, const Direction& ed1, double r1a, double r1b,
                     unsigned char mult);

  void InsertVertex(const int& id, const Point2d& param);

  vector_vec3d compute_cp_normals_by_footpoints() const;

  inline void SetQuiet(bool val)
  {
    quiet = val;
    std::list<TsplinePatch*>::iterator it;
    for(it=patchlist.begin(); it!=patchlist.end(); it++)
      (*it)->SetQuiet(val);
  }

  inline bool IsQuiet() const
  {
    return quiet;
  }

  inline TsplinePatch* GetPatch(const int& id)
  {
    std::list<TsplinePatch*>::iterator it;
    for(it=patchlist.begin(); it!=patchlist.end(); it++)
      if((*it)->GetID() == id)
        return (*it);
    return NULL;
  }

  /** @brief returns non-dublicate controlpoints of Tspline */
  std::vector<tspline::Tspline::Vertex_iterator> get_controlpoints();
  std::vector<tspline::Tspline::Vertex_const_iterator> get_controlpoints() const;

  /** @brief returns controlpoints of Tspline including dublicates (i.e. linked)*/
  std::vector<tspline::Tspline::Vertex_iterator> get_vertices();
  std::vector<tspline::Tspline::Vertex_const_iterator> get_vertices() const;

  size_t number_of_controlpoints() const; // vertices without dublicates (i.e. linked)
  size_t number_of_vertices() const; // all vertices

  Eigen::Vector3d center() const;

  bool empty() const { return patchlist.empty(); }



  //std::list<TsplinePatch*>::const_iterator GetPatch(const size_t& id) const;

};

}


#endif
