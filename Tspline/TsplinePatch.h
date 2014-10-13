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


#ifndef _TSPLINE_TSPLINE_PATCH_H_
#define _TSPLINE_TSPLINE_PATCH_H_

#include "Tspline.h"
#include "Transition.h"
#include "TsplineMultiPatch.h"

namespace tspline
{

class TsplinePatch : public Tspline
{
protected:
  using Tspline::insert_vertex;
  using Tspline::shoot_left;
  using Tspline::shoot_right;
  using Tspline::shoot_up;
  using Tspline::shoot_down;

  struct LinkSort
  {
    LinkSort(double _d, Direction _dir){ d=_d; dir=_dir; }
    double d;
    Direction dir;
    bool operator<(const LinkSort& a) const { return d < a.d; }
  };

  bool quiet;

public:
  std::vector<NorthTransition*> m_north;
  std::vector<EastTransition*> m_east;
  std::vector<SouthTransition*> m_south;
  std::vector<WestTransition*> m_west;

  int id;

protected:
  Vertex_iterator insert_vertex_no_transition(const double &s, const double &t);
  virtual Tspline::Vertex_iterator insert_vertex(Halfedge_iterator &e, Point2d &param);

private:
  TsplinePatch(const TsplinePatch& p) {}
  void operator*(const TsplinePatch& p) {}

public:
  TsplinePatch(const size_t &_id, const Tspline& tsp) : Tspline(tsp), quiet(true), id(_id) { }
  TsplinePatch(const size_t &_id) : quiet(true), id(_id) { }
  ~TsplinePatch();

  bool IsNorth(double s) const;
  bool IsEast(double t) const;
  bool IsSouth(double s) const;
  bool IsWest(double t) const;

  Transition* GetNorth(double s);
  Transition* GetEast(double t);
  Transition* GetSouth(double s);
  Transition* GetWest(double t);

  const Transition* GetNorthConst(double s) const;
  const Transition* GetEastConst(double t) const;
  const Transition* GetSouthConst(double s) const;
  const Transition* GetWestConst(double t) const;

  bool HasTransition(Direction source_dir, int target_id);
  bool HasTransition(int target_id);

  Transition* AddTransition(Direction source_dir, TsplinePatch* target, Direction target_dir,
                            unsigned char mult);
  Transition* AddTransition(Direction d0, double r0a, double r0b,
                            TsplinePatch* p1, Direction d1, double r1a, double r1b,
                            unsigned char mult);

  //  virtual void evaluate_basis(const double &s, const double &t, const TVertex &vext, double &B) const;
  virtual void evaluate_basis(const double &s, const double &t, const TVertex &vext, CPbasis &basis) const;

  virtual void evaluate_basis_sum(const double &s, const double &t, BasisSum &basis,
                                  std::vector<const TVertex*>& links,
                                  std::vector<tspline::CPbasis> &cpbasis) const;
  virtual void evaluate_basis_sum(const double &s, const double &t, BasisSum &basis, const Transition *trans,
                                  std::vector<const TVertex*>& links,
                                  std::vector<tspline::CPbasis> &cpbasis) const;

  virtual Eigen::Vector3d evaluate(const double &s, const double &t) const;
  virtual Eigen::Vector3d evaluate(const double &s, const double &t, Eigen::Vector3d &ts, Eigen::Vector3d &tt) const;

  virtual Eigen::Vector3d evaluate(const double &s, const double &t,
                                   std::vector<tspline::CPbasis> &cpbasis) const;
  virtual Eigen::Vector3d evaluate(const double &s, const double &t, Eigen::Vector3d &ts, Eigen::Vector3d &tt,
                                   std::vector<tspline::CPbasis> &cpbasis) const;

  inline const int& GetID() const
  {
    return id;
  }

  inline void SetQuiet(bool val)
  {
    quiet = val;
  }

  friend class TsplineMultiPatch;
  friend class Transition;

//public:
//  virtual void update_knot_vectors();

protected:
  virtual void shoot_left (Vertex_iterator &v, const unsigned &n, std::vector<double> &knots);
  virtual void shoot_right (Vertex_iterator &v, const unsigned &n, std::vector<double> &knots);
  virtual void shoot_up (Vertex_iterator &v, const unsigned &n, std::vector<double> &knots);
  virtual void shoot_down (Vertex_iterator &v, const unsigned &n, std::vector<double> &knots);

};

}


#endif
