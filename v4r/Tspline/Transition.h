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


#ifndef _TSPLINE_TRANSITION_H_
#define _TSPLINE_TRANSITION_H_

// This file defines transitions between T-spline patches
// It implements coordinate conversion between neighboring domains, multi-patch traversal,
// range computation and inter-transition basis function evaluation

#include "TsplinePatch.h"

namespace tspline
{

class TsplinePatch;

/** @brief  directed transition between T-spline patches
 *          for coordinate conversion, transition creation, T-spline multi-patch traversal */
class Transition
{
protected:
  Direction source_dir;       /// direction (i.e. side of border) of source domain
  Direction target_dir;       /// direction (i.e. side of border) of source domain
  double source_range[2];     /// range of transition on the border of the source domain
  double target_range[2];     /// range of transition on the border of the target domain
  unsigned char multiplicity; /// parametric knot multiplicity over the transition (0 = smooth, 1 = sharp)

public:
  TsplinePatch* source; /// pointer to source patch
  TsplinePatch* target; /// pointer to target patch
  Transition* twin;     /// twin transition (pointing in the other direction)

public:
  Transition() : multiplicity(0), source(NULL), target(NULL)
  {
    source_range[0] = 0.0;
    source_range[1] = 1.0;
    target_range[0] = 0.0;
    target_range[1] = 1.0;
  }

  Transition(TsplinePatch* p0, Direction d0,
             TsplinePatch* p1, Direction d1,
             unsigned char m)
  {
    if(p0==NULL || p1==NULL)
      throw std::runtime_error("[Transition::Transition] Error, patch not valid");

    source = p0;
    source_dir = d0;
    target = p1;
    target_dir = d1;
    source_range[0] = 0.0;
    source_range[1] = 1.0;
    target_range[0] = 0.0;
    target_range[1] = 1.0;
    multiplicity = m;
  }

  Transition(TsplinePatch* p0, Direction d0, double r0a, double r0b,
             TsplinePatch* p1, Direction d1, double r1a, double r1b,
             unsigned char m)
  {
    if(p0==NULL || p1==NULL)
      throw std::runtime_error("[Transition::Transition] Error, patch not valid");

    source = p0;
    source_dir = d0;
    target = p1;
    target_dir = d1;
    source_range[0] = r0a;
    source_range[1] = r0b;
    target_range[0] = r1a;
    target_range[1] = r1b;
    multiplicity = m;
  }
  virtual ~Transition() {}

  const Direction& GetSourceDir() const { return source_dir; }
  const Direction& GetTargetDir() const { return target_dir; }
  void GetSourceRange(double &a, double &b) const { a=source_range[0]; b=source_range[1]; }
  void GetTargetRange(double &a, double &b) const { a=target_range[0]; b=target_range[1]; }
  const unsigned char& GetMultiplicity() const { return multiplicity; }

  bool InSourceRange(const double &p) const;
  bool InTargetRange(const double &p) const;
  virtual void shoot(const Point2d& p0, unsigned n, std::vector<double> &knots);

  virtual Point2d ConvertParam(const Point2d& p0) const = 0;
  virtual Point2d RotateParamInverse(const Point2d& p0) const = 0;

  friend class TsplinePatch;
  friend class TsplineMultiPatch;
};

/** @brief  north specific functions
 *          overrides functions for coordinate conversion */
class NorthTransition : public Transition
{
public:
  NorthTransition() : Transition() {}
  NorthTransition(TsplinePatch* p0, Direction d0,
                  TsplinePatch* p1, Direction d1,
                  unsigned char m) : Transition(p0, d0, p1, d1, m){}
  NorthTransition(TsplinePatch* p0, Direction d0, double r0a, double r0b,
                  TsplinePatch* p1, Direction d1, double r1a, double r1b,
                  unsigned char m)
    : Transition(p0, d0, r0a, r0b, p1, d1, r1a, r1b, m){}
  virtual Point2d ConvertParam(const Point2d& p0) const;
  virtual Point2d RotateParamInverse(const Point2d& p0) const;
};

/** @brief  east specific functions
 *          overrides functions for coordinate conversion */
class EastTransition : public Transition
{
public:
  EastTransition() : Transition() {}
  EastTransition(TsplinePatch* p0, Direction d0,
                 TsplinePatch* p1, Direction d1,
                 unsigned char m) : Transition(p0, d0, p1, d1, m){}
  EastTransition(TsplinePatch* p0, Direction d0, double r0a, double r0b,
                 TsplinePatch* p1, Direction d1, double r1a, double r1b,
                 unsigned char m)
    : Transition(p0, d0, r0a, r0b, p1, d1, r1a, r1b, m){}
  virtual Point2d ConvertParam(const Point2d& p0) const;
  virtual Point2d RotateParamInverse(const Point2d& p) const;
};

/** @brief  south specific functions
 *          overrides functions for coordinate conversion */
class SouthTransition : public Transition
{
public:
  SouthTransition() : Transition() {}
  SouthTransition(TsplinePatch* p0, Direction d0,
                  TsplinePatch* p1, Direction d1,
                  unsigned char m) : Transition(p0, d0, p1, d1, m){}
  SouthTransition(TsplinePatch* p0, Direction d0, double r0a, double r0b,
                  TsplinePatch* p1, Direction d1, double r1a, double r1b,
                  unsigned char m)
    : Transition(p0, d0, r0a, r0b, p1, d1, r1a, r1b, m){}
  virtual Point2d ConvertParam(const Point2d& p0) const;
  virtual Point2d RotateParamInverse(const Point2d& p) const;
};

/** @brief  west specific functions
 *          overrides functions for coordinate conversion */
class WestTransition : public Transition
{
public:
  WestTransition() : Transition() {}
  WestTransition(TsplinePatch* p0, Direction d0,
                 TsplinePatch* p1, Direction d1,
                 unsigned char m) : Transition(p0, d0, p1, d1, m){}
  WestTransition(TsplinePatch* p0, Direction d0, double r0a, double r0b,
                 TsplinePatch* p1, Direction d1, double r1a, double r1b,
                 unsigned char m)
    : Transition(p0, d0, r0a, r0b, p1, d1, r1a, r1b, m){}
  virtual Point2d ConvertParam(const Point2d& p0) const;
  virtual Point2d RotateParamInverse(const Point2d& p) const;
};

}


#endif
