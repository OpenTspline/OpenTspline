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

#ifndef _TSPLINE_TSPLINE_TYPES_H_
#define _TSPLINE_TSPLINE_TYPES_H_

// This file defines the basic objects for a T-spline like control-points (CP),
// edges (halfedes) and faces, refinable basis and blending functions and more
// efficient structures for basis function values and their derivatives.

#include "Math.hpp"
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Arr_naive_point_location.h>

namespace tspline {

/** @brief additional data for control point
 *         physical position, parametric position, id, parametric knot vectors, linkage */
class TVertex
{
protected:
  Point4d cp; /// control point position: 3D coordinates + weight (1 by default)
  std::vector<TVertex*> cp_linked; /// linked vertices

public:
  int id;                 /// unique id
  int patch_id;           /// patch id, this vertex belongs to
  Point2d param;          /// parametric position of this vertex
  std::vector<double> s;  /// knot vector
  std::vector<double> t;  /// knot vector
  bool updated;           /// flag used during T-spline parameter update
  bool is_primary;        /// indicates if this is the primary vertex when linked (i.e. the vertex with the smallest id)

  TVertex();

  void SetCP(const Point4d& p);           /// set control point and weight
  void SetCP(const Point3d &p);           /// set control point (weight = 1)
  void SetCP(const Eigen::Vector4d& p);   /// set control point and weight
  void SetCP(const Eigen::Vector3d& p);   /// set control point (weight = 1)
  const Point4d &GetCP() const;           /// get control point and weight
  void GetCP(Eigen::Vector4d &p) const;   /// get control point and weight
  void GetCP(Eigen::Vector3d &p) const;   /// get control point

  void LinkCP(TVertex *vext);             /// link CP, re-link vext to this CP, compute dominant CP
  void UnLink();                          /// remove linkage, re-compute dominant CP if necessary
  bool IsLinked() const;                  /// return !cp_linked.empty()
  bool IsLinked(int pid, int vid) const;  /// check if this CP is linked to a specific CP

  bool IsIn(std::vector<const TVertex*>& links) const;      /// helper function checking if this CP is within a list of CPs

  bool IsLinked(std::vector<const TVertex*>& links) const;  /// helper function checking if this CP is linked to any of the listed CPs
  size_t GetNumLinks() const;                               /// return cp_linked.size();
  void GetLinks(std::vector<int> &_patch_ids) const;        /// return patches this CP is linked with
  void GetLinks(std::vector<int> &_patch_ids, std::vector<int> &_vertex_ids) const; /// return patches with CPs this CP is linked with
  int GetPrimaryID() const; /// get id of primary CP this CP is linked to
  void GetPrimaryID(int &_patch_id, int &_vertex_id) const; /// get id (patch and CP) of primary CP this CP is linked to

  void PrintKnotVectors() const;

protected:
  TVertex* ComputeDominant(int pid, int vid); /// computes the vertex with the lowest pid,vid and sets it dominant
};

/** @brief additional data for control point
 *         parametric edge length
*/
struct THalfedge
{
  double d;     /// parametric distance of halfedge
  THalfedge() :
    d(0.0)
  {
  }
};
struct TFace
{
  double error_sqr;     /// squared error of face:
                        ///     1) point fitting: sum of squared euclidean distance
                        ///     2) photometric: squared difference of intensities or NCC of the face-patch
  unsigned point_count; /// number of points (fitting) or pixels in the reference view (photometric)
  TFace() :
    error_sqr(0.0), point_count(0)
  {
  }
};

typedef CGAL::Arr_extended_dcel<Traits_2, TVertex, THalfedge, TFace> Dcel;
typedef CGAL::Arrangement_2<Traits_2, Dcel> Arrangement_2;
typedef CGAL::Arr_naive_point_location<Arrangement_2> NaivePointLocation;
typedef std::pair<Arrangement_2::Vertex_iterator,Arrangement_2::Vertex_iterator> VertexPair;

/** @brief storage for basis functions and derivatives (boosts efficiency when T-mesh is not changing) */
struct CPbasis
{
  Arrangement_2::Vertex_const_iterator vit; /// respective control point
  double b;                                 /// basis function value
  double bs;                                /// derivative of basis function value in s-direction
  double bt;                                /// derivative of basis function value in t-direction
  CPbasis() :
    b(0.0), bs(0.0), bt(0.0)
  {
  }
};

/** @brief sum of basis functions and derivatives */
struct BasisSum
{
  Eigen::Vector3d cp, ts, tt;
  double B, Bs, Bt;

  BasisSum() : cp(0.0, 0.0, 0.0),
    ts(0.0, 0.0, 0.0),
    tt(0.0, 0.0, 0.0),
    B(0.0), Bs(0.0), Bt(0.0)
  {}
};

/** @brief relates domains of multiple T-spline patches to each other
 *         multi-patch transitions are defined by the side of the borders the domains connect (i.e. their directions)
*/
enum Direction
{
  NORTH = 0,
  EAST = 1,
  SOUTH = 2,
  WEST = 3
};

// ******************************************************
// The following structures are for congruent refinement
// [1] Sederberg et al. "T-spline Simplification and Local Refinement", Transactions on Graphics, 2004

/** @brief  Scale (i.e. weight) of blending functions during knot insertion*/
struct Scale
{
  double c;
  int i;
  int j;
};

class Tspline;
class BlendingFunction;

/** @brief  enables recursively collecting blending functions that must be split during refinement
 *          BasisFunction are top level in this tree and are copied from the T-spline mesh
 * */
class BasisFunction
{
private:
  BasisFunction(const BasisFunction& b) {}    // disable copies
  void operator=(const BasisFunction& b) {}   // disable copies

protected:
  Tspline* tsp;
  bool is_split;
  std::vector<double> s;
  std::vector<double> t;

  BlendingFunction* B[2];

  Arrangement_2::Vertex_iterator vit;

  BasisFunction(Tspline* _tsp)
    : tsp(_tsp), is_split(false)
  {
    B[0] = NULL;
    B[1] = NULL;
  }

public:
  BasisFunction(Tspline* _tsp,
                const std::vector<double>& _s, const std::vector<double>& _t,
                Arrangement_2::Vertex_iterator _vit)
    : tsp(_tsp), is_split(false)
  {
    s = _s;
    t = _t;
    vit = _vit;
    B[0] = NULL;
    B[1] = NULL;
  }

  virtual ~BasisFunction();

  int id();
  bool split() { return is_split; }

  /** @brief  splits a blending function according to [1] */
  static void split(const std::vector<double> &N, const double &k,
                    std::vector<double> &N0, std::vector<double> &N1,
                    double& c0, double& c1);

  /** @brief  splits a blending at parametric position k */
  bool split_s(const double& k);

  /** @brief  splits a blending at parametric position k */
  bool split_t(const double& k);

  /** @brief  recursively computes violations of type 1 (see [1]) of splitting tree */
  virtual bool compute_violation_1();

  /** @brief  recursively computes violations of type 2 (see [1]) of splitting tree */
  virtual bool compute_violation_2();

  /** @brief  return type of this class (basis function as opposed to blending function) */
  virtual bool is_basis_function() { return true; }

  /** @brief  recursively computes scale and returns it */
  void get_scale(std::vector<Scale>& g_scale);

  friend class BlendingFunction;

};

Scale operator*(const Scale& s1, const Scale& s2);
std::vector<BasisFunction*>::iterator find(std::vector<BasisFunction *> &f, int id);

/** @brief enables recursively collecting blending functions that must be split during refinement */
class BlendingFunction : public BasisFunction
{
protected:
  Scale scale;
  BasisFunction* parent;

  /** @brief  recursively computes scale and returns it */
  void get_scale(std::vector<Scale>& g_scale, Scale parent_scale);

public:
  BlendingFunction(Tspline* _tsp, BasisFunction* _parent)
    : BasisFunction(_tsp), parent(_parent)
  { }

  BlendingFunction(Tspline* _tsp, BasisFunction* _parent,
                   const std::vector<double>& _s, const std::vector<double>& _t,
                   Arrangement_2::Vertex_iterator _vit)
    : BasisFunction(_tsp, s, t, vit), parent(_parent)
  { }
  virtual ~BlendingFunction() {}

  /** @brief  recursively computes violations of type 1 (see [1]) of splitting tree */
  virtual bool compute_violation_1();

  /** @brief  recursively computes violations of type 2 (see [1]) of splitting tree */
  virtual bool compute_violation_2();

  /** @brief  return type of this class (blending function as opposed to basis function) */
  virtual bool is_basis_function() { return false; }

  /** @brief  recursively computes scale and returns it */
  void get_scale(std::vector<Scale>& g_scale);

  friend class BasisFunction;

};

} // namespace tspline

#endif
