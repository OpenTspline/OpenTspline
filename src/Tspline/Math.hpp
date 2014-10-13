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

// This file contains mathematical definitions derifed from CGAL and Eigen
// It implements some helper functions like arithmetic mean, principal component
// analysis (PCA) and comparison of floating point values and vectors

#ifndef _TSPLINE_MATH_H_
#define _TSPLINE_MATH_H_

#include <vector>
#include <stdio.h>
#include <stdexcept>
#include <limits>

#undef Success
#include <Eigen/Eigen>

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_extended_dcel.h>

#define SIZE_T_MAX std::numeric_limits<std::size_t>::max()

namespace tspline
{

// derive algebraic objects from CGAL types
typedef CGAL::Cartesian<double> Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel> Traits_2;
typedef Kernel::Point_2 Point2d;
typedef Kernel::Point_3 Point3d;
typedef Kernel::Vector_2 Vector2d;
typedef Kernel::Vector_3 Vector3d;
typedef Kernel::Ray_3 Ray;
typedef Kernel::Direction_3 Direction3d;
typedef Traits_2::X_monotone_curve_2 Segment2;

// numerical limits
static double epsilon = 10.0 * std::numeric_limits<float>::epsilon(); // float because tgModel using OpenGL is float
static double digits = 1e10;
static double div_digits = 1e-10;

/** @brief extends Point3d by the weight entry for control points  (weight is 1.0 by default) */
class Point4d : public Point3d
{
protected:
  double weight;

public:
  Point4d() : Point3d(), weight(1.0) { }
  Point4d(const Point3d& a, const double& w=1.0) : Point3d(a), weight(w) { }
  Point4d(const double& x, const double& y, const double& z, const double& w) :
    Point3d(x,y,z), weight(w) { }
  Point4d(const double& x, const double& y, const double& z) :
    Point3d(x,y,z), weight(1.0) { }

  void operator=(const Point3d& a)
  {
    *this = Point4d(a);
  }

  double w() const { return weight; }

};

/** @brief aligned Eigen::Vector classes */
typedef std::vector<Eigen::Vector4d, Eigen::aligned_allocator<Eigen::Vector4d> > vector_vec4d;
typedef std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > vector_vec3d;
typedef std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d> > vector_vec2d;

/** @brief compute the arithmetic mean of a vector of Vector3d
 *  @param data input data */
inline Eigen::Vector3d compute_mean (const vector_vec3d &data)
{
  Eigen::Vector3d u (0.0, 0.0, 0.0);

  unsigned s = unsigned (data.size ());
  double ds = 1.0 / s;

  for (unsigned i = 0; i < s; i++)
    u += (data[i] * ds);

  return u;
}

/** @brief compute the principal components of a vector of Vector3d
 *  @param data in: input data
 *  @param mean out: the arithmetic mean of the data
 *  @param eigenvectors out: the eigenvectors sorted descending by their eigenvalues
 *  @param eigenvalues out: descending eigenvalues */
inline void pca (const vector_vec3d &data, Eigen::Vector3d &mean, Eigen::Matrix3d &eigenvectors,
                 Eigen::Vector3d &eigenvalues)
{
  if (data.empty ())
    throw std::runtime_error ("[Math::pca] Error, data is empty\n");

  mean = compute_mean (data);

  unsigned s = unsigned (data.size ());

  Eigen::MatrixXd Q (3, s);

  for (unsigned i = 0; i < s; i++)
    Q.col (i) << (data[i] - mean);

  Eigen::Matrix3d C = Q * Q.transpose ();

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver (C);
  if (eigensolver.info () != Eigen::Success)
    throw std::runtime_error ("[Math::pca] Can not find eigenvalues.\n");

  // reverse ordering of eigenvalues and eigenvectors
  for (int i = 0; i < 3; ++i)
  {
    eigenvalues (i) = eigensolver.eigenvalues () (2 - i);
    if (i == 2)
      eigenvectors.col (2) = eigenvectors.col (0).cross (eigenvectors.col (1));
    else
      eigenvectors.col (i) = eigensolver.eigenvectors ().col (2 - i);
  }
}

// comparison functions for floating point data types

// ==
static bool equal(const double& a, const double& b)
{
  if(std::abs<double>(a-b) < epsilon)
    return true;
  else
    return false;
}

// >
static bool greater(const double &a, const double &b)
{
  if(equal(a,b)) // equal
    return false;
  else if(a < b)
    return false;
  else
    return true;
}

// <
static bool smaller(const double &a, const double &b)
{
  if(equal(a,b)) // equal
    return false;
  else if(a > b)
    return false;
  else
    return true;
}

// >=
static bool gequal(const double &a, const double &b)
{
  if(equal(a,b)) // equal
    return true;
  else if(a > b)
    return true;
  else
    return false;
}

// <=
static bool sequal(const double &a, const double &b)
{
  if(equal(a,b))
    return true;
  else if(a < b)
    return true;
  else
    return false;
}

// ==
static bool equal( const Point2d& a, const Point2d& b)
{
  if(equal(a.x(), b.x()) && equal(a.y(), b.y()))
    return true;
  else
    return false;
}

// ==
static bool equal( const Point3d& a, const Point3d& b)
{
  if(equal(a.x(), b.x()) && equal(a.y(), b.y()) && equal(a.z(), b.z()))
    return true;
  else
    return false;
}

// ==
static bool equal( const Eigen::Vector3d& a, const Eigen::Vector3d& b)
{
  if(equal(a(0), b(0)) && equal(a(1), b(1)) && equal(a(2),b(2)))
    return true;
  else
    return false;
}

/** @brief adjust a double value to a grid of minimal resolution */
static double adjust(const double &a)
{
  return a;
//  return double(float(a));
  double b = round(a * digits) * div_digits;
//  if(equal(b,0.0))
//    b = 0.0;
  return b;
}

/** @brief adjust a Point2D value to a grid of minimal resolution */
static Point2d adjust(const Point2d &p)
{
  return Point2d(adjust(p.x()), adjust(p.y()));
}

/** @brief adjust a Point3D value to a grid of minimal resolution */
static Point3d adjust(const Point3d &p)
{
  return Point3d(adjust(p.x()), adjust(p.y()), adjust(p.z()));
}

/** @brief dot product of CGAL::Vector3d (NOT Eigen)*/
static double dot( const Vector3d& a, const Vector3d& b )
{
  return ( a.x()*b.x() + a.y()*b.y() + a.z()*b.z() );
}

/** @brief L2 norm of CGAL::Vector3d (NOT Eigen) */
static double norm( const Vector3d& a )
{
  return sqrt(a.squared_length());
}

/** @brief angle between two CGAL::Vector3d (NOT Eigen) */
static double angle( const Vector3d& a, const Vector3d& b )
{
  return acos( dot(a,b) / (norm(a) * norm(b)) );
}

}

#endif
