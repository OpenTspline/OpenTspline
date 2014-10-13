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

#ifndef _TSPLINE_FILE_H_
#define _TSPLINE_FILE_H_

// This file implement the file IO for T-splines and multi-patch T-splines

#include "Tspline.h"
#include "TsplineMultiPatch.h"

namespace tspline{

/** @brief  T-spline file IO
 *          for saving and loading T-spline files.
 *          naming convention: .tsp for T-spline files
 *                             .tsmp/.tmp for T-spline multi-patch files (possibly containing extraordinary vertices)
 */
class File
{
private:
  /** @brief helper structure for storing linkage information */
  struct Link
  {
    size_t patch_id;
    int vertex_id;
    std::vector<int> patch_ids;
    std::vector<int> vertex_ids;
  };

  /** @brief helper structure for storing transition information */
  struct Trans
  {
    size_t source_patch_id;
    Direction source_dir;
    size_t target_patch_id;
    Direction target_dir;
    unsigned char multiplicity;
    double source_range[2];
    double target_range[2];
    Trans(){}
    Trans(size_t _source, Direction _s_dir, double r0a, double r0b,
          size_t _target, Direction _t_dir, double r1a, double r1b,
          unsigned char _m) :
      source_patch_id(_source), source_dir(_s_dir),
      target_patch_id(_target), target_dir(_t_dir),
      multiplicity(_m)
    {
      source_range[0] = r0a;
      source_range[1] = r0b;
      target_range[0] = r1a;
      target_range[1] = r1b;
    }
  };

  template<typename T>
  static size_t read_vector (std::vector<T> &vec, FILE* pFile);
  static size_t read_string(std::string &str, FILE* pFile);
  static size_t read_point(Point2d &p, FILE* pFile);
  static size_t read_point(Point3d &p, FILE* pFile);
  static size_t read_point(Point4d &p, FILE* pFile);

  static size_t read_vertex(Tspline &tsp, FILE* pFile, Link &link);
  static size_t read_edge(Tspline &tsp, FILE* pFile);
  static size_t read_tspline(Tspline &tsp, FILE *pFile, std::vector<Link> &links);

  static size_t read_tspline_patch(TsplinePatch &tsp, FILE* pFile);
  static size_t read_tspline_multipatch(TsplineMultiPatch &tsp, FILE* pFile);


  template<typename T>
  static void write_vector(const std::vector<T> &vec, FILE* pFile);
  static void write_string(const std::string& str, FILE* pFile);
  static void write_point(const Point2d &p, FILE *pFile);
  static void write_point(const Point3d &p, FILE *pFile);
  static void write_point(const Point4d &p, FILE *pFile);

  static void write_vertex(Tspline::Vertex_const_iterator &vit, FILE* pFile);
  static void write_edge(Tspline::Edge_const_iterator &eit, FILE* pFile);
  static void write_tspline(const tspline::Tspline &tsp, FILE* pFile);

  static void write_tspline_patch(const TsplinePatch &tsp, FILE* pFile, std::vector<Trans> &transitions);
  static void write_tspline_multipatch(TsplineMultiPatch &tsmp, FILE* pFile);

public:

  /** @brief save a T-spline to a file */
  static bool Save(const tspline::Tspline &tsp, const std::string& filename);

  /** @brief load a T-spline from a file */
  static bool Load(tspline::Tspline &tsp, const std::string& filename);

  /** @brief save a multi-patch T-spline to a file */
  static bool Save(TsplineMultiPatch &tsmp, const std::string& filename);

  /** @brief load a multi-patch T-spline from a file */
  static bool Load(tspline::TsplineMultiPatch &tsp, const std::string& filename);

};

}

#endif
