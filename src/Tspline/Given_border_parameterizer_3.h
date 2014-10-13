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


#ifndef _TSPLINE_CGAL_GIVENBORDERPARAMETERIZER_3_H_
#define _TSPLINE_CGAL_GIVENBORDERPARAMETERIZER_3_H_

// This file implements a border parameterizer when the border is already given
// It implements function handles necessary for CGAL parameterization of submeshes
// embedded on a polyhedron

#include <CGAL/surface_mesh_parameterization_assertions.h>
#include <CGAL/Parameterizer_traits_3.h>

#include <cfloat>
#include <climits>
#include <vector>

namespace CGAL {


template<class ParameterizationMesh_3>
class Given_border_parameterizer_3
{
public:
  /// Export ParameterizationMesh_3 template parameter.
  typedef ParameterizationMesh_3          Adaptor;

private:
  typedef typename Adaptor::NT            NT;
  typedef typename Adaptor::Point_2       Point_2;
  typedef typename Adaptor::Point_3       Point_3;
  typedef typename Adaptor::Vector_2      Vector_2;
  typedef typename Adaptor::Vector_3      Vector_3;
  typedef typename Adaptor::Facet         Facet;
  typedef typename Adaptor::Facet_handle  Facet_handle;
  typedef typename Adaptor::Facet_const_handle  Facet_const_handle;
  typedef typename Adaptor::Facet_iterator Facet_iterator;
  typedef typename Adaptor::Facet_const_iterator  Facet_const_iterator;
  typedef typename Adaptor::Vertex        Vertex;
  typedef typename Adaptor::Vertex_handle Vertex_handle;
  typedef typename Adaptor::Vertex_const_handle  Vertex_const_handle;
  typedef typename Adaptor::Vertex_iterator Vertex_iterator;
  typedef typename Adaptor::Vertex_const_iterator  Vertex_const_iterator;
  typedef typename Adaptor::Border_vertex_iterator  Border_vertex_iterator;
  typedef typename Adaptor::Border_vertex_const_iterator  Border_vertex_const_iterator;
  typedef typename Adaptor::Vertex_around_facet_circulator  Vertex_around_facet_circulator;
  typedef typename Adaptor::Vertex_around_facet_const_circulator  Vertex_around_facet_const_circulator;
  typedef typename Adaptor::Vertex_around_vertex_circulator  Vertex_around_vertex_circulator;
  typedef typename Adaptor::Vertex_around_vertex_const_circulator  Vertex_around_vertex_const_circulator;

  typedef typename std::vector<double>    Offset_map;


public:
  typename Parameterizer_traits_3<Adaptor>::Error_code  parameterize_border(Adaptor& mesh)
  {
    return Parameterizer_traits_3<Adaptor>::OK;
  }

  bool  is_border_convex () { return true; }
};


} //namespace CGAL

#endif //TSPLINE_CGAL_GIVENBORDERPARAMETERIZER_3_H
