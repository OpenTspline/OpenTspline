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

#ifndef _TSPLINE_TSPLINE_H_
#define _TSPLINE_TSPLINE_H_

#include "Types.h"
#include "Cox.h"

namespace tspline {

/** @brief  T-spline main class
 *          implements all necessary T-spline functions
 *          TODO: clean up, access
 */
class Tspline : public Arrangement_2
{
private:
  void update_params(Vertex_iterator &v);

public:
  bool quiet;

  unsigned degree;
  bool clamped;
  Point2d param_min;
  Point2d param_max;
  Point2d point_min;
  Point2d point_max;

  std::string texture;

  Tspline() :
    quiet(true), degree(3), clamped(true)
  {
  }

  void update_params();
  void update_knot_vectors(); // applies Rule 1 (Sederberg 2003)
  void update_vertex_ids();

  //********************
  // refinement methods

  /** @brief split faces which exceed an average error threshold and cover more than 'min_points' (TFace.error_square / TFace.point_count) */
  unsigned refine_average_error(double avg_error, unsigned min_points);

  /** @brief  split a face
   *          assumes that only rectangular faces exist in grid
   *          insert edges if it violates rule 1
   *          update knot vectors of T-spline (exhaustive T-mesh traversal) */
  void refine(Face_handle f, bool insert_edges=true, bool update_knots=true);

  /** @brief  split a face without changing the physical geometry
   *          assumes that only rectangular faces exist in grid */
  void refine_congruent(Face_handle f);

  /** @brief  split all faces that cover more than 'max_num_points' (TFace.point_count) */
  unsigned refine_point_count(unsigned max_num_points=8);
  /** @brief  split a faces if it covers more than 'max_num_points' (TFace.point_count) */
  bool refine_point_count(Face_handle f, unsigned max_num_points=8);

  /** @brief  split all faces with a chord length (in physical domain) greater the threshold */
  unsigned refine_chord(double chord_length_threshold);
  /** @brief  split a face if chord length (in physical domain) exceeds the threshold */
  bool refine_chord(Face_handle f, double chord_length_threshold);

  /** @brief  insert vertex, insert missing edges and update knot vectors (exhaustive T-mesh traversal) */
  Vertex_iterator insert_vertex(const double &s, const double &t,
                                bool insert_edges=true, bool update_knots=true);

  /** @brief  split a face horizontal, inserting missing edges and updating knot vectors (not geometry-preserving) */
  Halfedge_iterator split_horizontal(Face_iterator f);
  /** @brief  split a face vertical, inserting missing edges and updating knot vectors (not geometry-preserving) */
  Halfedge_iterator split_vertical(Face_iterator f);

  /** @brief  split face horizontal, insert missing edges, make congruent and update knot vectors (geometry-preserving) */
  Halfedge_iterator split_horizontal_congruent(Face_iterator &f, double y);
  /** @brief  split face vertical, insert missing edges, make congruent and update knot vectors (geometry-preserving)  */
  Halfedge_iterator split_vertical_congruent(Face_iterator &f, double x);
  /** @brief  split edge, insert missing edges, make congruent and update knot vectors (geometry-preserving)  */
  Vertex_iterator split_edge_congruent(Halfedge_iterator &h, Point2d param);

  //*****************
  // removal methods
  void remove_vertex(Vertex_handle vh,
                     bool remove_isolated=true,
                     bool insert_edges=true,
                     bool update_knots = true,
                     bool remove_I_junctions=true);
  size_t remove_redundant_vertices(bool Ijunctions=true, bool loose=true, bool isolated=true);

  void reduce_vertices(double min_cage_curvature=0.04);
  bool reduce_vertices(double fitting_threshold, double sampling_resolution);

  //********************
  // evaluation methods
  /** @brief evaluate basis function at parameter position (s,t) */
  virtual void evaluate_basis(const double &s, const double &t, const TVertex &vext, double &B) const;
  /** @brief evaluate basis function and derivatives at parameter position (s,t) */
  virtual void evaluate_basis(const double &s, const double &t, const TVertex &vext, double &B, double &Bs, double &Bt) const;

  /** @brief evaluate basis function at parameter position (s,t) */
  virtual Eigen::Vector3d evaluate(const double &s, const double &t) const;
  /** @brief evaluate basis function and derivatives at parameter position (s,t) */
  virtual Eigen::Vector3d evaluate(const double &s, const double &t, Eigen::Vector3d &ts, Eigen::Vector3d &tt) const;
  /** @brief evaluate basis function at parameter position (s,t) and stores basis function values in cpbasis */
  virtual Eigen::Vector3d evaluate(const double &s, const double &t,
                                   std::vector<tspline::CPbasis> &cpbasis) const;
  /** @brief evaluate basis function and derivatives at parameter position (s,t) and stores basis function values in cpbasis */
  virtual Eigen::Vector3d evaluate(const double &s, const double &t,
                                   Eigen::Vector3d &ts, Eigen::Vector3d &tt,
                                   std::vector<tspline::CPbasis> &cpbasis) const;

  /** @brief find closest point on T-spline surface (hint must be sufficiently close to result)
   *  @param point in: the reference point
   *  @param hint  in: parametric coordinates close to the closest point on the T-spline surface
   *  @param iterations in: number of maximum iterations
   *  @param accuracy in: termination criterion
   *  @return closest point */
  Eigen::Vector2d closest_point(const Eigen::Vector3d &point, const Eigen::Vector2d &hint, unsigned iterations = 200,
                                double accuracy = 1e-4) const;
  /** @brief find closest point on T-spline surface returning derivatives at closest point */
  Eigen::Vector2d closest_point(const Eigen::Vector3d &point, const Eigen::Vector2d &hint, Eigen::Vector3d &p,
                                Eigen::Vector3d &ts, Eigen::Vector3d &tt, unsigned iterations, double accuracy) const;

  /** @brief approximate normal of a control point using the T-mesh geometry in euclidean space */
  vector_vec3d compute_cp_normals_by_cage();
  /** @brief approximate the normal of a control point using the footpoint of a control point */
  vector_vec3d compute_cp_normals_by_footpoints();

  // inline functions
  inline void set_quiet(bool v) { quiet = v; }

  // ##############################################################
  // operators (implemented in Tspline_operators.cpp)

  // T-spline operations
  CGAL::Object locate_param(const double &s, const double &t);
  CGAL::Object locate_param(Point2d param);
  CGAL::Object locate_point(const Point2d &p);
  Eigen::Vector3d center() const; // arithmetic mean of control points
  bool is_rational() const; // all weigths = 1
  bool is_NURBS() const; // regular control mesh
  virtual size_t number_of_controlpoints() const; // number of control points that are primary
  void get_bounding_box(Eigen::Vector3d& bb_min, Eigen::Vector3d& bb_max) const; // of control points

  void get_NURBS_dimension(unsigned& width, unsigned& height); // regular control grid resolution
  Vertex_iterator get_NURBS_vertex(unsigned col, unsigned row); // regular control grid indexing

  // Face operations
  Face_iterator locate_face(const double &s, const double &t); // param domain
  Eigen::Vector2d face_center(Face_const_handle fit) const; // param domain
  Eigen::Vector2d face_center(Face_const_handle fit, Eigen::Vector3d &point) const; // param domain
  void face_bounding_box(Face_const_handle f, Point2d &bbmin, Point2d &bbmax) const; // param domain
  void face_bounding_box(Face_const_handle f, Eigen::Vector4d &bb) const; // param domain
  void face_bounding_box_by_point(Face_const_handle f, Point2d &bbmin, Point2d &bbmax) const; // CGAL arr domain (point)
  void face_bounding_box_by_point(Face_const_handle f, Eigen::Vector4d &bb) const; // CGAL arr domain (point)
  bool face_is_collapsed(Face_const_handle f) const; // param domain
  bool face_is_corner(Face_const_handle f) const; // CGAL arr domain (point)
  bool face_is_bottom_left_corner(Face_const_handle f) const; // CGAL arr domain (point)
  void face_dimensions_squared(Face_handle f, double& ds, double& dt, const Eigen::Vector4d &bb); // physical domain
  void face_dimensions_squared(Face_handle f, double& ds, double& dt); // physical domain

  Face_iterator face(Face_const_iterator fc);

  // Vertex operations
  Vertex_iterator get_vertex(const int &id);
  bool is_boundary(Vertex_const_iterator v) const;
  bool is_boundary(Vertex_const_iterator v, Halfedge_const_iterator &he) const;
  bool is_corner(Vertex_const_iterator v) const;
  bool is_Xjunction(Vertex_const_iterator v) const;
  bool is_Tjunction(Vertex_const_iterator v) const;
  bool is_Ljunction(Vertex_const_iterator v) const;
  bool is_Ijunction(Vertex_const_iterator v) const;
  bool is_bottom_left_vertex(Vertex_const_iterator v) const;


  Vertex_iterator get_closest_vertex_by_point(const Point2d &point, double &d_sqr); // CGAL arr domain (point)
  Halfedge_iterator get_closest_edge_by_point(const Point2d &point, double &d_sqr); // CGAL arr domain (point)
  Face_const_iterator get_face_by_point(const Point2d &point); // CGAL arr domain (point)

  Vertex_iterator get_closest_vertex(const Point_2 &param); // param domain

  Eigen::Vector4d grid_normal_curvature(Vertex_iterator v) const; // param domain

  double edge_length_by_point(Halfedge_const_iterator h); // CGAL arr domain (point)
  Halfedge_handle get_shortest_edge(); // param domain

public:
  Halfedge_iterator get_left_halfedge_at_param(Face_iterator &f, const double &y, double &x);
  Halfedge_iterator get_right_halfedge_at_param(Face_iterator &f, const double &y, double &x);
  Halfedge_iterator get_top_halfedge_at_param(Face_iterator &f, const double &x, double &y);
  Halfedge_iterator get_bottom_halfedge_at_param(Face_iterator &f, const double &x, double &y);

  Halfedge_iterator get_left_halfedge_at_point(Face_iterator f, const Point2d &p);
  Halfedge_iterator get_right_halfedge_at_point(Face_iterator f, const Point2d &p);
  Halfedge_iterator get_top_halfedge_at_point(Face_iterator f, const Point2d &p);
  Halfedge_iterator get_bottom_halfedge_at_point(Face_iterator f, const Point2d &p);

  virtual void shoot_left(Vertex_iterator &v, const unsigned &n, std::vector<double> &knots);
  virtual void shoot_right(Vertex_iterator &v, const unsigned &n, std::vector<double> &knots);
  virtual void shoot_up(Vertex_iterator &v, const unsigned &n, std::vector<double> &knots);
  virtual void shoot_down(Vertex_iterator &v, const unsigned &n, std::vector<double> &knots);

  void shoot_left(Halfedge_iterator& he, const double& y_curr, const unsigned &n, std::vector<double> &knots);
  void shoot_right(Halfedge_iterator& he, const double& y_curr, const unsigned &n, std::vector<double> &knots);
  void shoot_up(Halfedge_iterator& he, const double& x_curr, const unsigned &n, std::vector<double> &knots);
  void shoot_down(Halfedge_iterator& he, const double& x_curr, const unsigned &n, std::vector<double> &knots);

  Halfedge_iterator get_right_halfedge(Vertex_iterator &v);
  Halfedge_iterator get_left_halfedge(Vertex_iterator &v);
  Halfedge_iterator get_top_halfedge(Vertex_iterator &v);
  Halfedge_iterator get_bottom_halfedge(Vertex_iterator &v);

  //  std::vector<Halfedge_iterator> get_right_halfedges(Face_iterator &f);
  //  std::vector<Halfedge_iterator> get_left_halfedges(Face_iterator &f);
  //  std::vector<Halfedge_iterator> get_top_halfedges(Face_iterator &f);
  //  std::vector<Halfedge_iterator> get_bottom_halfedges(Face_iterator &f);

  void get_halfedges(Face_iterator &f, Point_2 &param, Halfedge_iterator &right, Halfedge_iterator &left,
                     Halfedge_iterator &top, Halfedge_iterator &bottom);

  bool is_connectable(Vertex_iterator v0, Vertex_iterator v1);

public:

  /** @brief inserts missing edges of a face (when two opposing points are not connected)
   *  @brief L junctions are priorized */
  Halfedge_iterator insert_missing_edges(Face_iterator fit);

  /** @brief inserts missing edges (when two opposing points are not connected)
   *  @brief L junctions are priorized
   *  @brief applies Rule 2 (Sederberg 2003) */
  void insert_missing_edges();

  /** @brief resolves L junctions by inserting edges
   *  @brief possibly splitting opposing edges, introducing new vertices, if neccessary */
  void insert_edges_at_L_junctions();

  void compute_knot_vectors(Vertex_iterator vit);
  void compute_knot_vectors(Vertex_iterator vit,
                            std::vector<double>& s, std::vector<double>& t);

protected:
  bool compute_violation_1(std::vector<BasisFunction*> &blend, int num_vertices_prev);
  bool compute_violation_2(std::vector<BasisFunction*> &blend);
  void scale_2_matrix(const std::vector<Scale>& scales, Eigen::MatrixXd& M, int num_vertices_prev);
  void update_control_points(const std::vector<Scale>& scales, int num_vertices_prev, Eigen::MatrixXd& M);
  void make_congruent(int num_vertices_prev);
  void make_congruent(int num_vertices_prev, Eigen::MatrixXd& M);

  virtual Arrangement_2::Vertex_iterator insert_vertex(Halfedge_iterator e, Point2d param);
  Arrangement_2::Vertex_iterator insert_vertex(Face_iterator f, Point2d param);

  Arrangement_2::Halfedge_iterator split_horizontal(Face_iterator f, double y);
  Arrangement_2::Halfedge_iterator split_vertical(Face_iterator f, double x);

  friend class BasisFunction;
  friend class BlendingFunction;

};

} // namespace tspline

#endif
