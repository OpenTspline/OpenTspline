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

#include "TsplineCreator.h"

using namespace tspline;


//void TsplineCreator::CreateUniformPlaneXY(tspline::Tspline &tsp,
//                                          double x0, double y0, double z0,
//                                          double width, double height,
//                                          unsigned segX, unsigned segY)
//{
//  double dsegX(width / segX);
//  double dsegY(height / segY);
//  double x1 = x0+width;
//  double y1 = y0+height;

//  // init vertices
//  for(unsigned j=0; j<segY+1; j++)
//  {
//    double y = y0+j*dsegY;
//    for(unsigned i=0; i<segX+1; i++)
//    {
//      double x = x0+i*dsegX;
//      CGAL::insert_point(tsp, Point2d( x, y ));
//    }
//  }

//  // init grid
//  for(unsigned i=0; i<segX+1; i++)
//  {
//    double x = x0 + i * dsegX;
//    CGAL::insert( tsp, Segment2(Point2d (x, y0), Point2d (x, y1)) );
//  }
//  for(unsigned j=0; j<segY+1; j++)
//  {
//    double y = y0 + j * dsegY;
//    CGAL::insert( tsp, Segment2(Point2d (x0, y), Point2d (x1, y)) );
//  }

//  //  double stepX(width / segX);
//  //  double stepY(height / segY);
//  //  double epsX = stepX * 0.1;
//  //  double epsY = stepY * 0.1;
//  //  double x1 = x0+width;
//  //  double y1 = y0+height;

//  //  for(double j=y0; j<y1+epsY; j+=stepY)
//  //    for(double i=x0; i<x1+epsX; i+=stepX)
//  //      CGAL::insert_point(tsp, Point2d(i,j));

//  //  // init grid
//  //  for(double i=x0; i<x1+epsX; i+=stepX)
//  //    CGAL::insert( tsp, Segment2(Point2d (i, y0), Point2d (i, y1)) );

//  //  for(double j=y0; j<y1+epsY; j+=stepY)
//  //    CGAL::insert( tsp, Segment2(Point2d (x0, j), Point2d (x1, j)) );

//  // init edge distances
//  tspline::Tspline::Halfedge_iterator eit;
//  for (eit = tsp.halfedges_begin (); eit != tsp.halfedges_end (); eit++)
//  {
//    tspline::Tspline::Vertex_iterator v1 = eit->source();
//    tspline::Tspline::Vertex_iterator v2 = eit->target();
//    const Point2d& p1 = v1->point();
//    const Point2d& p2 = v2->point();

//    tspline::THalfedge& ext = eit->data();
//    ext.d = std::sqrt( (p1 - p2).squared_length());

//    if(equal(p1.y(),p2.y())) // horizontal
//    {
//      if(tsp.get_left_halfedge(v1) ==tsp.halfedges_end() ||
//         tsp.get_left_halfedge(v2) ==tsp.halfedges_end() ||
//         tsp.get_right_halfedge(v1)==tsp.halfedges_end() ||
//         tsp.get_right_halfedge(v2)==tsp.halfedges_end() )
//      {
//        ext.d = 0.0;
//      }else{
//        ext.d = ext.d / (segX-2) * segX;
//      }
//    }
//    if(equal(p1.x(),p2.x())) // vertical
//    {
//      if(tsp.get_top_halfedge(v1) ==tsp.halfedges_end() ||
//         tsp.get_top_halfedge(v2) ==tsp.halfedges_end() ||
//         tsp.get_bottom_halfedge(v1)==tsp.halfedges_end() ||
//         tsp.get_bottom_halfedge(v2)==tsp.halfedges_end() )
//      {
//        ext.d = 0.0;
//      }else{
//        ext.d = ext.d / (segY-2) * segY;
//      }
//    }

//    ext.d = adjust(ext.d);
//    //    for(size_t d=0; d<tsp.degree-2; d++)
//    //    {

//    //      size_t d1 = d+1;
//    //      if ( (p1.y() == d && p2.y() == d1) || (p1.y() == d1 && p2.y() == d))
//    //        e.d = 0;
//    //      if ( (p1.x() == d && p2.x() == d1) || (p1.x() == d1 && p2.x() == d))
//    //        e.d = 0;
//    //      if ( (p1.y() == dim-d1 && p2.y() == dim-d) || (p1.y() == dim-d && p2.y() == dim-d1))
//    //        e.d = 0;
//    //      if ( (p1.x() == dim-d1 && p2.x() == dim-d) || (p1.x() == dim-d && p2.x() == dim-d1))
//    //        e.d = 0;
//    //    }
//  }
//  tsp.clamped = false;

//  // init control points
//  int vid(0);
//  tspline::Tspline::Vertex_iterator vit;
//  for (vit = tsp.vertices_begin(); vit != tsp.vertices_end(); vit++)
//  {
//    Point2d& p = vit->point();
//    p = adjust(p);

//    tspline::TVertex vext = vit->data();
//    vext.SetCP(Point3d(p.x(), p.y(), z0));
//    vext.id = vid++;
//    vit->set_data(vext);
//  }


//  // update parameter space and knot vectors
//  tsp.update_params ();
//  tsp.update_knot_vectors ();
//}

void TsplineCreator::CreateClampedPlaneXY(tspline::Tspline &tsp,
                                          double x0, double y0, double z0,
                                          double width, double height,
                                          unsigned segX, unsigned segY)
{
  double dsegX(width / segX);
  double dsegY(height / segY);
  double dgridX(width/(segX-2));
  double dgridY(height/(segY-2));

  // init vertices
  int vid(0);
  for(unsigned j=0; j<segY+1; j++)
  {
    double ycp = y0+j*dsegY;
    double y = y0 + j * dgridY;
    for(unsigned i=0; i<segX+1; i++)
    {
      double xcp = x0+i*dsegX;
      double x = x0 + i * dgridX;
      Tspline::Vertex_handle v = CGAL::insert_point(tsp, Point2d( x, y ));
      v->data().id = vid++;
      v->data().SetCP(Point3d(xcp,ycp,z0));
    }
  }

  // init grid
  for(unsigned i=0; i<segX+1; i++)
  {
    double x = x0 + i * dgridX;
    CGAL::insert( tsp, Segment2(Point2d (x, y0), Point2d (x, y0+height+2*dgridY)) );
  }
  for(unsigned j=0; j<segY+1; j++)
  {
    double y = y0 + j * dgridY;
    CGAL::insert( tsp, Segment2(Point2d (x0, y), Point2d (x0+width+2*dgridX, y)) );
  }

  // init edge distances
  tspline::Tspline::Halfedge_iterator eit;
  for (eit = tsp.halfedges_begin(); eit != tsp.halfedges_end(); eit++)
  {
    tspline::Tspline::Vertex_iterator v1 = eit->source();
    tspline::Tspline::Vertex_iterator v2 = eit->target();
    const Point2d& p1 = v1->point();
    const Point2d& p2 = v2->point();

    if(equal(p1.y(),p2.y())) // horizontal
    {
      if(tsp.get_left_halfedge(v1)==tsp.halfedges_end() || // left border
         tsp.get_left_halfedge(v2)==tsp.halfedges_end() ||
         tsp.get_right_halfedge(v1)==tsp.halfedges_end()||  // right border
         tsp.get_right_halfedge(v2)==tsp.halfedges_end())
        eit->data().d = 0.0;
      else
        eit->data().d = dgridX;
    }
    else if(equal(p1.x(),p2.x())) // vertical
    {
      if(tsp.get_bottom_halfedge(v1)==tsp.halfedges_end() || // bottom border
         tsp.get_bottom_halfedge(v2)==tsp.halfedges_end() ||
         tsp.get_top_halfedge(v1)==tsp.halfedges_end() ||    // top border
         tsp.get_top_halfedge(v2)==tsp.halfedges_end())
        eit->data().d = 0.0;
      else
        eit->data().d = dgridY;
    }
    else
      throw std::runtime_error("[TsplineCreator::CreatePlaneXY] Error, invalid edge.");
  }

  // update parameter space of tsp
  tsp.update_params();
  tsp.update_knot_vectors();
  tsp.clamped = false;
}

void TsplineCreator::CreatePlaneXY(tspline::Tspline &tsp,
                                   double x0, double y0, double z0,
                                   double width, double height,
                                   unsigned segX, unsigned segY)
{
  double dsegX(width / segX);
  double dsegY(height / segY);
  double x1 = x0+width;
  double y1 = y0+height;
  // init vertices
  int vid(0);
  for(unsigned j=0; j<segY+1; j++)
  {
    double y = y0+j*dsegY;
    for(unsigned i=0; i<segX+1; i++)
    {
      double x = x0+i*dsegX;
      Tspline::Vertex_handle v = CGAL::insert_point(tsp, Point2d( x, y ));
      v->data().id = vid++;
      v->data().SetCP(Point3d(x,y,z0));
    }
  }
  // init grid
  for(unsigned i=0; i<segX+1; i++)
  {
    double x = x0 + i * dsegX;
    CGAL::insert( tsp, Segment2(Point2d (x, y0), Point2d (x, y1)) );
  }
  for(unsigned j=0; j<segY+1; j++)
  {
    double y = y0 + j * dsegY;
    CGAL::insert( tsp, Segment2(Point2d (x0, y), Point2d (x1, y)) );
  }
  // init edge distances
  tspline::Tspline::Halfedge_iterator eit;
  for (eit = tsp.halfedges_begin(); eit != tsp.halfedges_end(); eit++)
  {
    const Point2d& p1 = eit->source()->point();
    const Point2d& p2 = eit->target()->point();
    eit->data().d = std::sqrt((p1 - p2).squared_length());
  }
  // update parameter space of tsp
  tsp.update_params();
  tsp.update_knot_vectors();
}

void TsplineCreator::CreatePlaneYZ(tspline::Tspline &tsp,
                                   double x0, double y0, double z0,
                                   double width, double height,
                                   unsigned segY, unsigned segZ)
{
  double dsegY(width / segY);
  double dsegZ(height / segZ);
  double y1 = y0+width;
  double z1 = z0+width;

  // init vertices
  int vid(0);
  for(unsigned j=0; j<segZ+1; j++)
  {
    double z = z0+j*dsegZ;
    for(unsigned i=0; i<segY+1; i++)
    {
      double y = y0+i*dsegY;
      Tspline::Vertex_handle v = CGAL::insert_point(tsp, Point2d( y, z ));
      v->data().id = vid++;
      v->data().SetCP(Point3d(x0, y, z));
    }
  }

  // init grid
  for(unsigned i=0; i<segY+1; i++)
  {
    double y = y0 + i * dsegY;
    CGAL::insert( tsp, Segment2(Point2d (y, z0), Point2d (y, z1)) );
  }
  for(unsigned j=0; j<segZ+1; j++)
  {
    double z = z0 + j * dsegZ;
    CGAL::insert( tsp, Segment2(Point2d (y0, z), Point2d (y1, z)) );
  }

  // init edge distances
  tspline::Tspline::Halfedge_iterator eit;
  for (eit = tsp.halfedges_begin(); eit != tsp.halfedges_end(); eit++)
  {
    const Point2d& p1 = eit->source()->point();
    const Point2d& p2 = eit->target()->point();
    eit->data().d = std::sqrt((p1 - p2).squared_length());
  }

  // update parameter space of tsp
  tsp.update_params();
  tsp.update_knot_vectors();
}

void TsplineCreator::CreatePlaneXZ(tspline::Tspline &tsp,
                                   double x0, double y0, double z0,
                                   double width, double height,
                                   unsigned segX, unsigned segZ)
{
  double dsegX(width / segX);
  double dsegZ(height / segZ);
  double x1=x0+width;
  double z1=z0+height;

  // init vertices
  int vid(0);
  for(unsigned j=0; j<segZ+1; j++)
  {
    double z = z0+j*dsegZ;
    for(unsigned i=0; i<segX+1; i++)
    {
      double x = x0+i*dsegX;
      Tspline::Vertex_handle v = CGAL::insert_point(tsp, Point2d( x, z ));
      v->data().id = vid++;
      v->data().SetCP(Point3d(x,y0,z));
    }
  }

  // init grid
  for(unsigned i=0; i<segX+1; i++)
  {
    double x = x0 + i * dsegX;
    CGAL::insert( tsp, Segment2(Point2d (x, z0), Point2d (x, z1)) );
  }
  for(unsigned j=0; j<segZ+1; j++)
  {
    double z = z0 + j * dsegZ;
    CGAL::insert( tsp, Segment2(Point2d (x0, z), Point2d (x1, z)) );
  }

  // init edge distances
  tspline::Tspline::Halfedge_iterator eit;
  for (eit = tsp.halfedges_begin(); eit != tsp.halfedges_end(); eit++)
  {
    const Point2d& p1 = eit->source()->point();
    const Point2d& p2 = eit->target()->point();
    eit->data().d = std::sqrt((p1 - p2).squared_length());
  }

  // update parameter space of tsp
  tsp.update_params();
  tsp.update_knot_vectors();
}

void TsplineCreator::CreatePlaneXY(tspline::Tspline &tsp,
                                   std::vector<double> paramS,
                                   std::vector<double> paramT,
                                   std::vector<Point4d> controlpoints)
{
  double x0 = paramS.front();
  double x1 = paramS.back();
  double y0 = paramT.front();
  double y1 = paramT.back();

  bool use_cps(false);
  if(controlpoints.size() == paramS.size() * paramT.size())
    use_cps = true;

  // init vertices
  int vid(0);
  for(unsigned j=0; j<paramT.size(); j++)
  {
    double y = paramT[j];
    for(unsigned i=0; i<paramS.size(); i++)
    {
      double x = paramS[i];
      Tspline::Vertex_handle v = CGAL::insert_point(tsp, Point2d( x, y ));
      v->data().id = vid++;
      if(use_cps)
        v->data().SetCP(controlpoints[j*paramS.size()+i]);
      else
        v->data().SetCP(Point3d(x,y,0));
    }
  }

  // init grid
  for(unsigned i=0; i<paramS.size(); i++)
  {
    double x = paramS[i];
    CGAL::insert( tsp, Segment2(Point2d (x, y0), Point2d (x, y1)) );
  }
  for(unsigned j=0; j<paramT.size(); j++)
  {
    double y = paramT[j];
    CGAL::insert( tsp, Segment2(Point2d (x0, y), Point2d (x1, y)) );
  }

  // init edge distances
  tspline::Tspline::Halfedge_iterator eit;
  for (eit = tsp.halfedges_begin(); eit != tsp.halfedges_end(); eit++)
  {
    const Point2d& p1 = eit->source()->point();
    const Point2d& p2 = eit->target()->point();
    eit->data().d = std::sqrt((p1 - p2).squared_length());
  }

  // update parameter space of tsp
  tsp.update_params();
  tsp.update_knot_vectors();
}

void TsplineCreator::CreateBox(tspline::TsplineMultiPatch &tsmp,
                               unsigned res, unsigned char multiplicity)
{
  std::vector<int> id;
  for(size_t i=0; i<6; i++)
    id.push_back(tsmp.AddPatch(res, res));

  std::list<TsplinePatch*>::iterator it;
  for(it=tsmp.patchlist.begin(); it!=tsmp.patchlist.end(); it++)
  {
    TsplinePatch* tsp = (*it);

    for(Tspline::Vertex_iterator vit=tsp->vertices_begin(); vit!=tsp->vertices_end(); vit++)
    {
      TVertex &vext = vit->data();
      Eigen::Vector3d cp;

      if(tsp->GetID() == id[0])
        cp = Eigen::Vector3d(0.5-vext.param.x(), vext.param.y()-0.5, -0.5);

      else if(tsp->GetID() == id[1])
        cp = Eigen::Vector3d(0.5-vext.param.x(), 0.5, vext.param.y()-0.5);

      else if(tsp->GetID() == id[2])
        cp = Eigen::Vector3d(0.5-vext.param.x(), 0.5-vext.param.y(), 0.5);

      else if(tsp->GetID() == id[3])
        cp = Eigen::Vector3d(0.5-vext.param.x(), -0.5, 0.5-vext.param.y());

      else if(tsp->GetID() == id[4])
        cp = Eigen::Vector3d(-0.5, vext.param.x()-0.5, 0.5-vext.param.y());

      else if(tsp->GetID() == id[5])
        cp = Eigen::Vector3d(0.5, 0.5-vext.param.x(), 0.5-vext.param.y());

      else
        throw std::runtime_error("[TsplineCreator::CreateBox] Error, id does not exist.");

      vext.SetCP(cp);
    }
  }

  tsmp.AddTransition(id[0], NORTH, id[1], SOUTH, multiplicity);
  tsmp.AddTransition(id[1], NORTH, id[2], SOUTH, multiplicity);
  tsmp.AddTransition(id[2], NORTH, id[3], SOUTH, multiplicity);
  tsmp.AddTransition(id[3], NORTH, id[0], SOUTH, multiplicity);

  tsmp.AddTransition(id[4], NORTH, id[0], EAST, multiplicity);
  tsmp.AddTransition(id[4], EAST,  id[1], EAST, multiplicity);
  tsmp.AddTransition(id[4], SOUTH, id[2], EAST, multiplicity);
  tsmp.AddTransition(id[4], WEST,  id[3], EAST, multiplicity);

  tsmp.AddTransition(id[5], NORTH, id[0], WEST, multiplicity);
  tsmp.AddTransition(id[5], WEST,  id[1], WEST, multiplicity);
  tsmp.AddTransition(id[5], SOUTH, id[2], WEST, multiplicity);
  tsmp.AddTransition(id[5], EAST,  id[3], WEST, multiplicity);
}

void TsplineCreator::CreateSphere(tspline::TsplineMultiPatch &tsmp,
                                  unsigned res)
{
  CreateBox(tsmp, res);

  std::list<TsplinePatch*>::iterator it;
  for(it=tsmp.patchlist.begin(); it!=tsmp.patchlist.end(); it++)
  {
    TsplinePatch* tsp = (*it);

    for(Tspline::Vertex_iterator vit=tsp->vertices_begin(); vit!=tsp->vertices_end(); vit++)
    {
      TVertex &vext = vit->data();
      Eigen::Vector3d cp;
      vext.GetCP(cp);

      double norm = cp.norm();
      if(!equal(norm,0.0))
        cp /= norm;

      vext.SetCP(cp);
    }
  }
}
