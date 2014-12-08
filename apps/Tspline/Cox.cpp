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

/** @brief Testing basis function evaluation */

#include "v4r/TomGine/tgTomGineThread.h"
#include "Tspline/Cox.h"
#include <iostream>

using namespace tspline;
using namespace std;

TomGine::tgTomGineThread viewer (800, 600, "NURSS Basis Functions");

struct Spline
{
  std::vector<double> knots;
  int degree;
  std::vector<double> N, Nd;
  double evaluate (const double &xi)
  {
    cox (xi, degree, knots, N);
    return N[knots.size () * degree];
  }
  double evaluate (const double &xi, double &dxi)
  {
    cox (xi, degree, knots, N);
    coxder (degree, knots, N, Nd);
    dxi = Nd[knots.size () * degree];
    return N[knots.size () * degree];
  }
};

struct Point{
  double x,y;
  Point():x(0.0),y(0.0){}
  Point(double _x, double _y) : x(_x), y(_y){}
  Point operator*(const float& f) { return Point(x*f,y*f); }
  void operator+=(const Point& p) { x+=p.x; y+=p.y; }
};

int main ()
{
  Spline tmp;
  tmp.degree = 3;

  std::vector<Spline> splines (7, tmp);
  std::vector<Point> cps(7);

  int a (-2);
  int b (2);

  for (int i = a; i <= b; i++)
  {
    splines[0].knots.push_back (i - 3);
    splines[1].knots.push_back (i - 2);
    splines[2].knots.push_back (i - 1);
    splines[3].knots.push_back (i + 0);
    splines[4].knots.push_back (i + 1);
    splines[5].knots.push_back (i + 2);
    splines[6].knots.push_back (i + 3);
  }

  for (size_t j = 0; j < splines.size (); j++)
  {
    printf ("spline[%lu] ", j);
    for (size_t i = 0; i < splines[j].knots.size (); i++)
    {
      double &v = splines[j].knots[i];
      if (v < a)
        v = a;
      else if (v > b)
        v = b;
      printf ("%.1f ", v);
    }
    printf ("\n");
  }

//  for (size_t i = 0; i <= cps.size(); i++)
//  {
//    double xi = a + double(b-a) * i / (cps.size()-1);
//    viewer.AddLine3D (xi, 0.0, 0.0, xi, 1.0, 0.0, 255, 255, 255, 1.0);
//  }

  // knots
  for (int i = a; i <= b; i++)
  {
    double xi = double(i);
    viewer.AddLine3D (xi, 0.0, 0.0, xi, 1.0, 0.0, 255, 255, 255, 1.0);
  }

  for(size_t i=0; i<cps.size(); i++)
  {
    double x = a + double(b-a) * i / (cps.size()-1);
    Point &p = cps[i];
    p = Point(x, 1.1); //+i%2);
    viewer.AddPoint3D(p.x, p.y, 0.0, 0,255,0, 10.0);
    if(i>0)
      viewer.AddLine3D(cps[i-1].x,cps[i-1].y,0.0, p.x,p.y,0.0, 0,255,0,1.0);
  }

  for (double xi = a; xi <= b + 1e-6; xi += 0.02)
  {
    double V (0);
    Point c;
    //size_t j = 0;
    for (size_t j = 0; j < splines.size (); j++)
    {
      double dxi;
      double v = splines[j].evaluate (xi, dxi);
      viewer.AddPoint3D (xi, v, 0.0, 255, 0, 0, 3.0);
      //      viewer.AddPoint3D (xi, dxi, 0.0, 0, 255, 0, 5.0); // derivative
      V += v;
      c += (cps[j]*v);
    }
    viewer.AddPoint3D (xi, V, 0.0, 255,255,255, 5.0);
    viewer.AddPoint3D (c.x, c.y, 0.0, 0,255,0, 3.0);
  }

  TomGine::tgCamera cam = viewer.GetCamera();
  TomGine::mat4 e;
  e.identity();
  cam.SetExtrinsic(e);
  cam.TranslateF(-4.0f);
  viewer.SetCamera(cam);

  viewer.SetRotationCenter(0.5*(a+b), 1.0, 0.0);
  viewer.Update ();
  viewer.WaitForEvent (TomGine::TMGL_Press, TomGine::TMGL_Space);
  return 0;
}
