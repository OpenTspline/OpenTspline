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

#ifndef _TSPLINE_TSPLINE_CREATOR_H_
#define _TSPLINE_TSPLINE_CREATOR_H_

// This file defines helper functions for creation certain T-splines
// It implements planes, boxes and cylinders

#include "Tspline.h"
#include "TsplineMultiPatch.h"

namespace tspline{

class TsplineCreator
{
public:
  static void CreateUniformPlaneXY(tspline::Tspline &tsp,
                                  double x0, double y0, double z0,
                                  double width, double height,
                                  unsigned segX, unsigned segY);

  static void CreatePlaneXY(tspline::Tspline &tsp,
                            double x0, double y0, double z0,
                            double width, double height,
                            unsigned segX, unsigned segY);

  static void CreatePlaneYZ(tspline::Tspline &tsp,
                            double x0, double y0, double z0,
                            double width, double height,
                            unsigned segY, unsigned segZ);

  static void CreatePlaneXZ(tspline::Tspline &tsp,
                            double x0, double y0, double z0,
                            double width, double height,
                            unsigned segX, unsigned segZ);

  static void CreatePlaneXY(tspline::Tspline &tsp,
                            std::vector<double> paramS,
                            std::vector<double> paramT,
                            std::vector<Point4d> controlpoints = std::vector<Point4d>());

  static void CreateBox(tspline::TsplineMultiPatch &tsmp, unsigned res, unsigned char multiplicity=0);

  static void CreateSphere(tspline::TsplineMultiPatch &tsmp, unsigned res);



};


}

#endif
