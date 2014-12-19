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

// [1] Isogeometric Analysis, Toward Integration of CAD and FEA; J.Austin Cottrell, Thomas J.R. Hughes, Yuri Bazilevs

#include "Math.hpp"
#include "Cox.h"
#include <math.h>

using namespace tspline;

namespace tspline
{
  void cox (const double &xi, const unsigned &degree, const std::vector<double> &knots, std::vector<double> &N)
  {
    size_t nknots = knots.size();
    N.assign (nknots * (degree + 1), 0.0);

    for (unsigned p = 0; p < (degree+1); p++)
    { // loop from lower degree to higher -> unwrapped recursion

      for (unsigned s = 0; s < (nknots-1); s++)
      { // evaluate the basis N for each knotspan s

        if (p == 0)
        {
          // Equation (2.1) in [1]
          if(gequal(xi, knots[s]))
            if(smaller(xi, knots[s+1]))
              N[s] = 1.0;

          // check for interpolating knot on the right side [hacky]
          if(equal(xi,knots[nknots-1]))
          {
            if(s+1+degree < nknots)
              if(equal(knots[s+1], knots[s+1+degree]))
                N[s] = 1.0;

            if(s+degree < nknots)
              if(equal(knots[s+1], knots[s+degree]))
                N[s] = 1.0;
          }

        }
        else
        {
          // Equation (2.2)
          double A(0.0);
          double B(0.0);

          if( !equal(knots[s],knots[s+p]) )
            A = (xi - knots[s]) / (knots[s+p] - knots[s]);

          if( !equal(knots[s+1],knots[s+p+1]) )
            B = (knots[s+p+1] - xi) / (knots[s+p+1] - knots[s+1]);

          if(A < 0.0)
            A = 0.0;
          if(B < 0.0)
            B = 0.0;

          N[s + p*nknots] = A * N[s + (p-1)*nknots] + B * N[(s+1) + (p-1)*nknots];

        }
      }
    }
  }

  void coxder (const unsigned &degree, const std::vector<double> &knots, const std::vector<double> &N,
               std::vector<double> &Nd)
  {
    size_t nknots = knots.size();
    Nd.assign (nknots * (degree + 1), 0.0);
    unsigned p = degree;

    for (unsigned s = 0; s < nknots - 1; s++)
    {
      // Equation (2.12)
      double A(0.0);
      double B(0.0);

      if(!equal (knots[s+p], knots[s]))
        A = p / (knots[s+p] - knots[s]);

      if(!equal (knots[s+p+1], knots[s+1]))
        B = p / (knots[s+p+1] - knots[s+1]);

      if (A < 0.0)
        A = 0.0;
      if (B < 0.0)
        B = 0.0;

      Nd[s + p*nknots] = A * N[s + (p-1)*nknots] - B * N[(s+1) + (p-1)*nknots];

    }
  }

}
