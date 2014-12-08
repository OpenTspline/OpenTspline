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

#include "Tspline/TsplineMultiPatch.h"
#include "Tspline/TsplineCreator.h"
#include "Tspline/Utils.hpp"

#include "Tspline/File.h"

using namespace tspline;
using namespace TomGine;

tgTomGineThread viewer (800, 600, "T-spline Box", true);

int main ()
{
  viewer.SetInputSpeeds(0.2f);
  viewer.SetClearColor(1.0);

  TomGine::tgRenderModel meshBox;

  TsplineMultiPatch tsBox;
  TsplineCreator::CreateBox(tsBox, 4);
  convertTspline2tgModel(tsBox, meshBox, 32, 32);
  convertTsplineControl2tgModel(tsBox, meshBox);

//  tgTsplineMultiPatch mpBox(8,8);
//  TsplineCreator::CreateBox(mpBox, 4, 0);
//  meshBox = mpBox.UpdateModel();
//  mpBox.GetControlGrid(meshBox);

  meshBox.SetColor(0.6f, 0.4f, 0.2f);
  meshBox.m_line_width = 4.0;
  viewer.AddModel3D(meshBox);
  viewer.Update();

//    File::Save(tsBox, "box.tmp");

  viewer.Update ();
  viewer.WaitForEvent (TMGL_Press, TMGL_Escape);
  return (0);
}
