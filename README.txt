Copyright (c) 2014, Thomas Mörwald, Vienna University of Technology.
License available in the "License.txt" file.

Author: Thomas Moerwald (Vienna University of Technology)
Email: moerwald@acin.tuwien.ac.at

Tspline geometry toolkit (Version 0.1)

Modules:
	- B-spline Basis function evaluation (cox-de-boor)
	- T-spline evaluation, control point insertion/removal, face splitting
	- T-spline merge/split
	- Multi-patch T-splines enabling extraordinary vertices
	- Rectangular domain creation
	- T-spline creation (plane, box, sphere)

Dependencies:
  Required:
  	- CGAL: Computational Geometry Algorithms Library (libcgal-dev)
  	- Eigen: Template library for linear algebra (libeigen3-dev)
  Optional:
    - OpenMP: Open Multi-Processing library (www.openmp.com)
    - TomGine 4.1: Visualization and rendering engine (http://users.acin.tuwien.ac.at/tmoerwald/?site=4)

TomGine control:
  [left mousebutton] rotate view
  [right mousebutton] pan view
  [middle mousebutton/scroll] zoom view
  [d] show dots on/off
  [l] show lines on/off
  [m] show mesh on/off   
  [o] print out current view matrix and camera matrix   
  [w] wireframe on/off
  [z] reset view
  [F11] save screenshot to file
  [q/Esc] Quit
  

Installation:
  mkdir build
  cd build
  cmake ..
  make

  * Note that apps require TomGine 4.1 for visualization (http://users.acin.tuwien.ac.at/tmoerwald/?site=4)
