/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRBSpline2D.h"
#include "GoTools/lrsplines2D/BSplineUniLR.h"
#include "GoTools/lrsplines2D/Mesh2DUtils.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"

#include <iostream>
#include <fstream>
#include <string.h>

using std::vector;

using namespace Go;

int main(int argc, char *argv[])
{
  if (argc != 3) {
    std::cout << "Usage: lrspline_in (.g2) lrspline_out.g2 " << std::endl;
    return -1;
  }

  std::ifstream filein(argv[1]);
  std::ofstream fileout(argv[2]);

  // Create the default factory
  GoTools::init();
  Registrator<LRSplineSurface> r293;

  // Read input surface
  ObjectHeader header;
  try {
    header.read(filein);
  }
  catch (...)
    {
      std::cerr << "Exiting" << std::endl;
      exit(-1);
    }
  shared_ptr<GeomObject> geom_obj(Factory::createObject(header.classType()));
  geom_obj->read(filein);
  
  shared_ptr<ParamSurface> sf = dynamic_pointer_cast<ParamSurface, GeomObject>(geom_obj);
  if (!sf.get())
    {
      std::cerr << "Input file contains no surface" << std::endl;
      exit(-1);
    }

  shared_ptr<LRSplineSurface> lrsf = 
    dynamic_pointer_cast<LRSplineSurface, ParamSurface>(sf);
  if (!lrsf.get())
    {
      shared_ptr<SplineSurface> splsf = 
	dynamic_pointer_cast<SplineSurface, ParamSurface>(sf);
      if (splsf.get())
	lrsf = shared_ptr<LRSplineSurface>(new LRSplineSurface(splsf.get(), 1.0e-6));
    }
  if (!lrsf.get())
    {
      std::cerr << "Input file contains no spline surface" << std::endl;
      exit(-1);
    }
    
  int dim = lrsf->dimension();
  const Mesh2D mesh = lrsf->mesh();
  int end1 = mesh.numDistinctKnots(XFIXED);
  int end2 = mesh.numDistinctKnots(YFIXED);

  int degree1 = lrsf->degree(XFIXED);
  int degree2 = lrsf->degree(YFIXED);
  vector<shared_ptr<BSplineUniLR> > uni1;
  vector<shared_ptr<BSplineUniLR> > uni2;
  std::map<LRSplineSurface::BSKey, shared_ptr<LRBSpline2D> > bspl2d;
  BSplineUniLR *u1, *u2;
  size_t kv;
  Point cg(dim);
  double wgt = 1.0, gamma = 1.0;
  for (int kj=0; kj<end2; ++kj)
    for (int ki=0; ki<end1; ++ki)
      {
	vector<vector<int> > kvec_u, kvec_v;
	Mesh2DUtils::derive_Bspline_mesh(mesh, degree1, degree2,
					 ki, kj, kvec_u, kvec_v);

	for (size_t kr=0; kr<kvec_u.size(); ++kr)
	  {
	    shared_ptr<BSplineUniLR> curr_u(new BSplineUniLR(1, degree1,
							     kvec_u[kr].begin(),
							     &mesh));
	    for (kv=0; kv<uni1.size(); ++kv)
	      if (uni1[kv].get() == curr_u.get())
		break;
	    if (kv == uni1.size())
	      uni1.push_back(curr_u);
	    u1 = uni1[kv].get();

	    shared_ptr<BSplineUniLR> curr_v(new BSplineUniLR(2, degree2,
							     kvec_v[kr].begin(),
							     &mesh));
	    for (kv=0; kv<uni2.size(); ++kv)
	      if (uni2[kv].get() == curr_v.get())
		break;
	    if (kv == uni2.size())
	      uni2.push_back(curr_v);
	    u2 = uni2[kv].get();

	    shared_ptr<LRBSpline2D> bspl(new LRBSpline2D(cg, wgt, u1, u2,
							 gamma));

	    LRSplineSurface::BSKey key = 
	      LRSplineSurface::generate_key(*bspl, mesh);
	    bspl2d.insert(std::make_pair(key, bspl));
	  }
	int stop_break = 1;
      }

  // Write to file
  fileout << "293 1 0 0" << std::endl;
  fileout << "0 1e-05" << std::endl;
  fileout << mesh;
  fileout << bspl2d.size() << std::endl;
  for (auto it=bspl2d.begin(); it!=bspl2d.end(); ++it)
    it->second->write(fileout);

  return 0;
}
