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

#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/GoIntersections.h"
#include <fstream>
#include <algorithm>


using namespace Go;
using std::vector;
using std::pair;


int main(int argc, char** argv)
{
  if (argc != 10)
    {
      std::cout << "Input parameters: infile (spline curve, g2), lower left (2 floats), upper right (2 floats), number of lines, par dir (0/1), tolerance, outfile" << std:: endl;
      exit(1);
    }

  std::ifstream incv(argv[1]);
  double ll1 = atof(argv[2]);
  double ll2 = atof(argv[3]);
  double ur1 = atof(argv[4]);
  double ur2 = atof(argv[5]);
  int nmb = atoi(argv[6]);
  int pdir = atoi(argv[7]);
  double tol = atof(argv[8]);
  std::ofstream outcvs(argv[9]);

  if (ur1 <= ll1 && ur2 <= ll2)
    {
      std::cout << "Empty or negative interval" << std::endl;
      exit(1);
    }
  
  // Read curve
  ObjectHeader head;
  head.read(incv);

  shared_ptr<SplineCurve> crv1(new SplineCurve());
  crv1->read(incv);

  Point pos1(ll1, ll2);
  Point pos2 = pos1;
  Point ur(ur1, ur2);
  pos2[1-pdir] = ur[1-pdir];
  Point del(0.0, 0.0);
  del[pdir] = (ur[pdir] - pos1[pdir])/(double)(nmb-1);

  for (int ki=0; ki<nmb; ++ki, pos1+=del, pos2+=del)
    {
      shared_ptr<SplineCurve> crv2(new SplineCurve(pos1, pos2));

      vector<pair<double, double> > intpar;
      vector<int> pretop;
      vector<pair<pair<double,double>, pair<double,double> > > intcrvs;
      intersect2Dcurves(crv1.get(), crv2.get(), tol, intpar, pretop, intcrvs);

      for (size_t kj=0; kj<intpar.size(); ++kj)
	{
	  Point pos3 = crv1->ParamCurve::point(intpar[kj].first);
	  shared_ptr<SplineCurve> crv3(new SplineCurve(pos1, pos3));

	  crv3->writeStandardHeader(outcvs);
	  crv3->write(outcvs);
	}
    }
}

  
