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

#include "GoTools/lrsplines2D/Element2DAccuracyHistory.h"
#include "GoTools/lrsplines2D/Element2DAccuracyInfo.h"
#include "GoTools/geometry/LineCloud.h"
#include "GoTools/geometry/PointCloud.h"
#include <fstream>
#include <iostream>

using namespace Go;
using std::vector;

void Element2DAccuracyHistory::checkAccuracyChange(int level, char* filename,
						   char* filename2, double frac)
{
  if (level >= elementhist_.size()-1)
    return;   // No comparisement possible

  int colour[4][3] = {{0, 0, 255}, {0, 255, 0}, {100, 55, 100}, {255, 0, 0}};
  vector<vector<double> > elem_bd(4);
  vector<double> elem_in;
  vector<double> elem_out;
  double eps = 1.0e-6;

  int nmb_none = 0, nmb_under = 0, nmb_div = 0;;
  std::ofstream os(filename);
  for (size_t ki=0; ki<elementhist_[level].size(); ++ki)
    {
      double umin, umax, vmin, vmax;
      double maxerr, averr;
      int nmb_pts, nmb_out;
      elementhist_[level][ki]->getElementDomain(umin, umax, vmin, vmax);
      elementhist_[level][ki]->getAccuracyInfo(nmb_pts, nmb_out, maxerr, averr);
      vector<Element2DAccuracyInfo*> next = 
	elementhist_[level][ki]->getNextElements();

      os << "Element " << ki+1 << ", domain: [" << umin << "," << umax;
      os << "] x [" << vmin << "," << vmax <<"]" << std::endl;
      os << nmb_pts << ", " << nmb_out << ", " << maxerr << ", " << averr;
      os << ", " << (int)next.size() << ": ";

      Point mid(0.5*(umin+umax), 0.5*(vmin+vmax), 0.0);
      if (nmb_out == 0)
	elem_in.insert(elem_in.end(), mid.begin(), mid.end());
      else
	elem_out.insert(elem_out.end(), mid.begin(), mid.end());

      if (next.size() > 1)
	nmb_div++;
      if (nmb_pts == 0 && next.size() > 1)
	nmb_none++;
      if (nmb_pts > 0 && next.size() > 1 && nmb_out == 0)
	nmb_under++;

      double maxerr2 = 0.0;
      double averr2 = 0.0;
      double nmb_pts2 = 0.0;
      double nmb_out2 = 0.0;
      vector<double> curr_bd;
      if (next.size() > 1)
	{
	  curr_bd.resize(next.size() == 2 ? 6 : 12);

	  double umin2, umax2, vmin2, vmax2;
	  next[0]->getElementDomain(umin2, umax2, vmin2, vmax2);
	  if (next.size() == 2)
	    {
	      if (umax2-umin2 < umax-umin-eps)
		{
		  if (umin2-umin > umax-umax2)
		    curr_bd[0] = curr_bd[3] = umin2;
		  else
		    curr_bd[0] = curr_bd[3] = umax2;
		  curr_bd[1] = vmin;
		  curr_bd[4] = vmax;
		  curr_bd[2] = curr_bd[5] = 0.0;
		}
	      else
	      if (vmax2-vmin2 < vmax-vmin-eps)
		{
		  if (vmin2-vmin > vmax-vmax2)
		    curr_bd[1] = curr_bd[4] = vmin2;
		  else
		    curr_bd[1] = curr_bd[4] = vmax2;
		  curr_bd[0] = umin;
		  curr_bd[3] = umax;
		  curr_bd[2] = curr_bd[5] = 0.0;
		}
	    }
	  else
	    {
	      if (umin2-umin > umax-umax2)
		curr_bd[0] = curr_bd[3] = umin2;
	      else
		curr_bd[0] = curr_bd[3] = umax2;
	      curr_bd[1] = vmin;
	      curr_bd[4] = vmax;
	      if (vmin2-vmin > vmax-vmax2)
		curr_bd[7] = curr_bd[10] = vmin2;
	      else
		curr_bd[7] = curr_bd[10] = vmax2;
	      curr_bd[6] = umin;
	      curr_bd[9] = umax;
	      curr_bd[2] = curr_bd[5] = curr_bd[8] = curr_bd[11] = 0.0;
	    }
	}

      for (size_t kj=0; kj<next.size(); ++kj)
	{
	  double maxerr3, averr3;
	  int nmb_pts3, nmb_out3;
	  next[kj]->getAccuracyInfo(nmb_pts3, nmb_out3, maxerr3, averr3);
	  maxerr2 = std::max(maxerr2, maxerr3);
	  averr2 += nmb_pts3*averr3;
	  nmb_pts2 += nmb_pts3;
	  nmb_out2 += nmb_out3;
	  os << nmb_pts3 << ", ";
	}
      if (nmb_pts2 > 0)
	averr2 /= (double)nmb_pts2;

      if (next.size() > 1)
	{
	  if (nmb_pts == 0)
	    elem_bd[1].insert(elem_bd[1].end(), curr_bd.begin(), curr_bd.end());
	  else if (nmb_out == 0)
	    elem_bd[3].insert(elem_bd[3].end(), curr_bd.begin(), curr_bd.end());
	  else if (maxerr2 > maxerr || averr2 > frac*averr || nmb_out2 >= nmb_out)
	    elem_bd[2].insert(elem_bd[2].end(), curr_bd.begin(), curr_bd.end());
	  else
	    elem_bd[0].insert(elem_bd[0].end(), curr_bd.begin(), curr_bd.end());
	}

      os << std::endl;
      os << nmb_pts2 << ", " << nmb_out2;
      os << ", " << maxerr2 << ", " << averr2 << std::endl << std::endl;
    }

  int nmb2 = (level+1 < (int)elementhist_.size()) ? elementhist_[level+1].size() : 0;
  std::cout << "Number of elements: " << elementhist_[level].size();
  std::cout << ", " << nmb2 << std::endl;
  std::cout << "Number of elements divided: " << nmb_div << std::endl;
  std::cout << "Divided elements without points: " << nmb_none << std::endl;
  std::cout << "Divided elements inside tolerance: " << nmb_under << std::endl;
  
  std::ofstream os2(filename2);
  os2 << "400 1 0 4 0 200 200 255" << std::endl;
  PointCloud3D pos_in(elem_in.begin(), elem_in.size()/3);
  pos_in.write(os2);
  
  os2 << "400 1 0 4 200 0 200 255" << std::endl;
  PointCloud3D pos_out(elem_out.begin(), elem_out.size()/3);
  pos_out.write(os2);
  for (size_t ki=0; ki<elem_bd.size(); ++ki)
    {
      if (elem_bd[ki].size() == 0)
	continue;
      LineCloud bds(elem_bd[ki].begin(), elem_bd[ki].size()/6);
      os2 << "410 1 0 4 " << colour[ki][0] << " " << colour[ki][1];
      os2 << " " << colour[ki][2] << " 255" << std::endl;
      bds.write(os2);
    }
  
}
