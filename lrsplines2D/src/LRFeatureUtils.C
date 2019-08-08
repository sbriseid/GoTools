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


#include "GoTools/lrsplines2D/LRFeatureUtils.h"
#include "GoTools/geometry/RectDomain.h"

#include <fstream>

using namespace Go;
using std::vector;
using std::cout;
using std::endl;
using std::ifstream;

//#define DEBUG


//===========================================================================
void  LRFeatureUtils::writeCellInfo(const LRSplineSurface& srf, 
				    double tol, int ncell,
				    std::ostream &out)
//===========================================================================
{
  RectDomain dom = srf.containingDomain();
  double u1 = dom.umin();
  double u2 = dom.umax();
  double v1 = dom.vmin();
  double v2 = dom.vmax();
  double udel = (u2 - u1)/(double)ncell;
  double vdel = (v2 - v1)/(double)ncell;
  int ix = (srf.dimension() == 3) ? 4 : 2;

  int nc2 = ncell*ncell;
  vector<int> nmb_pts(3*nc2, 0);
  vector<double> cellinfo(6*nc2, 0.0);
  int *npt = &nmb_pts[0];
  int *nout_over = npt+nc2;
  int *nout_under = nout_over+nc2;
  double *avdist = &cellinfo[0];
  double *max_over = avdist+nc2;
  double *max_under = max_over+nc2;
  double *minheight = max_under+nc2;
  double *maxheight = minheight+nc2;
  double *nel = maxheight+nc2;

  // Adjust default
  for (int kr=0; kr<nc2; ++kr)
    {
      max_over[kr] = std::numeric_limits<double>::lowest();
      max_under[kr] = std::numeric_limits<double>::max();
      maxheight[kr] = std::numeric_limits<double>::lowest();
      minheight[kr] = std::numeric_limits<double>::max();
    }

  // Traverse elements and relate point and element information to
  // cell
  for (LRSplineSurface::ElementMap::const_iterator it=srf.elementsBegin();
       it != srf.elementsEnd(); ++it)
    {
      // Identify overlapping cells
      double uel1 = it->second->umin();
      double uel2 = it->second->umax();
      double vel1 = it->second->vmin();
      double vel2 = it->second->vmax();
      int i1 = (int)((uel1 - u1)/udel);
      int i2 = (int)((uel2 - u1)/udel);
      i1 = std::max(0, std::min(ncell-1, i1));
      i2 = std::max(0, std::min(ncell-1, i2));
      int j1 = (int)((vel1 - v1)/vdel);
      int j2 = (int)((vel2 - v1)/vdel);
      j1 = std::max(0, std::min(ncell-1, j1));
      j2 = std::max(0, std::min(ncell-1, j2));

      // Compute element fraction in cell
      for (int ki=i1; ki<=i2; ++ki)
	{
	  double ufrac = (std::min(uel2, u1+(ki+1)*udel) -
			  std::max(uel1, u1+ki*udel))/(uel2 - uel1);
	  for (int kj=j1; kj<=j2; ++kj)
	    {
	      double vfrac = (std::min(vel2, v1+(kj+1)*vdel) -
			      std::max(vel1, v1+kj*vdel))/(vel2 - vel1);
	      nel[kj*ncell+ki] += (ufrac*vfrac);
	    }
	}

      // Fetch associated data points
      vector<double>& points = it->second->getDataPoints();
      int nmb = it->second->nmbDataPoints();
      int del = it->second->getNmbValPrPoint();
      for (int kr=0; kr<nmb; ++kr)
	{
	  // Identify cell
	  i1 = (int)((points[kr*del] - u1)/udel);
	  i1 = std::max(0, std::min(ncell-1, i1));
	  j1 = (int)((points[kr*del+1] - v1)/vdel);
	  j1 = std::max(0, std::min(ncell-1, j1));
	  int ixcell = j1*ncell + i1;
	  avdist[ixcell] += points[kr*del+ix+1];
	  max_over[ixcell] = std::max(max_over[ixcell], points[kr*del+ix+1]);
	  max_under[ixcell] = std::min(max_under[ixcell], points[kr*del+ix+1]);
	  maxheight[ixcell] = std::max(maxheight[ixcell], points[kr*del+ix]);
	  minheight[ixcell] = std::min(minheight[ixcell], points[kr*del+ix]);
	  npt[ixcell]++;
	  if (points[kr*del+ix+1] > tol)
	    nout_over[ixcell]++;
	  if (points[kr*del+ix+1] < -tol)
	    nout_under[ixcell]++;
	}
    }
  
  for (int kr=0; kr<nc2; ++kr)
    {
      if (npt[kr] > 0)
	avdist[kr] /= (double)npt[kr];
    }

  // Write to file
  out << ncell << "  " << ncell << "  " << "8" << std::endl;
  for (int kr=0; kr<nc2; ++kr)
    {
      if (npt[kr] == 0)
	out << "0.0  0.0  0.0  ";
      else
	{
	  out << (double)(nout_under[kr]+nout_over[kr])/(double)npt[kr] << "  "; 
	  out << (double)nout_under[kr]/(double)npt[kr] << "  ";
	  out << (double)nout_over[kr]/(double)npt[kr] << "  ";
	}
      out << avdist[kr] << "  " ;
      if (npt[kr] == 0)
	out << "0.0  0.0  0.0";
      else
	{
	  out << avdist[kr]/std::max(fabs(max_under[kr]), fabs(max_over[kr]));
	  double max_diff = std::max(max_over[kr], 0.0) -
	    std::min(max_under[kr], 0.0);
	  double height_diff = maxheight[kr] - minheight[kr];
	  out <<"  " << max_diff << "  " << height_diff;
	}
      out << "  " << nel[kr] << std::endl;
    }
}

