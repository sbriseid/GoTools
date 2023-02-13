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

#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/creators/ApproxCurve.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/utils/Array.h"
#include <iostream>
#include <fstream>

using namespace Go;

using std::cout;
using std::endl;
using std::vector;
using std::string;

int compare_u_par(const void* el1, const void* el2)
{
  if (((double*)el1)[0] < ((double*)el2)[0])
    return -1;
  else if (((double*)el1)[0] > ((double*)el2)[0])
    return 1;
  else
    return 0;
}

int compare_v_par(const void* el1, const void* el2)
{
  if (((double*)el1)[1] < ((double*)el2)[1])
    return -1;
  else if (((double*)el1)[1] > ((double*)el2)[1])
    return 1;
  else
    return 0;
}

int main(int argc, char* argv[])
{
  if (argc != 4)
    {
      cout << "Usage: input point cloud (g2), number of cells in one direction, output point cloud" << std::endl;
      exit(0);
    }

  std::ifstream filein(argv[1]);
  int nmb_div = atoi(argv[2]);
  std::ofstream fileout(argv[3]);

  ObjectHeader header;
  PointCloud3D points;
  try {
    header.read(filein);
    points.read(filein);
  }
  catch (...)
    {
      std::cout << "ERROR: Not a valid point file" << std::endl;
      return -1;
    }

  BoundingBox box = points.boundingBox();
  double umin = box.low()[0];
  double umax = box.high()[0];
  double vmin = box.low()[1];
  double vmax = box.high()[1];
  double l1 = umax - umin;
  double l2 = vmax - vmin;

  int nm = nmb_div*nmb_div;
  double dom = (umax-umin)*(vmax-vmin);
  double c1 = std::pow((double)nm/dom, 1.0/2.0);
  int div1, div2;
  div1 = (int)(c1*(umax-umin));
  ++div1;
  div2 = (int)(c1*(vmax-vmin));
  ++div2;
  double udel = (umax - umin)/(double)(div1);
  double vdel = (vmax - vmin)/(double)(div2);

  std::ofstream of("division_lines.g2");
  of << "410 1 0 4 255 0 0 255" << std::endl;
  of << div1+div2+2 << std::endl;
  for (int ki=0; ki<=div1; ++ki)
    {
      Point p1(umin+ki*udel, vmin, 0.0);
      Point p2(umin+ki*udel, vmax, 0.0);
      of << p1 << " " << p2 << std::endl;
    }
  for (int ki=0; ki<=div2; ++ki)
    {
      Point p1(umin, vmin+ki*vdel, 0.0);
      Point p2(umax, vmin+ki*vdel, 0.0);
      of << p1 << " " << p2 << std::endl;
    }

  int nmb_pts = points.numPoints();
  vector<double> param(2*nmb_pts);

  for (int ki=0; ki<nmb_pts; ++ki)
    {
      Vector3D curr = points.point(ki);
      param[2*ki] = curr[0];
      param[2*ki+1] = curr[1];
    }

  std::ofstream ofpar("parpoints_init.g2");
  int nmbpar = param.size()/2;
  ofpar << "400 1 0 0" << std::endl;
  ofpar << nmbpar << std::endl;
  for (int ka=0; ka<nmbpar; ++ka)
    ofpar << param[2*ka] << " " << param[2*ka+1] << "  0.0" << std::endl;

  // Define raster
  vector<vector<int> > raster(div2);
  for (size_t kr=0; kr<raster.size(); ++kr)
    {
      raster[kr].resize(div1, 0);
    }

  qsort(&param[0], nmb_pts, 2*sizeof(double), compare_v_par);
  int pp0, pp1;
  int ka, kb;
  double upar, vpar;
  for (vpar=vmin+vdel, pp0=0, kb=0; kb<div2; ++kb, vpar+=vdel)
    {
      for (pp1=pp0; pp1<2*nmb_pts && param[pp1+1]<=vpar; pp1+=2);
      qsort(&param[pp0], (pp1-pp0)/2, 2*sizeof(double), compare_u_par);

      int pp2, pp3;
      for (upar=umin+udel, pp2=pp0, ka=0; ka<div1; ++ka, upar+=udel)
	{
	  for (pp3=pp2; pp3<pp1 && param[pp3]<=upar; pp3+=2);
	  raster[kb][ka] = (pp3-pp2)/2;
	  pp2 = pp3;
	}
      pp0 = pp1;
    }

  // Count number of empty cells
  int nmb_zero = 0;
  for (kb=0; kb<div2; ++kb)
    for (ka=0; ka<div1; ++ka)
      if (raster[kb][ka] == 0)
	++nmb_zero;

  // Find start corner
  int corner[4][2] = {{0, 0}, {div1-1, 0}, {0, div2-1}, {div1-1, div2-1}};
  int sgn1 = 1, sgn2 = 1;
  int ki, kj, kc, kd;
  for (kc=0, kj=0, sgn2=1; kj<2; ++kj,sgn2*=-1)
    for (sgn1=1, ki=0; ki<2; ++ki,  ++kc, sgn1*=-1)
      {
	for (ka=corner[kc][0], kd=0; kd<div1; ++kd, ka+=sgn1)
	  if (raster[corner[kc][1]][ka] > 0)
	    break;
	for (kb=corner[kc][1], kd=0; kd<div2; ++kd, kb+=sgn2)
	  if (raster[kb][corner[kc][0]] > 0)
	    break;
	if (fabs(corner[kc][1]-kb) < fabs(corner[kc][1]-ka))
	  corner[kc][1] = kb;
	else
	  corner[kc][0] = ka;
      }
  

  vector<vector<double> > bdhit(4);
  int i1, i2, i3, i4, ix;
  for (kc=0, ix=0; kc<2; ++kc, ix=div2-1)
    for (i1=i2=0; i1<div1; i1=i2)
      {
	for (; i1<div1; ++i1)
	  if (raster[ix][i1] > 0)
	    break;
	for (i2=i1; i2<div1; ++i2)
	  if (raster[ix][i2] == 0)
	    break;
	if (i2 > i1)
	  {
	    bdhit[2*kc].push_back(0.5*(i1+i2));
	    bdhit[2*kc].push_back(0.5*(i2-i1-1));
	  }
      }

  for (kc=0, ix=0; kc<2; ++kc, ix=div1-1)
    for (i1=i2=0; i1<div2; i1=i2)
      {
	for (; i1<div2; ++i1)
	  if (raster[i1][ix] > 0)
	    break;
	for (i2=i1; i2<div2; ++i2)
	  if (raster[i2][ix] == 0)
	    break;
	if (i2 > i1)
	  {
	    bdhit[2*kc+1].push_back(0.5*(i1+i2));
	    bdhit[2*kc+1].push_back(0.5*(i2-i1-1));
	  }
      }

  vector<Point> cand_par;
  for (size_t kr=0; kr<bdhit.size(); ++kr)
    {
      i1 = (kr==0 || kr==2) ? 1 : 0;
      i2 = 1 - i1;
      double t1 = (kr == 0 || kr == 1) ? 0.0 :
	((kr == 2) ? (double)(div2) : (double)(div1));
      for (size_t kh=0; kh<bdhit[kr].size(); kh+=2)
	{
	  Point curr(2);
	  curr[i1] = t1;
	  curr[i2] = bdhit[kr][kh];
	  cand_par.push_back(curr);
	}
    }

  int ix1=-1, ix2=-1;
  double maxd = 0.0;
  for (size_t kr=0; kr<cand_par.size(); ++kr)
    for (size_t kh=kr+1; kh<cand_par.size(); ++kh)
      {
	double dd = cand_par[kr].dist(cand_par[kh]);
	if (dd > maxd)
	  {
	    ix1 = (int)kr;
	    ix2 = (int)kh;
	    maxd = dd;
	  }
      }

  vector<double> ptpar;
  ptpar.push_back(cand_par[ix1][0]);
  ptpar.push_back(cand_par[ix1][1]);

  double curr[2], curr1[2], curr2[2];
  for (kc=0; kc<(int)ptpar.size(); kc+=2)
    {
      std::ofstream ofc("midpol.g2");
      ofc << "410 1 0 4 0 0 0 255" << std::endl;
      ofc << ptpar.size()/2 << std::endl;
      ofc << umin+udel*ptpar[0] << " " << vmin+vdel*ptpar[1] << "  0.0" << std::endl;
      for (size_t kr=2; kr<ptpar.size(); kr+=2)
	{
	  ofc << umin+udel*ptpar[kr] << " " << vmin+vdel*ptpar[kr+1] << "  0.0" << std::endl;
	  ofc << umin+udel*ptpar[kr] << " " << vmin+vdel*ptpar[kr+1] << "  0.0" << std::endl;
	}
      ofc << umin+udel*cand_par[ix2][0] << " " << vmin+vdel*cand_par[ix2][1] << " 0.0" << std::endl;

      curr[0] = curr1[0] = curr2[0] = ptpar[kc];
      curr[1] = curr1[1] = curr2[1] = ptpar[kc+1];
      int sgn1 = (cand_par[ix2][0] > curr[0]) ? 1 : -1;
      int sgn2 = (cand_par[ix2][1] > curr[1]) ? 1 : -1;
      double fac1 = (curr[0] == 0.0) ? 0.5 : 1.0;
      double fac2 = (curr[1] == 0.0) ? 0.5 : 1.0;
      curr[0] += fac1*sgn1;
      curr[1] += fac2*sgn2;
      curr1[0] += fac1*sgn1;
      curr2[1] += fac2*sgn2;
      Point prev(ptpar[kc], ptpar[kc+1]);
      if (prev.dist(cand_par[ix2]) <= 1.0)
	break;
      Point curr_2(curr[0],curr[1]);
      if ((cand_par[ix2]-prev)*(cand_par[ix2]-curr_2) < 0.0)
	break;
      curr[0] = std::min(curr[0], div1-0.5);
      curr1[0] = std::min(curr1[0], div1-0.5);
      curr[1] = std::min(curr[1], div2-0.5);
      curr2[1] = std::min(curr2[1], div2-0.5);

      int j1 = (int)curr1[0];
      int j2 = (int)curr2[1];
      for (i1=(int)curr1[1]; i1<div2; ++i1)
	if (raster[i1][j1] > 0)
	  break;
      if (i1 == div2)
	i1=(int)curr1[1];
      for (; i1<div2; ++i1)
	if (raster[i1][j1] == 0)
	  break;
      for (i2=(int)curr1[1]; i2>=0; --i2)
	if (raster[i2][j1] > 0)
	  break;
      if (i2 < 0)
	i2=(int)curr1[1];
      for (; i2>=0; --i2)
	if (raster[i2][j1] == 0)
	  break;
      curr1[1] = 0.5*(double)(i1+i2+1);

      for (i3=(int)curr2[0]; i3<div1; ++i3)
	if (raster[j2][i3] > 0)
	  break;
      if (i3 == div1)
	i3=(int)curr2[0];
      for (; i3<div1; ++i3)
	if (raster[j2][i3] == 0)
	  break;
      for (i4=(int)curr2[0]; i4>=0; --i4)
	if (raster[j2][i4] > 0)
	  break;
      if (i4 < 0)
	i4=(int)curr2[0]; 
      for (; i4>=0; --i4)
	if (raster[j2][i4] == 0)
	  break;
      curr2[0] = 0.5*(double)(i3+i4+1);
      // if (fabs(curr[1]-ptpar[kc+1]) > fabs(curr[0]-ptpar[kc]) /*&&
      // 								(kc == 0 || fabs(cand_par[ix2][0]-curr[0]) < fabs(cand_par[ix2][0]-ptpar[kc-2]))*/)
      // 	curr[1] = ptpar[kc+1] + fac2*sgn2;
      // else
      // 	curr[0] = ptpar[kc] + fac1*sgn1;

      if (i1-i2 < i3-i4)
	{
	  ptpar.push_back(curr1[0]);
	  ptpar.push_back(curr1[1]);
	}
      else
	{
	  ptpar.push_back(curr2[0]);
	  ptpar.push_back(curr2[1]);
	}
	
      
      int stop_break0 = 1;      
    }
  Point last(ptpar[ptpar.size()-2], ptpar[ptpar.size()-1]);
  if (last.dist(cand_par[ix2]) > 0.0)
    {
      ptpar.push_back(cand_par[ix2][0]);
      ptpar.push_back(cand_par[ix2][1]);
    }

  // Translate to real coordinates
  vector<double> ptpar2(ptpar.size());
  for (size_t kr=0; kr<ptpar.size(); kr+=2)
    {
      ptpar2[kr] = umin+udel*ptpar[kr];
      ptpar2[kr+1] = vmin+vdel*ptpar[kr+1];
    }

  std::ofstream ofp("parpt.g2");
  ofp << "400 1 0 4 200 55 0 255" << std::endl;
  ofp << ptpar2.size()/2 << std::endl;
  for (size_t kr=0; kr<ptpar2.size(); kr+=2)
    ofp << ptpar2[kr] << " " << ptpar2[kr+1] << " 0.0" << std::endl;

  vector<double> par2(ptpar.size()/2);
  par2[0] = 0.0;
  for (size_t kr=2; kr<ptpar2.size(); kr+=2)
    {
      double dd = Utils::distance_squared(&ptpar2[kr-2], &ptpar2[kr], &ptpar2[kr]);
      par2[kr/2] = par2[(kr-2)/2] + sqrt(dd);
    }

  int in = 4, ik = 4;
  int dim = 2;
  double eps = 0.1;
  ApproxCurve approx(ptpar2, par2, dim, eps, in, ik);
  
  double maxdist, avdist;
  shared_ptr<SplineCurve> midcv = approx.getApproxCurve(maxdist, avdist, 1);

  std::ofstream ofmid("midcv.g2");
  midcv->writeStandardHeader(ofmid);
  midcv->write(ofmid);

  // Reparameterize
  vector<double> param2(param.size());
  double tmin = midcv->startparam();
  double tmax = midcv->endparam();
  double umin2 = std::numeric_limits<double>::max();
  double umax2 = std::numeric_limits<double>::lowest();
  double vmin2 = std::numeric_limits<double>::max();
  double vmax2 = std::numeric_limits<double>::lowest();
  for (size_t kr=0; kr<param.size(); kr+=2)
    {
      Point currpar(param[kr], param[kr+1]);
      double tpar, dist;
      Point close;
      midcv->closestPoint(currpar, tmin, tmax, tpar, close, dist);

      vector<Point> der(2);
      midcv->point(der, tpar, 1);
      Point vec = close - currpar;
      Point vec2(-vec[1], vec[0]);
      if (vec2*der[1] < 0.0)
	dist *= -1.0;
      param2[kr] = tpar;
      param2[kr+1] = dist;
      umin2 = std::min(umin2, tpar);
      umax2 = std::max(umax2, tpar);
      vmin2 = std::min(vmin2, dist);
      vmax2 = std::max(vmax2, dist);
      int stop_break0 = 1;
    }
  
  std::ofstream ofpar2("parpoints_repar.g2");
  nmbpar = param2.size()/2;
  ofpar2 << "400 1 0 0" << std::endl;
  ofpar2 << nmbpar << std::endl;
  for (int ka=0; ka<nmbpar; ++ka)
    ofpar2 << param2[2*ka] << " " << param2[2*ka+1] << "  0.0" << std::endl;

  double dom2 = (umax2-umin2)*(vmax2-vmin2);
  double c2 = std::pow((double)nm/dom2, 1.0/2.0);
  int div3, div4;
  div3 = (int)(c2*(umax2-umin2));
  ++div3;
  div4 = (int)(c2*(vmax2-vmin2));
  ++div4;
  double udel2 = (umax2 - umin2)/(double)(div3);
  double vdel2 = (vmax2 - vmin2)/(double)(div4);
  
  // Define raster
  vector<vector<int> > raster2(div4);
  for (size_t kr=0; kr<raster2.size(); ++kr)
    {
      raster2[kr].resize(div3, 0);
    }

  qsort(&param2[0], nmb_pts, 2*sizeof(double), compare_v_par);
  for (vpar=vmin2+vdel2, pp0=0, kb=0; kb<div4; ++kb, vpar+=vdel2)
    {
      for (pp1=pp0; pp1<2*nmb_pts && param2[pp1+1]<=vpar; pp1+=2);
      qsort(&param2[pp0], (pp1-pp0)/2, 2*sizeof(double), compare_u_par);

      int pp2, pp3;
      for (upar=umin2+udel2, pp2=pp0, ka=0; ka<div3; ++ka, upar+=udel2)
	{
	  for (pp3=pp2; pp3<pp1 && param2[pp3]<=upar; pp3+=2);
	  raster2[kb][ka] = (pp3-pp2)/2;
	  pp2 = pp3;
	}
      pp0 = pp1;
    }

  // Count number of empty cells
  int nmb_zero2 = 0;
  for (kb=0; kb<div4; ++kb)
    for (ka=0; ka<div3; ++ka)
      if (raster2[kb][ka] == 0)
	++nmb_zero2;

  std::ofstream of2("division_lines2.g2");
  of2 << "410 1 0 4 255 0 0 255" << std::endl;
  of2 << div3+div4+2 << std::endl;
  for (int ki=0; ki<=div3; ++ki)
    {
      Point p1(umin2+ki*udel2, vmin2, 0.0);
      Point p2(umin2+ki*udel2, vmax2, 0.0);
      of2 << p1 << " " << p2 << std::endl;
    }
  for (int ki=0; ki<=div4; ++ki)
    {
      Point p1(umin2, vmin2+ki*vdel2, 0.0);
      Point p2(umax2, vmin2+ki*vdel2, 0.0);
      of2 << p1 << " " << p2 << std::endl;
    }

  int stop_break = 1;
}

