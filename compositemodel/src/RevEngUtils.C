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

#include "GoTools/compositemodel/RevEngUtils.h"
#include "GoTools/utils/MatrixXD.h"
#include "GoTools/creators/SmoothSurf.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SISLconversion.h"
#include "sislP.h"
#include "newmat.h"
#include "newmatap.h"
#include <fstream>

using namespace Go;
using std::vector;

typedef MatrixXD<double, 3> Matrix3D;

//===========================================================================
void RevEngUtils::principalAnalysis(Point& curr, vector<Point>& points, 
				    double lambda[3], double eigenvec[3][3])
//===========================================================================
{
  // Compute covariance matrix
  int numpt = (int)points.size();
  double comat[3][3];
  Point mean = curr;
  for (int kr=0; kr<numpt; ++kr)
    mean += points[kr];
  mean /= (double)(numpt+1);

  vector<double> wgt(numpt+1);
  double div = (double)((numpt+1)*(numpt+1));
  wgt[0] = (mean.dist(curr))/div;
  for (int kr=0; kr<numpt; ++kr)
    {
      double dr2 = mean.dist2(points[kr]);
      wgt[kr+1] = exp(-dr2/div);
    }
  
  for (int ki=0; ki<3; ++ki)
    for (int kj=0; kj<3; ++kj)
      {
	double tmp = 0.0;
	tmp += wgt[0]*(curr[ki] - mean[ki])*(curr[kj] - mean[kj]);
	for (int kr=0; kr<numpt; ++kr)
	  {
	    tmp += wgt[kr+1]*(points[kr][ki] - mean[ki])*(points[kr][kj] - mean[kj]);
	  }
	comat[ki][kj] = tmp/(double)(numpt);
      }
  
  // Compute singular values
  NEWMAT::Matrix nmat;
  nmat.ReSize(3, 3);
  for (int ki = 0; ki < 3; ++ki) {
    for (int kj = 0; kj < 3; ++kj) {
      nmat.element(ki, kj) = comat[ki][kj];
    }
  }
      
  static NEWMAT::DiagonalMatrix diag;
  static NEWMAT::Matrix V;
  try {
    NEWMAT::SVD(nmat, diag, nmat, V);
  } catch(...) {
    std::cout << "Exception in SVD" << std::endl;
    return;
  }
  
  // Singular values
  for (int ki=0; ki<3; ++ki)
    {
      lambda[ki] = diag.element(ki, ki);
      for (int kj=0; kj<3; ++kj)
	eigenvec[ki][kj] = V.element(kj, ki);
    }
}

//===========================================================================
void RevEngUtils::computeMonge(Point& curr, std::vector<Point>& points,
			       Point& vec1, Point& vec2, Point& normal, Point& mincvec,
			       double& minc, Point& maxcvec, double& maxc,
			       double& currdist, double& avdist)
//===========================================================================
{
  // Transform points to coordinate system given by vec1 (x-axis) and vec2 (y-axis)
  Matrix3D mat1, mat2, rotmat;
  Vector3D vec1_2(vec1[0], vec1[1], vec1[2]);
  Vector3D vec2_2(vec2[0], vec2[1], vec2[2]);
  Vector3D xaxis(1, 0, 0);
  Vector3D zaxis(0, 0, 1);
  mat1.setToRotation(vec1_2, xaxis);
  Vector3D v1 = mat1*vec1_2;
  Vector3D vec2_3 = mat1*vec2_2;
  mat2.setToRotation(vec2_3, zaxis);
  Vector3D v2 = mat2*vec2_3;
  rotmat = mat2*mat1;
  //rotmat.identity();

  // Perform rotation and sort parameter values and z-value
  int nmbpts = (int)points.size() + 1;
  vector<double> par(2*nmbpts);
  vector<double> zval(nmbpts);
  double xmin, xmax, ymin, ymax;
  xmin = xmax = par[0] = curr[0];
  ymin = ymax = par[1] = curr[1];
  zval[0] = curr[2];
  for (int ki=1; ki<nmbpts; ++ki)
    {
      Point dv = points[ki-1] - curr;
      Vector3D dv2(dv[0], dv[1], dv[2]);
      Vector3D dvrot = rotmat*dv2;
      //Vector3D dvrot = mat2*dvrot0;
      par[2*ki] = curr[0] + dvrot[0];
      par[2*ki+1] = curr[1] + dvrot[1];
      zval[ki] = curr[2] + dvrot[2];
      xmin = std::min(xmin, par[2*ki]);
      xmax = std::max(xmax, par[2*ki]);
      ymin = std::min(ymin, par[2*ki+1]);
      ymax = std::max(ymax, par[2*ki+1]);
    }

  // Approximate z-component by biquadratic Bezier function in x and y
  // First make quadratic B-spline function with coefficients equal to zero
  int order = 3;
  vector<double> knots1(2*order), knots2(2*order);
  for (int kj=0; kj<order; ++kj)
    {
      knots1[kj] = xmin;
      knots1[order+kj] = xmax;
      knots2[kj] = ymin;
      knots2[order+kj] = ymax;
    }
  vector<double> coefs(order*order, 0.0);
  shared_ptr<SplineSurface> bez(new SplineSurface(order, order, order, order, 
						  &knots1[0], &knots2[0], &coefs[0], 1));

  // Approximate
  SmoothSurf approx;
  vector<int> coef_known(nmbpts, 0);
  int seem[2];
  seem[0] = seem[1] = 0;
  vector<double> ptwgt(nmbpts, 1.0);
  approx.attach(bez, seem, &coef_known[0]);

  double wgt1 = 0.0, wgt2 = 0.0001, wgt3 = 0.0;
  double approxwgt = 1.0 - wgt1 - wgt2 - wgt3;
  approx.setOptimize(wgt1, wgt2, wgt3);
  approx.setLeastSquares(zval, par, ptwgt, approxwgt);

  shared_ptr<SplineSurface> mongesf;
  approx.equationSolve(mongesf);

  vector<double> coefs2(3*order*order);
  std::vector<double>::iterator cf = mongesf->coefs_begin();
  for (int ka=0; ka<order; ++ka)
    {
      double vpar = mongesf->basis_v().grevilleParameter(ka);
      for (int kb=0; kb<order; ++kb, ++cf)
	{
	  double upar = mongesf->basis_u().grevilleParameter(kb);
	  coefs2[(ka*order+kb)*3] = upar;
	  coefs2[(ka*order+kb)*3+1] = vpar;
	  coefs2[(ka*order+kb)*3+2] = *cf;
	}
    }
  shared_ptr<SplineSurface> tmp(new SplineSurface(order, order, order, order, 
						  &knots1[0], &knots2[0], &coefs2[0], 3));
  std::ofstream of("approx_sf.g2");
  tmp->writeStandardHeader(of);
  tmp->write(of);
  of << "400 1 0 4 0 255 0 255" << std::endl;
  of << 1 << std::endl;
  of << curr << std::endl;
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << nmbpts << std::endl;
  for (int ka=0; ka<nmbpts; ++ka)
    {
      Point tmppt(par[2*ka], par[2*ka+1], zval[ka]);
      of << tmppt << std::endl;
    }
  
  
  // Compute surface normal in curr
  vector<Point> der(3);
  mongesf->point(der, par[0], par[1], 1);
  Vector3D norm(-der[1][0], -der[2][0], 1.0);
  norm.normalize();

  // Accuracy of approximation
  currdist = fabs(zval[0] - der[0][0]);
  avdist = currdist;
  for (int ki=1; ki<nmbpts; ++ki)
    {
      Point pos;
      mongesf->point(pos, par[2*ki], par[2*ki+1]);
      avdist += fabs(zval[ki] - pos[0]);
    }
  avdist /= (double)nmbpts;
  
  // Compute principal curvatures in curr
  shared_ptr<SISLSurf> sislsf(GoSurf2SISL(*mongesf, false));
  int left1 = 0, left2 = 0;
  int stat = 0;
  double k1, k2;
  double d1[2], d2[2];
  s2542(sislsf.get(), 0, 0, 0, &par[0], &left1, &left2, &k1, &k2, d1, d2, &stat);
  Vector3D du(1.0, 0.0, der[1][0]);
  Vector3D dv(0.0, 1.0, der[2][0]);
  Vector3D cvec1 = d1[0]*du + d1[1]*dv;
  Vector3D cvec2 = d2[0]*du + d2[1]*dv;
  minc = k1;
  maxc = k2;
  
  // Transform results to original coordinate system
  Matrix3D mat3, mat4, rotmat2;
  mat4.setToRotation(zaxis, vec2_3);
  mat3.setToRotation(xaxis, vec1_2);
  rotmat2 = mat3*mat4;
  //rotmat2.identity();
  //Vector3D norm0 = mat4*norm;
  Vector3D norm2 = rotmat2*norm;
  normal = Point(norm2[0], norm2[1], norm2[2]);
  
  Vector3D cvec3 = rotmat2*cvec1;
  mincvec = Point(cvec3[0], cvec3[1],cvec3[2]); 
  Vector3D cvec4 = rotmat2*cvec2;
  maxcvec = Point(cvec4[0], cvec4[1],cvec4[2]); 

  int stop_break = 1;
}
