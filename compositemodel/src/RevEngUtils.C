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
#include "GoTools/compositemodel/ImplicitApprox.h"
#include "GoTools/utils/MatrixXD.h"
#include "GoTools/utils/LUDecomp.h"
#include "GoTools/creators/SmoothSurf.h"
#include "GoTools/creators/ApproxSurf.h"
#include "GoTools/creators/SmoothCurve.h"
#include "GoTools/creators/ApproxCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/Cone.h"
#include "GoTools/geometry/Sphere.h"
#include "GoTools/geometry/Torus.h"
#include "GoTools/geometry/Plane.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/geometry/GeometryTools.h"
#include "sislP.h"
#include "newmat.h"
#include "newmatap.h"
#include <fstream>

using namespace Go;
using std::vector;
using std::pair;

typedef MatrixXD<double, 3> Matrix3D;

//===========================================================================
void RevEngUtils::principalAnalysis(std::vector<RevEngPoint*>& points, 
				    double lambda[3], double eigenvec[3][3])
//===========================================================================
{
  if (points.size() < 5)
    return;
  vector<Point> remaining(points.size()-1);
  Vector3D xyz = points[0]->getPoint();
  Point curr(xyz[0], xyz[1], xyz[2]);
  for (size_t ki=1; ki<points.size(); ++ki)
    {
      xyz = points[ki]->getPoint();
      Point curr2(xyz[0], xyz[1], xyz[2]);
      remaining[ki-1] = curr2;
    }
  principalAnalysis(curr, remaining, lambda, eigenvec);
}

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
    //std::cout << "Exception in SVD" << std::endl;
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
void RevEngUtils::TaubinCurvature(Point curr, std::vector<Point>& points,
				  Point& tvec, Point& normal, Point& mincvec,
				  double& minc, Point& maxcvec, double& maxc)
//===========================================================================
{
  // Define matrix
  double mat[3][3];
  for (int ka=0; ka<3; ++ka)
    for (int kb=0; kb<3; ++kb)
      mat[ka][kb] = 0.0;

  double wijsum = 0.0;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      Point vec = points[ki] - curr;
      double len = vec.length();
      wijsum += (1.0/len);
    }
  
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      Point vec = points[ki] - curr;
      double len2 = vec.length2();
      double kij = 2.0*(normal*vec)/len2;
      double wij = 1.0/sqrt(len2);
      wij /= wijsum;

      vec -= (vec*normal)*normal;
      double fac = wij*kij;
      for (int ka=0; ka<3; ++ka)
	for (int kb=0; kb<3; ++kb)
	  mat[ka][kb] += fac*vec[ka]*vec[kb];

      // double phi = tvec.angle(vec);
      // double phicos2 = cos(phi);
      // phicos2 = phicos2*phicos2;
      // double phisin2 = sin(phi);
      // phisin2 = phisin2*phisin2;
      // A += wij*phicos2*phicos2;
      // B += wij*phisin2*phicos2;
      // C += wij*phisin2*phisin2;
    }
  
  // Compute singular values
  NEWMAT::Matrix nmat;
  nmat.ReSize(3, 3);
  for (int ka = 0; ka < 3; ++ka) {
    for (int kb = 0; kb < 3; ++kb) {
      nmat.element(ka, kb) = mat[ka][kb];
    }
  }
      
  static NEWMAT::DiagonalMatrix diag;
  static NEWMAT::Matrix V;
  try {
    NEWMAT::SVD(nmat, diag, nmat, V);
  } catch(...) {
    //std::cout << "Exception in SVD" << std::endl;
    return;
  }
  
  // Singular values
  double lambda1 = diag.element(0, 0);
  double lambda2 = diag.element(1, 1);
  Point cvec1 = Point(V.element(0, 0), V.element(0, 1), V.element(0, 2));
  Point cvec2 = Point(V.element(1, 0), V.element(1, 1), V.element(1, 2));

  double A=0.0, B=0.0, C=0.0;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      Point vec = points[ki] - curr;
      double len2 = vec.length2();
      double wij = 1.0/sqrt(len2);
      wij /= wijsum;
      vec -= (vec*normal)*normal;
      double phi = cvec1.angle(vec);
      double phicos2 = cos(phi);
      phicos2 = phicos2*phicos2;
      double phisin2 = sin(phi);
      phisin2 = phisin2*phisin2;
      A += wij*phicos2*phicos2;
      B += wij*phisin2*phicos2;
      C += wij*phisin2*phisin2;
    }
  
  double div = B*B - A*C;
  if (fabs(div) < 1.0e-10)
    {
      minc = maxc = 0.0;
      mincvec = maxcvec = Point(0.0, 0.0, 0.0);
      return;
    }
  
  double k1 = (B*lambda2 - C*lambda1)/div;
  double k2 = (B*lambda1 - A*lambda2)/div;
  if (fabs(k1) < fabs(k2))
    {
      minc = k1;
      mincvec = cvec1;
      maxc = k2;
      maxcvec = cvec2;
    }
  else
    {
      minc = k2;
      mincvec = cvec2;
      maxc = k1;
      maxcvec = cvec1;
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
  par[0] = curr[0];
  par[1] = curr[1];
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
    }

  // Approximate z-component by biquadratic Bezier function in x and y
  int order = 3;
  shared_ptr<SplineSurface> mongesf = RevEngUtils::surfApprox(zval, 1, par, order,
							      order, order, order);

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
						  mongesf->basis_u().begin(),
						  mongesf->basis_v().begin(), &coefs2[0], 3));
  int writesurface = 0;
  if (writesurface)
    {
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

  // Vector3D origin(par[0], par[1], zval[0]);
  // of << "410 1 0 4 0 0 0 255" << std::endl;
  // of << "1" << std::endl;
  // of << origin << " " << origin+norm << std::endl;

  // of << "410 1 0 4 0 55 155 255" << std::endl;
  // of << "1" << std::endl;
  // of << origin << " " << origin+cvec1 << std::endl;

  
  // of << "410 1 0 4 155 55 0 255" << std::endl;
  // of << "1" << std::endl;
  // of << origin << " " << origin+cvec2 << std::endl;

  
  
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

//===========================================================================
shared_ptr<SplineSurface>
RevEngUtils::surfApprox(vector<double>& data, int dim, vector<double>& param,
			int order1, int order2, int ncoef1, int ncoef2,
			bool close1, bool close2,
			int max_iter, double tol, double& maxd, double& avd, 
			int& num_out, vector<double>& parvals, double belt_frac)
//===========================================================================
{
  // Create initial spline space
  double umin = std::numeric_limits<double>::max();
  double umax = std::numeric_limits<double>::lowest();
  double vmin = std::numeric_limits<double>::max();
  double vmax = std::numeric_limits<double>::lowest();
  for (size_t ki=0; ki<param.size(); ki+=2)
    {
      umin = std::min(umin, param[ki]);
      umax = std::max(umax, param[ki]);
      vmin = std::min(vmin, param[ki+1]);
      vmax = std::max(vmax, param[ki+1]);
    }
  double udel = (umax - umin);
  double vdel = (vmax - vmin);
  if (!close1)
    {
      umin -= belt_frac*udel;
      umax += belt_frac*udel;
    }
  if (!close2)
    {
      vmin -= belt_frac*vdel;
      vmax += belt_frac*vdel;
    }
  udel = (umax - umin)/(double)(ncoef1 - order1 + 1);
  vdel = (vmax - vmin)/(double)(ncoef2 - order2 + 1);
  vector<double> et1(order1+ncoef1);
  vector<double> et2(order2+ncoef2);
  vector<double> coef(ncoef1*ncoef2*dim, 0.0);
  for (int ka=0; ka<order1; ++ka)
    {
      et1[ka] = umin;
      et1[ka+ncoef1] = umax;
    }
  for (int ka=0; ka<order2; ++ka)
    {
      et2[ka] = vmin;
      et2[ka+ncoef2] = vmax;
    }
  for (int ka=0; ka<ncoef1-order1; ++ka)
    et1[order1+ka] = umin + (ka+1)*udel; 
  for (int ka=0; ka<ncoef2-order2; ++ka)
    et2[order2+ka] = vmin + (ka+1)*vdel; 

  
 shared_ptr<SplineSurface> surf(new SplineSurface(ncoef1, ncoef2, order1,
						  order2, &et1[0], 
						  &et2[0], &coef[0], dim));
 shared_ptr<SplineSurface> surf1(new SplineSurface(ncoef1, ncoef2, order1,
						   order2, &et1[0], 
						   &et2[0], &coef[0], dim));

 // Approximate
  SmoothSurf approx2;
  vector<int> coef_known(ncoef1*ncoef2, 0);
  int seem[2];
  seem[0] = close1 ? 1 : 0;
  seem[1] = close2 ? 1 : 0;
  int nmbpts = (int)data.size()/dim;
  vector<double> ptwgt(nmbpts, 1.0);
  approx2.attach(surf1, seem, &coef_known[0]);

  //double wgt1 = 0.0, wgt2 = 0.001, wgt3 = 0.001;
  double wgt1 = 0.0, wgt2 = 0.001, wgt3 = 0.001;
  double approxwgt = 1.0 - wgt1 - wgt2 - wgt3;
  approx2.setOptimize(wgt1, wgt2, wgt3);
  approx2.setLeastSquares(data, param, ptwgt, approxwgt);
  shared_ptr<SplineSurface> surf3;
  approx2.equationSolve(surf3);

  std::ofstream of("first_approx.g2");
  surf3->writeStandardHeader(of);
  surf3->write(of);

  ApproxSurf approx(surf3, data, param, dim, tol);
  //ApproxSurf approx(surf3, data, param, dim, tol, 0, false, true, 0, true);
  approx.setMBA(true);
  approx.setFixBoundary(false);
  double acc_frac = 0.6;
  approx.setAccuracyCrit(1, acc_frac);
 shared_ptr<SplineSurface> surf2 = approx.getApproxSurf(maxd, avd, num_out, max_iter);
 if (surf2.get())
   parvals = approx.getParvals();
 return surf2;
}

//===========================================================================
shared_ptr<SplineSurface> RevEngUtils::surfApprox(vector<double>& data, int dim,
						  vector<double>& param, int order1,
						  int order2, int nmb_coef1, int nmb_coef2,
						  double del)
//===========================================================================
{
  // Define spline space
  double umin, umax, vmin, vmax;
  umin = umax = param[0];
  vmin = vmax = param[1];
  for (size_t kj=2; kj<param.size(); kj+=2)
    {
      umin = std::min(umin, param[kj]);
      umax = std::max(umax, param[kj]);
      vmin = std::min(vmin, param[kj+1]);
      vmax = std::max(vmax, param[kj+1]);
    }
  umin -= del;
  umax += del;
  vmin -= del;
  vmax += del;

  return surfApprox(data, dim, param, order1, order2, nmb_coef1, nmb_coef2,
		    umin, umax, vmin, vmax);
}
  
//===========================================================================
shared_ptr<SplineSurface> RevEngUtils::surfApprox(vector<double>& data, int dim,
						  vector<double>& param, int order1,
						  int order2, int nmb_coef1, int nmb_coef2,
						  double umin, double umax,
						  double vmin, double vmax)
//===========================================================================
{
  double udel = (umax - umin)/(double)(nmb_coef1-order1+1);
  double vdel = (vmax - vmin)/(double)(nmb_coef2-order2+1);
  
  vector<double> knots1(order1+nmb_coef1), knots2(order2+nmb_coef2);
  for (int ka=0; ka<order1; ++ka)
    {
      knots1[ka] = umin;
      knots1[nmb_coef1+ka] = umax;
    }
  for (int ka=order1; ka<nmb_coef1; ++ka)
    knots1[ka] = umin + (ka-order1+1)*udel;
  
  for (int ka=0; ka<order2; ++ka)
    {
      knots2[ka] = vmin;
      knots2[nmb_coef2+ka] = vmax;
    }
  for (int ka=order2; ka<nmb_coef2; ++ka)
    knots2[ka] = vmin + (ka-order2+1)*vdel;
  

  vector<double> coefs((nmb_coef1+order1)*(nmb_coef2+order2)*dim, 0.0);
  shared_ptr<SplineSurface> bez(new SplineSurface(nmb_coef1, nmb_coef2, order1, order2, 
						  &knots1[0], &knots2[0], &coefs[0], dim));

  // Approximate
  SmoothSurf approx;
  vector<int> coef_known((nmb_coef1+order1)*(nmb_coef2+order2), 0);
  int seem[2];
  seem[0] = seem[1] = 0;
  int nmbpts = (int)data.size()/dim;
  vector<double> ptwgt(nmbpts, 1.0);
  approx.attach(bez, seem, &coef_known[0]);

  //double wgt1 = 0.0, wgt2 = 0.001, wgt3 = 0.001;
  double wgt1 = 0.0, wgt2 = 0.05, wgt3 = 0.05;
  double approxwgt = 1.0 - wgt1 - wgt2 - wgt3;
  approx.setOptimize(wgt1, wgt2, wgt3);
  approx.setLeastSquares(data, param, ptwgt, approxwgt);

  shared_ptr<SplineSurface> surf;
  approx.equationSolve(surf);
  return surf;
 }

//===========================================================================
bool RevEngUtils::parameterizeOnPrimary(vector<RevEngPoint*>& points,
					shared_ptr<ParamSurface> surf,
					vector<double>& data, vector<double>& param,
					int& inner1, int& inner2, bool& close1,
					bool& close2)
//===========================================================================
{
  double eps = 1.0e-2;
  double angtol = 0.1;

  // Make sure to use an untrimmed surface
  shared_ptr<ParamSurface> sf = surf;
  shared_ptr<BoundedSurface> bdsf = dynamic_pointer_cast<BoundedSurface,ParamSurface>(sf);
  if (bdsf.get())
    sf = bdsf->underlyingSurface();
  shared_ptr<Cylinder> cyl = dynamic_pointer_cast<Cylinder,ParamSurface>(sf);
  shared_ptr<Cone> cone = dynamic_pointer_cast<Cone,ParamSurface>(sf);
  shared_ptr<Sphere> sph = dynamic_pointer_cast<Sphere,ParamSurface>(sf);
  shared_ptr<Torus> torus = dynamic_pointer_cast<Torus,ParamSurface>(sf);
  close1 = (cyl.get() || cone.get() || sph.get() || torus.get());
  close2 = (sph.get() || torus.get());
  RectDomain dom = sf->containingDomain();
  double u1 = dom.umin();
  double u2 = dom.umax();
  double udel = u2 - u1;
  double v1 = dom.vmin();
  double v2 = dom.vmax();
  double vdel = v2 - v1;

  // Parameterize
  int dim = surf->dimension();
  param.resize(2*points.size());
  data.reserve(dim*points.size());
  double umin = std::numeric_limits<double>::max();
  double umax = std::numeric_limits<double>::lowest();
  double vmin = std::numeric_limits<double>::max();
  double vmax = std::numeric_limits<double>::lowest();
  double sfac = 0.5;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      Vector3D xyz = points[ki]->getPoint();
      Point pos(xyz[0], xyz[1], xyz[2]);
      double upar, vpar, dist;
      Point close;
      sf->closestPoint(pos, upar, vpar, close, dist, eps);

      if (sf->isBounded() && dist > 0.1) //eps) TEST
	{
	  // Check boundary
	  double uparbd, vparbd, distbd;
	  Point closebd;
	  double seed[2];
	  seed[0] = upar;
	  seed[1] = vpar;
	  RectDomain dom = sf->containingDomain();
	  sf->closestBoundaryPoint(pos, uparbd, vparbd, closebd, distbd, eps, &dom, seed);
	  if (fabs(distbd - dist) < 0.1*dist)
	    {
	      // Angular check to see if a true closest point is found
	      Point norm;
	      sf->normal(norm, upar, vpar);
	      Point vec = pos - close;
	      double ang = norm.angle(vec);
	      ang = std::min(ang, fabs(M_PI - ang));
	      if (ang > angtol)
		return false;
	    }
	}

      // Check with seam
      // if (ki > 0 && close1 && ((upar-u1 < umin-upar && u2-umax < umin-u1) ||
      // 			       (u2-upar < upar-umax && umin-u1 < u2-umax)))
      // if (ki > 0 && close1 && (upar < umin || upar > umax) &&
      // 	  (fabs(umin-upar+udel) < fabs(upar-umin) || fabs(upar+udel-umax) < fabs(umax-upar)))
	  if (ki > 0 && close1 && (upar < umin || upar > umax) &&
	      std::min(fabs(umin-upar+udel),fabs(upar+udel-umax)) < std::min(fabs(upar-umin),fabs(upar-umax)))
	{
	  // if (ki >= 2)
	  //   {
	  //     // Compare distances
	  //     double d1 = points[ki-2]->pntDist(points[ki-1]);
	  //     double d2 = points[ki-1]->pntDist(points[ki]);
	  //     Point p1(param[2*(ki-2)], param[2*(ki-2)+1]);
	  //     Point p2(param[2*(ki-1)], param[2*(ki-1)+1]);
	  //     Point p3(upar,vpar);
	  //     double dp1 = p1.dist(p2);
	  //     double dp2 = p2.dist(p3);
	  //     if (dp1/dp2 < sfac*d1/d2)
	  // 	{
		  if (upar-u1 < u2-upar)
		    upar += (u2-u1);
		  else
		    upar -= (u2-u1);
	    // 	}
									 
	    //   int stop_break1 = 1;
	    // }
	}
      // if (ki > 0 && close2 && ((vpar-v1 < vmin-vpar && v2-vmax < vmin-v1) ||
      // 			       (v2-vpar < vpar-vmax && vmin-v1 < v2-vmax)))
	  // if (ki > 0 && close2 && (vpar < vmin || vpar > vmax) &&
	  //     (fabs(vmin-vpar+vdel) < fabs(vpar-vmin) || fabs(vpar+vdel-vmax) < fabs(vpar-vmax)))
	  if (ki > 0 && close2 && (vpar < vmin || vpar > vmax) &&
	      std::min(fabs(vmin-vpar+vdel),fabs(vpar+vdel-vmax)) < std::min(fabs(vpar-vmin),fabs(vpar-vmax)))
	{
	  // if (ki >= 2)
	  //   {
	  //     // Compare distances
	  //     double d1 = points[ki-2]->pntDist(points[ki-1]);
	  //     double d2 = points[ki-1]->pntDist(points[ki]);
	  //     Point p1(param[2*(ki-2)], param[2*(ki-2)+1]);
	  //     Point p2(param[2*(ki-1)], param[2*(ki-1)+1]);
	  //     Point p3(upar,vpar);
	  //     double dp1 = p1.dist(p2);
	  //     double dp2 = p2.dist(p3);
	  //     if (dp1/dp2 < sfac*d1/d2)
	  // 	{
		  if (vpar-v1 < v2-vpar)
		    vpar += (v2-v1);
		  else
		    vpar -= (v2-v1);
	    // 	}
									 
	    //   int stop_break2 = 1;
	    // }
	}
      param[2*ki] = upar;
      param[2*ki+1] = vpar;
      umin = std::min(umin, upar);
      umax = std::max(umax, upar);
      vmin = std::min(vmin, vpar);
      vmax = std::max(vmax, vpar);
      data.insert(data.end(), pos.begin(), pos.end());
    }

  if (close1 && umax-umin < 2*M_PI-eps)
    close1 = false;
  if (close2 && vmax-vmin < 2*M_PI-eps)
    close2 = false;
  
  // Reparameterize for surfaces with circular properties
  inner1 = inner2 = 0;
  if (cyl.get())
    {
      double rad = cyl->getRadius();
      for (size_t ki=0; ki<param.size(); ki+=2)
	param[ki] *= rad;
      inner1 = 2.0*(umax - umin)/M_PI;
    }

  if (cone.get())
    {
      double rad1 = cone->radius(0.0, vmin);
      double rad2 = cone->radius(0.0, vmax);
      double rad = 0.5*(rad1 + rad2);
      for (size_t ki=0; ki<param.size(); ki+=2)
	param[ki] *= rad;
      inner1 = 2.0*(umax - umin)/M_PI;
    }
  
  if (sph.get())
    {
      double rad = sph->getRadius();
      for (size_t ki=0; ki<param.size(); ++ki)
	param[ki] *= rad;
      inner1 = 2.0*(umax - umin)/M_PI;
      inner2 = 2.0*(vmax - vmin)/M_PI;
    }

  if (torus.get())
    {
      double rad1 = torus->getMajorRadius();
      double rad2 = torus->getMinorRadius();
      for (size_t ki=0; ki<param.size(); ki+=2)
	{
	  param[ki] *= rad1;
	  param[ki+1] *= rad2;
	}
      inner1 = 2.0*(umax - umin)/M_PI;
      inner2 = 2.0*(vmax - vmin)/M_PI;
    }
  return true;
}

//===========================================================================
void RevEngUtils::parameterizeWithPlane(vector<RevEngPoint*>& pnts,
					const BoundingBox& bbox,
					const Point& vec1, const Point& vec2,
					vector<double>& data, vector<double>& param)
//===========================================================================
{
  vector<Point> pos(pnts.size());
  for (size_t ki=0; ki<pnts.size(); ++ki)
    {
      Vector3D xyz = pnts[ki]->getPoint();
      pos[ki] = Point(xyz[0], xyz[1], xyz[2]);
    }

  parameterizeWithPlane(pos, bbox, vec1, vec2, data, param);
}

//===========================================================================
void RevEngUtils::parameterizeWithPlane(vector<Point>& pnts, const BoundingBox& bbox,
					const Point& vec1, const Point& vec2,
					vector<double>& data, vector<double>& param)
//===========================================================================
{
  double eps = 1.0e-6;
  Point mid = 0.5*(bbox.low() + bbox.high());
  int dim = mid.dimension();
  double diag = bbox.low().dist(bbox.high());
  int order = 2;
  double et[4];
  et[0] = et[1] = -diag;
  et[2] = et[3] = diag;
  vector<double> coefs;
  int sgn1, sgn2;
  int ka, kb;
  for (kb=0, sgn2=-1; kb<2; ++kb, sgn2=1)
    for (ka=0, sgn1=-1; ka<2; ++ka, sgn1=1)
      {
	Point pos = mid+sgn1*diag*vec1+sgn2*diag*vec2;
	coefs.insert(coefs.end(), pos.begin(), pos.end());
      }

  shared_ptr<SplineSurface> surf(new SplineSurface(order, order, order, order, &et[0], 
						   &et[0], coefs.begin(), dim));
  std::ofstream of("parplane.g2");
  surf->writeStandardHeader(of);
  surf->write(of);
  
  param.resize(2*pnts.size());
  data.reserve(dim*pnts.size());
  for (size_t ki=0; ki<pnts.size(); ++ki)
    {
      double upar, vpar, dist;
      Point close;
      surf->closestPoint(pnts[ki], upar, vpar, close, dist, eps);
      param[2*ki] = upar;
      param[2*ki+1] = vpar;
      data.insert(data.end(), pnts[ki].begin(), pnts[ki].end());
    }
}

//===========================================================================
void RevEngUtils::computeAxis(vector<pair<vector<RevEngPoint*>::iterator,
			      vector<RevEngPoint*>::iterator> >& points,
			      Point& axis, Point& Cx, Point& Cy)
//===========================================================================
{
  double Cmat[3][3];
  for (int ka=0; ka<3; ++ka)
    for (int kb=0; kb<3; ++kb)
      Cmat[ka][kb] = 0.0;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      vector<RevEngPoint*>::iterator start = points[ki].first;
      vector<RevEngPoint*>::iterator end = points[ki].second;
      for (int ka=0; ka<3; ++ka)
	for (int kb=0; kb<3; ++kb)
	  {
	    for (auto it=start; it!=end; ++it)
	      {
		RevEngPoint *pt = *it;
		Point norm1 = pt->getMongeNormal();
		Point norm2 = pt->getTriangNormal();
		Point norm = norm1; //0.5*(norm1 + norm2);
		Cmat[ka][kb] += norm[ka]*norm[kb];
		//Cmat[ka][kb] += norm2[ka]*norm2[kb];
	      }
	  }
    }
  
  // Compute singular values
  NEWMAT::Matrix nmat;
  nmat.ReSize(3, 3);
  for (int ka = 0; ka < 3; ++ka) {
    for (int kb = 0; kb < 3; ++kb) {
      nmat.element(ka, kb) = Cmat[ka][kb];
    }
  }
      
  static NEWMAT::DiagonalMatrix diag;
  static NEWMAT::Matrix V;
  try {
    NEWMAT::SVD(nmat, diag, nmat, V);
  } catch(...) {
    //std::cout << "Exception in SVD" << std::endl;
    exit(-1);
  }
  Cx = Point(V.element(0,0), V.element(1,0), V.element(2,0));
  Cy = Point(V.element(0,1), V.element(1,1), V.element(2,1));
  axis = Point(V.element(0,2), V.element(1,2), V.element(2,2));

}


//===========================================================================
void RevEngUtils::computeAxis(vector<Point>& points,
			      Point& axis, Point& Cx, Point& Cy)
//===========================================================================
{
  double Cmat[3][3];
  for (int ka=0; ka<3; ++ka)
    for (int kb=0; kb<3; ++kb)
      {
	Cmat[ka][kb] = 0.0;
	for (size_t ki=0; ki<points.size(); ++ki)
	  {
	    Cmat[ka][kb] += points[ki][ka]*points[ki][kb];
	  }
      }
  
  // Compute singular values
  NEWMAT::Matrix nmat;
  nmat.ReSize(3, 3);
  for (int ka = 0; ka < 3; ++ka) {
    for (int kb = 0; kb < 3; ++kb) {
      nmat.element(ka, kb) = Cmat[ka][kb];
    }
  }
      
  static NEWMAT::DiagonalMatrix diag;
  static NEWMAT::Matrix V;
  try {
    NEWMAT::SVD(nmat, diag, nmat, V);
  } catch(...) {
    //std::cout << "Exception in SVD" << std::endl;
    exit(-1);
  }
  Cx = Point(V.element(0,0), V.element(1,0), V.element(2,0));
  Cy = Point(V.element(0,1), V.element(1,1), V.element(2,1));
  axis = Point(V.element(0,2), V.element(1,2), V.element(2,2));

}


//===========================================================================
void RevEngUtils::coneAxis(vector<pair<vector<RevEngPoint*>::iterator,
			      vector<RevEngPoint*>::iterator> >& points,
			      Point& axis, Point& Cx, Point& Cy)
//===========================================================================
{
  size_t numpt = 0;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      numpt += (points[ki].second - points[ki].first);
    }
  double wgt = 1.0/(double)numpt;

  Point mid(3);
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      vector<RevEngPoint*>::iterator start = points[ki].first;
      vector<RevEngPoint*>::iterator end = points[ki].second;
      for (auto it=start; it!=end; ++it)
	{
	  RevEngPoint *pt = *it;
	  Point norm = pt->getMongeNormal();
	  mid += wgt*norm;
	}
    }

  double Cmat[3][3];
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      vector<RevEngPoint*>::iterator start = points[ki].first;
      vector<RevEngPoint*>::iterator end = points[ki].second;
      for (int ka=0; ka<3; ++ka)
	for (int kb=0; kb<3; ++kb)
	  {
	    Cmat[ka][kb] = 0.0;
	    for (auto it=start; it!=end; ++it)
	      {
		RevEngPoint *pt = *it;
		Point norm = pt->getMongeNormal();
		Point vec = norm - mid;
		Cmat[ka][kb] += vec[ka]*vec[kb];
	      }
	  }
    }
  
  // Compute singular values
  NEWMAT::Matrix nmat;
  nmat.ReSize(3, 3);
  for (int ka = 0; ka < 3; ++ka) {
    for (int kb = 0; kb < 3; ++kb) {
      nmat.element(ka, kb) = Cmat[ka][kb];
    }
  }
      
  static NEWMAT::DiagonalMatrix diag;
  static NEWMAT::Matrix V;
  try {
    NEWMAT::SVD(nmat, diag, nmat, V);
  } catch(...) {
    //std::cout << "Exception in SVD" << std::endl;
    exit(-1);
  }
  Cx = Point(V.element(0,0), V.element(1,0), V.element(2,0));
  Cy = Point(V.element(0,1), V.element(1,1), V.element(2,1));
  axis = Point(V.element(0,2), V.element(1,2), V.element(2,2));

}

//===========================================================================
void RevEngUtils::coneApex(vector<pair<vector<RevEngPoint*>::iterator,
			      vector<RevEngPoint*>::iterator> >& points,
			   Point axis, Point& apex, double& phi)
//===========================================================================
{
  double Mmat[3][3], Mi[3][3];
  double bvec[3], bi[3];
  for (int ka=0; ka<3; ++ka)
    {
      bvec[ka] = 0.0;
      for (int kb=0; kb<3; ++kb)
	Mmat[ka][kb] = 0.0;
    }

  int nmb = 0;
  for (size_t ki=0; ki<points.size(); ++ki)
    nmb += (int)(points[ki].second - points[ki].first);

  vector<Point> dird;
  vector<Point> pp;
  dird.reserve(nmb);
  pp.reserve(nmb);
  double wg = 1.0/(double)nmb;
    for (size_t ki=0; ki<points.size(); ++ki)
    {
      vector<RevEngPoint*>::iterator start = points[ki].first;
      vector<RevEngPoint*>::iterator end = points[ki].second;
      for (auto it=start; it!=end; ++it)
	{
	  RevEngPoint *pt = *it;
	  Point norm = pt->getMongeNormal();
	  Point tmp = norm.cross(axis);
	  Point di = tmp.cross(norm);
	  di.normalize_checked();
	  dird.push_back(di);
	  Vector3D xyz = pt->getPoint();
	  Point pos(xyz[0], xyz[1], xyz[2]);
	  pp.push_back(pos);

	  for (int ka=0; ka<3; ++ka)
	    for (int kb=0; kb<3; ++kb)
	      {
		if (ka == kb)
		  continue;
		Mi[ka][kb] = -di[ka]*di[kb];
	      }
	  Mi[0][0] = di[1]*di[1] + di[2]*di[2];
	  Mi[1][1] = di[0]*di[0] + di[2]*di[2];
	  Mi[2][2] = di[0]*di[0] + di[1]*di[1];

	  bi[0] = pos[0]*di[1]*di[1] - pos[1]*di[0]*di[1] - pos[2]*di[0]*di[2] + pos[0]*di[2]*di[2];
	  bi[1] = pos[1]*di[2]*di[2] - pos[2]*di[1]*di[2] - pos[0]*di[1]*di[0] + pos[1]*di[0]*di[0];
	  bi[2] = pos[2]*di[0]*di[0] - pos[0]*di[2]*di[0] - pos[1]*di[2]*di[1] + pos[2]*di[1]*di[1];
	  
	  for (int ka=0; ka<3; ++ka)
	    {
	      bvec[ka] += wg*bi[ka];
	      for (int kb=0; kb<3; ++kb)
		Mmat[ka][kb] += wg*Mi[ka][kb];
	    }
	}
    }

    std::ofstream of("directions.g2");
    of << "410 1 0 4 155 200 0 255" << std::endl;
    of << dird.size() << std::endl;
    for (size_t ki=0; ki<dird.size(); ++ki)
      of << pp[ki] << " " << pp[ki]+dird[ki] << std::endl;
    
    double det = 0.0;
    int sgn = 1;
    int ka, kb, kc;
    double ax=0.0, ay=0.0, az=0.0;
    // for (ka=0; ka<3; ++ka, sgn*=-1)
    //   {
    // 	kb = (ka+1)%3;
    // 	kc = (kb+1)%3;
    // 	det += sgn*Mmat[0][ka]*(Mmat[1][kb]*Mmat[2][kc]-Mmat[2][kb]*Mmat[1][kc]);
    // 	ax += sgn*bvec[ka]*(Mmat[1][kb]*Mmat[2][kc]-Mmat[2][kb]*Mmat[1][kc]);
    // 	ay += sgn*Mmat[0][ka]*(bvec[kb]*Mmat[2][kc]-Mmat[2][kb]*bvec[kc]);
    // 	az += sgn*Mmat[0][ka]*(Mmat[1][kb]*bvec[kc]-bvec[kb]*Mmat[1][kc]);
    //   }
    // apex = Point(ax/det, ay/det, az/det);

    double det2 = Mmat[0][0]*(Mmat[1][1]*Mmat[2][2] - Mmat[1][2]*Mmat[2][1]) -
      Mmat[0][1]*(Mmat[1][0]*Mmat[2][2] - Mmat[1][2]*Mmat[2][0]) +
      Mmat[0][2]*(Mmat[1][0]*Mmat[2][1] - Mmat[1][1]*Mmat[2][0]);
    double ax2 = bvec[0]*(Mmat[1][1]*Mmat[2][2] - Mmat[1][2]*Mmat[2][1]) -
      bvec[1]*(Mmat[1][0]*Mmat[2][2] - Mmat[1][2]*Mmat[2][0]) +
      bvec[2]*(Mmat[1][0]*Mmat[2][1] - Mmat[1][1]*Mmat[2][0]);
    double ay2 = Mmat[0][0]*(bvec[1]*Mmat[2][2] - bvec[2]*Mmat[2][1]) -
      Mmat[0][1]*(bvec[0]*Mmat[2][2] - bvec[2]*Mmat[2][0]) +
      Mmat[0][2]*(bvec[0]*Mmat[2][1] - bvec[1]*Mmat[2][0]);
    double az2 = Mmat[0][0]*(Mmat[1][1]*bvec[2] - Mmat[1][2]*bvec[1]) -
      Mmat[0][1]*(Mmat[1][0]*bvec[2] - Mmat[1][2]*bvec[0]) +
      Mmat[0][2]*(Mmat[1][0]*bvec[1] - Mmat[1][1]*bvec[0]);
    apex = (fabs(det2) < 1.0e-6) ? Point(0.0, 0.0, 0.0) : Point(ax2/det2, ay2/det2, az2/det2);
    
    // // std::cout << det << " " << det2 << std::endl;
    // for (int ka=0; ka<3; ++ka)
    //   {
    // 	double tmp = 0.0;
    // 	for (kb=0; kb<3; ++kb)
    // 	  tmp += Mmat[ka][kb]*apex[kb];
    // 	std::cout << tmp << " " << bvec[ka] << std::endl;
    //   }

    double nom=0.0, denom=0.0;
    for (size_t ki=0; ki<points.size(); ++ki)
    {
      vector<RevEngPoint*>::iterator start = points[ki].first;
      vector<RevEngPoint*>::iterator end = points[ki].second;
      for (auto it=start; it!=end; ++it)
	{
	  RevEngPoint *pt = *it;
	  Vector3D xyz = pt->getPoint();
	  Point pos(xyz[0], xyz[1], xyz[2]);
	  Point tmp1 = pos - apex;
	  Point tmp2 = tmp1.cross(axis);
	  nom += tmp2.length();
	  denom += tmp1*axis;
	}
    }
    
    double tanphi = nom/denom;
    phi = atan(tanphi);
}

//===========================================================================
void
RevEngUtils::computeSphereProp(vector<pair<vector<RevEngPoint*>::iterator,
			       vector<RevEngPoint*>::iterator> >& points,
			       Point& centre, double& radius)
//===========================================================================
{
  double Amat[4][4];
  double bvec[4];
  // vector<vector<double> > Amat(4, vector<double>(4,0.0));
  // vector<double> bvec(4, 0.0);
  size_t numpt = 0;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      numpt += (points[ki].second - points[ki].first);
    }
  
  vector<vector<double> > A1(numpt);
  vector<double> b1(numpt);
  for (size_t kj=0; kj<numpt; ++kj)
    A1[kj].resize(4);
      
  for (int ka=0; ka<4; ++ka)
    {
      for (int kb=0; kb<4; ++kb)
  	Amat[ka][kb] = 0.0;
      bvec[ka] = 0.0;
    }

  double wgt = 1.0/(double)numpt;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      vector<RevEngPoint*>::iterator start = points[ki].first;
      vector<RevEngPoint*>::iterator end = points[ki].second;
      size_t kj = 0;
      for (auto it=start; it!=end; ++it, ++kj)
	{
	  Vector3D curr = (*it)->getPoint();
	  A1[kj][0] = 2*curr[0];
	  A1[kj][1] = 2*curr[1];
	  A1[kj][2] = 2*curr[2];
	  A1[kj][3] = 1.0;
	  b1[kj] = curr.length2();
	}
    }

  for (int ka=0; ka<4; ++ka)
    {
      for (int kb=0; kb<4; ++kb)
	{
	  for (size_t kr=0; kr<numpt; ++kr)
	    Amat[ka][kb] += wgt*A1[kr][ka]*A1[kr][kb];
	}
      for (size_t kr=0; kr<numpt; ++kr)
	bvec[ka] += wgt*A1[kr][ka]*b1[kr];
    }

  // double detA = 0.0;
  // double bx[3];
  // bx[0] = bx[1] = bx[2] = 0.0;
  // int sgn1 = 1, sgn2 = 1;
  // for (int kb=0; kb<4; ++kb, sgn1*=(-1))
  //   {
  //     for (int kc=0; kc<4; ++kc)
  // 	{
  // 	  if (kc == kb)
  // 	    continue;
  // 	  int ka1 = (kb == 0);
  // 	  int ka2 = 3 - (kb == 3);
  // 	  detA += sgn1*Amat[0][kb]*
  // 	    (sgn2*Amat[1][kc]*(Amat[2][ka1]*Amat[3][ka2]-
  // 			       Amat[3][ka1]*Amat[2][ka2]));
  // 	  bx[0] += sgn1*bvec[kb]*
  // 	    (sgn2*Amat[1][kc]*(Amat[2][ka1]*Amat[3][ka2]-
  // 			       Amat[3][ka1]*Amat[2][ka2]));
  // 	  bx[1] += sgn1*Amat[0][kb]*
  // 	    (sgn2*bvec[kc]*(Amat[2][ka1]*Amat[3][ka2]-
  // 			    Amat[3][ka1]*Amat[2][ka2]));
  // 	  bx[2] += sgn1*Amat[0][kv]*
  // 	    (sgn2*Amat[1][kc]*(bvec[ka1]Amat[3][ka2]-
  // 			       bvec[ka2]*Amat[3][ka1]));
  // 	  bx[3] += sgn1*Amat[0][kb]*
  // 	    sgn2*Amat[1][kc]*((Amat[2][ka1]*bvec[ka2]-Amat[1][ka2]*bvec[ka1]);
  // 	  sgn2 += -1;
  // 	}
  //   }
  // double sx = bx[0]/detA;
  // double sy = bx[1]/detA;
  // double sz = bx[2]/detA;
  // double r2 = bx[3]/detA;
  LUsolveSystem(Amat, 4, &bvec[0]);
  double sx = bvec[0];
  double sy = bvec[1];
  double sz = bvec[2];
  double r2 = bvec[3];

  centre = Point(sx,sy,sz);

  radius = (r2 + sx*sx + sy*sy + sz*sz < 0.0) ? 0.0 : sqrt(r2 + sx*sx + sy*sy + sz*sz);
 }

//===========================================================================
void RevEngUtils::computeCylPosRadius(vector<pair<vector<RevEngPoint*>::iterator,
				      vector<RevEngPoint*>::iterator> >& points,
				      Point& low, Point& high,
				      Point& axis, Point& Cx, Point& Cy,
				      Point& pos, double& radius)
//===========================================================================
{
  double Amat[3][3];
  double bvec[3];
  size_t numpt = 0;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      numpt += (points[ki].second - points[ki].first);
    }
  Point mid = 0.5*(low+high);
  
  vector<vector<double> > A1(numpt);
  vector<double> b1(numpt);
  for (size_t kj=0; kj<numpt; ++kj)
    A1[kj].resize(3);
      
  for (int ka=0; ka<3; ++ka)
    {
      for (int kb=0; kb<3; ++kb)
	Amat[ka][kb] = 0.0;
      bvec[ka] = 0.0;
    }

  for (size_t ki=0; ki<points.size(); ++ki)
    {
      vector<RevEngPoint*>::iterator start = points[ki].first;
      vector<RevEngPoint*>::iterator end = points[ki].second;
      size_t kj = 0;
      for (auto it=start; it!=end; ++it, ++kj)
	{
	  Vector3D curr = (*it)->getPoint();
	  Point curr2(curr[0], curr[1], curr[2]);
	  curr2 -= mid;
	  double pxy[3];
	  double px = curr2*Cx;
	  double py = curr2*Cy;
	  pxy[0] = 2*px;
	  pxy[1] = 2*py;
	  pxy[2] = 1.0;
	  double plen2 = px*px + py*py;
	  A1[kj][0] = 2*px;
	  A1[kj][1] = 2*py;
	  A1[kj][2] = 1.0;
	  b1[kj] = plen2;
	  for (int ka=0; ka<3; ++ka)
	    {
	      for (int kb=0; kb<3; ++kb)
		Amat[ka][kb] += pxy[ka]*pxy[kb];
	      bvec[ka] += pxy[ka]*plen2;
	    }
	}
    }

  double Amat2[3][3];
  double bvec2[3];
  for (int ka=0; ka<3; ++ka)
    {
      for (int kb=0; kb<3; ++kb)
	Amat2[ka][kb] = 0.0;
      bvec2[ka] = 0.0;
    }
  for (int ka=0; ka<3; ++ka)
    {
      for (int kb=0; kb<3; ++kb)
	{
	  for (size_t kr=0; kr<numpt; ++kr)
	    Amat2[ka][kb] += A1[kr][ka]*A1[kr][kb];
	}
      for (size_t kr=0; kr<numpt; ++kr)
	bvec2[ka] += A1[kr][ka]*b1[kr];
    }

  double detA = 0.0;
  double bx[3];
  bx[0] = bx[1] = bx[2] = 0.0;
  int sgn = 1;
  for (int kb=0; kb<3; ++kb, sgn*=(-1))
    {
      int ka1 = (kb == 0);
      int ka2 = 2 - (kb == 2);
      detA += sgn*Amat[0][kb]*(Amat[1][ka1]*Amat[2][ka2]-Amat[2][ka1]*Amat[1][ka2]);
      bx[0] += sgn*bvec[kb]*(Amat[1][ka1]*Amat[2][ka2]-Amat[2][ka1]*Amat[1][ka2]);
      bx[1] += sgn*Amat[0][kb]*(bvec[ka1]*Amat[2][ka2]-bvec[ka2]*Amat[2][ka1]);
      bx[2] += sgn*Amat[0][kb]*(Amat[1][ka1]*bvec[ka2]-Amat[1][ka2]*bvec[ka1]);
    }
  double sx = bx[0]/detA;
  double sy = bx[1]/detA;
  double r2 = bx[2]/detA;

  Point pos2 = sx*Cx + sy*Cy;
  pos2 += mid;

  radius = (r2 + sx*sx + sy*sy < 0.0) ? 0.0 : sqrt(r2 + sx*sx + sy*sy);
  double len = low.dist(high);
  Point vec = mid - pos2;
  Point ax = axis;
  ax.normalize();
  pos = pos2 + (vec*ax)*ax;
 }

 

//===========================================================================
void RevEngUtils::computeCircPosRadius(vector<Point>& points,
				      const Point& axis, const Point& Cx, 
				      const Point& Cy, Point& pos, double& radius)
//===========================================================================
{
  double Amat[3][3];
  double bvec[3];
  size_t numpt = points.size();
  double wgt = 1.0/(double)numpt;
  
  vector<vector<double> > A1(numpt);
  vector<double> b1(numpt);
  for (size_t kj=0; kj<numpt; ++kj)
    A1[kj].resize(3);
      
  for (int ka=0; ka<3; ++ka)
    {
      for (int kb=0; kb<3; ++kb)
	Amat[ka][kb] = 0.0;
      bvec[ka] = 0.0;
    }

  Point mid(0.0, 0.0, 0.0);
  for (size_t ki=0; ki<points.size(); ++ki)
      mid += wgt*points[ki];
  
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      double pxy[3];
      double px = (points[ki]-mid)*Cx;
      double py = (points[ki]-mid)*Cy;
      pxy[0] = 2*px;
      pxy[1] = 2*py;
      pxy[2] = 1.0;
      double plen2 = px*px + py*py;
      A1[ki][0] = 2*px;
      A1[ki][1] = 2*py;
      A1[ki][2] = 1.0;
      b1[ki] = plen2;
      for (int ka=0; ka<3; ++ka)
	{
	  for (int kb=0; kb<3; ++kb)
	    Amat[ka][kb] += pxy[ka]*pxy[kb];
	  bvec[ka] += pxy[ka]*plen2;
	}
    }

  double Amat2[3][3];
  double bvec2[3];
  for (int ka=0; ka<3; ++ka)
    {
      for (int kb=0; kb<3; ++kb)
	Amat2[ka][kb] = 0.0;
      bvec2[ka] = 0.0;
    }
  for (int ka=0; ka<3; ++ka)
    {
      for (int kb=0; kb<3; ++kb)
	{
	  for (size_t kr=0; kr<numpt; ++kr)
	    Amat2[ka][kb] += A1[kr][ka]*A1[kr][kb];
	}
      for (size_t kr=0; kr<numpt; ++kr)
	bvec2[ka] += A1[kr][ka]*b1[kr];
    }

  double detA = 0.0;
  double bx[3];
  bx[0] = bx[1] = bx[2] = 0.0;
  int sgn = 1;
  for (int kb=0; kb<3; ++kb, sgn*=(-1))
    {
      int ka1 = (kb == 0);
      int ka2 = 2 - (kb == 2);
      detA += sgn*Amat[0][kb]*(Amat[1][ka1]*Amat[2][ka2]-Amat[2][ka1]*Amat[1][ka2]);
      bx[0] += sgn*bvec[kb]*(Amat[1][ka1]*Amat[2][ka2]-Amat[2][ka1]*Amat[1][ka2]);
      bx[1] += sgn*Amat[0][kb]*(bvec[ka1]*Amat[2][ka2]-bvec[ka2]*Amat[2][ka1]);
      bx[2] += sgn*Amat[0][kb]*(Amat[1][ka1]*bvec[ka2]-Amat[1][ka2]*bvec[ka1]);
    }
  if (fabs(detA) < 1.0e-12 || std::isnan(detA))
    THROW("Circle computation fail");
  //std::cout << "Circposradius, detA:" << detA << std::endl;
    
  double sx = bx[0]/detA;
  double sy = bx[1]/detA;
  double r2 = bx[2]/detA;

  pos = sx*Cx + sy*Cy;
  pos += mid;
  pos -= ((pos-mid)*axis)*axis;

  radius = (r2 + sx*sx + sy*sy < 0.0) ? 0.0 : sqrt(r2 + sx*sx + sy*sy);
 }



//===========================================================================
void RevEngUtils::computeRadius(vector<Point>& points, Point& axis, 
				Point& Cx, Point& Cy, double& radius)
//===========================================================================
{
  double Amat[3][3];
  double bvec[3];
  size_t numpt = points.size();
  double wgt = 1.0/numpt;
  
  vector<vector<double> > A1(numpt);
  vector<double> b1(numpt);
  for (size_t kj=0; kj<numpt; ++kj)
    A1[kj].resize(3);
      
  for (int ka=0; ka<3; ++ka)
    {
      for (int kb=0; kb<3; ++kb)
	Amat[ka][kb] = 0.0;
      bvec[ka] = 0.0;
    }

  Point mid(0.0, 0.0, 0.0);
  for (size_t ki=0; ki<points.size(); ++ki)
      mid += wgt*points[ki];
  
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      double pxy[3];
      double px = (points[ki]-mid)*Cx;
      double py = (points[ki]-mid)*Cy;
      pxy[0] = 2*px;
      pxy[1] = 2*py;
      pxy[2] = 1.0;
      double plen2 = px*px + py*py;
      A1[ki][0] = 2*px;
      A1[ki][1] = 2*py;
      A1[ki][2] = 1.0;
      b1[ki] = plen2;
      for (int ka=0; ka<3; ++ka)
	{
	  for (int kb=0; kb<3; ++kb)
	    Amat[ka][kb] += pxy[ka]*pxy[kb];
	  bvec[ka] += pxy[ka]*plen2;
	}
    }

  double Amat2[3][3];
  double bvec2[3];
  for (int ka=0; ka<3; ++ka)
    {
      for (int kb=0; kb<3; ++kb)
	Amat2[ka][kb] = 0.0;
      bvec2[ka] = 0.0;
    }
  for (int ka=0; ka<3; ++ka)
    {
      for (int kb=0; kb<3; ++kb)
	{
	  for (size_t kr=0; kr<numpt; ++kr)
	    Amat2[ka][kb] += A1[kr][ka]*A1[kr][kb];
	}
      for (size_t kr=0; kr<numpt; ++kr)
	bvec2[ka] += A1[kr][ka]*b1[kr];
    }

  double detA = 0.0;
  double bx[3];
  bx[0] = bx[1] = bx[2] = 0.0;
  int sgn = 1;
  for (int kb=0; kb<3; ++kb, sgn*=(-1))
    {
      int ka1 = (kb == 0);
      int ka2 = 2 - (kb == 2);
      detA += sgn*Amat[0][kb]*(Amat[1][ka1]*Amat[2][ka2]-Amat[2][ka1]*Amat[1][ka2]);
      bx[0] += sgn*bvec[kb]*(Amat[1][ka1]*Amat[2][ka2]-Amat[2][ka1]*Amat[1][ka2]);
      bx[1] += sgn*Amat[0][kb]*(bvec[ka1]*Amat[2][ka2]-bvec[ka2]*Amat[2][ka1]);
      bx[2] += sgn*Amat[0][kb]*(Amat[1][ka1]*bvec[ka2]-Amat[1][ka2]*bvec[ka1]);
    }
  if (fabs(detA) < 1.0e-12)
    THROW("Circle with infinite radius");
  
  double sx = bx[0]/detA;
  double sy = bx[1]/detA;
  double r2 = bx[2]/detA;

  radius = (r2 + sx*sx + sy*sy < 0.0) ? 0.0 : sqrt(r2 + sx*sx + sy*sy);
 }

//===========================================================================
void RevEngUtils::computePlane(vector<pair<vector<RevEngPoint*>::iterator,
			       vector<RevEngPoint*>::iterator> >& points,
			       Point& pos, Point& norm)
//===========================================================================
{
  double Cmat[4][4];
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      vector<RevEngPoint*>::iterator start = points[ki].first;
      vector<RevEngPoint*>::iterator end = points[ki].second;
      for (int ka=0; ka<4; ++ka)
	for (int kb=0; kb<4; ++kb)
	  {
	    Cmat[ka][kb] = 0.0;
	    for (auto it=start; it!=end; ++it)
	      {
		RevEngPoint *pt = *it;
		Vector3D xyz = pt->getPoint();
		double tmp[4] = {xyz[0], xyz[1], xyz[2], 1};
		Cmat[ka][kb] += tmp[ka]*tmp[kb];
	      }
	  }
    }
  // Compute singular values
  NEWMAT::Matrix nmat;
  nmat.ReSize(4, 4);
  for (int ka = 0; ka < 4; ++ka) {
    for (int kb = 0; kb < 4; ++kb) {
      nmat.element(ka, kb) = Cmat[ka][kb];
    }
  }
      
  static NEWMAT::DiagonalMatrix diag;
  static NEWMAT::Matrix V;
  try {
    NEWMAT::SVD(nmat, diag, nmat, V);
  } catch(...) {
    //std::cout << "Exception in SVD" << std::endl;
    exit(-1);
  }

  double sigma[4];
  double coefs[4];
  int ixv = 3;
  for (int ka=0; ka<4; ++ka)
    {
      sigma[ka] = diag.element(ka,ka);
      coefs[ka] = V.element(ka, ixv);
    }


  int num = 0;
  double maxd =0.0, avd = 0.0;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      vector<RevEngPoint*>::iterator start = points[ki].first;
      vector<RevEngPoint*>::iterator end = points[ki].second;
      for (auto it=start; it!=end; ++it)
	{
	  Vector3D xyz = (*it)->getPoint();
	  double dist = coefs[0]*xyz[0] + coefs[1]*xyz[1] +
	    coefs[2]*xyz[2] + coefs[3];
	  maxd = std::max(maxd, fabs(dist));
	  avd += dist;
	  int stop_break = 1;
	}
    }
  avd /= (double)num;

  int ix = -1;
  double mm = 0.0;
  for (int ka=0; ka<3; ++ka)
    if (fabs(coefs[ka]) > mm)
      {
	mm = fabs(coefs[ka]);
	ix = ka;
      }

  double t1 = 0.0, t2 = 0.0, t3 = 0.0;
  Vector3D pt1 = (*points[0].first)->getPoint();
  vector<RevEngPoint*>::iterator it = points[0].second;
  it--;
  Vector3D pt2 = (*it)->getPoint();
  for (int ka=0; ka<3; ++ka)
    {
      if (ka == ix)
	continue;
      t1 += coefs[ka]*pos[ka];
      t2 += coefs[ka]*pt1[ka];
      t3 += coefs[ka]*pt2[ka];
    }
  pos[ix] = -(t1 + coefs[3])/coefs[ix];
  pt1[ix] = -(t2 + coefs[3])/coefs[ix];
  pt2[ix] = -(t3 + coefs[3])/coefs[ix];
  norm = Point(coefs[0], coefs[1], coefs[2]);
  norm.normalize();

  Point pt1_2(pt1[0], pt1[1], pt1[2]);
  Point pt2_2(pt2[0], pt2[1], pt2[2]);
  Point vec1 = pt1_2 - pos;
  Point vec2 = pt2_2 - pos;
  Point norm2 = vec1.cross(vec2);
  norm2.normalize();

  int stop_break = 1;
}

//===========================================================================
void RevEngUtils::projectToPlane(vector<RevEngPoint*>& points,
				 Point& axis, Point& mid, std::vector<Point>& projected,
				 double& maxdist, double& avdist)
//===========================================================================
{
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > points2;
  points2.push_back(std::make_pair(points.begin(), points.end()));
  projectToPlane(points2, axis, mid, projected, maxdist, avdist);
}

//===========================================================================
void RevEngUtils::projectToPlane(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
				 std::vector<RevEngPoint*>::iterator> >& points,
				 Point& axis, Point& mid, std::vector<Point>& projected,
				 double& maxdist, double& avdist)
//===========================================================================
{
  maxdist = 0.0;
  avdist = 0.0;
  int nmb = 0;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      vector<RevEngPoint*>::iterator start = points[ki].first;
      vector<RevEngPoint*>::iterator end = points[ki].second;
      for (auto it=start; it!=end; ++it)
	{
	  RevEngPoint *pt = *it;
	  Vector3D pnt = pt->getPoint();
	  Point curr(pnt[0], pnt[1], pnt[2]);
	  Point curr2 = curr - mid;
	  curr2 -= ((curr2*axis)*axis);
	  curr2 += mid;
	  projected.push_back(curr2);
	  double dist = curr.dist(curr2);
	  maxdist = std::max(maxdist, dist);
	  avdist += dist;
	  nmb++;
	}
    }
  avdist /= (double)nmb;
}

//===========================================================================
void RevEngUtils::rotateToPlane(vector<pair<vector<RevEngPoint*>::iterator,
				vector<RevEngPoint*>::iterator> >& points,
				Point& xvec, Point& axis, Point& mid,
				vector<Point>& rotated)
//===========================================================================
{
  Point yvec = xvec.cross(axis);
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      vector<RevEngPoint*>::iterator start = points[ki].first;
      vector<RevEngPoint*>::iterator end = points[ki].second;
      for (auto it=start; it!=end; ++it)
	{
	  RevEngPoint *pt = *it;
	  Vector3D pnt = pt->getPoint();
	  Point curr(pnt[0], pnt[1], pnt[2]);
	  curr -= mid;
	  pnt = Vector3D(curr[0], curr[1], curr[2]);
	  curr -= (curr*axis)*axis;
	  double angle = curr.angle(xvec);
	  if (curr*yvec < 0.0)
	    angle *= -1.0;
	  Matrix3D mat;
	  mat.setToRotation(angle, axis[0], axis[1], axis[2]);
	  Vector3D pnt2 = mat*pnt;
	  rotated.push_back(mid + Point(pnt2[0], pnt2[1], pnt2[2]));
	}
    }
 }

//===========================================================================
void RevEngUtils::rotateToPlane(vector<Point>& points,
				Point& xvec, Point& axis, Point& mid,
				vector<Point>& rotated)
//===========================================================================
{
  rotated.resize(points.size());
  Point yvec = xvec.cross(axis);
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      Point curr = points[ki];
      curr -= mid;
      Vector3D pnt(curr[0], curr[1], curr[2]);
      curr -= (curr*axis)*axis;
      double angle = curr.angle(xvec);
      if (curr*yvec < 0.0)
	angle *= -1.0;
      Matrix3D mat;
      mat.setToRotation(angle, axis[0], axis[1], axis[2]);
      Vector3D pnt2 = mat*pnt;
      rotated[ki] = mid + Point(pnt2[0], pnt2[1], pnt2[2]);
    }
}

//===========================================================================
void RevEngUtils::distToSurf(vector<RevEngPoint*>::iterator start,
			     vector<RevEngPoint*>::iterator end,
			     shared_ptr<ParamSurface> surf, double tol,
			     double& maxdist, double& avdist, int& num_inside,
			     vector<RevEngPoint*>& in, vector<RevEngPoint*>& out,
			     vector<double>& parvals,
			     vector<pair<double,double> >& distang,
			     double angtol)
//===========================================================================
{
  double eps = 1.0e-6;
  maxdist = avdist = 0.0;
  num_inside = 0;
  int num = 0;
  double *seed = 0;
  double seed2[2];
  Point prev;
  double fac = 100.0;
  for (auto it=start; it!=end; ++it)
    {
      Vector3D xyz = (*it)->getPoint();
      Point pnt(xyz[0], xyz[1], xyz[2]);
      if (prev.dimension() == pnt.dimension() && prev.dist(pnt) < fac*tol)
	seed = seed2;
      
      double upar, vpar, dist;
      Point close;
      Point norm1, norm2;
      surf->closestPoint(pnt, upar, vpar, close, dist, eps, 0, seed);
      parvals.push_back(upar);
      parvals.push_back(vpar);
      surf->normal(norm1, upar, vpar);
      norm2 = (*it)->getMongeNormal();
      maxdist = std::max(maxdist, dist);
      avdist += dist;
      double ang = norm1.angle(norm2);
      ang = std::min(M_PI-ang, ang);
      distang.push_back(std::make_pair(dist, ang));
      if (dist <= tol && (angtol < 0.0 || ang <= angtol))
	{
	  in.push_back(*it);
	  ++num_inside;
	}
	else
	  {
	  out.push_back(*it);
	  int stop_break = 1;
	}
      seed2[0] = upar;
      seed2[1] = vpar;
      prev = pnt;
      ++num;
    }
  avdist /= (double)num;
}

//===========================================================================
void RevEngUtils::distToSurf(vector<Point>& points,
			     shared_ptr<ParamSurface> surf, double tol,
			     double& maxdist, double& avdist, int& num_inside,
			     vector<double>& distance)
//===========================================================================
{
  double eps = 1.0e-6;
  maxdist = avdist = 0.0;
  num_inside = 0;
  int num = 0;
  distance.resize(points.size());
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      double upar, vpar, dist;
      Point close;
      surf->closestPoint(points[ki], upar, vpar, close, dist, eps);
      maxdist = std::max(maxdist, dist);
      avdist += dist;
      distance[ki] = dist;
      if (dist <= tol)
	++num_inside;
      else
	{
	  int stop_break = 1;
	}
      ++num;
    }
  avdist /= (double)num;
}

//===========================================================================
void RevEngUtils::distToCurve(vector<Point>& points,
			     shared_ptr<ParamCurve> curve, double tol,
			     double& maxdist, double& avdist, int& num_inside)
//===========================================================================
{
  double eps = 1.0e-6;
  maxdist = avdist = 0.0;
  num_inside = 0;
  int num = 0;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      double tpar, dist;
      Point close;
      curve->closestPoint(points[ki], curve->startparam(), curve->endparam(),
			  tpar, close, dist);
      maxdist = std::max(maxdist, dist);
      avdist += dist;
      if (dist <= tol)
	++num_inside;
      else
	{
	  int stop_break = 1;
	}
      ++num;
    }
  avdist /= (double)num;
}

//===========================================================================
shared_ptr<ParamSurface> RevEngUtils::doMergePlanes(vector<pair<vector<RevEngPoint*>::iterator,
						    vector<RevEngPoint*>::iterator> > points,
						    const BoundingBox& bbox,
						    vector<int>& nmbpts,
						    bool set_bound)
//===========================================================================
{
  int totnmb = 0;
  for (size_t kh=0; kh<nmbpts.size(); ++kh)
    totnmb += nmbpts[kh];
  
  Point pos(0.0, 0.0, 0.0);
  Point norm(0.0, 0.0, 0.0);
  vector<RevEngPoint*> all_pts;
  all_pts.reserve(totnmb);
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      double wgt = 1.0/(double)totnmb;

      for (auto it=points[ki].first; it!=points[ki].second; ++it)
	{
	  Point curr = (*it)->getMongeNormal();
	  Vector3D xyz = (*it)->getPoint();
	  pos +=  wgt*Point(xyz[0], xyz[1], xyz[2]);
	  norm += wgt*curr;
	  all_pts.push_back(*it);
	}
    }
  
  // Perform approximation with combined point set
  shared_ptr<ImplicitApprox> impl(new ImplicitApprox());
  impl->approx(points, 1);
  Point pos2, normal2;
  impl->projectPoint(pos, norm, pos2, normal2);
  // std::ofstream outviz("implsf_merge.g2");
  // impl->visualize(all_pts, outviz);
 
  shared_ptr<Plane> surf(new Plane(pos2, normal2));
  Point low = bbox.low();
  Point high = bbox.high();
  if (set_bound)
    {
      double len = low.dist(high);
      surf->setParameterBounds(-0.5*len, -0.5*len, 0.5*len, 0.5*len);
    }

  return surf;
}

//===========================================================================
shared_ptr<ParamSurface> RevEngUtils::doMergeCylinders(vector<pair<vector<RevEngPoint*>::iterator,
						    vector<RevEngPoint*>::iterator> > points,
						    const BoundingBox& bbox,
						    vector<int>& nmbpts,
						    bool set_bound)
//===========================================================================
{
  // Estimate cylinder axis
  Point axis, Cx, Cy;
  RevEngUtils::computeAxis(points, axis, Cx, Cy);

  // Estimate radius and point on axis
  double rad;
  Point pnt;
  Point low = bbox.low();
  Point high = bbox.high();
  RevEngUtils::computeCylPosRadius(points, low, high,
				   axis, Cx, Cy, pnt, rad);
  shared_ptr<Cylinder> surf(new Cylinder(rad, pnt, axis, Cy));
  if (set_bound)
    {
      double len = low.dist(high);
      surf->setParamBoundsV(-len, len);
    }

  return surf;
}

//===========================================================================
shared_ptr<ParamSurface> RevEngUtils::doMergeSpheres(vector<pair<vector<RevEngPoint*>::iterator,
						    vector<RevEngPoint*>::iterator> > points,
						    const BoundingBox& bbox,
						     vector<int>& nmbpts, Point& normal)
//===========================================================================
{
  Point centre;
  double radius;
  try {
    RevEngUtils::computeSphereProp(points, centre, radius);
  }
  catch (...)
    {
      shared_ptr<Sphere> dummy;
      return dummy;
    }

  int totnmb = 0;
  for (size_t kh=0; kh<nmbpts.size(); ++kh)
    totnmb += nmbpts[kh];
  vector<Point> pnts;
  pnts.reserve(totnmb);
  for (size_t ki=0; ki<points.size(); ++ki)
    for (auto it=points[ki].first; it!=points[ki].second; ++it)
      {
	Vector3D xyz = (*it)->getPoint();
	pnts.push_back(Point(xyz[0], xyz[1], xyz[2]));
      }

  double eigenvec[3][3];
  double lambda[3];
  Point eigen1, eigen2, eigen3;
  RevEngUtils::principalAnalysis(pnts[0], pnts, lambda, eigenvec);
  Point z_axis = Point(eigenvec[1][0], eigenvec[1][1], eigenvec[1][2]);
  Point x_axis = normal.cross(z_axis);

  shared_ptr<Sphere> sph(new Sphere(radius, centre, z_axis, normal));

  return sph;
}

//===========================================================================
shared_ptr<ParamSurface> RevEngUtils::doMergeTorus(vector<pair<vector<RevEngPoint*>::iterator,
						   vector<RevEngPoint*>::iterator> > points,
						   const BoundingBox& bbox,
						   vector<int>& nmbpts)
//===========================================================================
{
  shared_ptr<ParamSurface> dummy;
  
  int totnmb = 0;
  for (size_t kh=0; kh<nmbpts.size(); ++kh)
    totnmb += nmbpts[kh];
  
  // Compute mean curvature and initial point in plane
  double k2mean = 0.0;
  double wgt = 1.0/(double)totnmb;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      for (auto it=points[ki].first; it!=points[ki].second; ++it)
	{
	  double kmax = (*it)->maxPrincipalCurvature();
	  k2mean += wgt*kmax;
	}
    }
  double rd = 1.0/k2mean;
  
  vector<Point> centr(totnmb);
  Point mid(0.0, 0.0, 0.0);
  size_t kr = 0;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      for (auto it=points[ki].first; it!=points[ki].second; ++it, ++kr)
	{
	  Point norm = (*it)->getMongeNormal();
	  Vector3D xyz = (*it)->getPoint();
	  Point xyz2(xyz[0], xyz[1], xyz[2]);
	  centr[kr] = xyz2 + rd*norm;
	  mid += wgt*centr[kr];
	}
    }
  
  shared_ptr<ImplicitApprox> impl(new ImplicitApprox());
  impl->approxPoints(centr, 1);

  double val;
  Point grad;
  impl->evaluate(mid, val, grad);
  grad.normalize_checked();
  Point pos, normal;
  impl->projectPoint(mid, grad, pos, normal);
  double eps1 = 1.0e-8;
  if (normal.length() < eps1)
    return dummy;
  
  Point Cx = centr[0] - mid;
  Cx -= (Cx*normal)*normal;
  Cx.normalize();
  Point Cy = Cx.cross(normal);
  
  double rad;
  Point pnt;
  RevEngUtils::computeCircPosRadius(centr, normal, Cx, Cy, pnt, rad);
  pnt -= ((pnt - pos)*normal)*normal;

  vector<Point> rotated;
  RevEngUtils::rotateToPlane(points, Cx, normal, pnt, rotated);
  Point cpos;
  double crad;
  RevEngUtils::computeCircPosRadius(rotated, Cy, Cx, normal, cpos, crad);
  Point cvec = cpos - pnt;
  double R1 = (cvec - (cvec*normal)*normal).length();
  double R2 = (cvec*normal)*normal.length();
 
  shared_ptr<Torus> surf(new Torus(R1, crad, pnt+R2*normal, normal, Cy));

  return surf;
}


//===========================================================================
shared_ptr<SplineCurve> RevEngUtils::midCurve(shared_ptr<SplineCurve>& cv1,
					      shared_ptr<SplineCurve>& cv2)
//===========================================================================
{
  shared_ptr<SplineCurve> spl1(cv1->clone());
  shared_ptr<SplineCurve> spl2(cv2->clone());

  // Check orientation
  Point pt1 = spl1->ParamCurve::point(spl1->startparam());
  Point pt2 = spl1->ParamCurve::point(spl1->endparam());
  Point pt3 = spl2->ParamCurve::point(spl2->startparam());
  Point pt4 = spl2->ParamCurve::point(spl2->endparam());
  double len1 = pt1.dist(pt3);
  double len2 = pt1.dist(pt4);
  if (len2 < len1)
    spl2->reverseParameterDirection();

  // Ensure same spline room
  spl2->setParameterInterval(spl1->startparam(), spl1->endparam());

  double tol = 1.0e-4;
  vector<shared_ptr<SplineCurve> > curves(2);
  curves[0] = spl1;
  curves[1] = spl2;
  GeometryTools::unifyCurveSplineSpace(curves, tol);

  shared_ptr<SplineCurve> midcv = GeometryTools::curveSum(*curves[0], 0.5,
							  *curves[1], 0.5);
  return midcv;
}


//===========================================================================
shared_ptr<SplineCurve> RevEngUtils::createCurve(vector<RevEngPoint*>& points,
						 int degree, double tol, int maxiter)
//===========================================================================
{
  shared_ptr<SplineCurve> cv;
  if (points.size() < 2)
    return cv;
  
  // Parameterize curves and fetch data points
  vector<double> param(points.size(), 0.0);
  vector<double> pts;
  pts.reserve(3*points.size());
  Vector3D prev = points[0]->getPoint();
  double tmp = 0.0;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      Vector3D xyz = points[ki]->getPoint();
      pts.insert(pts.end(), xyz.begin(), xyz.end());
      param[ki] = tmp + prev.dist(xyz);
      prev = xyz;
      tmp = param[ki];
    }

  double smoothwgt = 0.1;
  ApproxCurve approx(pts, param, 3, tol, degree+1, degree+1);
  approx.setSmooth(smoothwgt);

  double maxdist, avdist;
  cv = approx.getApproxCurve(maxdist, avdist, maxiter);

  return cv;
}
