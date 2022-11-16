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

#include "GoTools/compositemodel/ImplicitApprox.h"
#include "GoTools/compositemodel/RevEngPoint.h"
#include "GoTools/implicitization/ImplicitizePointCloudAlgo.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/LineCloud.h"
#include "GoTools/geometry/SweepSurfaceCreator.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/ftPlane.h"
#include "GoTools/compositemodel/ftCurve.h"
#include "GoTools/implicitization/BernsteinPoly.h"
#include "sisl.h"
#include <fstream>

using namespace Go;
using std::vector;

//===========================================================================
ImplicitApprox::ImplicitApprox()
  : eps_(1.0e-12)
//===========================================================================
{
}

//===========================================================================
ImplicitApprox::~ImplicitApprox()
//===========================================================================
{
}

//===========================================================================
void ImplicitApprox::approx(vector<RevEngPoint*> points, int degree)
//===========================================================================
{
  // Extract xyz values
  vector<Vector3D> xyz(points.size());
  for (size_t ki=0; ki<points.size(); ++ki)
    xyz[ki] = points[ki]->getPoint();

  PointCloud3D pointset(xyz);

  // Implicitize
  degree_ = degree;
  ImplicitizePointCloudAlgo implicitize(pointset, degree);
  implicitize.perform();
  
  // Get result
  implicitize.getResultData(implicit_, bc_, sigma_min_);
  
  // Differentiate
  Vector4D bdir1(1.0, 0.0, 0.0, 0.0);
  Vector4D bdir2(0.0, 1.0, 0.0, 0.0);
  Vector4D bdir3(0.0, 0.0, 1.0, 0.0);
  Vector4D bdir4(0.0, 0.0, 0.0, 1.0);
  implicit_.deriv(1, bdir1, deriv1_);
  implicit_.deriv(1, bdir2, deriv2_);
  implicit_.deriv(1, bdir3, deriv3_);
  implicit_.deriv(1, bdir4, deriv4_);

}

//===========================================================================
void ImplicitApprox::approxPoints(vector<Point> points, int degree)
//===========================================================================
{
  // Extract xyz values
  vector<Vector3D> xyz(points.size());
  for (size_t ki=0; ki<points.size(); ++ki)
    xyz[ki] = Vector3D(points[ki].begin());

  PointCloud3D pointset(xyz);

  // Implicitize
  degree_ = degree;
  ImplicitizePointCloudAlgo implicitize(pointset, degree);
  implicitize.perform();
  
  // Get result
  implicitize.getResultData(implicit_, bc_, sigma_min_);
  
  // Differentiate
  Vector4D bdir1(1.0, 0.0, 0.0, 0.0);
  Vector4D bdir2(0.0, 1.0, 0.0, 0.0);
  Vector4D bdir3(0.0, 0.0, 1.0, 0.0);
  Vector4D bdir4(0.0, 0.0, 0.0, 1.0);
  implicit_.deriv(1, bdir1, deriv1_);
  implicit_.deriv(1, bdir2, deriv2_);
  implicit_.deriv(1, bdir3, deriv3_);
  implicit_.deriv(1, bdir4, deriv4_);

}

//===========================================================================
void ImplicitApprox::approx(vector<pair<vector<RevEngPoint*>::iterator,
			    vector<RevEngPoint*>::iterator> >& points,
			    int degree)
//===========================================================================
{
  // Extract xyz values
  vector<Vector3D> xyz;
  for (size_t kj=0; kj<points.size(); ++kj)
    {
      for (auto it=points[kj].first; it!=points[kj].second; ++it)
	{
	  Vector3D curr = (*it)->getPoint();
	  xyz.push_back(curr);
	}
    }
  PointCloud3D pointset(xyz);

  // Implicitize
  degree_ = degree;
  ImplicitizePointCloudAlgo implicitize(pointset, degree);
  implicitize.perform();
  
  // Get result
  implicitize.getResultData(implicit_, bc_, sigma_min_);
  
  // Differentiate
  Vector4D bdir1(1.0, 0.0, 0.0, 0.0);
  Vector4D bdir2(0.0, 1.0, 0.0, 0.0);
  Vector4D bdir3(0.0, 0.0, 1.0, 0.0);
  Vector4D bdir4(0.0, 0.0, 0.0, 1.0);
  implicit_.deriv(1, bdir1, deriv1_);
  implicit_.deriv(1, bdir2, deriv2_);
  implicit_.deriv(1, bdir3, deriv3_);
  implicit_.deriv(1, bdir4, deriv4_);

}

//===========================================================================
void ImplicitApprox::evaluate(Point& pt, double& val, Point& grad)
//===========================================================================
{
  Vector3D xyz(pt[0], pt[1], pt[2]);
  Vector4D bary = bc_.cartToBary(xyz);
  val = implicit_(bary);
  double d1 = deriv1_(bary);
  double d2 = deriv2_(bary);
  double d3 = deriv3_(bary);
  double d4 = deriv4_(bary);
  Vector4D dv(d1,d2,d3,d4);
  Vector4D bary2 = bary+dv;
  Vector3D pt2 = bc_.baryToCart(bary2);
  Vector3D grad2 = pt2 - xyz;
  grad = Point(grad2[0], grad2[1], grad2[2]);
}

//===========================================================================
double ImplicitApprox::estimateDist(RevEngPoint* pt)
//===========================================================================
{
  Vector3D xyz = pt->getPoint();
  Vector4D bary = bc_.cartToBary(xyz);
  double dist0 = implicit_(bary);
  double d1 = deriv1_(bary);
  double d2 = deriv2_(bary);
  double d3 = deriv3_(bary);
  double d4 = deriv4_(bary);
  Vector4D dv(d1,d2,d3,d4);
  Vector4D bary2 = bary+dv;
  Vector3D pt2 = bc_.baryToCart(bary2);
  Vector3D grad = pt2 - xyz;
  double len = grad.length();
  double dist = (len > eps_) ? dist0/len : dist0;

  Point norm = pt->getMongeNormal();
  Vector3D norm2(norm[0], norm[1], norm[2]);
  norm2 *= 100;
  Vector3D xyz2 = xyz + norm2;
  Vector3D xyz3 = xyz - norm2;
  Vector4D bary3 = bc_.cartToBary(xyz2);
  Vector4D bary4 = bc_.cartToBary(xyz3);
  BernsteinPoly line = implicit_.pickLine(bary3, bary4);

  // Compute zeroes of bernstein polynomial
  // First make sisl curve
  int ik = degree_ + 1;
  vector<double> et(2*ik, 0.0);  // Knot vector of line curve
  for (int ki=0; ki<ik; ++ki)
    et[ik+ki] = 1.0;
  vector<double> ecoef(line.coefsBegin(), line.coefsEnd());
  SISLCurve *qc = newCurve(ik, ik, &et[0], &ecoef[0], 1, 1, 1);
  double zero = 0.0;
  
  // Intersect
  double eps = 1.0e-6;
  int kstat = 0;
  int kcrv=0, kpt=0;
  double *epar = 0;
  SISLIntcurve **intcv = 0;
  if (qc)
    s1871(qc, &zero, 1, eps, &kpt, &epar, &kcrv, &intcv, &kstat);
  if (qc)
    freeCurve(qc);
  
  // Compute cartesian points and curves associated with intersections
  double dd = std::numeric_limits<double>::max();
  for (int kr=0; kr<kpt; ++kr)
    {
      Vector4D barypt = (1.0 - epar[kr])*bary3 + epar[kr]*bary4;
      Vector3D pos = bc_.baryToCart(barypt);
      double dd2 = xyz.dist(pos);
      if (dd2 < dd)
	dd = dd2;
    }
  return dd; //dist;
}

//===========================================================================
void ImplicitApprox::projectPoint(Point point, Point dir,
				  Point& projpos, Point& normal)
//===========================================================================
{
  double len = 100.0;
  dir.normalize();
  
  Point xdir(1.0, 0.0, 0.0);
  Point ydir(0.0, 1.0, 0.0);
  Point zdir(0.0, 0.0, 1.0);
  double a1 = xdir.angle(dir);
  double a2 = ydir.angle(dir);
  double a3 = zdir.angle(dir);
  Point dir2;
  if (a1 > std::min(a2, a3))
    dir2 = xdir;
  else if (a2 > a3)
    dir2 = ydir;
  else
    dir2 = zdir;
  Point dir3 = dir%dir2;
  dir2 = dir%dir3;
  dir2.normalize();
  dir3.normalize();
  Point points[3];
  points[0] = point;
  points[1] = point + dir2;
  points[2] = point + dir3;

  Vector3D proj[3];
  for (int ka=0; ka<3; ++ka)
    {
      Vector3D xyz(points[ka].begin());
      Point p1 = points[ka] - len*dir;
      Point p2 = points[ka] + len*dir;

      Vector3D cart1(p1.begin());
      Vector3D cart2(p2.begin());
      Vector4D bary1 = bc_.cartToBary(Vector3D(cart1[0], cart1[1], cart1[2]));
      Vector4D bary2 = bc_.cartToBary(Vector3D(cart2[0], cart2[1], cart2[2]));

      // Pick line
      BernsteinPoly line = implicit_.pickLine(bary1, bary2);

      // Compute zeroes of bernstein polynomial
      // First make sisl curve
      int ik = degree_ + 1;
      vector<double> et(2*ik, 0.0);  // Knot vector of line curve
      for (int ki=0; ki<ik; ++ki)
	et[ik+ki] = 1.0;
      vector<double> ecoef(line.coefsBegin(), line.coefsEnd());
      SISLCurve *qc = newCurve(ik, ik, &et[0], &ecoef[0], 1, 1, 1);
      double zero = 0.0;

      // Intersect
      double eps = 1.0e-6;
      int kstat = 0;
      int kcrv=0, kpt=0;
      double *epar = 0;
      SISLIntcurve **intcv = 0;
      if (qc)
	s1871(qc, &zero, 1, eps, &kpt, &epar, &kcrv, &intcv, &kstat);
      if (qc)
	freeCurve(qc);

      // Compute cartesian points and curves associated with intersections
      double mindist = std::numeric_limits<double>::max();
      for (int kr=0; kr<kpt; ++kr)
	{
	  Vector4D barypt = (1.0 - epar[kr])*bary1 + epar[kr]*bary2;
		
	  Vector3D pos = bc_.baryToCart(barypt);
	  double dist = pos.dist(xyz);
	  if (dist < mindist)
	    {
	      mindist = dist;
	      proj[ka] = pos;
	    }
	}
    }
  projpos = Point(proj[0][0], proj[0][1], proj[0][2]);
  Point pt2(proj[1][0], proj[1][1], proj[1][2]);
  Point pt3(proj[2][0], proj[2][1], proj[2][2]);

  Point vec1 = pt2 - projpos;
  Point vec2 = pt3 - projpos;
  normal = vec1.cross(vec2);
  normal.normalize_checked();
}

//===========================================================================
void ImplicitApprox::visualize(vector<RevEngPoint*> points, std::ostream& os)
//===========================================================================
{
  // View direction
  Point dir = points[0]->getMongeNormal();
  dir.normalize();

  BoundingBox bb(3);
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      Vector3D xyz = points[ki]->getPoint();
      bb.addUnionWith(Point(xyz[0], xyz[1], xyz[2]));
    }
  Point low = bb.low();
  Point high = bb.high();
  Point bmid = 0.5*(low + high);
  Point diag = high - low;
  double diaglen = diag.length();

  double gap = 1.0e-6;
  Point xdir(1.0, 0.0, 0.0);
  Point ydir(0.0, 1.0, 0.0);
  Point zdir(0.0, 0.0, 1.0);
  CompositeModelFactory factory(gap, gap, 10.0*gap, 0.01, 0.05);
  shared_ptr<SurfaceModel> boxmod(factory.createFromBox(low-2.0*diag, xdir, ydir, 
							5*diag[0], 5*diag[1], 5*diag[2]));
    
    // Find the coordinate direction with the largest angle with the view direction
    double a1 = xdir.angle(dir);
    double a2 = ydir.angle(dir);
    double a3 = zdir.angle(dir);
    Point dir2;
    if (a1 > std::min(a2, a3))
      dir2 = xdir;
    else if (a2 > a3)
      dir2 = ydir;
    else
      dir2 = zdir;
    Point dir3 = dir%dir2;
    dir2 = dir%dir3;
    if (dir2*(high-low) < 0.0)
      dir2 *= -1.0;
    if (dir3*(high-low) < 0.0)
      dir3 *= -1.0;
    dir2.normalize();
    dir3.normalize();
    double len = low.dist(high);
    int nmb_sample = 100;
    shared_ptr<SplineCurve> cv1(new SplineCurve(bmid-len*dir2, 0.0, bmid+len*dir2, 1.0));
    shared_ptr<SplineCurve> cv2(new SplineCurve(bmid-len*dir3, 0.0, bmid+len*dir3, 1.0));
    SweepSurfaceCreator sweep;
    shared_ptr<SplineSurface> ssf(sweep.linearSweptSurface(*cv1, *cv2, bmid));
    double del = 1.0/(double)(nmb_sample-1);
    double p1, p2;
    int ki, kj, kr;
    int ik = degree_ + 1;
    vector<double> et(2*ik, 0.0);  // Knot vector of line curve
    for (int ki=0; ki<ik; ++ki)
      et[ik+ki] = 1.0;

    vector<double> sfpoints;
    vector<double> vecs;
    vector<double> linesegs;
    vector<double> der;
    vector<double> der2;
    vector<double> lineder;
    // Evaluate line
    vector<double> tmpline;
    for (kj=0, p2=0.0; kj<nmb_sample; ++kj, p2+=del)
      {
	for (ki=0, p1=0.0; ki<nmb_sample; ++ki, p1+=del)
	  {
	    // Compute barysentric coordinates of end points of line
	    // First cartesian
	    Point sfpos = ssf->ParamSurface::point(p1,p2);
	    Point cart1 = sfpos + len*dir;
	    Point cart2 = sfpos - len*dir;
	    tmpline.insert(tmpline.end(), cart1.begin(), cart1.end());
	    tmpline.insert(tmpline.end(), cart2.begin(), cart2.end());

	    Vector4D bary1 = bc_.cartToBary(Vector3D(cart1[0], cart1[1], cart1[2]));
	    Vector4D bary2 = bc_.cartToBary(Vector3D(cart2[0], cart2[1], cart2[2]));

	    Vector3D tp1 = bc_.baryToCart(bary1);
	    Vector3D tp2 = bc_.baryToCart(bary2);
	    
	    // Pick line
	    BernsteinPoly line = implicit_.pickLine(bary1, bary2);

	    // Compute zeroes of bernstein polynomial
	    // First make sisl curve
	    vector<double> ecoef(line.coefsBegin(), line.coefsEnd());
	    SISLCurve *qc = newCurve(ik, ik, &et[0], &ecoef[0], 1, 1, 1);
	    double zero = 0.0;

	    // Intersect
	    double eps = 1.0e-6;
	    int kstat = 0;
	    int kcrv=0, kpt=0;
	    double *epar = 0;
	    SISLIntcurve **intcv = 0;
	    if (qc)
	      s1871(qc, &zero, 1, eps, &kpt, &epar, &kcrv, &intcv, &kstat);
	    if (qc)
	      freeCurve(qc);

	    // Compute cartesian points and curves associated with intersections
	    for (kr=0; kr<kpt; ++kr)
	      {
		Vector4D barypt = (1.0 - epar[kr])*bary1 + epar[kr]*bary2;
		int kb;
		for (kb=0; kb<4; ++kb)
		  if (barypt[kb] < -0.001 || barypt[kb] > 1.001)
		    break;
		if (kb < 4)
		  continue;
		
		Vector3D pos = bc_.baryToCart(barypt);
		sfpoints.insert(sfpoints.end(), pos.begin(), pos.end());

	      }
	  }
      }
    
    // Output
    if (sfpoints.size() > 0)
      {
	PointCloud3D ptcloud(&sfpoints[0], sfpoints.size()/3);
	os << "400 1 0 4 255 0 0 255" << std::endl;
	ptcloud.write(os);
      }
}

//===========================================================================
void ImplicitApprox::visualize(vector<Point> points, Point& dir, std::ostream& os)
//===========================================================================
{
  BoundingBox bb(3);
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      bb.addUnionWith(points[ki]);
    }
  Point low = bb.low();
  Point high = bb.high();
  Point bmid = 0.5*(low + high);
  Point diag = high - low;
  double diaglen = diag.length();

  double gap = 1.0e-6;
  Point xdir(1.0, 0.0, 0.0);
  Point ydir(0.0, 1.0, 0.0);
  Point zdir(0.0, 0.0, 1.0);
  CompositeModelFactory factory(gap, gap, 10.0*gap, 0.01, 0.05);
  shared_ptr<SurfaceModel> boxmod(factory.createFromBox(low-2.0*diag, xdir, ydir, 
							5*diag[0], 5*diag[1], 5*diag[2]));
    
    // Find the coordinate direction with the largest angle with the view direction
    double a1 = xdir.angle(dir);
    double a2 = ydir.angle(dir);
    double a3 = zdir.angle(dir);
    Point dir2;
    if (a1 > std::min(a2, a3))
      dir2 = xdir;
    else if (a2 > a3)
      dir2 = ydir;
    else
      dir2 = zdir;
    Point dir3 = dir%dir2;
    dir2 = dir%dir3;
    if (dir2*(high-low) < 0.0)
      dir2 *= -1.0;
    if (dir3*(high-low) < 0.0)
      dir3 *= -1.0;
    dir2.normalize();
    dir3.normalize();
    double len = low.dist(high);
    int nmb_sample = 100;
    shared_ptr<SplineCurve> cv1(new SplineCurve(bmid-len*dir2, 0.0, bmid+len*dir2, 1.0));
    shared_ptr<SplineCurve> cv2(new SplineCurve(bmid-len*dir3, 0.0, bmid+len*dir3, 1.0));
    SweepSurfaceCreator sweep;
    shared_ptr<SplineSurface> ssf(sweep.linearSweptSurface(*cv1, *cv2, bmid));
    double del = 1.0/(double)(nmb_sample-1);
    double p1, p2;
    int ki, kj, kr;
    int ik = degree_ + 1;
    vector<double> et(2*ik, 0.0);  // Knot vector of line curve
    for (int ki=0; ki<ik; ++ki)
      et[ik+ki] = 1.0;

    vector<double> sfpoints;
    vector<double> vecs;
    vector<double> linesegs;
    vector<double> der;
    vector<double> der2;
    vector<double> lineder;
    // Evaluate line
    vector<double> tmpline;
    for (kj=0, p2=0.0; kj<nmb_sample; ++kj, p2+=del)
      {
	for (ki=0, p1=0.0; ki<nmb_sample; ++ki, p1+=del)
	  {
	    // Compute barysentric coordinates of end points of line
	    // First cartesian
	    Point sfpos = ssf->ParamSurface::point(p1,p2);
	    Point cart1 = sfpos + len*dir;
	    Point cart2 = sfpos - len*dir;
	    tmpline.insert(tmpline.end(), cart1.begin(), cart1.end());
	    tmpline.insert(tmpline.end(), cart2.begin(), cart2.end());

	    Vector4D bary1 = bc_.cartToBary(Vector3D(cart1[0], cart1[1], cart1[2]));
	    Vector4D bary2 = bc_.cartToBary(Vector3D(cart2[0], cart2[1], cart2[2]));

	    Vector3D tp1 = bc_.baryToCart(bary1);
	    Vector3D tp2 = bc_.baryToCart(bary2);
	    
	    // Pick line
	    BernsteinPoly line = implicit_.pickLine(bary1, bary2);

	    // Compute zeroes of bernstein polynomial
	    // First make sisl curve
	    vector<double> ecoef(line.coefsBegin(), line.coefsEnd());
	    SISLCurve *qc = newCurve(ik, ik, &et[0], &ecoef[0], 1, 1, 1);
	    double zero = 0.0;

	    // Intersect
	    double eps = 1.0e-6;
	    int kstat = 0;
	    int kcrv=0, kpt=0;
	    double *epar = 0;
	    SISLIntcurve **intcv = 0;
	    if (qc)
	      s1871(qc, &zero, 1, eps, &kpt, &epar, &kcrv, &intcv, &kstat);
	    if (qc)
	      freeCurve(qc);

	    // Compute cartesian points and curves associated with intersections
	    for (kr=0; kr<kpt; ++kr)
	      {
		Vector4D barypt = (1.0 - epar[kr])*bary1 + epar[kr]*bary2;
		int kb;
		for (kb=0; kb<4; ++kb)
		  if (barypt[kb] < -0.001 || barypt[kb] > 1.001)
		    break;
		if (kb < 4)
		  continue;
		
		Vector3D pos = bc_.baryToCart(barypt);
		sfpoints.insert(sfpoints.end(), pos.begin(), pos.end());

	      }
	  }
      }
    
    // Output
    if (sfpoints.size() > 0)
      {
	PointCloud3D ptcloud(&sfpoints[0], sfpoints.size()/3);
	os << "400 1 0 4 255 0 0 255" << std::endl;
	ptcloud.write(os);
      }
}
