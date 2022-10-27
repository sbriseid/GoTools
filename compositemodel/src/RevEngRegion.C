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

#include "GoTools/compositemodel/RevEngRegion.h"
#include "GoTools/compositemodel/RevEngPoint.h"
#include "GoTools/compositemodel/RevEngUtils.h"
#include "GoTools/compositemodel/HedgeSurface.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/Cone.h"
#include "GoTools/geometry/Sphere.h"
#include "GoTools/geometry/Torus.h"
#include "GoTools/geometry/Circle.h"
#include "GoTools/geometry/PointCloud.h"
//#include "GoTools/utils/Array.h"
//#include "GoTools/utils/Point.h"
#include "newmat.h"
#include "newmatap.h"
#include <vector>
#include <fstream>

using namespace Go;
using std::vector;
  
//===========================================================================
RevEngRegion::RevEngRegion()
//===========================================================================
  : classification_type_(CLASSIFICATION_UNDEF), associated_sf_(0),
    mink1_(0.0), maxk1_(0.0), mink2_(0.0), maxk2_(0.0), visited_(false)
{
  bbox_ = BoundingBox(3);
}

//===========================================================================
RevEngRegion::RevEngRegion(int classification_type)
//===========================================================================
  : classification_type_(classification_type), associated_sf_(0),
    mink1_(0.0), maxk1_(0.0), mink2_(0.0), maxk2_(0.0), visited_(false)
{
  bbox_ = BoundingBox(3);
}

//===========================================================================
RevEngRegion::RevEngRegion(int classification_type,
			   vector<RevEngPoint*> & points)
//===========================================================================
  : group_points_(points), classification_type_(classification_type),
    associated_sf_(0), visited_(false)
{
  for (size_t kj=0; kj<group_points_.size(); ++kj)
    group_points_[kj]->setRegion(this);
  
  // Bounding box and principal curvature summary
  maxk2_ = std::numeric_limits<double>::lowest();
  mink2_ = std::numeric_limits<double>::max();
  maxk1_ = std::numeric_limits<double>::lowest();
  mink1_ = std::numeric_limits<double>::max();
  bbox_ = BoundingBox(3);
  for  (size_t kj=0; kj<group_points_.size(); ++kj)
    {
      double k1 = group_points_[kj]->minPrincipalCurvature();
      double k2 = group_points_[kj]->maxPrincipalCurvature();
      mink1_ = std::min(mink1_, fabs(k1));
      maxk1_ = std::max(maxk1_, fabs(k1));
      mink2_ = std::min(mink2_, fabs(k2));
      maxk2_ = std::max(maxk2_, fabs(k2));
      Vector3D point = group_points_[kj]->getPoint();
      Point point2(point[0], point[1], point[2]);
      bbox_.addUnionWith(point2);
    }
  
  normalcone_ = DirectionCone(group_points_[0]->getPCANormal());
  for  (size_t kj=1; kj<group_points_.size(); ++kj)
    {
      normalcone_.addUnionWith(group_points_[kj]->getPCANormal());
    }
}

//===========================================================================
void RevEngRegion::collect(RevEngPoint *pt)
//===========================================================================
{
  if (pt->hasRegion() && pt->region() != this)
    return;  // Cannot grow
  group_points_.push_back(pt);
  if (classification_type_ == CLASSIFICATION_UNDEF)
    return; // Cannot grow
  int type = pt->surfaceClassification(classification_type_);
  double meancurv = pt->meanCurvature();
  if (type == C1_UNDEF)  // SI_UNDEF == C1_UNDEF
    return; // Cannot grow

  vector<RevEngPoint*> grouped;
  grouped.push_back(pt);
  Vector3D xyz = pt->getPoint();
  Point xyz2(xyz[0], xyz[1], xyz[2]);
  bbox_.addUnionWith(xyz2);
  if (normalcone_.dimension() == 0)
    normalcone_ = DirectionCone(pt->getPCANormal());
  else
    normalcone_.addUnionWith(pt->getPCANormal());
  bool planar = planartype();
  for (size_t kj=0; kj<grouped.size(); ++kj)
    {
      vector<ftSamplePoint*> next = grouped[kj]->getNeighbours();
      for (size_t ki=0; ki<next.size(); ++ki)
	{
	  RevEngPoint *curr = dynamic_cast<RevEngPoint*>(next[ki]);
	  if (!curr)
	    continue;  // Should not happen
	  if (curr->hasRegion())
	    continue;  // Already belonging to a segment
	  if (curr->isEdge())
	    continue;  // An edge point
	  int type2 = curr->surfaceClassification(classification_type_);
	  if (type2 != type)
	    continue;   // Different classification
	  if ((!planar) && meancurv*curr->meanCurvature() < 0.0)
	    continue;  // Not compatible with one primary surface

	  // Add to region
	  curr->setRegion(this);
	  group_points_.push_back(curr);
	  xyz = curr->getPoint();
	  xyz2 = Point(xyz[0], xyz[1], xyz[2]);
	  bbox_.addUnionWith(xyz2);
	  normalcone_.addUnionWith(curr->getPCANormal());

	  // Continue growing from this point
	  grouped.push_back(curr);
	}
    }

  // Principal curvature summary
  maxk2_ = std::numeric_limits<double>::lowest();
  mink2_ = std::numeric_limits<double>::max();
  maxk1_ = std::numeric_limits<double>::lowest();
  mink1_ = std::numeric_limits<double>::max();
  for  (size_t kj=0; kj<group_points_.size(); ++kj)
    {
      double k1 = group_points_[kj]->minPrincipalCurvature();
      double k2 = group_points_[kj]->maxPrincipalCurvature();
      mink1_ = std::min(mink1_, fabs(k1));
      maxk1_ = std::max(maxk1_, fabs(k1));
      mink2_ = std::min(mink2_, fabs(k2));
      maxk2_ = std::max(maxk2_, fabs(k2));
    }
}

//===========================================================================
RevEngPoint* RevEngRegion::seedPoint()
//===========================================================================
{
  int min_next = std::max(10, (int)group_points_.size()/100);
  double rfac = 3;
  double min_out = std::numeric_limits<double>::max();
  int min_ix = -1;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      double local_len = group_points_[ki]->getMeanEdgLen();
      double radius = 2.0*rfac*local_len;
      vector<RevEngPoint*> nearpts;
      group_points_[ki]->fetchClosePoints2(radius, min_next, 5*min_next, nearpts);

      // Count deviant points
      int deviant = 0;
      for (size_t kj=0; kj<nearpts.size(); ++kj)
      {
	if (nearpts[kj]->region() != this)
	  ++deviant;
      }

      if (deviant == 0)
	return group_points_[ki];   // Seed point found

      double curr_dev = (double)deviant/(double)nearpts.size();
      if (curr_dev < min_out)
	{
	  min_out = curr_dev;
	  min_ix = (int)ki;
	}
    }
  return (min_ix >= 0) ? group_points_[min_ix] : 0;
}

//===========================================================================
void RevEngRegion::growLocal(RevEngPoint* seed, double tol, double radius,
			     int min_close, vector<RevEngPoint*>& outpts)
//===========================================================================
{
  // Fetch nearby points belonging to the same region
  vector<RevEngPoint*> nearpts;
  seed->fetchClosePoints2(radius, min_close, 5*min_close, nearpts, this);
  nearpts.insert(nearpts.begin(), seed);
  for (size_t ki=0; ki<nearpts.size(); )
    {
      if (nearpts[ki]->getPointDistance() > 1.5*nearpts[ki]->getAveragePointDistance())
	nearpts.erase(nearpts.begin()+ki);
      else
	++ki;
    }

  // Approximate with implicit algebraic surface
  int degree = 2;
  impl_ = shared_ptr<ImplicitApprox>(new ImplicitApprox());
  impl_->approx(nearpts, degree);
  std::ofstream outviz("implsf.g2");
  impl_->visualize(nearpts, outviz);

  std::ofstream no("noproject.g2");
  double large = 1.0e3;
  double maxdist = 0.0, avdist = 0.0;
  int nodist = 0;
  for (size_t ki=0; ki<nearpts.size(); ++ki)
    {
      double dist = fabs(impl_->estimateDist(nearpts[ki]));
      if (dist > large)
	{
	  no << "400 1 0 4 0 0 255 255" << std::endl;
	  no << "1" << std::endl;
	  no << nearpts[ki]->getPoint() << std::endl;
	}
      else
	{
	  maxdist = std::max(maxdist, dist);
	  avdist += dist;
	  ++nodist;
	}
    }
  avdist /= (double)nodist;
  double tol2 = 1.5*maxdist; //std::min(tol, 10.0*maxdist);
  std::ofstream ofnear("nearpts.g2");
  ofnear << "400 1 0 4 0 255 0 255" << std::endl;
  ofnear << nearpts.size() << std::endl;
  for (size_t ki=0; ki<nearpts.size(); ++ki)
    ofnear << nearpts[ki]->getPoint() << std::endl;
   
  // Check accuracy
  std::ofstream of("grow_local.g2");
  vector<RevEngPoint*> out, in;
  vector<RevEngPoint*> core;
  core.push_back(seed);
  seed->setVisited();
  int max_new_out = std::max(50, (int)group_points_.size()/100);
  int prev_out = 0;
  for (size_t ki=0; ki<core.size(); ++ki)
    {
      vector<ftSamplePoint*> next = core[ki]->getNeighbours();
      for (size_t kj=0; kj<next.size(); ++kj)
	{
	  RevEngPoint *pt = dynamic_cast<RevEngPoint*>(next[kj]);
	  if (pt->visited())
	    continue;
	  if (pt->region() != this)
	    continue;    // Can consider storing the point for later access
	  pt->setVisited();
	  core.push_back(pt);
	  
	  double dist = fabs(impl_->estimateDist(pt));
	  if (dist > large)
	    {
	      no << "400 1 0 4 0 0 255 255" << std::endl;
	      no << "1" << std::endl;
	      no << pt->getPoint() << std::endl;
	    }
 	  if (dist <= tol2)
	    {
	      in.push_back(pt);
	    }
	  else
	    out.push_back(pt);
	}

      if ((int)out.size() > prev_out + max_new_out)
	{
	  // Update implicit
	  vector<RevEngPoint*> curr_pts(in.begin(), in.end());
	  curr_pts.insert(curr_pts.end(), out.begin(), out.end());
	  shared_ptr<ImplicitApprox> impl2(new ImplicitApprox());
	  impl2->approx(curr_pts, degree);
	  std::ofstream outviz2("implsf2.g2");
	  impl2->visualize(curr_pts, outviz2);
	  std::ofstream ofnear2("nearpts2.g2");
	  ofnear2 << "400 1 0 4 0 255 0 255" << std::endl;
	  ofnear2 << curr_pts.size() << std::endl;
	  for (size_t ki=0; ki<curr_pts.size(); ++ki)
	    ofnear2 << curr_pts[ki]->getPoint() << std::endl;
 	  
	  vector<RevEngPoint*> out2, in2;
	  for (size_t kr=0; kr<curr_pts.size(); ++kr)
	    {
	      double dist = fabs(impl2->estimateDist(curr_pts[kr]));
	      if (dist <= tol2)
		{
		  in2.push_back(curr_pts[kr]);
		}
	      else
		out2.push_back(curr_pts[kr]);
	    }
	  
	  if (out2.size() < out.size())
	    {
	      impl_ = impl2;
	      in = in2;
	      out = out2;
	    }
	  else
	    break;
	  prev_out = (int)out.size();
	}
    }

  // Remaining points
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(group_points_[ki]);
      double dist = fabs(impl_->estimateDist(pt));
      if (dist <= tol2)
	{
	  in.push_back(pt);
	}
      else
	out.push_back(pt);
    }
  
  if (in.size() > 0)
    {
      of << "400 1 0 0" << std::endl;
      of << in.size() << std::endl;
      for (size_t kr=0; kr<in.size(); ++kr)
	of << in[kr]->getPoint() << std::endl;
      of << std::endl;
    }
  
  if (out.size() > 0)
    {
      of << "400 1 0 0" << std::endl;
      of << out.size() << std::endl;
      for (size_t kr=0; kr<out.size(); ++kr)
	of << out[kr]->getPoint() << std::endl;
    }
  
  of << "400 1 0 0" << std::endl;
  of << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    of << group_points_[kr]->getPoint() << std::endl;
  of << std::endl;

  int stop_break = 1;
  
}

//===========================================================================
bool RevEngRegion::extractPlane(double tol, int min_pt,
				vector<shared_ptr<HedgeSurface> >& hedgesfs,
				vector<HedgeSurface*>& prevsfs,
				std::ostream& fileout)
//===========================================================================
{
  std::ofstream of("curr_region.g2");
  of << "400 1 0 4 100 0 155 255" << std::endl;
  of << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    of << group_points_[kr]->getPoint() << std::endl;
  of << "410 1 0 4 200 55 0 255" << std::endl;
  of << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      Vector3D xyz = group_points_[kr]->getPoint();
      Point xyz2(xyz[0], xyz[1], xyz[2]);
      Point norm = group_points_[kr]->getMongeNormal();
      of << xyz2 << " " << xyz2 + norm << std::endl;
    }
  
  std::ofstream of2("curr_normals.g2");
  of2 << "400 1 0 4 100  0 155 255" << std::endl;
  of2 << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(group_points_[kr]);
      Point norm = pt->getMongeNormal();
      of2 << norm << std::endl;
    }
  Sphere sph(1.0, Point(0.0, 0.0, 0.0), Point(0.0, 0.0, 1.0),
	     Point(1.0, 0.0, 0.0));
  sph.writeStandardHeader(of2);
  sph.write(of2);

  bool found = false;
  int min_nmb = 20;
  if (group_points_.size() < min_nmb)
    return false;
  Point normal1 = normalcone_.centre();
  Point normal(0.0, 0.0, 0.0);
  Point pos(0.0, 0.0, 0.0);
  double wgt = 1.0/(double)(group_points_.size());
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Point curr = group_points_[ki]->getPCANormal();
      Vector3D xyz = group_points_[ki]->getPoint();
      normal += wgt*curr;
      pos += wgt*Point(xyz[0], xyz[1], xyz[2]);
    }

  double angtol = 0.1;
  vector<RevEngPoint*> in, out;
  analysePlaneProperties(normal, angtol, in, out);
  
  shared_ptr<Plane> surf1(new Plane(pos, normal1));
  shared_ptr<Plane> surf(new Plane(pos, normal));

  Point pos2 = pos;
  Point normal2;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(group_points_.begin(), group_points_.end()));
  RevEngUtils::computePlane(group, pos2, normal2);
  shared_ptr<Plane> surf2(new Plane(pos2, normal2));
  
  impl_ = shared_ptr<ImplicitApprox>(new ImplicitApprox());
  impl_->approx(group_points_, 1);
  Point pos3, normal3;
  impl_->projectPoint(pos, normal, pos3, normal3);

  vector<RevEngPoint*> in3, out3;
  analysePlaneProperties(normal3, angtol, in3, out3);
  
  shared_ptr<Plane> surf3(new Plane(pos3, normal3));

  // Check accuracy
  double maxdist, avdist;
  double maxdist1, avdist1;
  double maxdist2, avdist2;
  double maxdist3, avdist3;
  int num_inside, num_inside1, num_inside2, num_inside3;
  vector<RevEngPoint*> inpt, outpt;
  // RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
  // 			  surf, tol, maxdist, avdist, num_inside);
  // RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
  // 			  surf1, tol, maxdist1, avdist1, num_inside1);
  // RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
  // 			  surf2, tol, maxdist2, avdist2, num_inside2);
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			  surf3, tol, maxdist3, avdist3, num_inside3, inpt, outpt);
  std::ofstream ofd("in_out.g2");
  ofd << "400 1 0 4 155 50 50 255" << std::endl;
  ofd << inpt.size() << std::endl;
  for (size_t kr=0; kr<inpt.size(); ++kr)
    ofd << inpt[kr]->getPoint() << std::endl;
  ofd << "400 1 0 4 50 155 50 255" << std::endl;
  ofd << outpt.size() << std::endl;
  for (size_t kr=0; kr<outpt.size(); ++kr)
    ofd << outpt[kr]->getPoint() << std::endl;

  Point low = bbox_.low();
  Point high = bbox_.high();
  double len = low.dist(high);
  surf3->setParameterBounds(-0.5*len, -0.5*len, 0.5*len, 0.5*len);
  std::ofstream plane("curr_plane.g2");
  surf3->writeStandardHeader(plane);
  surf3->write(plane);

  int num = (int)group_points_.size();
  if (num_inside3 > min_pt && num_inside3 > num/2 && avdist3 <= tol) 
    {
      found = true;
      setAccuracy(maxdist3, avdist3, num_inside3);
      
      std::cout << "Plane. N1: " << num << ", N2: " << num_inside3 << ", max: " << maxdist3 << ", av: " << avdist3 << std::endl;

      shared_ptr<HedgeSurface> hedge(new HedgeSurface(surf3, this));
      for (size_t kh=0; kh<associated_sf_.size(); ++kh)
	prevsfs.push_back(associated_sf_[kh]);
      associated_sf_.clear();
      associated_sf_.push_back(hedge.get());
      hedgesfs.push_back(hedge);
	  
      surf3->writeStandardHeader(fileout);
      surf3->write(fileout);
      //impl_->visualize(group_points_, fileout);
    }

  return found;
}

//===========================================================================
void RevEngRegion::analysePlaneProperties(Point avnorm, double angtol,
					  vector<RevEngPoint*>& in,
					  vector<RevEngPoint*> out)
//===========================================================================
{
  vector<double> midang(group_points_.size());
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Point curr = group_points_[ki]->getPCANormal();
      double ang = curr.angle(avnorm);
      midang[ki] = ang;
      if (ang <= angtol)
	in.push_back(group_points_[ki]);
      else
	out.push_back(group_points_[ki]);
    }

  std::sort(midang.begin(), midang.end());
  int stop_break = 1;
}

//===========================================================================
void RevEngRegion::analyseNormals(double tol, Point& normal, Point& centre,
				  double& radius) //double& beta)
//===========================================================================
{
  std::ofstream of2("curr_normals.g2");
  of2 << "400 1 0 4 100  0 155 255" << std::endl;
  of2 << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(group_points_[kr]);
      Point norm = pt->getMongeNormal();
      of2 << norm << std::endl;
    }
  Sphere sph(1.0, Point(0.0, 0.0, 0.0), Point(0.0, 0.0, 1.0),
	     Point(1.0, 0.0, 0.0));
  sph.writeStandardHeader(of2);
  sph.write(of2);

  Point axis;
  Point Cx;
  Point Cy;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(group_points_.begin(), group_points_.end()));
  RevEngUtils::computeAxis(group, axis, Cx, Cy);

  Point pnt(0.0, 0.0, 0.0);
  vector<Point> vec(group_points_.size());
  double wgt = 1.0/(double)(group_points_.size());
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      vec[ki] = group_points_[ki]->getPCANormal();
      pnt += wgt*vec[ki];
    }
  Point tmpnorm;
  double maxang = 0.0;
  size_t ix = 0;
  for (size_t kj=1; kj<vec.size(); ++kj)
    {
      double ang = vec[0].angle(vec[kj]);
      ang = std::min(ang, M_PI-ang);
      if (ang > maxang)
	{
	  maxang = ang;
	  ix = kj;
	}
    }
  tmpnorm = vec[0].cross(vec[ix]);
  tmpnorm.normalize();

  ImplicitApprox impl;
  impl.approxPoints(vec, 1);

  //Point pos; //, normal;
  impl.projectPoint(vec[0], tmpnorm, pnt, normal);
  shared_ptr<Plane> surf(new Plane(pnt, axis));

  Point origo(0.0, 0.0,0.0);
  of2 << "410 1 0 4 255 0 0 255" << std::endl;
  of2 << "1" << std::endl;
  of2 << origo << " " << axis << std::endl;
  
  double maxdist, avdist;
  int num_inside;
  int num = (int)vec.size();
  RevEngUtils::distToSurf(vec, surf, tol, maxdist, avdist, num_inside);

  //double radius;
  RevEngUtils::computeRadius(vec, axis, Cx, Cy, radius);
  centre = pnt;
  // double r2 = std::min(radius, 1.0);
  // double phi = acos(r2);
  // double delta = r2*tan(phi);
  // double alpha = (delta > 1.0e-10) ? atan((1.0-r2)/delta) : M_PI;
  // beta = M_PI - alpha;
  int stop_break = 1;

}

//===========================================================================
bool RevEngRegion::extractCylinder(double tol, int min_pt, double mean_edge_len,
				   vector<shared_ptr<HedgeSurface> >& hedgesfs,
				   vector<HedgeSurface*>& prevsfs,
				   std::ostream& fileout)
//===========================================================================
{
  std::ofstream of("curr_region.g2");
  of << "400 1 0 4 100 0 155 255" << std::endl;
  of << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    of << group_points_[kr]->getPoint() << std::endl;
  of << "410 1 0 4 200 55 0 255" << std::endl;
  of << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      Vector3D xyz = group_points_[kr]->getPoint();
      Point xyz2(xyz[0], xyz[1], xyz[2]);
      Point norm = group_points_[kr]->getPCANormal();
      of << xyz2 << " " << xyz2 + norm << std::endl;
    }
  
  of << "410 1 0 4 0 100  155 255" << std::endl;
  of << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      Vector3D xyz = group_points_[kr]->getPoint();
      Point xyz2(xyz[0], xyz[1], xyz[2]);
      Point vec = group_points_[kr]->minCurvatureVec();
      of << xyz2 << " " << xyz2 + vec << std::endl;
    }

      Point vec0 = group_points_[0]->minCurvatureVec();
Point avvec(0.0, 0.0, 0.0);
double wgt = 1.0/(double)group_points_.size();
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      Point vec = group_points_[kr]->minCurvatureVec();
if (vec*vec0 < 0.0)
vec *= -1;
vec *= wgt;
avvec += vec;
}

double angtol= 0.1;
vector<RevEngPoint*> in, out;
analyseCylinderProperties(avvec, angtol, in, out);
  // of << "410 1 0 4 100 0  155 255" << std::endl;
  // of << group_points_.size() << std::endl;
  // for (size_t kr=0; kr<group_points_.size(); ++kr)
  //   {
  //     vector<Point> nearpts;
  //     double local_len = group_points_[kr]->getMeanEdgLen();
  //     double radius = 3.0*(local_len + mean_edge_len);
  //     double min_next = 10;
  //     Point curr = group_points_[kr]->fetchClosePoints(radius, min_next, nearpts);
  //     Point normal, mincvec, maxcvec;
  //     double minc, maxc;
  //     double currdist, avdist;
  //     Point eigen1 = group_points_[kr]->getPCAEigen1();
  //     Point eigen3 = group_points_[kr]->getPCANormal();
  //     RevEngUtils::computeMonge(curr, nearpts, eigen1, eigen3,
  // 				normal, mincvec, minc,
  // 				maxcvec, maxc, currdist, avdist);

  //     Vector3D xyz = group_points_[kr]->getPoint();
  //     Point xyz2(xyz[0], xyz[1], xyz[2]);
  //     of << xyz2 << " " << xyz2 + mincvec << std::endl;
  //   }
  
  std::ofstream of2("curr_normals.g2");
  of2 << "400 1 0 4 100  0 155 255" << std::endl;
  of2 << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(group_points_[kr]);
      Point norm = pt->getPCANormal();
      of2 << norm << std::endl;
    }
  Sphere sph(1.0, Point(0.0, 0.0, 0.0), Point(0.0, 0.0, 1.0),
	     Point(1.0, 0.0, 0.0));
  sph.writeStandardHeader(of2);
  sph.write(of2);

  of2 << "400 1 0 4 0 100  155 255" << std::endl;
  of2 << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(group_points_[kr]);
      Point vec = pt->minCurvatureVec();
      of2 << vec << std::endl;
    }
  of2 << std::endl;
  
  
  // // Associate points with radius based on density on Gauss sphere
  //   extendWithGaussRad();
  
  // of2 << "400 1 0 4 100 155 0 255" << std::endl;
  // of2 << group_points_.size() << std::endl;
  // for (size_t kr=0; kr<group_points_.size(); ++kr)
  //   {
  //     RevEngPoint *pt = dynamic_cast<RevEngPoint*>(group_points_[kr]);
  //     Point norm = pt->getMongeNormal();
  //     double Grad = group_points_[kr]->getGaussRad();
  //     of2 << Grad*norm << std::endl;
  //   }

  //double beta;
  // Point normalG, centreG;
  // double radG;
  // analyseNormals(tol, normalG, centreG, radG); //beta);
  
  // Cylinder orientation by covariance matrix of normal vectors
  bool found = false;
  int min_nmb = 20;
  double eps = 1.0e-6;
  if (group_points_.size() < min_nmb)
    return false;
  Point axis;
  Point Cx;
  Point Cy;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(group_points_.begin(), group_points_.end()));
  RevEngUtils::computeAxis(group, axis, Cx, Cy);
  Point origo(0.0, 0.0, 0.0);
  of2 << "410 1 0 4 255 0 0 255" << std::endl;
  of2 << "1" << std::endl;
  of2 << origo << " " << axis << std::endl;

 

  vector<RevEngPoint*> in2, out2;
  analyseCylinderProperties(axis, angtol, in2, out2);
  
  Point axis2;
  Point Cx2;
  Point Cy2;
  RevEngUtils::coneAxis(group, axis2, Cx2, Cy2);
  of2 << "410 1 0 4 55 200 0 255" << std::endl;
  of2 << "1" << std::endl;
  of2 << origo << " " << axis2 << std::endl;

  vector<RevEngPoint*> in3, out3;
  analyseCylinderProperties(axis2, angtol, in3, out3);
  
  Point apex;
  double phi;
  RevEngUtils::coneApex(group, axis2, apex, phi);
  
  double rad;
  Point pnt;
  Point low = bbox_.low();
  Point high = bbox_.high();
  double len = low.dist(high);
  RevEngUtils::computeCylPosRadius(group, low, high,
				   axis, Cx, Cy, pnt, rad);

  vector<Point> rotated;
  RevEngUtils::rotateToPlane(group, Cx, axis, pnt, rotated);
  std::ofstream of3("rotated_pts.g2");
  of3 << "400 1 0 4 255 0 0 255" << std::endl;
  of3 << rotated.size() << std::endl;
  for (size_t kr=0; kr<rotated.size(); ++kr)
    of3 << rotated[kr] << std::endl;
  of3 << "410 1 0 4 0 0 255 255" << std::endl;
  of3 << "1" << std::endl;
  of3 << pnt-0.5*len*axis << " " << pnt+0.5*len*axis << std::endl;
  of3 << "410 1 0 4 0 255 0 255" << std::endl;
  of3 << "1" << std::endl;
  of3 << pnt-0.5*len*Cx << " " << pnt+0.5*len*Cx << std::endl;

  shared_ptr<Circle> circ(new Circle(rad, pnt, axis, Cx));
  vector<Point> projected;
  double maxdp, avdp;
  RevEngUtils::projectToPlane(group, axis, pnt, projected, maxdp, avdp);
  std::ofstream ofp3("projected_pts.g2");
  ofp3 << "400 1 0 4 255 0 0 255" << std::endl;
  ofp3 << projected.size() << std::endl;
  for (size_t kr=0; kr<projected.size(); ++kr)
    ofp3 << projected[kr] << std::endl;
  circ->writeStandardHeader(ofp3);
  circ->write(ofp3);
  
  vector<Point> rotated2;
  RevEngUtils::rotateToPlane(group, Cx2, axis2, apex, rotated2);
  std::ofstream of4("rotated_pts2.g2");
  of4 << "400 1 0 4 255 0 0 255" << std::endl;
  of4 << rotated2.size() << std::endl;
  for (size_t kr=0; kr<rotated2.size(); ++kr)
    of4 << rotated2[kr] << std::endl;
  of4 << "410 1 0 4 0 0 255 255" << std::endl;
  of4 << "1" << std::endl;
  of4 << apex-0.5*len*axis2 << " " << apex+0.5*len*axis2 << std::endl;
  of4 << "410 1 0 4 0 255 0 255" << std::endl;
  of4 << "1" << std::endl;
  of4 << apex-0.5*len*Cx2 << " " << apex+0.5*len*Cx2 << std::endl;
  // pluckerAxis();

  vector<Point> maxcvec(group_points_.size());
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(group_points_[kr]);
      maxcvec[kr] = pt->maxCurvatureVec();
    }
  Point axis3, Cx3, Cy3;
  RevEngUtils::computeAxis(maxcvec, axis3, Cx3, Cy3);
			
  
  shared_ptr<Cylinder> cyl(new Cylinder(rad, pnt, axis, Cy));
  cyl->setParamBoundsV(-0.5*len, 0.5*len);
  //PointCloud3D cloud(&xyzpoints[0], xyzpoints.size()/3);

  Point pos0 = apex + ((pnt - apex)*axis2)*axis2;
  double dd = apex.dist(pos0);
  double rad0 = dd*tan(phi);
  shared_ptr<Cone> cone(new Cone(rad0, pos0, axis2, Cx2, phi));
  cone->setParamBoundsV(-0.5*len, 0.5*len);
 // shared_ptr<Cone> cone2(new Cone(rad, pnt, -axis, -Cy, beta));

  std::ofstream ofs("sfs.g2");
  cyl->writeStandardHeader(ofs);
  cyl->write(ofs);
  cone->writeStandardHeader(ofs);
  cone->write(ofs);
  ofs <<"400 1 0 4 255 0 0 255" << std::endl;
  ofs << "1" << std::endl;
  ofs << pos0 << std::endl;
  
  // Check accuracy
  double maxd, avd, maxd2, avd2, maxd3, avd3;
  int num2, num3, num4;
  //shared_ptr<ParamSurface> surf = cyl;
  vector<RevEngPoint*> inpt, outpt, inpt2, outpt2;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			  cyl, tol, maxd, avd, num2, inpt, outpt);
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
  			  cone, tol, maxd2, avd2, num3, inpt2, outpt2);
  //  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
  // 			  cone2, tol, maxd3, avd3, num4);
  std::ofstream ofd("in_out.g2");
  ofd << "400 1 0 4 155 50 50 255" << std::endl;
  ofd << inpt.size() << std::endl;
  for (size_t kr=0; kr<inpt.size(); ++kr)
    ofd << inpt[kr]->getPoint() << std::endl;
  ofd << "400 1 0 4 50 155 50 255" << std::endl;
  ofd << outpt.size() << std::endl;
  for (size_t kr=0; kr<outpt.size(); ++kr)
    ofd << outpt[kr]->getPoint() << std::endl;


  int num = (int)group_points_.size();
  if (num2 > min_pt && num2 > num/2) //avd <= avlim && maxd <= maxlim)
    {
      found = true;
      setAccuracy(maxd, avd, num2);
      
      // Limit cylinder with respect to bounding box
      double gap = 1.0e-6;
      Point xdir(1.0, 0.0, 0.0);
      Point ydir(0.0, 1.0, 0.0);
      Point zdir(0.0, 0.0, 1.0);
      Point bbdiag = high - low;
      CompositeModelFactory factory(gap, gap, 10.0*gap, 0.01, 0.05);
      shared_ptr<SurfaceModel> boxmod(factory.createFromBox(low-2.0*bbdiag, xdir, ydir, 
							    5*bbdiag[0], 5*bbdiag[1],
							    5*bbdiag[2]));
      vector<shared_ptr<ParamSurface> > sfs;
      sfs.push_back(cyl);
      shared_ptr<SurfaceModel> cylmod(new SurfaceModel(gap, gap, 10.0*gap, 0.01,
						       0.05, sfs));
      vector<shared_ptr<SurfaceModel> > divcyl = cylmod->splitSurfaceModels(boxmod);
      
      std::cout << "Cylinder. N1: " << num << ", N2: " << num2 << ", max: " << maxd << ", av: " << avd << std::endl;
      associated_sf_.clear();
      for (int ka=0; ka<divcyl[0]->nmbEntities(); ++ka)
	{
	  shared_ptr<ParamSurface> cyl2 = divcyl[0]->getSurface(ka);
	  shared_ptr<HedgeSurface> hedge(new HedgeSurface(cyl2, this));
	  for (size_t kh=0; kh<associated_sf_.size(); ++kh)
	    prevsfs.push_back(associated_sf_[kh]);
	  associated_sf_.push_back(hedge.get());
	  hedgesfs.push_back(hedge);
	  
	  cyl2->writeStandardHeader(fileout);
	  cyl2->write(fileout);
	}
      // cloud2.writeStandardHeader(fileout);
      // cloud2.write(fileout);
      
      // if (dotest)
      // 	{
      // 	  double ang0 = axis0.angle(ax);
      // 	  ang0 = std::min(ang0, M_PI-ang0);
      // 	  Point vec0 = pnt - pos0;
      // 	  double dist0 = (vec0 - (vec0*ax)*ax).length();
      // 	  std::cout << "Angular diff: " << ang0 << ", center diff: " << dist0;
      // 	  std::cout << ", radius diff: " << fabs(rad-rad0) << std::endl;
      // 	}
    }
  int stop_break0 = 1;
  return found;
}

//===========================================================================
bool RevEngRegion::extractTorus(double tol, int min_pt, double mean_edge_len,
				vector<shared_ptr<HedgeSurface> >& hedgesfs,
				vector<HedgeSurface*>& prevsfs,
				std::ostream& fileout)
//===========================================================================
{
  bool found = false;
  
  std::ofstream of("curr_region.g2");
  of << "400 1 0 4 100 0 155 255" << std::endl;
  of << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    of << group_points_[kr]->getPoint() << std::endl;
  of << "410 1 0 4 200 55 0 255" << std::endl;
  of << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      Vector3D xyz = group_points_[kr]->getPoint();
      Point xyz2(xyz[0], xyz[1], xyz[2]);
      Point norm = group_points_[kr]->getPCANormal();
      of << xyz2 << " " << xyz2 + norm << std::endl;
    }
  
  of << "410 1 0 4 0 100  155 255" << std::endl;
  of << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      Vector3D xyz = group_points_[kr]->getPoint();
      Point xyz2(xyz[0], xyz[1], xyz[2]);
      Point vec = group_points_[kr]->minCurvatureVec();
      of << xyz2 << " " << xyz2 + vec << std::endl;
    }

  // Compute curvature variations
  double k1mean = 0.0, k2mean = 0.0;
  double k1_1 = 0.0, k1_2 = 0.0, k2_1 = 0.0, k2_2 = 0.0;
  double d1 = 0.0, d2 = 0.0;
  double wgt = 1.0/(double)group_points_.size();
  Point vec(0.0, 0.0, 0.0);
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      double kmin = group_points_[kr]->minPrincipalCurvature();
      double kmax = group_points_[kr]->maxPrincipalCurvature();
      d1 += kmin;
      d2 += kmax;
      k1_1 += kmin*kmin;
      k2_1 += kmax*kmax;
      k1_2 += fabs(kmin);
      k2_2 += fabs(kmax);
      k1mean += wgt*kmin;
      k2mean += wgt*kmax;

      Point norm = group_points_[kr]->getPCANormal();
      Point minvec = group_points_[kr]->minCurvatureVec();
      Point temp = norm.cross(minvec);
      temp.normalize_checked();
      vec += wgt*temp;
    }
  double vark1 = (group_points_.size()*k1_1)/fabs(d1) - fabs(k1_2);
  double vark2 = (group_points_.size()*k2_1)/fabs(d2) - fabs(k2_2);
  double rd1 = 1.0/k1mean; //(vark1 < vark2) ? 1.0/k1mean : 1.0/k2mean;
  double rd2 = 1.0/k2mean;

  vector<Point> centr1(group_points_.size());
  vector<Point> centr2(group_points_.size());
  Point mid1(0.0, 0.0, 0.0), mid2(0.0, 0.0, 0.0);
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      Vector3D xyz = group_points_[kr]->getPoint();
      Point xyz2(xyz[0], xyz[1], xyz[2]);
      Point norm = group_points_[kr]->getPCANormal();
      centr1[kr] = xyz2 + rd1*norm;
      centr2[kr] = xyz2 + rd2*norm;
      mid1 += wgt*centr1[kr];
      mid2 += wgt*centr2[kr];
    }
  
  // of << "400 1 0 4 155 100 0 255" << std::endl;
  // of << group_points_.size() << std::endl;
  // for (size_t kr=0; kr<group_points_.size(); ++kr)
  //   {
  //     of << centr1[kr] << std::endl;
  //   }
  of << "400 1 0 4 155 100 0 255" << std::endl;
  of << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      of << centr2[kr] << std::endl;
    }

  shared_ptr<ImplicitApprox> impl1(new ImplicitApprox());
  impl1->approxPoints(centr1, 1);
  shared_ptr<ImplicitApprox> impl2(new ImplicitApprox());
  impl2->approxPoints(centr2, 1);

  double val1, val2;
  Point grad1, grad2;
  impl1->evaluate(mid1, val1, grad1);
  impl2->evaluate(mid2, val2, grad2);
  grad1.normalize_checked();
  grad2.normalize_checked();
  
  Point pos1, normal1, pos2, normal2;
  impl1->projectPoint(mid1, grad1, pos1, normal1);
  impl2->projectPoint(mid2, grad2, pos2, normal2);
  double eps1 = 1.0e-8;
  if (normal1.length() < eps1 || normal2.length() < eps1)
    return false;

  std::ofstream ofimpl("impl_plane.g2");
  impl2->visualize(group_points_, ofimpl);
  
  Point low = bbox_.low();
  Point high = bbox_.high();
  double len = low.dist(high);
  // of << "400 1 0 4 0 0 255 255" << std::endl;
  // of << "1" << std::endl;
  // of << pos1 << std::endl;
  // of << "410 1 0 4 0 0 255 255" << std::endl;
  // of << "1" << std::endl;
  // of << pos1 << " " << pos1+0.5*len*normal1 << std::endl;
  
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << "1" << std::endl;
  of << pos2 << std::endl;
  of << "410 1 0 4 255 0 0 255" << std::endl;
  of << "1" << std::endl;
  of << pos2 << " " << pos2+0.5*len*normal2 << std::endl;

  Vector3D xyz = group_points_[0]->getPoint();
  Point xyz2(xyz[0], xyz[1], xyz[2]);
  Point Cx1 = centr1[0] - xyz2;
  Cx1 -= (Cx1*normal1)*normal1;
  Cx1.normalize();
  Point Cy1 = Cx1.cross(normal1);
  Point Cx2 = centr2[0] - xyz2;
  Cx2 -= (Cx2*normal2)*normal2;
  Cx2.normalize();
  Point Cy2 = Cx2.cross(normal2);
  
  double rad1, rad2;
  Point pnt1, pnt2;
  RevEngUtils::computeCircPosRadius(centr1, normal1, Cx1, Cy1, pnt1, rad1);
  pnt1 -= ((pnt1 - pos1)*normal1)*normal1;
  RevEngUtils::computeCircPosRadius(centr2, normal2, Cx2, Cy2, pnt2, rad2);
  pnt2 -= ((pnt2 - pos2)*normal2)*normal2;
  shared_ptr<Circle> circ1(new Circle(rad1, pnt1, normal1, Cx1));
  shared_ptr<Circle> circ2(new Circle(rad2, pnt2, normal2, Cx2));
  // circ1->writeStandardHeader(of);
  // circ1->write(of);
  circ2->writeStandardHeader(of);
  circ2->write(of);

  vector<Point> rotated;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(group_points_.begin(), group_points_.end()));
  RevEngUtils::rotateToPlane(group, Cx2, normal2, pnt2, rotated);
  std::ofstream of3("rotated_pts.g2");
  of3 << "400 1 0 4 255 0 0 255" << std::endl;
  of3 << rotated.size() << std::endl;
  for (size_t kr=0; kr<rotated.size(); ++kr)
    of3 << rotated[kr] << std::endl;
  of3 << "410 1 0 4 0 0 255 255" << std::endl;
  of3 << "1" << std::endl;
  of3 << pnt2-0.5*len*normal2 << " " << pnt2+0.5*len*normal2 << std::endl;
  of3 << "410 1 0 4 0 255 0 255" << std::endl;
  of3 << "1" << std::endl;
  of3 << pnt2-0.5*len*Cx2 << " " << pnt2+0.5*len*Cx2 << std::endl;

  vector<Point> projected;
  double maxdp, avdp;
  RevEngUtils::projectToPlane(group, normal2, pnt2, projected, maxdp, avdp);
  std::ofstream ofp3("projected_pts.g2");
  ofp3 << "400 1 0 4 255 0 0 255" << std::endl;
  ofp3 << projected.size() << std::endl;
  for (size_t kr=0; kr<projected.size(); ++kr)
    ofp3 << projected[kr] << std::endl;
  circ2->writeStandardHeader(ofp3);
  circ2->write(ofp3);

  shared_ptr<Torus> tor1(new Torus(rad1, fabs(rd1), pnt1, normal1, Cy1));
  shared_ptr<Torus> tor2(new Torus(rad2, fabs(rd2), pnt2, normal2, Cy2));
  // tor1->writeStandardHeader(of);
  // tor1->write(of);
  tor2->writeStandardHeader(of);
  tor2->write(of);
  
  // Check accuracy
  double maxd1, avd1, maxd2, avd2;
  int num1, num2;
  //shared_ptr<ParamSurface> surf = cyl;
  vector<RevEngPoint*> in1, out1,in2, out2;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			  tor1, tol, maxd1, avd1, num1, in1, out1);
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
  			  tor2, tol, maxd2, avd2, num2, in2, out2);
  std::ofstream ofd("in_out.g2");
  ofd << "400 1 0 4 155 50 50 255" << std::endl;
  ofd << in2.size() << std::endl;
  for (size_t kr=0; kr<in2.size(); ++kr)
    ofd << in2[kr]->getPoint() << std::endl;
  ofd << "400 1 0 4 50 155 50 255" << std::endl;
  ofd << out2.size() << std::endl;
  for (size_t kr=0; kr<out2.size(); ++kr)
    ofd << out2[kr]->getPoint() << std::endl;
  
  int num = (int)group_points_.size();
  if (num1 > min_pt && num1 > num/2 && num1 > num2) //avd <= avlim && maxd <= maxlim)
    {
      found = true;
      setAccuracy(maxd1, avd1, num1);
      
      // // Limit cylinder with respect to bounding box
      // double gap = 1.0e-6;
      // Point xdir(1.0, 0.0, 0.0);
      // Point ydir(0.0, 1.0, 0.0);
      // Point zdir(0.0, 0.0, 1.0);
      // Point bbdiag = high - low;
      // CompositeModelFactory factory(gap, gap, 10.0*gap, 0.01, 0.05);
      // shared_ptr<SurfaceModel> boxmod(factory.createFromBox(low-2.0*bbdiag, xdir, ydir, 
      // 							    5*bbdiag[0], 5*bbdiag[1],
      // 							    5*bbdiag[2]));
      // vector<shared_ptr<ParamSurface> > sfs;
      // sfs.push_back(tor1);
      // shared_ptr<SurfaceModel> tormod(new SurfaceModel(gap, gap, 10.0*gap, 0.01,
      // 						       0.05, sfs));
      // vector<shared_ptr<SurfaceModel> > divmod = tormod->splitSurfaceModels(boxmod);
      
      std::cout << "Torus 1. N1: " << num << ", N2: " << num1 << ", max: " << maxd1 << ", av: " << avd1 << std::endl;
      // associated_sf_.clear();
      // for (int ka=0; ka<divmod[0]->nmbEntities(); ++ka)
      // 	{
      //shared_ptr<ParamSurface> tor1_2 = divmod[0]->getSurface(ka);
	  shared_ptr<HedgeSurface> hedge(new HedgeSurface(tor1, this));
	  for (size_t kh=0; kh<associated_sf_.size(); ++kh)
	    prevsfs.push_back(associated_sf_[kh]);
	  associated_sf_.push_back(hedge.get());
	  hedgesfs.push_back(hedge);
	  
	  tor1->writeStandardHeader(fileout);
	  tor1->write(fileout);
	// }
	}
  else if (num2 > min_pt && num2 > num/2)
    {
      found = true;
      setAccuracy(maxd2, avd2, num2);
      
      // // Limit cylinder with respect to bounding box
      // double gap = 1.0e-6;
      // Point xdir(1.0, 0.0, 0.0);
      // Point ydir(0.0, 1.0, 0.0);
      // Point zdir(0.0, 0.0, 1.0);
      // Point bbdiag = high - low;
      // CompositeModelFactory factory(gap, gap, 10.0*gap, 0.01, 0.05);
      // shared_ptr<SurfaceModel> boxmod(factory.createFromBox(low-2.0*bbdiag, xdir, ydir, 
      // 							    5*bbdiag[0], 5*bbdiag[1],
      // 							    5*bbdiag[2]));
      // vector<shared_ptr<ParamSurface> > sfs;
      // sfs.push_back(tor2);
      // shared_ptr<SurfaceModel> tormod(new SurfaceModel(gap, gap, 10.0*gap, 0.01,
      // 						       0.05, sfs));
      // vector<shared_ptr<SurfaceModel> > divmod = tormod->splitSurfaceModels(boxmod);
      
      std::cout << "Torus 2. N1: " << num << ", N2: " << num2 << ", max: " << maxd2 << ", av: " << avd2 << std::endl;
      // associated_sf_.clear();
      // for (int ka=0; ka<divmod[0]->nmbEntities(); ++ka)
      // 	{
      // 	  shared_ptr<ParamSurface> tor2_2 = divmod[0]->getSurface(ka);
	  shared_ptr<HedgeSurface> hedge(new HedgeSurface(tor2, this));
	  for (size_t kh=0; kh<associated_sf_.size(); ++kh)
	    prevsfs.push_back(associated_sf_[kh]);
	  associated_sf_.push_back(hedge.get());
	  hedgesfs.push_back(hedge);
	  
	  tor2->writeStandardHeader(fileout);
	  tor2->write(fileout);
	// }
    }
  int stop_break = 1;
  return found;
}
 
//===========================================================================
void RevEngRegion::analyseCylinderProperties(Point avvec, double angtol,
					  vector<RevEngPoint*>& in,
					  vector<RevEngPoint*> out)
//===========================================================================
{
  vector<double> midang(group_points_.size());
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Point curr = group_points_[ki]->minCurvatureVec();
      double ang = curr.angle(avvec);
      ang = std::min(ang, M_PI-ang);
      midang[ki] = ang;
      if (ang <= angtol)
	in.push_back(group_points_[ki]);
      else
	out.push_back(group_points_[ki]);
    }

  std::sort(midang.begin(), midang.end());
  int stop_break = 1;
}

//===========================================================================
const Point& RevEngRegion::pluckerAxis()
//===========================================================================
{
  Point dummy(0.0, 0.0, 0.0);

  // Matrix in minimization expression
  double M[6][6];
  int ki, kj, kk;
  for (ki=0; ki<6; ++ki)
    for (kj=0; kj<6; ++kj)
      M[ki][kj] = 0.0;

  double del = (double)group_points_.size();
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(group_points_[kr]);
      Point norm = pt->getMongeNormal();
      Vector3D xyz = pt->getPoint();
      Point pos(xyz[0], xyz[1], xyz[2]);
      Point vec = norm.cross(pos);
      vec.normalize();
      for (int ka=0; ka<3; ++ka)
	{
	  for (int kb=0; kb<3; ++kb)
	    {
	      M[ka][kb] += norm[ka]*norm[kb]/del;
	      M[3+ka][kb] += norm[ka]*vec[kb]/del;
	      M[ka][3+kb] += norm[ka]*vec[kb]/del;
	      M[3+ka][3+kb] += vec[ka]*vec[kb]/del;
	    }
	}
    }

  // Compute coefficients in cubic function
  // Sub-determinants
  vector<vector<int> > ixs;
  for (ki=0; ki<6; ++ki)
    for (kj=ki+1; kj<6; ++kj)
      for (kk=kj+1; kk<6; ++kk)
	{
	  vector<int> currix(3);
	  currix[0] = ki;
	  currix[1] = kj;
	  currix[2] = kk;
	  ixs.push_back(currix);
	}

  vector<double> Sd(ixs.size());
  int sgn1, sgn2, sgn3;
  for (size_t kr=0; kr<Sd.size(); ++kr)
    {
      Sd[kr] = 0.0;
      for (ki=0, sgn1=1; ki<3; ++ki, sgn1*=-1)
	{
	  for (kj=0, sgn2=1; kj<3; ++kj)
	    {
	      if (kj == ki)
		continue;
	      for (kk=0; kk<3; ++kk)
		{
		  if (kk == ki || kk == kj)
		    continue;
		  double tmp =
		    sgn1*sgn2*M[3][ixs[kr][ki]]*M[4][ixs[kr][kj]]*M[5][ixs[kr][kk]];
		  Sd[kr] += tmp;
		}
	      sgn2 *= -1;
	    }
	}
    }


  // Coefficients in first version of polynomial
  sgn1 = sgn2 = sgn3 = 1;
  double L0=0.0, L1=0.0, L2=0.0, L3=0.0, L12=0.0, L13=0.0, L23=0.0, L123=0.0;
  for (ki=0; ki<6; ++ki)
    {
      for (kj=0; kj<6; ++kj)
  	{
  	  if (kj == ki)
  	    continue;
  	  for (kk=0; kk<6; ++kk)
  	    {
  	      if (kk == ki || kk == kj)
  		continue;
	      
	      // Find determinant of sub matrix
	      size_t kr;
	      for (kr=0; kr<Sd.size(); ++kr)
		{
		  int d1, d2, d3;
		  for (d1=0; d1<3; ++d1)
		    if (ixs[kr][d1] == ki)
		      break;
		  for (d2=0; d2<3; ++d2)
		    if (ixs[kr][d2] == kj)
		      break;
		  for (d3=0; d3<3; ++d3)
		    if (ixs[kr][d3] == kk)
		      break;
		  if (d1 == 3 && d2 == 3 && d3 == 3)
		    break;
		}
	      if (kr == Sd.size())
		return dummy;
	      
  	      double tmp = sgn1*sgn2*sgn3*Sd[kr];
  	      if (ki==0 && kj == 1 && kk == 2)
  		L123 += tmp;
  	      else if (ki == 0 && kj == 1)
  		L12 += M[2][kk]*tmp;
  	      else if (ki == 0 && kk == 2)
  		L13 += M[1][kj]*tmp;
  	      else if (kj == 1 && kk == 2)
  		L23 += M[0][ki]*tmp;
  	      else if (ki == 0)
  		L1 += M[1][kj]*M[2][kk]*tmp;
  	      else if (kj == 1)
  		L2 += M[0][ki]*M[2][kk]*tmp;
  	      else if (kk == 2)
  		L3 += M[1][kj]*M[2][kk]*tmp;
  	      else
  		L0 += M[0][ki]*M[1][kj]*M[2][kk]*tmp;
  	      sgn3 += -1;
  	    }
  	  sgn2 *= -1;
  	}
      sgn1 *= -1;
    }

  // Reorganize coefficients to correspond to (x^3, x^2, x, 1)
  double cf[4];
  cf[0] = -L123;
  cf[1] = (M[0][0] + M[1][1] + M[2][2])*L123 + L12 + L13 + L23;
  cf[2] = -(M[0][0]*M[1][1] + M[0][0]*M[2][2] + M[1][1]*M[2][2])*L123 -
    (M[0][0] + M[1][1])*L12 - (M[0][0] + M[2][2])*L13 - (M[1][1] + M[2][2])*L23;
  cf[3] = M[0][0]*M[1][1]*M[2][2]*L123 + M[0][0]*M[1][1]*L12 + M[0][0]*M[2][2]*L13 +
    M[1][1]*M[2][2]*L23 + M[0][0]*L1 + M[1][1]*L2 + M[2][2]*L3 + L0;

    std::ofstream ofl("lambda.g2");
    int nsample = 101;
    double tmin = -10.0, tmax = 10.0;
    double tdel = (tmax - tmin)/(double)(nsample-1);
    double tpar;
    ofl << "410 1 0 4 100 100 55 255" << std::endl;
    ofl << nsample-1 << std::endl;
    double val1 = cf[0]*tmin*tmin*tmin + cf[1]*tmin*tmin + cf[2]*tmin + cf[3];
for (ki=1, tpar=tmin+tdel; ki<nsample; ++ki, tpar+=tdel)
{
  ofl << tpar-tdel << " " << val1 << " 0.0 ";
double val2 = cf[0]*tpar*tpar*tpar + cf[1]*tpar*tpar + cf[2]*tpar + cf[3];
ofl << tpar << " " << val2 << " 0.0" << std::endl;
 val1 = val2;
}
    ofl << "410 1 0 4 0 0 255 255" << std::endl;
    ofl << 1 << std::endl;
ofl << tmin << " 0.0 0.0 " << tmax << " 0.0 0.0" << std::endl;
 ofl << "400 1 0 4 255 0 0 255" << std::endl;
 ofl << "1" << std::endl;
 ofl << "0.0 0.0 0.0" << std::endl;

NEWMAT::Matrix Mmat;
Mmat.ReSize(6,6);
for (int ka=0; ka<6; ++ka)
  for (int kb=0; kb<6; ++kb)
    Mmat.element(ka,kb) = M[ka][kb];

  static NEWMAT::DiagonalMatrix diag;
  static NEWMAT::Matrix V;
  try {
    NEWMAT::SVD(Mmat, diag, Mmat, V);
  } catch(...) {
    std::cout << "Exception in SVD" << std::endl;
    return dummy;
  }

double lambda[6];
double eigen[6][6];
for (int ka=0; ka<6; ++ka)
  {
    lambda[ka] = diag.element(ka,ka);
    for (int kb=0; kb<6; ++kb)
      eigen[ka][kb] = V.element(ka,kb);
  }

  return dummy;
}


//===========================================================================
void
RevEngRegion::splitFromSurfaceNormals(vector<RevEngPoint*>& smallrad,
				      vector<vector<RevEngPoint*> >& separate_group)
//===========================================================================
{
  // Associate radius to normals according to density in unit spehere
  extendWithGaussRad();
  extendWithGaussRad2();

  // Sort points according to normal radius and remove points with a small radius
  double radlim = 0.5;
  for (size_t kr=0; kr<group_points_.size(); )
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(group_points_[kr]);
      double Grad = group_points_[kr]->getGaussRad();
      if (Grad < radlim)
	{
	  group_points_[kr]->unsetRegion();
	  smallrad.push_back(group_points_[kr]);
	  group_points_.erase(group_points_.begin()+kr);
	}
      else
	++kr;
    }

  if (group_points_.size() == 0)
    {
      group_points_.insert(group_points_.end(), smallrad.begin(), smallrad.end());
      smallrad.clear();
    }
  
  std::ofstream of3("curr_region_split.g2");
  if (group_points_.size() > 0)
    {
      of3 << "400 1 0 4 155 0 100 255" << std::endl;
      of3 << group_points_.size() << std::endl;
      for (size_t kr=0; kr<group_points_.size(); ++kr)
	of3 << group_points_[kr]->getPoint() << std::endl;
    }
  if (smallrad.size() > 0)
    {
      of3 << "400 1 0 4 0 155 100 255" << std::endl;
      of3 << smallrad.size() << std::endl;
      for (size_t kr=0; kr<smallrad.size(); ++kr)
	of3 << smallrad[kr]->getPoint() << std::endl;
    }

  // Split remaining region points into disjunct groupe
  if (smallrad.size() > 0)
    splitRegion(separate_group);
  int stop_break = 1;
}

//===========================================================================
void RevEngRegion::splitRegion(vector<vector<RevEngPoint*> >& separate_groups)
//===========================================================================
{
  vector<vector<RevEngPoint*> > connected;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      if (group_points_[ki]->visited())
	continue;
      vector<RevEngPoint*> curr_group;
      group_points_[ki]->fetchConnected(this, (int)group_points_.size(), curr_group);
      connected.push_back(curr_group);
    }

  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      group_points_[ki]->unsetVisited();
    }

  std::ofstream of("curr_region_split2.g2");
  for (size_t kj=0; kj<connected.size(); ++kj)
    {
      of << "400 1 0 0" << std::endl;
      of << connected[kj].size() << std::endl;
      for (size_t kr=0; kr<connected[kj].size(); ++kr)
	of << connected[kj][kr]->getPoint() << std::endl;
    }
  
  group_points_.clear();
  group_points_ = connected[0];
  
  // Update bounding box and principal curvature summary
  maxk2_ = std::numeric_limits<double>::lowest();
  mink2_ = std::numeric_limits<double>::max();
  maxk1_ = std::numeric_limits<double>::lowest();
  mink1_ = std::numeric_limits<double>::max();
  bbox_ = BoundingBox(3);
  for  (size_t kj=0; kj<group_points_.size(); ++kj)
    {
      double k1 = group_points_[kj]->minPrincipalCurvature();
      double k2 = group_points_[kj]->maxPrincipalCurvature();
      mink1_ = std::min(mink1_, fabs(k1));
      maxk1_ = std::max(maxk1_, fabs(k1));
      mink2_ = std::min(mink2_, fabs(k2));
      maxk2_ = std::max(maxk2_, fabs(k2));
      Vector3D point = group_points_[kj]->getPoint();
      Point point2(point[0], point[1], point[2]);
      bbox_.addUnionWith(point2);
    }

  normalcone_ = DirectionCone(group_points_[0]->getPCANormal());
  for  (size_t kj=1; kj<group_points_.size(); ++kj)
    {
      normalcone_.addUnionWith(group_points_[kj]->getPCANormal());
    }
  

  for (size_t kj=1; kj<connected.size(); ++kj)
    separate_groups.push_back(connected[kj]);
  
}

//===========================================================================
void RevEngRegion::addPoint(RevEngPoint* point)
//===========================================================================
{
  group_points_.push_back(point);
  Vector3D point2 = point->getPoint();
  Point point3(point2[0], point2[1], point2[2]);
  bbox_.addUnionWith(point3);
  normalcone_.addUnionWith(point->getPCANormal());
}

//===========================================================================
void RevEngRegion::extendWithGaussRad()
//===========================================================================
{
  // Distribute points according to normals position in unit sphere
  double eps = 1.0e-6;
  int n1=81, n2=41;
  vector<vector<vector<RevEngPoint*> > > div(n2);
  for (int ka=0; ka<n2; ++ka)
    div[ka].resize(n1);
  double p1=-M_PI, p2=M_PI;
  double p3=-0.5*M_PI, p4=0.5*M_PI;
  double del1 = (p2-p1)/(double)(n1);
  double del2 = (p4-p3)/(double)(n2);

  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Point norm = group_points_[ki]->getMongeNormal();
      // Point vec(norm[0], norm[1]);
      // vec.normalize_checked();
      double phi2 = asin(std::max(-1.0, std::min(1.0, norm[2])));
      double h1 = cos(phi2);
      double phi1;
      if (norm[1] > 0.0)
      	phi1 = (h1 >= eps) ?
      	  acos(std::max(-1.0, std::min(1.0, norm[0]/h1))) : 0.0;
      else if (norm[0] > 0.0)
      	phi1 =  (h1 >= eps) ?
      	  asin(std::max(-1.0, std::min(1.0, norm[1]/h1))) : 0.0;
      else
      	{
      	  phi1 = (h1 >= eps) ?
      	    acos(std::max(-1.0, std::min(1.0, norm[0]/h1))) : 0.0;
      	  phi1 *= -1.0;
      	}
      // if (vec[1] > 0.0)
      // 	phi1 = (h1 >= eps) ?
      // 	  acos(std::max(-1.0, std::min(1.0, vec[0]))) : 0.0;
      // else if (vec[0] > 0.0)
      // 	phi1 =  (h1 >= eps) ?
      // 	  asin(std::max(-1.0, std::min(1.0, vec[1]))) : 0.0;
      // else
      // 	{
      // 	  phi1 = (h1 >= eps) ?
      // 	    std::max(-1.0, std::min(1.0, acos(vec[0]))) : 0.0;
      // 	  phi1 *= -1.0;
      // 	}

      double u1 = (phi1 - p1)/(p2 - p1);
      double u2 = (phi2 - p3)/(p4 - p3);
      double v1 = u1*n1;
      double v2 = u2*n2;
      int ix1 = (int)(v1);
      //ix1 += (v1 - (double)ix1 >=  0.5);
      ix1 = std::min(ix1, n1-1);
      int ix2 = (int)(v2);
      //ix2 += (v2 - (double)ix2 >= 0.5);
      ix2 = std::min(ix2, n2-1);
      div[ix2][ix1].push_back(group_points_[ki]);
    }

  // Statistics on distribution
  int minn = (int)group_points_.size();
  int maxn = 0;
  vector<int> nmbn(n1*n2);
  int kc=0;
  for (int kb=0; kb<n2; ++kb)
    for (int ka=0; ka<n1; ++ka, ++kc)
      {
	nmbn[kc] = (int)div[kb][ka].size();
	if (div[kb][ka].size() > 0)
	  minn = std::min(minn, (int)div[kb][ka].size());
	maxn = std::max(maxn, (int)div[kb][ka].size());
      }
  std::sort(nmbn.begin(), nmbn.end());
  
  // int level = 3;
  // vector<int> lim(level-1);
  double lim = std::min(0.9*minn + 0.1*maxn, 0.95*maxn);
  double minrad = 0.3;
  double maxrad = 1.0;
  //double rdel = (maxrad - minrad)/(double)(level-1);
  // lim[0] = (int)(0.9*minn + 0.1*maxn);
  // lim[1] = (int)(0.1*minn + 0.9*maxn);
  //  int ldel = (maxn - minn)/level;
  // for (int ka=1; ka<level; ++ka)
  //   lim[ka-1] = minn + ka*ldel;

  for (int kb=0; kb<n2; ++kb)
    for (int ka=0; ka<n1; ++ka)
      {
	int npt = (int)div[kb][ka].size();
	if (npt == 0)
	  continue;
	// size_t kj;
	// for (kj=0; kj<lim.size(); ++kj)
	//   if (lim[kj] > npt)
	//     break;

	std::ofstream of("curr_norm_group.g2");
	double Grad = (npt < lim) ? minrad : maxrad; //minrad + kj*rdel;
	of << "400 1 0 4 0 0 0 255" << std::endl;
	of << div[kb][ka].size() << std::endl;
	for (size_t ki=0; ki<div[kb][ka].size(); ++ki)
	  {
	    RevEngPoint *pt = div[kb][ka][ki];
	    pt->setGaussRad(Grad);
	    Point norm = pt->getMongeNormal();
	    of << Grad*norm << std::endl;
	  }
	int stop_break0 = 1;
      }

  std::ofstream of2("curr_normals2.g2");
  vector<Point> g1, g2;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(group_points_[kr]);
      Point norm = pt->getMongeNormal();
      double Grad = pt->getGaussRad();
      if (Grad < 0.5)
	g2.push_back(norm);
      else
	g1.push_back(norm);
    }
  of2 << "400 1 0 4 100  0 155 255" << std::endl;
  of2 << g1.size() << std::endl;
  for (size_t kr=0; kr<g1.size(); ++kr)
    of2 << g1[kr] << std::endl;

  of2 << "400 1 0 4 0 100  155 255" << std::endl;
  of2 << g2.size() << std::endl;
  for (size_t kr=0; kr<g2.size(); ++kr)
    of2 << g2[kr] << std::endl;

  Sphere sph(1.0, Point(0.0, 0.0, 0.0), Point(0.0, 0.0, 1.0),
	     Point(1.0, 0.0, 0.0));
  sph.writeStandardHeader(of2);
  sph.write(of2);
  int stop_break;
}


//===========================================================================
void RevEngRegion::extendWithGaussRad2()
//===========================================================================
{
  // Distribute points according to direction of minimu curvature position in unit sphere
  double eps = 1.0e-6;
  int n1=81, n2=41;
  vector<vector<vector<RevEngPoint*> > > div(n2);
  for (int ka=0; ka<n2; ++ka)
    div[ka].resize(n1);
  double p1=-M_PI, p2=M_PI;
  double p3=-0.5*M_PI, p4=0.5*M_PI;
  double del1 = (p2-p1)/(double)(n1);
  double del2 = (p4-p3)/(double)(n2);

  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Point cvec[2];
      cvec[0] = group_points_[ki]->minCurvatureVec();
      cvec[1] = -1*cvec[0];

      for (int ka=0; ka<2; ++ka)
	{
	  double phi2 = asin(std::max(-1.0, std::min(1.0, cvec[ka][2])));
	  double h1 = cos(phi2);
	  double phi1;
	  if (cvec[ka][1] > 0.0)
	    phi1 = (h1 >= eps) ?
	      acos(std::max(-1.0, std::min(1.0, cvec[ka][0]/h1))) : 0.0;
	  else if (cvec[ka][0] > 0.0)
	    phi1 =  (h1 >= eps) ?
	      asin(std::max(-1.0, std::min(1.0, cvec[ka][1]/h1))) : 0.0;
	  else
	    {
	      phi1 = (h1 >= eps) ?
		acos(std::max(-1.0, std::min(1.0, cvec[ka][0]/h1))) : 0.0;
	      phi1 *= -1.0;
	    }

	  double u1 = (phi1 - p1)/(p2 - p1);
	  double u2 = (phi2 - p3)/(p4 - p3);
	  double v1 = u1*n1;
	  double v2 = u2*n2;
	  int ix1 = (int)(v1);
	  //ix1 += (v1 - (double)ix1 >=  0.5);
	  ix1 = std::min(ix1, n1-1);
	  int ix2 = (int)(v2);
	  //ix2 += (v2 - (double)ix2 >= 0.5);
	  ix2 = std::min(ix2, n2-1);
	  div[ix2][ix1].push_back(group_points_[ki]);
	}
    }

  // Statistics on distribution
  int minn = (int)group_points_.size();
  int maxn = 0;
  vector<int> nmbn(n1*n2);
  int kc=0;
  for (int kb=0; kb<n2; ++kb)
    for (int ka=0; ka<n1; ++ka, ++kc)
      {
	nmbn[kc] = (int)div[kb][ka].size();
	if (div[kb][ka].size() > 0)
	  minn = std::min(minn, (int)div[kb][ka].size());
	maxn = std::max(maxn, (int)div[kb][ka].size());
      }
  std::sort(nmbn.begin(), nmbn.end());
  
  // int level = 3;
  // vector<int> lim(level-1);
  double lim = std::min(0.95*minn + 0.05*maxn, 0.9*maxn);
  double minrad = 0.2;
  double maxrad = 1.0;
  //double rdel = (maxrad - minrad)/(double)(level-1);
  // lim[0] = (int)(0.9*minn + 0.1*maxn);
  // lim[1] = (int)(0.1*minn + 0.9*maxn);
  //  int ldel = (maxn - minn)/level;
  // for (int ka=1; ka<level; ++ka)
  //   lim[ka-1] = minn + ka*ldel;

  for (int kb=0; kb<n2; ++kb)
    for (int ka=0; ka<n1; ++ka)
      {
	int npt = (int)div[kb][ka].size();
	if (npt == 0)
	  continue;
	// size_t kj;
	// for (kj=0; kj<lim.size(); ++kj)
	//   if (lim[kj] > npt)
	//     break;

	std::ofstream of("curr_cvec_group.g2");
	double Grad = (npt < lim) ? minrad : maxrad; //minrad + kj*rdel;
	of << "400 1 0 4 0 0 0 255" << std::endl;
	of << div[kb][ka].size() << std::endl;
	for (size_t ki=0; ki<div[kb][ka].size(); ++ki)
	  {
	    RevEngPoint *pt = div[kb][ka][ki];
	    double Grad0 = pt->getGaussRad();
	    if (Grad0 < maxrad)
	      continue;
	    pt->setGaussRad(Grad);
	    Point vec = pt->minCurvatureVec();
	    of << Grad*vec << std::endl;
	  }
	int stop_break0 = 1;
      }

  std::ofstream of2("curr_cvec.g2");
		      vector<Point> g1, g2, g3;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(group_points_[kr]);
      Point vec = pt->minCurvatureVec();
      double Grad = pt->getGaussRad();
      if (Grad < 0.25)
	g3.push_back(vec);
      else if (Grad < 0.5)
	g2.push_back(vec);
      else
	g1.push_back(vec);
    }
  of2 << "400 1 0 4 100  0 155 255" << std::endl;
  of2 << g1.size() << std::endl;
  for (size_t kr=0; kr<g1.size(); ++kr)
    of2 << g1[kr] << std::endl;

  of2 << "400 1 0 4 0 100  155 255" << std::endl;
  of2 << g2.size() << std::endl;
  for (size_t kr=0; kr<g2.size(); ++kr)
    of2 << g2[kr] << std::endl;

  of2 << "400 1 0 4 50 150  55 255" << std::endl;
  of2 << g3.size() << std::endl;
  for (size_t kr=0; kr<g3.size(); ++kr)
    of2 << g3[kr] << std::endl;

  Sphere sph(1.0, Point(0.0, 0.0, 0.0), Point(0.0, 0.0, 1.0),
	     Point(1.0, 0.0, 0.0));
  sph.writeStandardHeader(of2);
  sph.write(of2);
  int stop_break;
}


//===========================================================================
void RevEngRegion::growWithSurf(int max_nmb, double tol,
				vector<RevEngRegion*>& grown_regions,
				vector<HedgeSurface*>& adj_surfs)
//===========================================================================
{
  double eps = 1.0e-6;
  if (associated_sf_.size() == 0)
    return;  // No surface with which to check growt

  int sfcode;
  ClassType classtype = associated_sf_[0]->instanceType(sfcode);
  shared_ptr<ParamSurface> surf = associated_sf_[0]->surface();
  shared_ptr<BoundedSurface> bdsurf =
    dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf);
  if (bdsurf.get())
    surf = bdsurf->underlyingSurface();

  vector<RevEngPoint*> visited;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      vector<ftSamplePoint*> next = group_points_[ki]->getNeighbours();
      for (size_t kj=0; kj<next.size(); ++kj)
	{
	  RevEngPoint *curr = dynamic_cast<RevEngPoint*>(next[kj]);
	  if (curr->isEdge())
	    continue;  // Do not include points labelled as edge
	  
	  RevEngRegion *adj_reg = curr->region();
	  if (adj_reg)
	    {
	      if (adj_reg == this || adj_reg->visited())
		continue;
	      
	      adj_reg->setVisited(true);
	      
	      // Check criteria for growt
	      int num = adj_reg->numPoints();
	      if (/*adj_reg->hasSurface() ||*/ num > max_nmb)
		continue;

	      // std::ofstream of("curr_extend.g2");
	      // of << "400 1 0 4 0 255 0 255" << std::endl;
	      // of << group_points_.size() << std::endl;
	      // for (size_t kh=0; kh<group_points_.size(); ++kh)
	      // 	of << group_points_[kh]->getPoint() << std::endl;

	      // of << "400 1 0 4 255 0 0 255" << std::endl;
	      // of << num << std::endl;
	      // for (int ka=0; ka<num; ++ka)
	      // 	of << adj_reg->getPoint(ka)->getPoint() << std::endl;

	      // surf->writeStandardHeader(of);
	      // surf->write(of);
	      
	      if (!adj_reg->isCompatible(classtype, sfcode))
		continue;
	      
	      // Check accuracy
	      double maxd, avd;
	      int num_inside;
	      double maxd_init, avd_init;
	      int num_inside_init;
	      adj_reg->getAccuracy(maxd_init, avd_init, num_inside_init);
	      
	      vector<RevEngPoint*> in, out;
	      RevEngUtils::distToSurf(adj_reg->pointsBegin(), adj_reg->pointsEnd(),
				      surf, tol, maxd, avd, num_inside, in, out);
	      if ((adj_reg->hasSurface() && num_inside > num_inside_init && avd < avd_init) ||
		  ((!adj_reg->hasSurface()) && num_inside > num/2 && avd < tol))
		{
		  // Criteria must be updated
		  // Include adjacent region in present
		  bbox_.addUnionWith(adj_reg->boundingBox());
		  normalcone_.addUnionWith(adj_reg->getNormalCone());
		  maxdist_ = std::max(maxdist_, maxd);
		  double div = (double)((int)group_points_.size() + num);
		  avdist_ = (group_points_.size()*avdist_ + num*avd)/div;

		  double mink1, maxk1, mink2, maxk2;
		  adj_reg->getPrincipalCurvatureInfo(mink1, maxk1, mink2, maxk2);
		  mink1_ = std::min(mink1_, mink1);
		  maxk1_ = std::max(maxk1_, maxk1);
		  mink2_ = std::min(mink2_, mink2);
		  maxk2_ = std::max(maxk2_, maxk2);

		  for (auto it=adj_reg->pointsBegin(); it != adj_reg->pointsEnd(); ++it)
		    (*it)->setRegion(this);
		  group_points_.insert(group_points_.end(), adj_reg->pointsBegin(),
				       adj_reg->pointsEnd());
		  grown_regions.push_back(adj_reg);

		  int num_sf = adj_reg->numSurface();
		  for (int ka=0; ka<num_sf; ++ka)
		    adj_surfs.push_back(adj_reg->getSurface(ka));
		}
	    }
	  else
	    {
	      if (curr->visited())
		continue;
	      curr->setVisited();
	      visited.push_back(curr);
	      
	      // Check accuracy
	      Vector3D xyz = curr->getPoint();
	      Point pnt(xyz[0], xyz[1], xyz[2]);
	      double upar, vpar, dist;
	      Point close;
	      surf->closestPoint(pnt, upar, vpar, close, dist, eps);
	      if (dist <= tol)
		{
		  double k1 = curr->minPrincipalCurvature();
		  double k2 = curr->maxPrincipalCurvature();
		  mink1_ = std::min(mink1_, fabs(k1));
		  maxk1_ = std::max(maxk1_, fabs(k1));
		  mink2_ = std::min(mink2_, fabs(k2));
		  maxk2_ = std::max(maxk2_, fabs(k2));
		  bbox_.addUnionWith(pnt);
		  normalcone_.addUnionWith(curr->getPCANormal());
		  group_points_.push_back(curr);
 		}
	    }
	  
	}
    }
  for (size_t kj=0; kj<visited.size(); ++kj)
    visited[kj]->unsetVisited();
}


//===========================================================================
bool RevEngRegion::isCompatible(ClassType classtype, int sfcode)
//===========================================================================
{
  double anglim = 0.1;
  if (classtype == Class_Plane)
    {
      if (planartype() || (normalcone_.angle() <= anglim &&
			   (!normalcone_.greaterThanPi())))
	return true;
      else
	return false;
    }
  else
  if (classtype == Class_Cylinder)
    {
      return true;  // Must mike proper criteria
    }
    return false;  // To be extended
}
