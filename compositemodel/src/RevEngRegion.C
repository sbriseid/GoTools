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
#include "GoTools/geometry/Sphere.h"
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
    mink1_(0.0), maxk1_(0.0), mink2_(0.0), maxk2_(0.0)
{
  bbox_ = BoundingBox(3);
}

//===========================================================================
RevEngRegion::RevEngRegion(int classification_type)
//===========================================================================
  : classification_type_(classification_type), associated_sf_(0),
    mink1_(0.0), maxk1_(0.0), mink2_(0.0), maxk2_(0.0)
{
  bbox_ = BoundingBox(3);
}

//===========================================================================
RevEngRegion::RevEngRegion(int classification_type,
			   vector<RevEngPoint*> & points)
//===========================================================================
  : group_points_(points), classification_type_(classification_type),
    associated_sf_(0)
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
      group_points_[ki]->fetchClosePoints2(radius, min_next, nearpts);

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
  seed->fetchClosePoints2(radius, min_close, nearpts, this);
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

  shared_ptr<Plane> surf3(new Plane(pos3, normal3));
  
  // Check accuracy
  double maxdist, avdist;
  double maxdist1, avdist1;
  double maxdist2, avdist2;
  double maxdist3, avdist3;
  int num_inside, num_inside1, num_inside2, num_inside3;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			  surf, tol, maxdist, avdist, num_inside);
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			  surf1, tol, maxdist1, avdist1, num_inside1);
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			  surf2, tol, maxdist2, avdist2, num_inside2);
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			  surf3, tol, maxdist3, avdist3, num_inside3);

  Point low = bbox_.low();
  Point high = bbox_.high();
  double len = low.dist(high);
  surf->setParameterBounds(-len, -len, len, len);

  int num = (int)group_points_.size();
  if (num_inside3 > min_pt && num_inside3 > num/2 && avdist3 <= tol) 
    {
      found = true;
      setAccuracy(maxdist3, avdist3, num_inside3);
      
      std::cout << "Plane. N1: " << num << ", N2: " << num_inside3 << ", max: " << maxdist3 << ", av: " << avdist3 << std::endl;

      shared_ptr<HedgeSurface> hedge(new HedgeSurface(surf3, this));
      associated_sf_.push_back(hedge.get());
      hedgesfs.push_back(hedge);
	  
      surf->writeStandardHeader(fileout);
      surf->write(fileout);
      impl_->visualize(group_points_, fileout);
    }

  return found;
}


//===========================================================================
void RevEngRegion::analyseNormals(double tol)
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

  Point pos, normal;
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

  double radius;
  RevEngUtils::computeRadius(vec, axis, Cx, Cy, radius);
  int stop_break = 1;

}

//===========================================================================
bool RevEngRegion::extractCylinder(double tol, int min_pt,
				   vector<shared_ptr<HedgeSurface> >& hedgesfs,
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
  
  analyseNormals(tol);
  
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

  double rad;
  Point pnt;
  Point low = bbox_.low();
  Point high = bbox_.high();
  RevEngUtils::computeCylPosRadius(group, low, high,
				   axis, Cx, Cy, pnt, rad);
  // double Cmat[3][3];
  // for (int ka=0; ka<3; ++ka)
  //   for (int kb=0; kb<3; ++kb)
  //     {
  // 	Cmat[ka][kb] = 0.0;
  // 	for (size_t kj=0; kj<group_points_.size(); ++kj)
  // 	  {
  // 	    RevEngPoint *pt = dynamic_cast<RevEngPoint*>(group_points_[kj]);
  // 	    Point norm1 = pt->getMongeNormal();
  // 	    Point norm2 = pt->getPCANormal();
  // 	    Point norm = pt->getTriangNormal();
  // 	    Cmat[ka][kb] += norm[ka]*norm[kb];
  // 	  }
  //     }

  // // Compute singular values
  // NEWMAT::Matrix nmat;
  // nmat.ReSize(3, 3);
  // for (int ka = 0; ka < 3; ++ka) {
  //   for (int kb = 0; kb < 3; ++kb) {
  //     nmat.element(ka, kb) = Cmat[ka][kb];
  //   }
  // }
      
  // static NEWMAT::DiagonalMatrix diag;
  // static NEWMAT::Matrix V;
  // try {
  //   NEWMAT::SVD(nmat, diag, nmat, V);
  // } catch(...) {
  //   std::cout << "Exception in SVD" << std::endl;
  //   exit(-1);
  // }
  // Cx = Point(V.element(0,0), V.element(1,0), V.element(2,0));
  // Cy = Point(V.element(0,1), V.element(1,1), V.element(2,1));
  // axis = Point(V.element(0,2), V.element(1,2), V.element(2,2));

  // Center, radius and bounding box
  // double fac = 0.5;

  // // if (maxk2_ - mink2_ < fac*maxk2_)
  // //   return false;   // Not a likely cylinder
  // double Amat[3][3];
  // double bvec[3];
  // vector<vector<double> > A1(group_points_.size());
  // vector<double> b1(group_points_.size());
  // for (size_t kj=0; kj<group_points_.size(); ++kj)
  //   A1[kj].resize(3);
      
  // for (int ka=0; ka<3; ++ka)
  //   {
  //     for (int kb=0; kb<3; ++kb)
  // 	Amat[ka][kb] = 0.0;
  //     bvec[ka] = 0.0;
  //   }

  // Point low(std::numeric_limits<double>::max(),
  // 	    std::numeric_limits<double>::max(),
  // 	    std::numeric_limits<double>::max());
  // Point high(std::numeric_limits<double>::lowest(),
  // 	     std::numeric_limits<double>::lowest(),
  // 	     std::numeric_limits<double>::lowest());
  
  // vector<double> xyzpoints;
  // xyzpoints.reserve(group_points_.size()*3);
  // for (size_t kj=0; kj<group_points_.size(); ++kj)
  //   {
  //     Vector3D curr = group_points_[kj]->getPoint();
  //     Point curr2(curr[0], curr[1], curr[2]);
  //     xyzpoints.insert(xyzpoints.end(), curr2.begin(), curr2.end());
  //     for (int kb=0; kb<3; ++kb)
  // 	{
  // 	  low[kb] = std::min(low[kb],curr2[kb]);
  // 	  high[kb] = std::max(high[kb],curr2[kb]);
  // 	}
  //     double pxy[3];
  //     double px = curr2*Cx;
  //     double py = curr2*Cy;
  //     pxy[0] = 2*px;
  //     pxy[1] = 2*py;
  //     pxy[2] = 1.0;
  //     double plen2 = px*px + py*py;
  //     A1[kj][0] = 2*px;
  //     A1[kj][1] = 2*py;
  //     A1[kj][2] = 1.0;
  //     b1[kj] = plen2;
  //     for (int ka=0; ka<3; ++ka)
  // 	{
  // 	  for (int kb=0; kb<3; ++kb)
  // 	    Amat[ka][kb] += pxy[ka]*pxy[kb];
  // 	  bvec[ka] += pxy[ka]*plen2;
  // 	}
  //   }

  // double Amat2[3][3];
  // double bvec2[3];
  // for (int ka=0; ka<3; ++ka)
  //   {
  //     for (int kb=0; kb<3; ++kb)
  // 	Amat2[ka][kb] = 0.0;
  //     bvec2[ka] = 0.0;
  //   }
  // for (int ka=0; ka<3; ++ka)
  //   {
  //     for (int kb=0; kb<3; ++kb)
  // 	{
  // 	  for (size_t kr=0; kr<group_points_.size(); ++kr)
  // 	    Amat2[ka][kb] += A1[kr][ka]*A1[kr][kb];
  // 	}
  //     for (size_t kr=0; kr<group_points_.size(); ++kr)
  // 	bvec2[ka] += A1[kr][ka]*b1[kr];
  //   }

  // double detA = 0.0;
  // double bx[3];
  // bx[0] = bx[1] = bx[2] = 0.0;
  // int sgn = 1;
  // for (int kb=0; kb<3; ++kb, sgn*=(-1))
  //   {
  //     int ka1 = (kb == 0);
  //     int ka2 = 2 - (kb == 2);
  //     detA += sgn*Amat[0][kb]*(Amat[1][ka1]*Amat[2][ka2]-Amat[2][ka1]*Amat[1][ka2]);
  //     bx[0] += sgn*bvec[kb]*(Amat[1][ka1]*Amat[2][ka2]-Amat[2][ka1]*Amat[1][ka2]);
  //     bx[1] += sgn*Amat[0][kb]*(bvec[ka1]*Amat[2][ka2]-bvec[ka2]*Amat[2][ka1]);
  //     bx[2] += sgn*Amat[0][kb]*(Amat[1][ka1]*bvec[ka2]-Amat[1][ka2]*bvec[ka1]);
  //   }
  // double sx = bx[0]/detA;
  // double sy = bx[1]/detA;
  // double r2 = bx[2]/detA;

  // Point pos = sx*Cx + sy*Cy;
  // Point pos2 = sy*Cx + sx*Cy;
  //     //std::cout << "pos dist: " << pos.dist(pos2) << std::endl;
  // double rad = sqrt(r2 + sx*sx + sy*sy);
  // double len = low.dist(high);
  // Point mid = 0.5*(low + high);
  // Point vec = mid - pos;
  // Point ax = axis;
  // ax.normalize();
  // Point pnt = pos + (vec*ax)*ax;

  // Point ac(0.0, 0.0, 0.0);
  // double div = 0.0;
  // for (int ka=0; ka<(int)group_points_.size(); ++ka)
  //   {
  //     Point curr(&xyzpoints[3*ka], &xyzpoints[3*(ka+1)]);
  //     double ti = (curr - pnt)*axis;
  //     ac += (ti*(curr - pnt));
  //     div += (ti*ti);
  //   }
  // Point ax2 = (div > eps) ? ac/div : axis;

  // pluckerAxis();
  
  shared_ptr<Cylinder> cyl(new Cylinder(rad, pnt, axis, Cy));
  double len = low.dist(high);
  cyl->setParamBoundsV(-len, len);
  //PointCloud3D cloud(&xyzpoints[0], xyzpoints.size()/3);

  // Check accuracy
  double maxd, avd;
  int num2;
  //shared_ptr<ParamSurface> surf = cyl;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			  cyl, tol, maxd, avd, num2);
  // int num = cloud.numPoints();
  // double avd = 0.0;
  // double maxd = 0.0;
  // vector<double> xyz2;
  // for (int ka=0; ka<num; ++ka)
  //   {
  //     Point curr(&xyzpoints[3*ka], &xyzpoints[3*(ka+1)]);

  //     double upar, vpar, dist;
  //     Point close;
  //     cyl->closestPoint(curr, upar, vpar, close, dist, eps);
  //     maxd = std::max(maxd, dist);
  //     avd += dist;
  //     if (dist <= tol)
  // 	xyz2.insert(xyz2.end(), curr.begin(), curr.end());
  //   }
  // avd /= (double)num;
  // PointCloud3D cloud2(&xyz2[0], xyz2.size()/3);
  // int num2 = cloud2.numPoints();

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
      for (int ka=0; ka<divcyl[0]->nmbEntities(); ++ka)
	{
	  shared_ptr<ParamSurface> cyl2 = divcyl[0]->getSurface(ka);
	  shared_ptr<HedgeSurface> hedge(new HedgeSurface(cyl2, this));
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
NEWMAT::SVD(Mmat, diag, Mmat, V);

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
      group_points_[ki]->fetchConnected(this, curr_group);
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
  int n1=80, n2=40;
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
      double phi2 = asin(norm[2]);
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
	    std::max(-1.0, std::min(1.0, acos(norm[0]/h1))) : 0.0;
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

  // Statistics on distribution
  int minn = (int)group_points_.size();
  int maxn = 0;
  for (int kb=0; kb<n2; ++kb)
    for (int ka=0; ka<n1; ++ka)
      {
	if (div[kb][ka].size() > 0)
	  minn = std::min(minn, (int)div[kb][ka].size());
	maxn = std::max(maxn, (int)div[kb][ka].size());
}

// int level = 3;
// vector<int> lim(level-1);
double lim = 0.95*minn + 0.05*maxn;
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

	//std::ofstream of("curr_norm_group.g2");
	double Grad = (npt < lim) ? minrad : maxrad; //minrad + kj*rdel;
	// of << "400 1 0 4 0 0 0 255" << std::endl;
	// of << div[kb][ka].size() << std::endl;
	for (size_t ki=0; ki<div[kb][ka].size(); ++ki)
	  {
	    RevEngPoint *pt = div[kb][ka][ki];
	    pt->setGaussRad(Grad);
	    Point norm = pt->getMongeNormal();
	    //of << Grad*norm << std::endl;
	  }
	int stop_break0 = 1;
	     }

  std::ofstream of2("curr_normals2.g2");
  of2 << "400 1 0 4 100  0 155 255" << std::endl;
  of2 << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(group_points_[kr]);
      Point norm = pt->getMongeNormal();
      double Grad = pt->getGaussRad();
      of2 << Grad*norm << std::endl;
    }
  Sphere sph(1.0, Point(0.0, 0.0, 0.0), Point(0.0, 0.0, 1.0),
	     Point(1.0, 0.0, 0.0));
  sph.writeStandardHeader(of2);
  sph.write(of2);
  int stop_break;
}
