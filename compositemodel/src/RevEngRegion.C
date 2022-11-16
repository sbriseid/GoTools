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
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/Curvature.h"
#include "GoTools/creators/SmoothCurve.h"
#include "newmat.h"
#include "newmatap.h"
#include <vector>
#include <fstream>

using namespace Go;
using std::vector;
using std::set;
  
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
void RevEngRegion::updateRegion(double approx_tol, double anglim,
				vector<RevEngRegion*>& adapted_regions,
				vector<shared_ptr<RevEngRegion> >& outdiv_regions)
//===========================================================================
{
  if (group_points_.size() < 5)
    return;

  std::ofstream of("region_grow.g2");
  writeRegionInfo(of);

  double eps = 1.0e-6;
  shared_ptr<SplineSurface> surf = surfApprox(group_points_, bbox_);
  std::ofstream of2("approx_sf.g2");
  surf->writeStandardHeader(of2);
  surf->write(of2);
  
  // Check if the accuracy is sufficient as a basis for growth
      double avdist = 0.0;
  double wgt = 1.0/(double)group_points_.size();
  vector<RevEngPoint*> in, out;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Vector3D xyz = group_points_[ki]->getPoint();
      Point pos(xyz[0], xyz[1], xyz[2]);
      double upar, vpar, dist;
      Point close;
      surf->closestPoint(pos, upar, vpar, close, dist, eps);
      avdist += wgt*dist;
      if (dist < approx_tol)
	{
	  Point normal;
	  surf->normal(normal, upar, vpar);
	  double ang = normal.angle(group_points_[ki]->getPCANormal());
	  ang = std::min(ang, M_PI-ang);
	  if (ang < anglim)
	    in.push_back(group_points_[ki]);
	  else
	    out.push_back(group_points_[ki]);
	}
      else
	out.push_back(group_points_[ki]);
    }

  if (in.size() < group_points_.size()/2 || avdist > approx_tol)
    return;

  // Grow
  vector<RevEngPoint*> visited;
  for (int ka=0; ka<10; ++ka)
    {
      size_t in_size = in.size();
      for (size_t ki=0; ki<in.size(); ++ki)
	{
	  vector<ftSamplePoint*> next = in[ki]->getNeighbours();
	  for (size_t kj=0; kj<next.size(); ++kj)
	    {
	      RevEngPoint *curr = dynamic_cast<RevEngPoint*>(next[kj]);
	      if (curr->moved())
		continue;  // Already in an extended region
	      if (curr->isEdge())
		continue;  // Do not include points labelled as edge
	      if (curr->region() == this)
		continue;
	      if (curr->visited())
		continue;
	      curr->setVisited();
	      visited.push_back(curr);
	      Vector3D xyz = curr->getPoint();
	      Point pos(xyz[0], xyz[1], xyz[2]);
	      double upar, vpar, dist;
	      Point close;
	      surf->closestPoint(pos, upar, vpar, close, dist, eps);
	      if (dist < approx_tol)
		{
		  Point normal;
		  surf->normal(normal, upar, vpar);
		  double ang = normal.angle(curr->getPCANormal());
		  ang = std::min(ang, M_PI-ang);
		  if (ang < anglim)
		    in.push_back(curr);
		  else
		    out.push_back(curr);
		}
	      else
		out.push_back(curr);
	    }
	}
      std::ofstream of3("in_out.g2");
      of3 << "400 1 0 4 255 0 0 255" << std::endl;
      of3 << in.size() << std::endl;
      for (size_t kj=0; kj<in.size(); ++kj)
	of3 << in[kj]->getPoint() << std::endl;
      of3 << "400 1 0 4 0 255 0 255" << std::endl;
      of3 << out.size() << std::endl;
      for (size_t kj=0; kj<out.size(); ++kj)
	of3 << out[kj]->getPoint() << std::endl;

      if (in.size() == in_size)
	break;

      vector<RevEngPoint*> in_out(in.begin(), in.end());
      in_out.insert(in_out.end(), out.begin(), out.end());
      shared_ptr<SplineSurface> surf2 = surfApprox(in_out, bbox_);
      std::ofstream of4("approx_sf2.g2");
      surf2->writeStandardHeader(of4);
      surf2->write(of4);

      int nmb_in2 = 0;
      double avdist2 = 0.0;
      in.clear();
      out.clear();
      double wgt = 1.0/in_out.size();
      for (size_t ki=0; ki<in_out.size(); ++ki)
	{
	  Vector3D xyz = in_out[ki]->getPoint();
	  Point pos(xyz[0], xyz[1], xyz[2]);
	  double upar, vpar, dist;
	  Point close;
	  surf2->closestPoint(pos, upar, vpar, close, dist, eps);
	  avdist2 += wgt*dist;
	  if (dist < approx_tol)
	    {
	      Point normal;
	      surf2->normal(normal, upar, vpar);
	      double ang = normal.angle(in_out[ki]->getPCANormal());
	      ang = std::min(ang, M_PI-ang);
	      if (ang < anglim)
		in.push_back(in_out[ki]);
	      else
		out.push_back(in_out[ki]);
	    }
	  else
	    out.push_back(in_out[ki]);
	}

      if (in.size() < out.size() || avdist2 > approx_tol)
	{
	  break;
	}

      int stop_break0 = 1;
    }

  // Move in points
  set<RevEngRegion*> affected_reg;
  for (size_t ki=0; ki<in.size(); ++ki)
    {
      in[ki]->setMoved();
      
      if (in[ki]->region() == this)
	continue;

      RevEngRegion *other_reg = in[ki]->region();
      if (other_reg)
	{
	  other_reg->removePoint(in[ki]);
	  affected_reg.insert(other_reg);
	}
      in[ki]->setRegion(this);
      group_points_.push_back(in[ki]);
    }

  for (size_t ki=0; ki<visited.size(); ++ki)
    visited[ki]->unsetVisited();
      
  vector<RevEngRegion*> affected_reg2(affected_reg.begin(), affected_reg.end());
  for (size_t ki=0; ki<affected_reg2.size(); ++ki)
    if (affected_reg2[ki]->numPoints() == 0)
      adapted_regions.push_back(affected_reg2[ki]);
    else
      affected_reg2[ki]->updateInfo();

  updateInfo();
  int stop_break = 1;
  
}

//===========================================================================
shared_ptr<SplineSurface> RevEngRegion::surfApprox(vector<RevEngPoint*>& points,
						   const BoundingBox& bbox)
//===========================================================================
{
  // PCA analysis of entire group to orient plane for parametrerization
  Vector3D xyz = points[0]->getPoint();
  Point pos0(xyz[0], xyz[1], xyz[2]);
  vector<Point> pnts;
  pnts.reserve(points.size());
  double wgt = 1.0/(double)points.size();
  Point vec = wgt*points[0]->minCurvatureVec();
  for (size_t ki=1; ki<points.size(); ++ki)
    {
      xyz = points[ki]->getPoint();
      Point pos(xyz[0], xyz[1], xyz[2]);
      pnts.push_back(pos);
      Point vec0 = points[ki]->minCurvatureVec();
      if (vec0*vec< 0.0)
	vec0 *= -1;
      vec += wgt*vec0;
    }

  double lambda[3];
  double eigenvec[3][3];
  RevEngUtils::principalAnalysis(pos0, pnts, lambda, eigenvec);
  Point nn(eigenvec[2][0],eigenvec[2][1],eigenvec[2][2]);
  Point ydir = nn.cross(vec);
  ydir.normalize_checked();
  Point xdir = ydir.cross(nn);
  xdir.normalize_checked();

  // Parameterize points based on planar surface
  pnts.push_back(pos0);

  vector<double> data;
  vector<double> param;
  RevEngUtils::parameterizeWithPlane(pnts, bbox, xdir, ydir,
				     /*Point(eigenvec[0][0],eigenvec[0][1],eigenvec[0][2]),
				       Point(eigenvec[1][0],eigenvec[1][1],eigenvec[1][2]),*/
				     data, param);

  // Surface approximation
  int order = 3;
  double belt = 0.1*bbox.low().dist(bbox.high());
  shared_ptr<SplineSurface> surf = RevEngUtils::surfApprox(data, 3, param, order,
							   order, order, order, belt);
return surf;
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
  double gausscurv = pt->GaussCurvature();
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
	  if (classification_type_ == CLASSIFICATION_SHAPEINDEX &&
	      (!planar) && meancurv*curr->meanCurvature() < 0.0)
	    continue;  // Not compatible with one primary surface
	  double curv1 = meancurv*gausscurv;
	  double curv2 = curr->meanCurvature()*curr->GaussCurvature();
	  if (classification_type_ == CLASSIFICATION_SHAPEINDEX &&
	      (!planar) && curv1*curv2 < 0.0)
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
  writeRegionInfo(of);
  
  std::ofstream of2("curr_normals.g2");
  writeUnitSphereInfo(of2);

  bool found = false;
  int min_nmb = 20;
  if (group_points_.size() < min_nmb)
    return false;

  // double angtol = 0.1;
  // vector<RevEngPoint*> in, out;
  // analysePlaneProperties(normal, angtol, in, out);
  
  // shared_ptr<Plane> surf1(new Plane(pos, normal1));
  // shared_ptr<Plane> surf(new Plane(pos, normal));

  // Point pos2 = pos;
  // Point normal2;
  // vector<pair<vector<RevEngPoint*>::iterator,
  // 	      vector<RevEngPoint*>::iterator> > group;
  // group.push_back(std::make_pair(group_points_.begin(), group_points_.end()));
  // RevEngUtils::computePlane(group, pos2, normal2);
  // shared_ptr<Plane> surf2(new Plane(pos2, normal2));
  
  // vector<RevEngPoint*> in3, out3;
  // analysePlaneProperties(normal3, angtol, in3, out3);
  
  shared_ptr<Plane> surf3 = computePlane(); //(new Plane(pos3, normal3));

  // Check accuracy
  // double maxdist, avdist;
  // double maxdist1, avdist1;
  // double maxdist2, avdist2;
  double maxdist3, avdist3;
  int num_inside3; //, num_inside1, num_inside2, num_inside3;
  vector<RevEngPoint*> inpt, outpt;
  // RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
  // 			  surf, tol, maxdist, avdist, num_inside);
  // RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
  // 			  surf1, tol, maxdist1, avdist1, num_inside1);
  // RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
  // 			  surf2, tol, maxdist2, avdist2, num_inside2);
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			  surf3, tol, maxdist3, avdist3, num_inside3, inpt, outpt);
  std::ofstream ofd("in_out_plane.g2");
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
  surf3->setParameterBounds(-len, -len, len, len);
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
shared_ptr<Plane> RevEngRegion::computePlane()
//===========================================================================
{
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
  
  impl_ = shared_ptr<ImplicitApprox>(new ImplicitApprox());
  impl_->approx(group_points_, 1);
  Point pos3, normal3;
  impl_->projectPoint(pos, normal, pos3, normal3);
  
  shared_ptr<Plane> surf(new Plane(pos3, normal3));
  Point low = bbox_.low();
  Point high = bbox_.high();
  double len = low.dist(high);
  surf->setParameterBounds(-len, -len, len, len);
  return surf;
}
//===========================================================================
bool RevEngRegion::possiblePlane(double angtol, double inlim)
//===========================================================================
{
  if (planartype())
    return true;

  double wgt = 1.0/(double)group_points_.size();
  Point avnorm(0.0, 0.0, 0.0);
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Point curr = group_points_[ki]->getPCANormal();
      avnorm += wgt*curr;
    }

  int nmb_in = 0;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Point curr = group_points_[ki]->getPCANormal();
      double ang = curr.angle(avnorm);
      if (ang <= angtol)
	++nmb_in;
    }

  double frac = (double)nmb_in/(double)group_points_.size();
  return (frac >= inlim);
}

//===========================================================================
bool RevEngRegion::possibleCylinder(double angtol, double inlim)
//===========================================================================
{
  // if (cylindertype())
  //   return true;

  double wgt = 1.0/(double)group_points_.size();
  Point avvec = wgt*group_points_[0]->minCurvatureVec();
  for (size_t ki=1; ki<group_points_.size(); ++ki)
    {
      Point curr = group_points_[ki]->minCurvatureVec();
      if (curr*avvec < 0.0)
	curr *= -1;
      avvec += wgt*curr;
    }

  int nmb_in = 0;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Point curr = group_points_[ki]->minCurvatureVec();
      double ang = curr.angle(avvec);
      ang = std::min(ang, M_PI-ang);
      if (ang <= angtol)
	++nmb_in;
    }

  double frac = (double)nmb_in/(double)group_points_.size();
  return (frac >= inlim);
}

//===========================================================================
bool RevEngRegion::possibleCone(double tol, double inlim)
//===========================================================================
{
  return true;  // For the time being
}

//===========================================================================
bool RevEngRegion::possibleTorus(double tol, double inlim)
//===========================================================================
{
  // if (planartype() || cylindertype())
  //   return false;
  // Compute curvature variations
  double k1mean = 0.0, k2mean = 0.0;
  double k1_1 = 0.0, k1_2 = 0.0, k2_1 = 0.0, k2_2 = 0.0;
  double d1 = 0.0, d2 = 0.0;
  double wgt = 1.0/(double)group_points_.size();
  Point vec(0.0, 0.0, 0.0);

  double eps = 1.0e-3;
  vector<Vector3D> cneg, cpos, czero;
  double k1min = std::numeric_limits<double>::max();
  double k1max = std::numeric_limits<double>::lowest();
  double k2min = std::numeric_limits<double>::max();
  double k2max = std::numeric_limits<double>::lowest();
  int k1pos = 0, k1neg = 0, k2pos = 0, k2neg = 0;
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
      k1min = std::min(k1min, kmin);
      k1max = std::max(k1max, kmin);
      k2min = std::min(k2min, kmax);
      k2max = std::max(k2max, kmax);
      if (kmin < 0.0)
	k1neg++;
      else
	k1pos++;
      if (kmax < 0.0)
	k2neg++;
      else
	k2pos++;

      Point norm = group_points_[kr]->getPCANormal();
      Point minvec = group_points_[kr]->minCurvatureVec();
      Point temp = norm.cross(minvec);
      temp.normalize_checked();
      vec += wgt*temp;
      if (group_points_[kr]->GaussCurvature()*group_points_[kr]->meanCurvature() < -eps)
	cneg.push_back(group_points_[kr]->getPoint());
      else if (group_points_[kr]->GaussCurvature()*group_points_[kr]->meanCurvature() > eps)
	cpos.push_back(group_points_[kr]->getPoint());
      else
	czero.push_back(group_points_[kr]->getPoint());
    }
  double vark1 = (group_points_.size()*k1_1)/fabs(d1) - fabs(k1_2);
  double vark2 = (group_points_.size()*k2_1)/fabs(d2) - fabs(k2_2);

  double klim1 = 0.1*fabs(k1mean);
  double klim2 = 0.1*fabs(k2mean);
  int nmb1=0, nmb2=0;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      double kmin = group_points_[kr]->minPrincipalCurvature();
      double kmax = group_points_[kr]->maxPrincipalCurvature();
      if (fabs(kmin - k1mean) < klim1)
	++nmb1;
     if (fabs(kmax - k2mean) < klim2)
	++nmb2;
    }

  double frac1 = (double)nmb1/(double)group_points_.size();
  double frac2 = (double)nmb2/(double)group_points_.size();

  std::ofstream of("posnegcurv.g2");
  of << "400 1 0 4 0 100 155 255" << std::endl;
  of << cneg.size() << std::endl;
  for (size_t kr=0; kr<cneg.size(); ++kr)
    of << cneg[kr] << std::endl;
  of << "400 1 0 4 0 255 0 255" << std::endl;
  of << cpos.size() << std::endl;
  for (size_t kr=0; kr<cpos.size(); ++kr)
    of << cpos[kr] << std::endl;
  of << "400 1 0 4 200 55 0 255" << std::endl;
  of << czero.size() << std::endl;
  for (size_t kr=0; kr<czero.size(); ++kr)
    of << czero[kr] << std::endl;
  
  // Something with variation in curvature
  return true;  // For the time being
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
  writeRegionInfo(of);
  
  std::ofstream of2("curr_normals.g2");
  writeUnitSphereInfo(of2);

  bool found = false;
  int min_nmb = 20;
  double eps = 1.0e-6;
  if (group_points_.size() < min_nmb)
    return false;
  // Point vec0 = group_points_[0]->minCurvatureVec();
  // Point avvec(0.0, 0.0, 0.0);
  // double wgt = 1.0/(double)group_points_.size();
  // for (size_t kr=0; kr<group_points_.size(); ++kr)
  //   {
  //     Point vec = group_points_[kr]->minCurvatureVec();
  //     if (vec*vec0 < 0.0)
  // 	vec *= -1;
  //     vec *= wgt;
  //     avvec += vec;
  //   }

  // double angtol= 0.1;
  // vector<RevEngPoint*> in, out;
  // analyseCylinderProperties(avvec, angtol, in, out);
  
  // Point origo(0.0, 0.0, 0.0);
  // of2 << "410 1 0 4 255 0 0 255" << std::endl;
  // of2 << "1" << std::endl;
  // of2 << origo << " " << axis << std::endl;

 

  // vector<RevEngPoint*> in2, out2;
  // analyseCylinderProperties(axis, angtol, in2, out2);
  


  // vector<RevEngPoint*> in3, out3;
  // analyseCylinderProperties(axis2, angtol, in3, out3);
  

  // pluckerAxis();

  // vector<Point> maxcvec(group_points_.size());
  // for (size_t kr=0; kr<group_points_.size(); ++kr)
  //   {
  //     RevEngPoint *pt = dynamic_cast<RevEngPoint*>(group_points_[kr]);
  //     maxcvec[kr] = pt->maxCurvatureVec();
  //   }
  // Point axis3, Cx3, Cy3;
  // RevEngUtils::computeAxis(maxcvec, axis3, Cx3, Cy3);
			
  
  // shared_ptr<Cylinder> cyl(new Cylinder(rad, pnt, axis, Cy));
  // cyl->setParamBoundsV(-0.5*len, 0.5*len);
  //PointCloud3D cloud(&xyzpoints[0], xyzpoints.size()/3);

  // Point pos0 = apex + ((pnt - apex)*axis2)*axis2;
  // double dd = apex.dist(pos0);
  // double rad0 = dd*tan(phi);
  // shared_ptr<Cone> cone(new Cone(rad0, pos0, axis2, Cx2, phi));
  // cone->setParamBoundsV(-0.5*len, 0.5*len);
 // shared_ptr<Cone> cone2(new Cone(rad, pnt, -axis, -Cy, beta));

  shared_ptr<Cylinder> cyl = computeCylinder();
  std::ofstream ofs("cyl.g2");
  cyl->writeStandardHeader(ofs);
  cyl->write(ofs);
  // cone->writeStandardHeader(ofs);
  // cone->write(ofs);
  // ofs <<"400 1 0 4 255 0 0 255" << std::endl;
  // ofs << "1" << std::endl;
  // ofs << pos0 << std::endl;
  
  // Check accuracy
  double maxd, avd; //, maxd2, avd2, maxd3, avd3;
  int num2; //, num3, num4;
  //shared_ptr<ParamSurface> surf = cyl;
  vector<RevEngPoint*> inpt, outpt; //, inpt2, outpt2;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			  cyl, tol, maxd, avd, num2, inpt, outpt);
  // RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
  // 			  cone, tol, maxd2, avd2, num3, inpt2, outpt2);
  //  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
  // 			  cone2, tol, maxd3, avd3, num4);
  std::ofstream ofd("in_out_cyl.g2");
  ofd << "400 1 0 4 155 50 50 255" << std::endl;
  ofd << inpt.size() << std::endl;
  for (size_t kr=0; kr<inpt.size(); ++kr)
    ofd << inpt[kr]->getPoint() << std::endl;
  ofd << "400 1 0 4 50 155 50 255" << std::endl;
  ofd << outpt.size() << std::endl;
  for (size_t kr=0; kr<outpt.size(); ++kr)
    ofd << outpt[kr]->getPoint() << std::endl;


  int num = (int)group_points_.size();
  double maxd_init, avd_init;
  int num_init;
  getAccuracy(maxd_init, avd_init, num_init);
   if (num2 > min_pt && num2 > num/2) //avd <= avlim && maxd <= maxlim)
    {
     bool OK = true;
      if (associated_sf_.size() > 0)
	{
	  // Check with current approximating surface
	  if (num2 < num_init || avd > avd_init)
	    OK = false;
	}

      if (OK)
	{
	  found = true;
	  setAccuracy(maxd, avd, num2);
      
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
      // sfs.push_back(cyl);
      // shared_ptr<SurfaceModel> cylmod(new SurfaceModel(gap, gap, 10.0*gap, 0.01,
      // 						       0.05, sfs));
      // vector<shared_ptr<SurfaceModel> > divcyl = cylmod->splitSurfaceModels(boxmod);
      
	  std::cout << "Cylinder. N1: " << num << ", N2: " << num2 << ", max: " << maxd << ", av: " << avd << std::endl;
	  shared_ptr<HedgeSurface> hedge(new HedgeSurface(cyl, this));
	  for (size_t kh=0; kh<associated_sf_.size(); ++kh)
	    prevsfs.push_back(associated_sf_[kh]);
	  associated_sf_.clear();
	  associated_sf_.push_back(hedge.get());
	  hedgesfs.push_back(hedge);
      // for (int ka=0; ka<divcyl[0]->nmbEntities(); ++ka)
      // 	{
      // 	  shared_ptr<ParamSurface> cyl2 = divcyl[0]->getSurface(ka);
      // 	  shared_ptr<HedgeSurface> hedge(new HedgeSurface(cyl2, this));
      // 	  for (size_t kh=0; kh<associated_sf_.size(); ++kh)
      // 	    prevsfs.push_back(associated_sf_[kh]);
      // 	  associated_sf_.push_back(hedge.get());
      // 	  hedgesfs.push_back(hedge);
	  
	  cyl->writeStandardHeader(fileout);
	  cyl->write(fileout);
	}
    }
  int stop_break0 = 1;
  return found;
}

//===========================================================================
shared_ptr<Cylinder> RevEngRegion::computeCylinder()
//===========================================================================
{
  // Cylinder orientation by covariance matrix of normal vectors
  Point axis;
  Point Cx;
  Point Cy;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(group_points_.begin(), group_points_.end()));
  RevEngUtils::computeAxis(group, axis, Cx, Cy);
  
  std::ofstream of3("rotated_pts_cyl.g2");
  std::ofstream ofp3("projected_pts_cyl.g2");
  Point low = bbox_.low();
  Point high = bbox_.high();
  double len = low.dist(high);
  double rad;
  Point pnt;
// for (int ka=0; ka<2; ++ka)
// {
  RevEngUtils::computeCylPosRadius(group, low, high,
				   axis, Cx, Cy, pnt, rad);

  vector<Point> rotated;
  RevEngUtils::rotateToPlane(group, Cx, axis, pnt, rotated);
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
  ofp3 << "400 1 0 4 255 0 0 255" << std::endl;
  ofp3 << projected.size() << std::endl;
  for (size_t kr=0; kr<projected.size(); ++kr)
    ofp3 << projected[kr] << std::endl;
  circ->writeStandardHeader(ofp3);
  circ->write(ofp3);

  vector<Point> cpos;
  vector<double> crad;
  int nmb_split = 4;
  splitCylinderRad(pnt, axis, Cx, Cy, nmb_split, cpos, crad);
  Point axis2(0.0, 0.0, 0.0);
  for (int kb=1; kb<nmb_split; ++kb)
    {
      Point dir = cpos[kb] - cpos[kb-1];
      dir.normalize();
      axis2 += dir;
    }
  axis2.normalize();
  Point Cx2 = Cy.cross(axis2);
  Cx2.normalize();
  Point Cy2 = axis2.cross(Cx2);
  Cy2.normalize();
 //  if (ka == 0)
 //    {
 //      Cx = Cx2;
 //      Cy = Cy2;
 //      axis = axis2;
 //    }
 // }
  shared_ptr<Cylinder> cyl(new Cylinder(rad, pnt, axis, Cy));
  cyl->setParamBoundsV(-len, len);
  return cyl;
}

//===========================================================================
bool RevEngRegion::extractCone(double tol, int min_pt, double mean_edge_len,
			       vector<shared_ptr<HedgeSurface> >& hedgesfs,
			       vector<HedgeSurface*>& prevsfs,
			       std::ostream& fileout)
//===========================================================================
{
  std::ofstream of("curr_region.g2");
  writeRegionInfo(of);
  
  std::ofstream of2("curr_normals.g2");
  writeUnitSphereInfo(of2);

  bool found = false;
  int min_nmb = 20;
  double eps = 1.0e-6;
  if (group_points_.size() < min_nmb)
    return false;
  
  shared_ptr<Cone> cone = computeCone();
  std::ofstream ofs("cone.g2");
  cone->writeStandardHeader(ofs);
  cone->write(ofs);
  
  // Check accuracy
  double maxd, avd; 
  int num2; 
  vector<RevEngPoint*> inpt, outpt;
  double angtol = 0.1;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
  			  cone, tol, maxd, avd, num2, inpt, outpt, angtol);
  std::ofstream ofd("in_out_cone.g2");
  ofd << "400 1 0 4 155 50 50 255" << std::endl;
  ofd << inpt.size() << std::endl;
  for (size_t kr=0; kr<inpt.size(); ++kr)
    ofd << inpt[kr]->getPoint() << std::endl;
  ofd << "400 1 0 4 50 155 50 255" << std::endl;
  ofd << outpt.size() << std::endl;
  for (size_t kr=0; kr<outpt.size(); ++kr)
    ofd << outpt[kr]->getPoint() << std::endl;


  int num = (int)group_points_.size();
  double maxd_init, avd_init;
  int num_init;
  getAccuracy(maxd_init, avd_init, num_init);
   if (num2 > min_pt && num2 > num/2) //avd <= avlim && maxd <= maxlim)
    {
     bool OK = true;
      if (associated_sf_.size() > 0)
	{
	  // Check with current approximating surface
	  if (num2 < num_init || avd > avd_init)
	    OK = false;
	}

      if (OK)
	{
	  found = true;
	  setAccuracy(maxd, avd, num2);
	  std::cout << "Cone. N1: " << num << ", N2: " << num2 << ", max: " << maxd << ", av: " << avd << std::endl;
	  shared_ptr<HedgeSurface> hedge(new HedgeSurface(cone, this));
	  for (size_t kh=0; kh<associated_sf_.size(); ++kh)
	    prevsfs.push_back(associated_sf_[kh]);
	  associated_sf_.clear();
	  associated_sf_.push_back(hedge.get());
	  hedgesfs.push_back(hedge);
	  
	  cone->writeStandardHeader(fileout);
	  cone->write(fileout);
	}
    }
  int stop_break0 = 1;
  return found;
}


//===========================================================================
shared_ptr<Cone> RevEngRegion::computeCone()
//===========================================================================
{
  Point low = bbox_.low();
  Point high = bbox_.high();
  double len = low.dist(high);

  Point mid(0.0, 0.0, 0.0);
  double wgt = 1.0/(double)group_points_.size();
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      Vector3D xyz = group_points_[kr]->getPoint();
      Point pnt(xyz[0], xyz[1], xyz[2]);
      mid += wgt*pnt;
    }
  Point axis;
  Point Cx;
  Point Cy;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(group_points_.begin(), group_points_.end()));
  RevEngUtils::coneAxis(group, axis, Cx, Cy);

  Point apex;
  double phi;
  RevEngUtils::coneApex(group, axis, apex, phi);
  
  vector<Point> rotated;
  RevEngUtils::rotateToPlane(group, Cx, axis, apex, rotated);
  std::ofstream of("rotated_pts_cone.g2");
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << rotated.size() << std::endl;
  for (size_t kr=0; kr<rotated.size(); ++kr)
    of << rotated[kr] << std::endl;
  of << "410 1 0 4 0 0 255 255" << std::endl;
  of << "1" << std::endl;
  of << apex-0.5*len*axis << " " << apex+0.5*len*axis << std::endl;
  of << "410 1 0 4 0 255 0 255" << std::endl;
  of << "1" << std::endl;
  of << apex-0.5*len*Cx << " " << apex+0.5*len*Cx << std::endl;
  
  Point pos0 = apex + ((mid - apex)*axis)*axis;
  double dd = apex.dist(pos0);
  double rad0 = dd*tan(phi);
  shared_ptr<Cone> cone(new Cone(fabs(rad0), pos0, axis, Cx, phi));
  cone->setParamBoundsV(-0.2*len, 0.2*len);
  
  shared_ptr<Circle> circ(new Circle(rad0, pos0, axis, Cx));
  vector<Point> projected;
  double maxdp, avdp;
  RevEngUtils::projectToPlane(group, axis, pos0, projected, maxdp, avdp);
  std::ofstream ofp3("projected_pts_cone.g2");
  ofp3 << "400 1 0 4 255 0 0 255" << std::endl;
  ofp3 << projected.size() << std::endl;
  for (size_t kr=0; kr<projected.size(); ++kr)
    ofp3 << projected[kr] << std::endl;
  circ->writeStandardHeader(ofp3);
  circ->write(ofp3);

  return cone;
}

//===========================================================================
bool RevEngRegion::extractTorus(double tol, int min_pt, double mean_edge_len,
				vector<shared_ptr<HedgeSurface> >& hedgesfs,
				vector<HedgeSurface*>& prevsfs,
				std::ostream& fileout)
//===========================================================================
{
  bool found = false;
  int min_nmb = 20;
  if (group_points_.size() < min_nmb)
    return false;
  
  std::ofstream of("curr_region.g2");
  writeRegionInfo(of);

  shared_ptr<Torus> tor2;
  shared_ptr<Torus> tor1 = computeTorus(tor2);
  if (!tor1.get())
    return false;
  
  // Check accuracy
  double maxd1, avd1, maxd2, avd2;
  int num1, num2;
  //shared_ptr<ParamSurface> surf = cyl;
  vector<RevEngPoint*> in1, out1,in2, out2;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			  tor1, tol, maxd1, avd1, num1, in1, out1);
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
  			  tor2, tol, maxd2, avd2, num2, in2, out2);

  std::ofstream ofd("in_out_tor.g2");
  ofd << "400 1 0 4 155 50 50 255" << std::endl;
  ofd << in1.size() << std::endl;
  for (size_t kr=0; kr<in1.size(); ++kr)
    ofd << in1[kr]->getPoint() << std::endl;
  ofd << "400 1 0 4 50 155 50 255" << std::endl;
  ofd << out1.size() << std::endl;
  for (size_t kr=0; kr<out1.size(); ++kr)
    ofd << out1[kr]->getPoint() << std::endl;
  ofd << "400 1 0 4 155 50 50 255" << std::endl;
  ofd << in2.size() << std::endl;
  for (size_t kr=0; kr<in2.size(); ++kr)
    ofd << in2[kr]->getPoint() << std::endl;
  ofd << "400 1 0 4 50 155 50 255" << std::endl;
  ofd << out2.size() << std::endl;
  for (size_t kr=0; kr<out2.size(); ++kr)
    ofd << out2[kr]->getPoint() << std::endl;
  
  int num = (int)group_points_.size();
  double maxd_init, avd_init;
  int num_init;
  getAccuracy(maxd_init, avd_init, num_init);
  if (num1 > min_pt && num1 > num/2 && num1 > num2) //avd <= avlim && maxd <= maxlim)
    {
      bool OK = true;
      if (associated_sf_.size() > 0)
	{
	  // Check with current approximating surface
	  if (num1 < num_init || avd1 > avd_init)
	    OK = false;
	}

      if (OK)
	{
	  found = true;
	  setAccuracy(maxd1, avd1, num1);      
      
	  std::cout << "Torus 1. N1: " << num << ", N2: " << num1 << ", max: " << maxd1 << ", av: " << avd1 << std::endl;
	  // associated_sf_.clear();
	  // for (int ka=0; ka<divmod[0]->nmbEntities(); ++ka)
	  // 	{
	  //shared_ptr<ParamSurface> tor1_2 = divmod[0]->getSurface(ka);
	  shared_ptr<HedgeSurface> hedge(new HedgeSurface(tor1, this));
	  for (size_t kh=0; kh<associated_sf_.size(); ++kh)
	    prevsfs.push_back(associated_sf_[kh]);
	  associated_sf_.clear();
	  associated_sf_.push_back(hedge.get());
	  hedgesfs.push_back(hedge);
	  
	  tor1->writeStandardHeader(fileout);
	  tor1->write(fileout);
	}
	// }
	}
  else if (num2 > min_pt && num2 > num/2)
    {
      bool OK = true;
      if (associated_sf_.size() > 0)
	{
	  // Check with current approximating surface
	  if (num2 < num_init || avd2 > avd_init)
	    OK = false;
	}

      if (OK)
	{
	  found = true;
	  setAccuracy(maxd2, avd2, num2);
      
	  std::cout << "Torus 2. N1: " << num << ", N2: " << num2 << ", max: " << maxd2 << ", av: " << avd2 << std::endl;
	  // associated_sf_.clear();
	  // for (int ka=0; ka<divmod[0]->nmbEntities(); ++ka)
	  // 	{
	  // 	  shared_ptr<ParamSurface> tor2_2 = divmod[0]->getSurface(ka);
	  shared_ptr<HedgeSurface> hedge(new HedgeSurface(tor2, this));
	  for (size_t kh=0; kh<associated_sf_.size(); ++kh)
	    prevsfs.push_back(associated_sf_[kh]);
	  associated_sf_.clear();
	  associated_sf_.push_back(hedge.get());
	  hedgesfs.push_back(hedge);
	  
	  tor2->writeStandardHeader(fileout);
	  tor2->write(fileout);
	  // }
	}
    }
  int stop_break = 1;
  return found;
}
 
//===========================================================================
shared_ptr<Torus> RevEngRegion::computeTorus(shared_ptr<Torus>& torus2)
//===========================================================================
{
  std::ofstream of("torus_compute.g2");
  
  // Compute mean curvature and initial point in plane
  double k2mean = 0.0;
  double wgt = 1.0/(double)group_points_.size();
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      double kmax = group_points_[kr]->maxPrincipalCurvature();
      k2mean += wgt*kmax;
    }
  double rd = 1.0/k2mean;
  
  vector<Point> centr(group_points_.size());
  Point mid(0.0, 0.0, 0.0);
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      double kmax = group_points_[kr]->maxPrincipalCurvature();
      k2mean += wgt*kmax;

      Vector3D xyz = group_points_[kr]->getPoint();
      Point xyz2(xyz[0], xyz[1], xyz[2]);
      Point norm = group_points_[kr]->getPCANormal();
      centr[kr] = xyz2 + rd*norm;
      mid += wgt*centr[kr];
    }
  
  of << "400 1 0 4 155 100 0 255" << std::endl;
  of << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      of << centr[kr] << std::endl;
    }
  Point low = bbox_.low();
  Point high = bbox_.high();
  double len = low.dist(high);

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
    {
      shared_ptr<Torus> dummy;
      return dummy;
    }
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << "1" << std::endl;
  of << pos << std::endl;
  of << "410 1 0 4 255 0 0 255" << std::endl;
  of << "1" << std::endl;
  of << pos << " " << pos+0.5*len*normal << std::endl;


  std::ofstream ofimpl("impl_plane.g2");
  impl->visualize(group_points_, ofimpl);
  
  Vector3D xyz = group_points_[0]->getPoint();
  Point xyz2(xyz[0], xyz[1], xyz[2]);
  Point Cx = centr[0] - xyz2;
  Cx -= (Cx*normal)*normal;
  Cx.normalize();
  Point Cy = Cx.cross(normal);
  
  double rad;
  Point pnt;
  RevEngUtils::computeCircPosRadius(centr, normal, Cx, Cy, pnt, rad);
  pnt -= ((pnt - pos)*normal)*normal;
  shared_ptr<Circle> circ(new Circle(rad, pnt, normal, Cx));
  circ->writeStandardHeader(of);
  circ->write(of);

  vector<Point> rotated;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(group_points_.begin(), group_points_.end()));
  RevEngUtils::rotateToPlane(group, Cx, normal, pnt, rotated);
  std::ofstream of3("rotated_pts_tor.g2");
  of3 << "400 1 0 4 255 0 0 255" << std::endl;
  of3 << rotated.size() << std::endl;
  for (size_t kr=0; kr<rotated.size(); ++kr)
    of3 << rotated[kr] << std::endl;
  of3 << "410 1 0 4 0 0 255 255" << std::endl;
  of3 << "1" << std::endl;
  of3 << pnt-0.5*len*normal << " " << pnt+0.5*len*normal << std::endl;
  of3 << "410 1 0 4 0 255 0 255" << std::endl;
  of3 << "1" << std::endl;
  of3 << pnt-0.5*len*Cx << " " << pnt+0.5*len*Cx << std::endl;

  Point cpos;
  double crad;
  RevEngUtils::computeCircPosRadius(rotated, Cy, Cx, normal, cpos, crad);
  shared_ptr<Circle> circ2(new Circle(crad, cpos, Cy, Cx));
  circ2->writeStandardHeader(of3);
  circ2->write(of3);
  shared_ptr<SplineCurve> spl;
  curveApprox(rotated, circ2, spl);
  spl->writeStandardHeader(of3);
  spl->write(of3);
  
  vector<Point> projected;
  double maxdp, avdp;
  RevEngUtils::projectToPlane(group, normal, pnt, projected, maxdp, avdp);
  std::ofstream ofp3("projected_pts_tor.g2");
  ofp3 << "400 1 0 4 255 0 0 255" << std::endl;
  ofp3 << projected.size() << std::endl;
  for (size_t kr=0; kr<projected.size(); ++kr)
    ofp3 << projected[kr] << std::endl;
  circ->writeStandardHeader(ofp3);
  circ->write(ofp3);

  shared_ptr<Torus> tor1(new Torus(rad, fabs(rd), pnt, normal, Cy));
  tor1->writeStandardHeader(of);
  tor1->write(of);
  Point cvec = cpos - pnt;
  double R1 = (cvec - (cvec*normal)*normal).length();
  double R2 = (cvec*normal)*normal.length();
  shared_ptr<Torus> tor2(new Torus(R1, crad, pnt+R2*normal, normal, Cy));
  tor2->writeStandardHeader(of);
  tor2->write(of);
  std::cout << "Torus small radius: " << fabs(rd) << ", " << crad << std::endl;

  torus2 = tor2;
  return tor1;
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
  //extendWithGaussRad2();

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
void RevEngRegion::updateInfo()
//===========================================================================
{
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
void RevEngRegion::addPoint(RevEngPoint* point)
//===========================================================================
{
group_points_.push_back(point);
 Vector3D point2 = point->getPoint();
 Point point3(point2[0], point2[1], point2[2]);
 bbox_.addUnionWith(point3);
 normalcone_.addUnionWith(point->getPCANormal());
 double k1 = point->minPrincipalCurvature();
 double k2 = point->maxPrincipalCurvature();
 mink1_ = std::min(mink1_, fabs(k1));
 maxk1_ = std::max(maxk1_, fabs(k1));
 mink2_ = std::min(mink2_, fabs(k2));
 maxk2_ = std::max(maxk2_, fabs(k2));
}

//===========================================================================
void RevEngRegion::removePoint(RevEngPoint* point)
//===========================================================================
{
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    if (group_points_[ki] == point)
      {
	group_points_.erase(group_points_.begin()+ki);
	break;
      }
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
void RevEngRegion::setRegionAdjacency()
//===========================================================================
{
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      vector<ftSamplePoint*> next = group_points_[ki]->getNeighbours();
      for (size_t kj=0; kj<next.size(); ++kj)
	{
	  RevEngPoint *curr = dynamic_cast<RevEngPoint*>(next[kj]);
	  RevEngRegion *adj_reg = curr->region();
	  if (adj_reg == this)
	    continue;
	  if (adj_reg)
	    {
	      adj_reg->addAdjacentRegion(this);
	      addAdjacentRegion(adj_reg);
	    }
	}
    }
}

//===========================================================================
bool RevEngRegion::integrateInAdjacent(double mean_edge_len, int min_next,
				       int max_next, double tol, double angtol,
				       int max_nmb_outlier)
//===========================================================================
{
if (adjacent_regions_.size() != 1)
  return false;   // To be removed

if (group_points_.size() > max_next/2)
  return false;

// Fetch nearby points
 double local_len = group_points_[0]->getMeanEdgLen(10.0*mean_edge_len);
 double radius = 3.0*(local_len + mean_edge_len);
 radius = std::min(radius, 20.0*mean_edge_len);
 vector<RevEngPoint*> nearpts;
 group_points_[0]->fetchClosePoints2(radius, min_next,
				     max_next-(int)group_points_.size(),
				     nearpts, (*adjacent_regions_.begin()));
 size_t nearnmb = nearpts.size();
 if ((int)nearnmb <= max_nmb_outlier)
   {
     if ((int)group_points_.size() <= max_nmb_outlier)
       {
	 for (size_t ki=0; ki<group_points_.size(); ++ki)
	   {
	     group_points_[ki]->unsetRegion();
	     group_points_[ki]->setOutlier();
	   }
	 return true;
       }
     else
       return false;
   }

 nearpts.insert(nearpts.end(), group_points_.begin(), group_points_.end());
 BoundingBox bbox = bbox_;
 for (size_t ki=0; ki<nearnmb; ++ki)
   {
     Vector3D xyz = nearpts[ki]->getPoint();
     Point xyz2(xyz[0], xyz[1], xyz[2]);
     bbox.addUnionWith(xyz2);
   }
 shared_ptr<SplineSurface> surf = surfApprox(nearpts, bbox);

 // Check accuracy
 double maxd1, maxd2, avd1, avd2;
 int nmb_in1, nmb_in2;
 vector<RevEngPoint*> in1, in2, out1, out2;
 RevEngUtils::distToSurf(nearpts.begin(), nearpts.begin()+nearnmb, surf,
			 tol, maxd1, avd1, nmb_in1, in1, out1, angtol);
 RevEngUtils::distToSurf(nearpts.begin()+nearnmb, nearpts.end(), surf,
			 tol, maxd2, avd2, nmb_in2, in2, out2, angtol);

 if (avd2 < tol && //avd2 <= avd1 &&
     (double)nmb_in2/(double)group_points_.size() <= (double)nmb_in1/(double)nearnmb)
   {
     (*adjacent_regions_.begin())->addRegion(this);
     for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
       (*it)->removeAdjacentRegion(this);
     return true;
   }
 
 return false;
}
	   
//===========================================================================
void RevEngRegion::addRegion(RevEngRegion* reg, double maxd, double avd,
			     int num_inside)
//===========================================================================
{
  int num = reg->numPoints();
  bbox_.addUnionWith(reg->boundingBox());
  normalcone_.addUnionWith(reg->getNormalCone());
  if (num_inside >= 0)
    {
      maxdist_ = std::max(maxdist_, maxd);
      double div = (double)((int)group_points_.size() + num);
      avdist_ = (group_points_.size()*avdist_ + num*avd)/div;
      num_inside_ += num_inside;
    }

  double mink1, maxk1, mink2, maxk2;
  reg->getPrincipalCurvatureInfo(mink1, maxk1, mink2, maxk2);
  mink1_ = std::min(mink1_, mink1);
  maxk1_ = std::max(maxk1_, maxk1);
  mink2_ = std::min(mink2_, mink2);
  maxk2_ = std::max(maxk2_, maxk2);
  
  for (auto it=reg->pointsBegin(); it != reg->pointsEnd(); ++it)
    (*it)->setRegion(this);
  group_points_.insert(group_points_.end(), reg->pointsBegin(),
		       reg->pointsEnd());
  
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
		  addRegion(adj_reg, maxd, avd, num_inside);
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
		  maxdist_ = std::max(maxdist_, dist);
		  avdist_ = (group_points_.size()*avdist_ + dist)/(double)group_points_.size();
		  num_inside_ += 1;
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

  // Check if the surface should be updated
  shared_ptr<ParamSurface> updated, updated2;
  if (classtype == Class_Plane)
    {
      updated = computePlane();
    }
  else if (classtype == Class_Cylinder)
    {
      updated = computeCylinder();
    }
  else if (classtype == Class_Torus)
    {
      shared_ptr<Torus> torus2;
      updated = computeTorus(torus2);
      updated2 = torus2;
    }

  if (updated.get())
    {
      std::ofstream ofn("updated.g2");
      updated->writeStandardHeader(ofn);
      updated->write(ofn);
      if (updated2.get())
	{
	  updated2->writeStandardHeader(ofn);
	  updated2->write(ofn);
	}
    }
	

  shared_ptr<ParamSurface> replacesurf;
  if (updated.get())
    {
        double maxd, avd;
	int num_inside;
	vector<RevEngPoint*> inpt, outpt;
	RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
				updated, tol, maxd, avd, num_inside, inpt, outpt);

	if (updated2.get())
	  {
	    double maxd2, avd2;
	    int num_inside2;
	    vector<RevEngPoint*> inpt2, outpt2;
	    RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
				    updated2, tol, maxd2, avd2, num_inside2,
				    inpt2, outpt2);
	    if ((num_inside2 < num_inside || (num_inside2 == num_inside && avd2 < avd))
		&& (num_inside2 < num_inside_ || (num_inside2 == num_inside_ && avd2 < avdist_))
		&& avd2 < tol)
	      {
		replacesurf = updated2;
		setAccuracy(maxd2, avd2, num_inside2);
	      }
	    int stop_break2 = 1;
	  }
	if ((!replacesurf.get()) &&
	    (num_inside < num_inside_ || (num_inside == num_inside_ && avd < avdist_)) &&
	    avd < tol)
	  {
	    replacesurf = updated;
	    setAccuracy(maxd, avd, num_inside);
	  }
	int stop_break = 1;
    }

  if (replacesurf.get())
    {
      associated_sf_[0]->replaceSurf(replacesurf);
    }
}


//===========================================================================
bool RevEngRegion::isCompatible(ClassType classtype, int sfcode)
//===========================================================================
{
  return true;
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


//===========================================================================
void  RevEngRegion::splitCylinderRad(const Point& pos, const Point& axis,
				     const Point& Cx, const Point& Cy,
				     int nmb_split, vector<Point>& centr,
				     vector<double>& rad)
//===========================================================================
{
  centr.resize(nmb_split);
  rad.resize(nmb_split);
  shared_ptr<Line> line(new Line(pos, axis));
  vector<double> par(group_points_.size());
  double tmin = std::numeric_limits<double>::max();
  double tmax = std::numeric_limits<double>::lowest();
  double diag = bbox_.low().dist(bbox_.high());
  double tmin0 = -diag;
  double tmax0 = diag;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Vector3D xyz = group_points_[ki]->getPoint();
      Point curr(xyz[0], xyz[1], xyz[2]);
      double tpar, dist;
      Point close;
      line->closestPoint(curr, tmin0, tmax0, tpar, close, dist);
      par[ki] = tpar;
      tmin = std::min(tmin, tpar);
      tmax = std::max(tmax, tpar);
    }

  double tdel = (tmax - tmin)/(double)nmb_split;
  vector<Point> mid(nmb_split);
    double tpar = tmin + 0.5*tdel;
  for (int ka=0; ka<nmb_split; ++ka, tpar+=tdel)
    mid[ka] = line->ParamCurve::point(tpar);

  vector<vector<Point> > proj_pts(nmb_split);
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Vector3D pnt = group_points_[ki]->getPoint();
      size_t kj;
      for (kj=0, tpar=tmin; kj<mid.size(); ++kj, tpar+=tdel)
	if (par[ki] >= tpar && par[ki] < tpar+tdel)
	  break;
      kj = std::min(kj, mid.size()-1);
      Point curr(pnt[0], pnt[1], pnt[2]);
      Point curr2 = curr - mid[kj];
      curr2 -= ((curr2*axis)*axis);
      curr2 += mid[kj];
      proj_pts[kj].push_back(curr2);
    }

  std::ofstream of("split_circ.g2");
  for (size_t kj=0; kj<proj_pts.size(); ++kj)
    {
      RevEngUtils::computeCircPosRadius(proj_pts[kj], axis, Cx, Cy, centr[kj], rad[kj]);
      shared_ptr<Circle> circ(new Circle(rad[kj], centr[kj], axis, Cx));
      circ->writeStandardHeader(of);
      circ->write(of);
      of << "400 1 0 4 255 0 0 255" << std::endl;
      of << proj_pts[kj].size() << std::endl;
      for (size_t kr=0; kr<proj_pts[kj].size(); ++kr)
	of << proj_pts[kj][kr] << std::endl;
    }


  
  int stop_break = 1;
 }


//===========================================================================
void  RevEngRegion::curveApprox(vector<Point>& points,
				shared_ptr<Circle> circle, shared_ptr<SplineCurve>& curve)
//===========================================================================
{
  vector<double> pts;
  vector<double> param;
  double tmin = circle->startparam();
  double tmax = circle->endparam();
  double tmin2 = tmax;
  double tmax2 = tmin;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      double tpar, dist;
      Point close;
      circle->closestPoint(points[ki], tmin, tmax, tpar, close, dist);
      pts.insert(pts.end(), points[ki].begin(), points[ki].end());
      param.push_back(tpar);
      tmin2 = std::min(tmin2, tpar);
      tmax2 = std::max(tmax2, tpar);
    }

  int ik = 4;
  int in = 8;
  double tdel = (tmax2 - tmin2)/(double)(in - ik + 1);
  double et[12];
  for (int ka=0; ka<ik; ++ka)
    {
      et[ka] = tmin2;
      et[in+ka] = tmax2;
    }
  for (int ka=ik; ka<in; ++ka)
    et[ka] = tmin2 + (ka-ik+1)*tdel;

  vector<double> ecoef(3*in, 0.0);
  shared_ptr<SplineCurve> cv(new SplineCurve(in, ik, et, &ecoef[0], 3));

  SmoothCurve smooth(3);
  vector<int> cfn(in, 0);
  vector<double> wgts(param.size(), 1.0);
  smooth.attach(cv, &cfn[0]);

  smooth.setOptim(0.0, 0.001, 0.001);
  smooth.setLeastSquares(pts, param, wgts, 0.998);

  smooth.equationSolve(curve);
  int stop_break = 1;
}

void RevEngRegion::writeRegionInfo(std::ostream& of)
{
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
}


void RevEngRegion::writeUnitSphereInfo(std::ostream& of)
{
  of << "400 1 0 4 100  0 155 255" << std::endl;
  of << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(group_points_[kr]);
      Point norm = pt->getPCANormal();
      of << norm << std::endl;
    }
  Sphere sph(1.0, Point(0.0, 0.0, 0.0), Point(0.0, 0.0, 1.0),
	     Point(1.0, 0.0, 0.0));
  sph.writeStandardHeader(of);
  sph.write(of);

  of << "400 1 0 4 0 100  155 255" << std::endl;
  of << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(group_points_[kr]);
      Point vec = pt->minCurvatureVec();
      of << vec << std::endl;
    }
  of << std::endl;
}
