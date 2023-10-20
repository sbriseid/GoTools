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
#include "GoTools/compositemodel/RevEng.h"
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
#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/SweepSurfaceCreator.h"
#include "GoTools/geometry/GoIntersections.h"
#include "GoTools/creators/SmoothCurve.h"
#include "GoTools/creators/ApproxCurve.h"
#include "GoTools/creators/ApproxSurf.h"
// #include "GoTools/creators/SmoothSurf.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "newmat.h"
#include "newmatap.h"
#include <vector>
#include <fstream>

//#define DEBUG_JOIN
//#define DEBUG_UPDATE
//#define DEBUG_INTEGRATE
#define DEBUG_TORUSCONTEXT

using namespace Go;
using std::vector;
using std::set;

//===========================================================================
RevEngRegion::RevEngRegion(int edge_class_type)
//===========================================================================
  : classification_type_(CLASSIFICATION_UNDEF), edge_class_type_(edge_class_type),
    associated_sf_(0), mink1_(0.0), maxk1_(0.0), 
    mink2_(0.0), maxk2_(0.0), avH_(0.0), avK_(0.0), MAH_(0.0), MAK_(0.0),
    frac_norm_in_(0.0), maxdist_(0.0), avdist_(0.0), variance_(0.0), num_inside_(0), 
    prev_region_(0), maxdist_base_(0.0), avdist_base_(0.0), num_in_base_(0),
    maxdist_primary_(0.0), avdist_primary_(0.0), num_in_primary_(0),
    visited_(false)
{
  domain_[0] = domain_[1] = domain_[2] = domain_[3] = 0.0;
  bbox_ = BoundingBox(3);
}

//===========================================================================
RevEngRegion::RevEngRegion(int classification_type, int edge_class_type)
//===========================================================================
  : classification_type_(classification_type), edge_class_type_(edge_class_type),
    associated_sf_(0), mink1_(0.0), maxk1_(0.0), 
    mink2_(0.0), maxk2_(0.0), avH_(0.0), avK_(0.0), MAH_(0.0), MAK_(0.0),
    frac_norm_in_(0.0), maxdist_(0.0), avdist_(0.0), variance_(0.0), num_inside_(0), 
    prev_region_(0), maxdist_base_(0.0), avdist_base_(0.0), num_in_base_(0),
    maxdist_primary_(0.0), avdist_primary_(0.0), num_in_primary_(0),
    visited_(false)
{
  bbox_ = BoundingBox(3);
}

//===========================================================================
RevEngRegion::RevEngRegion(int classification_type,
			   int edge_class_type,
			   vector<RevEngPoint*> & points)
//===========================================================================
  : group_points_(points), classification_type_(classification_type),
    edge_class_type_(edge_class_type), associated_sf_(0),
    avH_(0.0), avK_(0.0), MAH_(0.0), MAK_(0.0),
    frac_norm_in_(0.0), maxdist_(0.0), avdist_(0.0), variance_(0.0), num_inside_(0), 
    prev_region_(0), maxdist_base_(0.0), avdist_base_(0.0), num_in_base_(0),
    maxdist_primary_(0.0), avdist_primary_(0.0), num_in_primary_(0),
    visited_(false)
{
  domain_[0] = domain_[1] = domain_[2] = domain_[3] = 0.0;

  for (size_t kj=0; kj<group_points_.size(); ++kj)
    group_points_[kj]->setRegion(this);
  
  // Bounding box and principal curvature summary
  maxk2_ = std::numeric_limits<double>::lowest();
  mink2_ = std::numeric_limits<double>::max();
  maxk1_ = std::numeric_limits<double>::lowest();
  mink1_ = std::numeric_limits<double>::max();
  bbox_ = BoundingBox(3);
  double fac = 1.0/(double)group_points_.size();
  for  (size_t kj=0; kj<group_points_.size(); ++kj)
    {
      double k1 = group_points_[kj]->minPrincipalCurvature();
      double k2 = group_points_[kj]->maxPrincipalCurvature();
      double H = group_points_[kj]->meanCurvature();
      double K = group_points_[kj]->GaussCurvature();
      mink1_ = std::min(mink1_, fabs(k1));
      maxk1_ = std::max(maxk1_, fabs(k1));
      mink2_ = std::min(mink2_, fabs(k2));
      maxk2_ = std::max(maxk2_, fabs(k2));
      avH_ += fac*H;
      avK_ += fac*K;
      MAH_ += fac*fabs(H);
      MAK_ += fac*fabs(K);
      Vector3D point = group_points_[kj]->getPoint();
      Point point2(point[0], point[1], point[2]);
      bbox_.addUnionWith(point2);
    }
  
  normalcone_ = DirectionCone(group_points_[0]->getMongeNormal());
  avnorm_ = Point(0.0, 0.0, 0.0);
  for  (size_t kj=1; kj<group_points_.size(); ++kj)
    {
      Point norm = group_points_[kj]->getMongeNormal();
      normalcone_.addUnionWith(norm);
      avnorm_ += fac*norm;
    }
}

//===========================================================================
RevEngRegion::~RevEngRegion()
//===========================================================================
{
  int stop_break = 1;
}

//===========================================================================
void RevEngRegion::setAccuracy(double maxdist, double avdist, int num_inside)
//===========================================================================
{
  maxdist_ = maxdist;
  avdist_ = avdist;
  num_inside_ = num_inside;

  variance_ = 0.0;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      double tmp = group_points_[ki]->getSurfaceDist()-avdist_;
      tmp = tmp*tmp;
      tmp /= (double)(group_points_.size()-1);
      variance_ += tmp;
    }
  
}

//===========================================================================
void RevEngRegion::joinRegions(double approx_tol, double anglim,
				vector<RevEngRegion*>& adapted_regions)
//===========================================================================
{
  if (group_points_.size() < 5)
    return;

#ifdef DEBUG_JOIN
  std::ofstream of("region_grow.g2");
  writeRegionInfo(of);
#endif

  //double eps = 1.0e-6;
  shared_ptr<ParamSurface> surf1;
  if (basesf_.get())
    surf1 = basesf_;
  else
    surf1 = surfApprox(group_points_, bbox_);
#ifdef DEBUG_JOIN
  std::ofstream of2("approx_sf.g2");
  surf1->writeStandardHeader(of2);
  surf1->write(of2);
#endif
  
  double avdist1, maxdist1;
  vector<RevEngPoint*> in1, out1;
  approximationAccuracy(group_points_, surf1, approx_tol, anglim, maxdist1, avdist1,
			in1, out1);
  if (in1.size() < group_points_.size()/2 || avdist1 > approx_tol)
    return;

  int kv=0;
  int numfac = 10;
  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end();)
    {
      ++kv;
      if ((*it)->hasSurface())
	{
	  ++it;
	  continue;  // Not the correct grow method
	}
      
      if ((*it)->visited())
	{
	  ++it;
	  continue;
	}
      (*it)->setVisited(true);

      int num = (int)group_points_.size();
      int num2 = (*it)->numPoints();
      if (num2 > numfac*num)
	{
	  ++it;
	  continue;
	}

      auto itnext = it;
      ++itnext;

#ifdef DEBUG_JOIN
      std::ofstream ofn("region_grow_cand.g2");
      (*it)->writeRegionInfo(ofn);
#endif
      
     vector<RevEngPoint*> points = (*it)->getPoints();

      // Check if the points can be approximated by the current surface
      double avdist2_0, maxdist2_0;
      vector<RevEngPoint*> in2_0, out2_0;
      approximationAccuracy(points, surf1, approx_tol, anglim, maxdist2_0, avdist2_0,
			    in2_0, out2_0);

      // Compute overall numbers
      int tot_num = (int)group_points_.size()+points.size();
      double frac = (double)(in1.size()+in2_0.size())/(double)(tot_num);
      double avd = (avdist1*(double)group_points_.size() +
		    avdist2_0*(double)points.size())/(double)tot_num;
      //if (in2_0.size() > points.size()/2 && avdist2_0 < approx_tol)
      if (frac > 0.5 && avd < approx_tol)
	{
	  // Include region in current
	  adapted_regions.push_back((*it));
	  for (auto it2=(*it)->adjacent_regions_.begin();
	       it2!=(*it)->adjacent_regions_.end(); ++it2)
	    {
	      if (*it2 != this)
		{
		  addAdjacentRegion(*it2);
		  (*it2)->addAdjacentRegion(this);
		  (*it2)->removeAdjacentRegion(*it);
		}
	    }
	  maxdist1 = std::max(maxdist1, maxdist2_0);
	  avdist1 = (num*avdist1 + num2*avdist2_0)/(double)(num+num2);
	  for (auto it3=(*it)->pointsBegin(); it3!=(*it)->pointsEnd(); ++it3)
	    (*it3)->addMove();
	  vector<pair<double,double> > dummy;
	  addRegion((*it), dummy);
	  removeAdjacentRegion(*it);
	  it = itnext;
	  continue;
	}
			     

      bool always_stop = true;
      if (num2 < 5 || num2 > num || always_stop)
	{
	  it = itnext;
	  continue;
	}
      
      shared_ptr<ParamSurface> surf2;

      if (basesf_.get() &&
	  ((!(*it)->basesf_.get()) ||
	   ((*it)->basesf_.get() &&
	    basesf_->instanceType() == (*it)->basesf_->instanceType())))
	{
	  vector<RevEngPoint*> all_pts;
	  all_pts.insert(all_pts.end(), group_points_.begin(), group_points_.end());
	  all_pts.insert(all_pts.end(), points.begin(), points.end());
	  if (basesf_->instanceType() == Class_Plane)
	    surf2 = computePlane(all_pts, avnorm_);
	  else if (basesf_->instanceType() == Class_Cylinder)
	    {
	      vector<vector<RevEngPoint*> > configs;
	      surf2 = computeCylinder(all_pts, approx_tol, configs);
	    }
	}
      else
	surf2 = surfApprox(points, (*it)->getBbox());

      if (!surf2.get())
	{
	  it = itnext;
	  continue;
	}
	
#ifdef DEBUG_JOIN
      std::ofstream of3("approx_sf2.g2");
      surf2->writeStandardHeader(of3);
      surf2->write(of3);
#endif
      
      double avdist2, maxdist2;
      vector<RevEngPoint*> in2, out2;
      approximationAccuracy(points, surf2, approx_tol, anglim,
			    maxdist2, avdist2, in2, out2);
      if (in2.size() < points.size()/2 && maxdist2 > approx_tol)
	{
	  it = itnext;
	  continue;
	}

      vector<RevEngPoint*> mergept(group_points_.begin(), group_points_.end());
      mergept.insert(mergept.end(), (*it)->pointsBegin(),
		     (*it)->pointsEnd());
      BoundingBox bbox = bbox_;
      bbox.addUnionWith((*it)->getBbox());

      shared_ptr<SplineSurface> surf3 = surfApprox(mergept, bbox);
#ifdef DEBUG_JOIN
      std::ofstream of4("approx_sf3.g2");
      surf3->writeStandardHeader(of4);
      surf3->write(of4);
#endif
      
      double avdist3_0, maxdist3_0;
      vector<RevEngPoint*> in3_0, out3_0;
      approximationAccuracy(points, surf3, approx_tol, anglim,
			    maxdist3_0, avdist3_0, in3_0, out3_0);
      if (in3_0.size() < points.size()/2 || avdist3_0 > approx_tol)
	{
	  it = itnext;
	  continue;
	}

      
      double avdist3, maxdist3;
      vector<RevEngPoint*> in3, out3;
      approximationAccuracy(mergept, surf3, approx_tol, anglim,
			    maxdist3, avdist3, in3, out3);
      
      if (/*maxdist3 <= 1.1*std::max(maxdist1, maxdist2) && */
	  avdist3 <= 1.1*std::max(avdist1, avdist2) && avdist3 < approx_tol &&
	  in3.size() > in1.size() && in3.size() > 3*(in1.size()+in2.size())/4)
	{
	  // Include region in current
	  adapted_regions.push_back((*it));
	  for (auto it2=(*it)->adjacent_regions_.begin();
	       it2!=(*it)->adjacent_regions_.end(); ++it2)
	    {
	      if (*it2 != this)
		{
		  addAdjacentRegion(*it2);
		  (*it2)->addAdjacentRegion(this);
		  (*it2)->removeAdjacentRegion(*it);
		}
	    }
	  maxdist1 = maxdist3;
	  avdist1 = avdist3;
	  for (auto it3=(*it)->pointsBegin(); it3!=(*it)->pointsEnd(); ++it3)
	    (*it3)->addMove();
	  vector<pair<double,double> > dummy;
	  addRegion((*it), dummy);
	  removeAdjacentRegion(*it);
	  in1 = in3;
	  out1 = out3;
	  surf1 = surf3;
	}

#ifdef DEBUG_JOIN
      std::ofstream ofl("region_grow2.g2");
      writeRegionInfo(ofl);
#endif
      int stop_break = 1;
      it = itnext;
    }

  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
    (*it)->setVisited(false);
      
}

  
//===========================================================================
void RevEngRegion::approximationAccuracy(vector<RevEngPoint*>& points,
					 shared_ptr<ParamSurface> surf,
					 double tol, double angtol,
					 double& maxd, double& avd,
					 vector<RevEngPoint*>& in,
					 vector<RevEngPoint*>& out)
//===========================================================================
{
  avd = maxd = 0.0;
  double eps = 1.0e-6;
  double wgt = 1.0/(double)points.size();
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      Vector3D xyz = points[ki]->getPoint();
      Point pos(xyz[0], xyz[1], xyz[2]);
      double upar, vpar, dist;
      Point close;
      surf->closestPoint(pos, upar, vpar, close, dist, eps);
      avd += wgt*dist;
			       maxd = std::max(maxd, dist);
      if (dist < tol)
	{
	  Point normal;
	  surf->normal(normal, upar, vpar);
	  double ang = normal.angle(points[ki]->getMongeNormal());
	  ang = std::min(ang, M_PI-ang);
	  if (ang < angtol)
	    in.push_back(points[ki]);
	  else
	    out.push_back(points[ki]);
	}
      else
	out.push_back(points[ki]);
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

#ifdef DEBUG_UPDATE
  std::ofstream of("region_grow.g2");
  writeRegionInfo(of);
#endif
  double eps = 1.0e-6;
  shared_ptr<SplineSurface> surf = surfApprox(group_points_, bbox_);
#ifdef DEBUG_UPDATE
  std::ofstream of2("approx_sf.g2");
  surf->writeStandardHeader(of2);
  surf->write(of2);
#endif  
  // Check if the accuracy is sufficient as a basis for growth
  double avdist = 0.0, maxdist = 0.0;
  vector<RevEngPoint*> in, out;
  approximationAccuracy(group_points_, surf, approx_tol, anglim, maxdist, avdist,
			in, out);
  if (in.size() < group_points_.size()/2 || avdist > approx_tol)
    return;

  setBaseSf(surf, maxdist, avdist, (int)in.size());
  
  // Grow
  vector<RevEngPoint*> visited;
  int numpt = (int)group_points_.size();
  int numfac = 10;
  // std::set<RevEngPoint*> tmpset0(in.begin(), in.end());
  // if (tmpset0.size() != in.size())
  //   std::cout << "Point number mismatch, init. " << tmpset0.size() << " " << in.size() << std::endl;
  size_t ki = 0;
  for (int ka=0; ka<10; ++ka)
    {
      size_t in_size = in.size();
      for (; ki<in.size(); ++ki)
	{
	  vector<ftSamplePoint*> next = in[ki]->getNeighbours();
	  for (size_t kj=0; kj<next.size(); ++kj)
	    {
	      RevEngPoint *curr = dynamic_cast<RevEngPoint*>(next[kj]);
	      if (curr->moved())
		continue;  // Already in an extended region
	      if (curr->isEdge(edge_class_type_))
		continue;  // Do not include points labelled as edge
	      if (curr->region() == this)
		continue;
	      if (curr->visited())
		continue;

	      // if (std::find(in.begin(), in.end(), curr) != in.end())
	      // 	{
	      // 	  std::cout << "Double point" << curr << std::endl;
	      // 	}
	      curr->setVisited();
	      visited.push_back(curr);
	      int numpt2 = 0;
	      if (curr->region())
		numpt2 = curr->region()->numPoints();
	      if (numpt2 > numfac*numpt)
		continue;
	      Vector3D xyz = curr->getPoint();
	      Point pos(xyz[0], xyz[1], xyz[2]);
	      double upar, vpar, dist;
	      Point close;
	      surf->closestPoint(pos, upar, vpar, close, dist, eps);
	      if (dist < approx_tol)
		{
		  Point normal;
		  surf->normal(normal, upar, vpar);
		  double ang = normal.angle(curr->getMongeNormal());
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
#ifdef DEBUG_UPDATE
      std::ofstream of3("in_out.g2");
      of3 << "400 1 0 4 255 0 0 255" << std::endl;
      of3 << in.size() << std::endl;
      for (size_t kj=0; kj<in.size(); ++kj)
	of3 << in[kj]->getPoint() << std::endl;
      of3 << "400 1 0 4 0 255 0 255" << std::endl;
      of3 << out.size() << std::endl;
      for (size_t kj=0; kj<out.size(); ++kj)
	of3 << out[kj]->getPoint() << std::endl;
#endif
      if (in.size() == in_size)
	break;

      vector<RevEngPoint*> in_out(in.begin(), in.end());
      in_out.insert(in_out.end(), out.begin(), out.end());
      shared_ptr<SplineSurface> surf2 = surfApprox(in_out, bbox_);
#ifdef DEBUG_UPDATE
      std::ofstream of4("approx_sf2.g2");
      surf2->writeStandardHeader(of4);
      surf2->write(of4);
#endif
      //int nmb_in2 = 0;
      double avdist2 = 0.0;
      double maxdist2 = 0.0;
      in.clear();
      out.clear();
      double wgt = 1.0/(double)in_out.size();
      for (size_t ki=0; ki<in_out.size(); ++ki)
	{
	  Vector3D xyz = in_out[ki]->getPoint();
	  Point pos(xyz[0], xyz[1], xyz[2]);
	  double upar, vpar, dist;
	  Point close;
	  surf2->closestPoint(pos, upar, vpar, close, dist, eps);
	  avdist2 += wgt*dist;
	  maxdist2 = std::max(maxdist2,dist);
	  if (dist < approx_tol)
	    {
	      Point normal;
	      surf2->normal(normal, upar, vpar);
	      double ang = normal.angle(in_out[ki]->getMongeNormal());
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

      setBaseSf(surf2, maxdist2, avdist2, (int)in.size());
      surf = surf2;
      // std::set<RevEngPoint*> tmpset(in.begin(), in.end());
      // if (tmpset.size() != in.size())
      // 	std::cout << "Point number mismatch. " << ka << " " << tmpset.size() << " " << in.size() << std::endl;
      int stop_break0 = 1;
    }

  //std::cout << "Ready to move points " << std::endl;
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
      in[ki]->addMove();
      group_points_.push_back(in[ki]);
    }

  for (size_t ki=0; ki<visited.size(); ++ki)
    visited[ki]->unsetVisited();
      
  vector<RevEngRegion*> affected_reg2(affected_reg.begin(), affected_reg.end());
  for (size_t ki=0; ki<affected_reg2.size(); ++ki)
    {
      if (affected_reg2[ki]->numPoints() == 0)
	{
	  for (auto it=affected_reg2[ki]->adjacent_regions_.begin();
	       it!=affected_reg2[ki]->adjacent_regions_.end(); ++it)
	    {
	      // if ((*it) != this)
	      //   {
	      // 	addAdjacentRegion(*it);
	      (*it)->removeAdjacentRegion(affected_reg2[ki]);
	      // }
	    }
	  // removeAdjacentRegion(affected_reg2[ki]);
	  adapted_regions.push_back(affected_reg2[ki]);
	}
      else
	{
	  vector<vector<RevEngPoint*> > separate_groups;
	  affected_reg2[ki]->splitRegion(separate_groups);
	  int classtype = affected_reg2[ki]->getClassificationType();
	  for (size_t ki=0; ki<separate_groups.size(); ++ki)
	    {
	      shared_ptr<RevEngRegion> reg(new RevEngRegion(classtype,
							    edge_class_type_,
							    separate_groups[ki]));
	      outdiv_regions.push_back(reg);
	    }
	  affected_reg2[ki]->updateInfo();
	}
    }

  updateInfo();
  int stop_break = 1;
  
}

//===========================================================================
bool RevEngRegion::segmentByPlanarContext(int min_point_in, double tol,
					  vector<RevEngRegion*>& adj_planar,
					  vector<vector<RevEngPoint*> >& added_groups)
//===========================================================================
{
  std::ofstream of("adj_planar.g2");
  for (size_t ki=0; ki<adj_planar.size(); ++ki)
    adj_planar[ki]->writeRegionInfo(of);

  double lim = 0.1;
  int num_pt = numPoints();
  for (size_t ki=0; ki<adj_planar.size(); ++ki)
    {
      if ((double)adj_planar[ki]->numPoints() > lim*num_pt)
	{
	  // Get seed points
	  vector<RevEngPoint*> adj_pts = extractNextToAdjacent(adj_planar[ki]);
	  
	  // Grow adjacent
	  adj_planar[ki]->growFromNeighbour(adj_pts, tol, this);
	}
      std::ofstream of2("updated_planar.g2");
      writeRegionInfo(of2);
      for (size_t ki=0; ki<adj_planar.size(); ++ki)
	adj_planar[ki]->writeRegionInfo(of2);

    }
  if (numPoints() != num_pt)
    {
      splitRegion(added_groups);
      updateInfo();
      return true;
    }
  return false;
}

//===========================================================================
bool RevEngRegion::segmentByDirectionContext(int min_point_in, double tol,
					     const Point& dir, double angtol, 
					     vector<vector<RevEngPoint*> >& added_groups)
//===========================================================================
{
  double pihalf = 0.5*M_PI;
  vector<vector<RevEngPoint*> > pnt_groups(3);
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Point norm = group_points_[ki]->getMongeNormal();
      double ang = dir.angle(norm);
      ang = std::min(ang, M_PI-ang);
      if (ang < angtol)
	pnt_groups[0].push_back(group_points_[ki]);
      else if (fabs(pihalf-ang) < angtol)
	pnt_groups[1].push_back(group_points_[ki]);
      else
	pnt_groups[2].push_back(group_points_[ki]);
    }

  int num_pnt_groups = (pnt_groups[0].size() > 0) + (pnt_groups[1].size() > 0) +
    (pnt_groups[2].size() > 0);
  if (num_pnt_groups <= 1)
    return false;

  std::ofstream of("point_groups.g2");
  for (int ka=0; ka<3; ++ka)
    {
      if (pnt_groups[ka].size() > 0)
	{
	  of << "400 1 0 0" << std::endl;
	  of << pnt_groups[ka].size() << std::endl;
	  for (size_t ki=0; ki<pnt_groups[ka].size(); ++ki)
	    of << pnt_groups[ka][ki]->getPoint() << std::endl;
	}
    }
  
  vector<vector<RevEngPoint*> > sep_groups;
  size_t g_ix = 0;
  for (int ka=0; ka<3; ++ka)
    {
      if (pnt_groups[ka].size() > 1)
	{
	  shared_ptr<RevEngRegion> reg(new RevEngRegion(classification_type_,
							edge_class_type_,
							pnt_groups[ka]));
	  vector<vector<RevEngPoint*> > curr_sep_groups;
	  reg->splitRegion(curr_sep_groups);
	  sep_groups.push_back(reg->getPoints());
	  if (curr_sep_groups.size() > 0)
	    sep_groups.insert(sep_groups.end(), curr_sep_groups.begin(),
			      curr_sep_groups.end());
	  for (; g_ix<sep_groups.size(); ++g_ix)
	    for (size_t kj=0; kj<sep_groups[g_ix].size(); ++kj)
	      sep_groups[g_ix][kj]->unsetRegion();
	}
      else if (pnt_groups[ka].size() == 1)
	pnt_groups[ka][0]->unsetRegion();
    }

  std::ofstream of2("sep_groups.g2");
  for (size_t kj=0; kj<sep_groups.size(); ++kj)
    {
      of2 << "400 1 0 0" << std::endl;
      of2 << sep_groups[kj].size() << std::endl;
      for (size_t ki=0; ki<sep_groups[kj].size(); ++ki)
	of2 << sep_groups[kj][ki]->getPoint() << std::endl;
    }

  int max_ix = -1;
  int num_pnts = 0;
  for (size_t ki=0; ki<sep_groups.size(); ++ki)
    {
      if ((int)sep_groups[ki].size() > num_pnts)
	{
	  num_pnts = (int)sep_groups[ki].size();
	  max_ix = (int)ki;
	}
    }

  for (size_t ki=0; ki<sep_groups.size(); ++ki)
    {
      if ((int)ki == max_ix)
	{
	  group_points_.clear();
	  group_points_ = sep_groups[max_ix];
	  for (size_t kj=0; kj<group_points_.size(); ++kj)
	    group_points_[kj]->setRegion(this);
	  updateInfo();
	}
      else
	{
	  for (size_t kj=0; kj<sep_groups[ki].size(); ++kj)
	    sep_groups[ki][kj]->unsetRegion();
	  added_groups.push_back(sep_groups[ki]);
	}
    }

  return (added_groups.size() > 0);
}

//===========================================================================
void RevEngRegion::getPCA(double lambda[3], Point& eigen1, Point& eigen2,
			  Point& eigen3)
//===========================================================================
{
  // PCA analysis of group points
  Vector3D xyz = group_points_[0]->getPoint();
  Point pos0(xyz[0], xyz[1], xyz[2]);
  vector<Point> pnts;
  pnts.reserve(group_points_.size());
  //double wgt = 1.0/(double)group_points_.size();
  for (size_t ki=1; ki<group_points_.size(); ++ki)
    {
      xyz = group_points_[ki]->getPoint();
      Point pos(xyz[0], xyz[1], xyz[2]);
      pnts.push_back(pos);
    }

  double eigenvec[3][3];
  RevEngUtils::principalAnalysis(pos0, pnts, lambda, eigenvec);
  eigen1 = Point(eigenvec[0][0], eigenvec[0][1], eigenvec[0][2]);
  eigen2 = Point(eigenvec[1][0], eigenvec[1][1], eigenvec[1][2]);
  eigen3 = Point(eigenvec[2][0], eigenvec[2][1], eigenvec[2][2]);
}


//===========================================================================
shared_ptr<SplineSurface> RevEngRegion::surfApprox(vector<RevEngPoint*>& points,
						   const BoundingBox& bbox)
//===========================================================================
{
  // PCA analysis of given points to orient plane for parametrerization
  Vector3D xyz = points[0]->getPoint();
  Point pos0(xyz[0], xyz[1], xyz[2]);
  vector<Point> pnts;
  pnts.reserve(points.size());
  Point vec(0.0, 0.0, 0.0);
  double wgt = 1.0/(double)(group_points_.size());
  for (size_t ki=1; ki<points.size(); ++ki)
    {
      Vector3D xyz = points[ki]->getPoint();
      Point pos(xyz[0], xyz[1], xyz[2]);
      pnts.push_back(pos);
      Point vec0 = points[ki]->minCurvatureVec();
      if (vec0*vec < 0.0)
	vec0 *= -1;
      vec += wgt*vec0;
    }
  double lambda[3];
  double eigenvec[3][3];
  RevEngUtils::principalAnalysis(pos0, pnts, lambda, eigenvec);
  Point eigen3(eigenvec[2][0], eigenvec[2][1], eigenvec[2][2]);
  Point ydir = eigen3.cross(vec);
  ydir.normalize_checked();
  Point xdir = ydir.cross(eigen3);
  xdir.normalize_checked();

  // Parameterize points based on planar surface
  pnts.push_back(pos0);
  vector<double> data;
  vector<double> param;
  RevEngUtils::parameterizeWithPlane(pnts, bbox, xdir, ydir,
				     data, param);

  // Surface approximation
  int order = 3;
  double belt = 0.1*bbox.low().dist(bbox.high());
  shared_ptr<SplineSurface> surf = RevEngUtils::surfApprox(data, 3, param, order,
							   order, order, order, belt);
return surf;
}

//===========================================================================
void RevEngRegion::collect(RevEngPoint *pt, RevEngRegion *prev)
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
    normalcone_ = DirectionCone(pt->getMongeNormal());
  else
    normalcone_.addUnionWith(pt->getMongeNormal());
  bool planar = planartype();
  for (size_t kj=0; kj<grouped.size(); ++kj)
    {
      vector<ftSamplePoint*> next = grouped[kj]->getNeighbours();
      for (size_t ki=0; ki<next.size(); ++ki)
	{
	  RevEngPoint *curr = dynamic_cast<RevEngPoint*>(next[ki]);
	  if (!curr)
	    continue;  // Should not happen
	  if (curr->isOutlier())
	    continue;
	  if (curr->hasRegion() && curr->region() != prev)
	    continue;  // Already belonging to a segment
	  if (curr->isEdge(edge_class_type_))
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
	  normalcone_.addUnionWith(curr->getMongeNormal());

	  // Continue growing from this point
	  grouped.push_back(curr);
	}
    }

  // Principal curvature summary
  updateInfo();
}

//===========================================================================
BoundingBox RevEngRegion::getParameterBox()
//===========================================================================
{
  BoundingBox parbox(2);
  if (associated_sf_.size() == 0)
    return parbox;   // No surface means that the points are not parameterized

  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Vector2D par = group_points_[ki]->getPar();
      Point par2(par[0], par[1]);
      parbox.addUnionWith(par2);
    }
  return parbox;
}


//===========================================================================
RevEngPoint* RevEngRegion::seedPointPlane(int min_next, double rfac, double angtol)
//===========================================================================
{
  double min_in = 0;
  int min_ix = -1;
  int multfac = 10;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      double local_len = group_points_[ki]->getMeanEdgLen();
      Point normal = group_points_[ki]->getMongeNormal();
      double radius = 2.0*rfac*local_len;
      vector<RevEngPoint*> nearpts;
      group_points_[ki]->fetchClosePoints2(radius, min_next, 5*min_next,
					   nearpts, this);

      // Count deviant points
      int deviant = 0;
      for (size_t kj=0; kj<nearpts.size(); ++kj)
      {
	Point curr_norm = nearpts[kj]->getMongeNormal();
	double ang = curr_norm.angle(normal);
	if (nearpts[kj]->region() != this || ang > angtol)
	  ++deviant;
      }

      if (deviant == 0)
	return group_points_[ki];   // Seed point found

      if ((int)nearpts.size() - deviant > min_in)
	{
	  min_in = (int)nearpts.size() - deviant;
	  min_ix = (int)ki;
	  if (min_in < multfac*min_next)
	    break;
	}
    }
  return (min_ix >= 0) ? group_points_[min_ix] : 0;
}

//===========================================================================
void RevEngRegion::growLocalPlane(double tol, vector<RevEngPoint*>& plane_pts,
				  shared_ptr<Plane>& plane_out)
//===========================================================================
{
  // Get seed point for local grow
  int min_next = std::max(10, (int)group_points_.size()/100);
  double rfac = 3;
  double angtol = 0.1;
  RevEngPoint* seed = seedPointPlane(min_next, rfac, angtol);
  if (!seed)
    return;
 
  std::ofstream of("seed.g2");
  of << "400 1 0 4 200 0 55 255" << std::endl;
  of << "1" << std::endl;
  of << seed->getPoint() << std::endl;

  // Fetch nearby points belonging to the same region
  double local_len = seed->getMeanEdgLen();
  double radius = 2.0*rfac*local_len;
  vector<RevEngPoint*> nearpts;
  seed->fetchClosePoints2(radius, min_next, 5*min_next, nearpts, this);
  nearpts.insert(nearpts.begin(), seed);
  if (nearpts.size() < min_next/2)
    return;

  // Approximate with plane
  shared_ptr<Plane> plane = computePlane(nearpts, avnorm_);

  double eps = 1.0e-6;
  double maxdist, avdist;
  int num_inside;
  vector<pair<double, double> > dist_ang;
  vector<double> parvals;
  vector<RevEngPoint*> inpt, outpt;
  RevEngUtils::distToSurf(nearpts.begin(), nearpts.end(),
			  plane, tol, maxdist, avdist, num_inside, inpt,
			  outpt, parvals, dist_ang, angtol);

  if (maxdist > tol)
    return;   // For the time being

  
  for (size_t ki=0; ki<nearpts.size(); ++ki)
    nearpts[ki]->setVisited();
  size_t prev_size = nearpts.size();
  while (true)
    {
      for (size_t ki=0; ki<nearpts.size(); ++ki)
	{
	  vector<ftSamplePoint*> next = nearpts[ki]->getNeighbours();
	  for (size_t kj=0; kj<next.size(); ++kj)
	    {
	      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(next[kj]);
	      if ((!pt->hasRegion()) || pt->region() != this)
		continue;
	      if (pt->visited())
		continue;
	      Vector3D xyz = pt->getPoint();

	      double upar, vpar, dist;
	      Point close;
	      plane->closestPoint(Point(xyz[0],xyz[1],xyz[2]), upar, vpar, close, dist, eps);
	      pt->setVisited();
	      if (dist < tol)
		nearpts.push_back(pt);
	    }
	}

      for (size_t ki=0; ki<group_points_.size(); ++ki)
	{
	  if (std::find(nearpts.begin(), nearpts.end(), group_points_[ki]) == nearpts.end())
	    group_points_[ki]->unsetVisited();
	}

      if (nearpts.size() == prev_size)
	break;
      
      shared_ptr<Plane> plane2 = computePlane(nearpts, avnorm_);

      double maxdist2, avdist2;
      int num_inside2;
      vector<pair<double, double> > dist_ang2;
      vector<double> parvals2;
     RevEngUtils::distToSurf(nearpts.begin(), nearpts.end(),
			      plane2, tol, maxdist2, avdist2, num_inside2, inpt,
			      outpt, parvals2, dist_ang2, angtol);
      if (maxdist2 > tol)
	break;
      
      prev_size = nearpts.size();
      plane = plane2;
      maxdist = maxdist2;
      avdist = avdist2;
      num_inside = num_inside2;
      parvals = parvals2;
      dist_ang = dist_ang2;
    }

  for (size_t ki=0; ki<group_points_.size(); ++ki)
    group_points_[ki]->unsetVisited();
  
  std::ofstream of2("plane_pts.g2");
  of2 << "400 1 0 4 100 155 0 255" << std::endl;
  of2 << nearpts.size() << std::endl;
  for (size_t ki=0; ki<nearpts.size(); ++ki)
    of2 << nearpts[ki]->getPoint() << std::endl;

  plane->writeStandardHeader(of);
  plane->write(of);

  plane_pts = nearpts;
  plane_out = plane;

  int stop_break = 1;

}

//===========================================================================
void RevEngRegion::segmentByPlaneGrow(double tol, int min_pt, 
				      vector<shared_ptr<HedgeSurface> >& hedgesfs,
				      vector<HedgeSurface*>& prevsfs,
				      vector<vector<RevEngPoint*> >& out_groups)
//===========================================================================
{
  // Extract group of planar points and associated plane
  int min_nmb_pts = std::max((int)group_points_.size()/100, 10);
  vector<RevEngPoint*> plane_pts;
  shared_ptr<Plane> plane;
  growLocalPlane(tol, plane_pts, plane);
  if ((int)plane_pts.size() < min_nmb_pts || plane_pts.size() == group_points_.size())
    return;

  // Fetch accuracy information
  double angtol = 0.1;
   double maxd, avd;
  int num_in;
  vector<pair<double, double> > dist_ang;
  vector<double> parvals;
  vector<RevEngPoint*> inpt, outpt;
  RevEngUtils::distToSurf(plane_pts.begin(), plane_pts.end(),
			  plane, tol, maxd, avd, num_in, inpt,
			  outpt, parvals, dist_ang, angtol);
  
  // Extract remaining points
  for (size_t ki=0; ki<plane_pts.size(); ++ki)
    plane_pts[ki]->setVisited();

  vector<RevEngPoint*> remaining_pts;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      if (!group_points_[ki]->visited())
	{
	  remaining_pts.push_back(group_points_[ki]);
	  removePoint(group_points_[ki]);
	}
    }

  for (size_t ki=0; ki<plane_pts.size(); ++ki)
    plane_pts[ki]->unsetVisited();

  if (hasSurface())
    {
      // No longer valid
      prevsfs.insert(prevsfs.end(), associated_sf_.begin(), associated_sf_.end());
      clearSurface();
    }
  
  if ((int)plane_pts.size() >= min_pt)
    {
      // Register current plane
      for (size_t ki=0; ki<plane_pts.size(); ++ki)
	{
	  plane_pts[ki]->setPar(Vector2D(parvals[2*ki],parvals[2*ki+1]));
	  plane_pts[ki]->setSurfaceDist(dist_ang[ki].first, dist_ang[ki].second);
	}
      setAccuracy(maxd, avd, num_in);
      shared_ptr<HedgeSurface> hedge(new HedgeSurface(plane, this));
      associated_sf_.push_back(hedge.get());
      hedgesfs.push_back(hedge);
    }

  // Store as base surface and primary surface
  setBaseSf(plane, maxd, avd, num_in);
  setPrimarySf(plane, maxd, avd, num_in);

  updateInfo();

  // Distribute remaining points into groups and create regions
  shared_ptr<RevEngRegion> reg(new RevEngRegion(classification_type_,
						 edge_class_type_,
						 remaining_pts));
  vector<vector<RevEngPoint*> > connected;
  reg->splitRegion(connected);
  out_groups.push_back(reg->getPoints());
  if (connected.size() > 0)
    out_groups.insert(out_groups.end(), connected.begin(), connected.end());
}

/// ===========================================================================
// void RevEngRegion::growLocal(RevEngPoint* seed, double tol, double radius,
// 			     int min_close, vector<RevEngPoint*>& outpts)
// ===========================================================================
// {
//   Fetch nearby points belonging to the same region
//   vector<RevEngPoint*> nearpts;
//   seed->fetchClosePoints2(radius, min_close, 5*min_close, nearpts, this);
//   nearpts.insert(nearpts.begin(), seed);
//   for (size_t ki=0; ki<nearpts.size(); )
//     {
//       if (nearpts[ki]->getPointDistance() > 1.5*nearpts[ki]->getAveragePointDistance())
// 	nearpts.erase(nearpts.begin()+ki);
//       else
// 	++ki;
//     }

//   Approximate with implicit algebraic surface
//   int degree = 2;
//   impl_ = shared_ptr<ImplicitApprox>(new ImplicitApprox());
//   impl_->approx(nearpts, degree);
//   std::ofstream outviz("implsf.g2");
//   impl_->visualize(nearpts, outviz);

//   std::ofstream no("noproject.g2");
//   double large = 1.0e3;
//   double maxdist = 0.0, avdist = 0.0;
//   int nodist = 0;
//   for (size_t ki=0; ki<nearpts.size(); ++ki)
//     {
//       double dist = fabs(impl_->estimateDist(nearpts[ki]));
//       if (dist > large)
// 	{
// 	  no << "400 1 0 4 0 0 255 255" << std::endl;
// 	  no << "1" << std::endl;
// 	  no << nearpts[ki]->getPoint() << std::endl;
// 	}
//       else
// 	{
// 	  maxdist = std::max(maxdist, dist);
// 	  avdist += dist;
// 	  ++nodist;
// 	}
//     }
//   avdist /= (double)nodist;
//   double tol2 = 1.5*maxdist; std::min(tol, 10.0*maxdist);
//   std::ofstream ofnear("nearpts.g2");
//   ofnear << "400 1 0 4 0 255 0 255" << std::endl;
//   ofnear << nearpts.size() << std::endl;
//   for (size_t ki=0; ki<nearpts.size(); ++ki)
//     ofnear << nearpts[ki]->getPoint() << std::endl;
   
//   Check accuracy
//   std::ofstream of("grow_local.g2");
//   vector<RevEngPoint*> out, in;
//   vector<RevEngPoint*> core;
//   core.push_back(seed);
//   seed->setVisited();
//   int max_new_out = std::max(50, (int)group_points_.size()/100);
//   int prev_out = 0;
//   for (size_t ki=0; ki<core.size(); ++ki)
//     {
//       vector<ftSamplePoint*> next = core[ki]->getNeighbours();
//       for (size_t kj=0; kj<next.size(); ++kj)
// 	{
// 	  RevEngPoint *pt = dynamic_cast<RevEngPoint*>(next[kj]);
// 	  if (pt->visited())
// 	    continue;
// 	  if (pt->region() != this)
// 	    continue;    Can consider storing the point for later access
// 	  pt->setVisited();
// 	  core.push_back(pt);
	  
// 	  double dist = fabs(impl_->estimateDist(pt));
// 	  if (dist > large)
// 	    {
// 	      no << "400 1 0 4 0 0 255 255" << std::endl;
// 	      no << "1" << std::endl;
// 	      no << pt->getPoint() << std::endl;
// 	    }
//  	  if (dist <= tol2)
// 	    {
// 	      in.push_back(pt);
// 	    }
// 	  else
// 	    out.push_back(pt);
// 	}

//       if ((int)out.size() > prev_out + max_new_out)
// 	{
// 	  Update implicit
// 	  vector<RevEngPoint*> curr_pts(in.begin(), in.end());
// 	  curr_pts.insert(curr_pts.end(), out.begin(), out.end());
// 	  shared_ptr<ImplicitApprox> impl2(new ImplicitApprox());
// 	  impl2->approx(curr_pts, degree);
// 	  std::ofstream outviz2("implsf2.g2");
// 	  impl2->visualize(curr_pts, outviz2);
// 	  std::ofstream ofnear2("nearpts2.g2");
// 	  ofnear2 << "400 1 0 4 0 255 0 255" << std::endl;
// 	  ofnear2 << curr_pts.size() << std::endl;
// 	  for (size_t ki=0; ki<curr_pts.size(); ++ki)
// 	    ofnear2 << curr_pts[ki]->getPoint() << std::endl;
 	  
// 	  vector<RevEngPoint*> out2, in2;
// 	  for (size_t kr=0; kr<curr_pts.size(); ++kr)
// 	    {
// 	      double dist = fabs(impl2->estimateDist(curr_pts[kr]));
// 	      if (dist <= tol2)
// 		{
// 		  in2.push_back(curr_pts[kr]);
// 		}
// 	      else
// 		out2.push_back(curr_pts[kr]);
// 	    }
	  
// 	  if (out2.size() < out.size())
// 	    {
// 	      impl_ = impl2;
// 	      in = in2;
// 	      out = out2;
// 	    }
// 	  else
// 	    break;
// 	  prev_out = (int)out.size();
// 	}
//     }

//   Remaining points
//   for (size_t ki=0; ki<group_points_.size(); ++ki)
//     {
//       RevEngPoint *pt = dynamic_cast<RevEngPoint*>(group_points_[ki]);
//       double dist = fabs(impl_->estimateDist(pt));
//       if (dist <= tol2)
// 	{
// 	  in.push_back(pt);
// 	}
//       else
// 	out.push_back(pt);
//     }
  
//   if (in.size() > 0)
//     {
//       of << "400 1 0 0" << std::endl;
//       of << in.size() << std::endl;
//       for (size_t kr=0; kr<in.size(); ++kr)
// 	of << in[kr]->getPoint() << std::endl;
//       of << std::endl;
//     }
  
//   if (out.size() > 0)
//     {
//       of << "400 1 0 0" << std::endl;
//       of << out.size() << std::endl;
//       for (size_t kr=0; kr<out.size(); ++kr)
// 	of << out[kr]->getPoint() << std::endl;
//     }
  
//   of << "400 1 0 0" << std::endl;
//   of << group_points_.size() << std::endl;
//   for (size_t kr=0; kr<group_points_.size(); ++kr)
//     of << group_points_[kr]->getPoint() << std::endl;
//   of << std::endl;

//   int stop_break = 1;
  
// }

//===========================================================================
bool RevEngRegion::extractPlane(double tol, int min_pt, double angtol,
				int prefer_elementary_,
				vector<shared_ptr<HedgeSurface> >& hedgesfs,
				vector<HedgeSurface*>& prevsfs,
				vector<vector<RevEngPoint*> >& out_groups)
//===========================================================================
{
  // std::ofstream of("curr_region.g2");
  // writeRegionInfo(of);
  
  // std::ofstream of2("curr_normals.g2");
  // writeUnitSphereInfo(of2);

  if (normalcone_.greaterThanPi())
    {
      std::cout << "Greater than pi" << std::endl;
      std::ofstream ofpi("pi_region.g2");
      writeRegionInfo(ofpi);
    }
  
  bool found = false;
  int min_nmb = 20;
  if ((int)group_points_.size() < min_nmb)
    return false;

  double angtol2 = 0.1; //*M_PI;
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

  vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*>  > adj_elem, adj_elem_base, adj_primary;
  getAdjacentElemInfo(adj_elem, adj_elem_base, adj_primary);
  
  shared_ptr<Plane> surf3 = computePlane(group_points_, avnorm_); //(new Plane(pos3, normal3));

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
  vector<pair<double, double> > dist_ang;
  vector<double> parvals;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			  surf3, tol, maxdist3, avdist3, num_inside3, inpt,
			  outpt, parvals, dist_ang, angtol);
  std::ofstream ofd("in_out_plane.g2");
  ofd << "400 1 0 4 155 50 50 255" << std::endl;
  ofd << inpt.size() << std::endl;
  for (size_t kr=0; kr<inpt.size(); ++kr)
    ofd << inpt[kr]->getPoint() << std::endl;
  ofd << "400 1 0 4 50 155 50 255" << std::endl;
  ofd << outpt.size() << std::endl;
  for (size_t kr=0; kr<outpt.size(); ++kr)
    ofd << outpt[kr]->getPoint() << std::endl;

  if (inpt.size() > group_points_.size()/3 && (int)inpt.size() > min_nmb &&
      (normalcone_.angle() > angtol2 ||
       normalcone_.centre().angle(surf3->getNormal()) > angtol2) &&
      (!(basesf_.get() && basesf_->instanceType() == Class_Cylinder)))
    {
      vector<RevEngPoint*> ang_points;
      identifyAngPoints(dist_ang, angtol2, ang_points);
      if (ang_points.size() < group_points_.size()/2)
	{
	  extractSpesPoints(ang_points, out_groups, true);
	  if (out_groups.size() > 0)
	    {
		shared_ptr<Plane> plane_in =
		  computePlane(group_points_, avnorm_);
	      vector<RevEngPoint*> inpt_in, outpt_in; //, inpt2, outpt2;
	      vector<pair<double, double> > dist_ang_in;
	      vector<double> parvals_in;
	      double maxd_in, avd_in;
	      int num2_in;
	      RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
				      plane_in, tol, maxd_in, avd_in, num2_in,
				      inpt_in, outpt_in, parvals_in, dist_ang_in, angtol);
  
	      std::ofstream ofd2("in_out_plane2.g2");
	      ofd2 << "400 1 0 4 155 50 50 255" << std::endl;
	      ofd2 << inpt_in.size() << std::endl;
	      for (size_t kr=0; kr<inpt_in.size(); ++kr)
		ofd2 << inpt_in[kr]->getPoint() << std::endl;
	      ofd2 << "400 1 0 4 50 155 50 255" << std::endl;
	      ofd2 << outpt_in.size() << std::endl;
	      for (size_t kr=0; kr<outpt_in.size(); ++kr)
		ofd2 << outpt_in[kr]->getPoint() << std::endl;

	      // if (num2_in > num_inside3 || avd_in < avdist3)
	      // 	{
	      std::swap(surf3, plane_in);
	      std::swap(num_inside3, num2_in);
	      std::swap(avdist3, avd_in);
	      std::swap(maxdist3, maxd_in);
	      std::swap(parvals, parvals_in);
	      std::swap(dist_ang, dist_ang_in);
	      // 	  std::cout << "Plane swap" << std::endl;
	      // 	  std::cout << group_points_.size() << " " << num2_in << " ";
	      // std::cout<< num_inside3 << " " << avd_in << " " << avdist3;
	      // 	  std::cout << " " << maxd_in << " " << maxdist3 << std::endl;
	      // }
	    }
	}
    }

  Point low = bbox_.low();
  Point high = bbox_.high();
  //double len = low.dist(high);
  //surf3->setParameterBounds(-len, -len, len, len);
  std::ofstream plane("curr_plane.g2");
  surf3->writeStandardHeader(plane);
  surf3->write(plane);

  int num = (int)group_points_.size();
  if (accuracyOK(min_pt, tol, num_inside3, avdist3))
    {
      found = true;
      for (size_t kh=0; kh<group_points_.size(); ++kh)
	{
	  group_points_[kh]->setPar(Vector2D(parvals[2*kh],parvals[2*kh+1]));
	  group_points_[kh]->setSurfaceDist(dist_ang[kh].first, dist_ang[kh].second);
	}
      setAccuracy(maxdist3, avdist3, num_inside3);
      
      std::cout << "Plane. N1: " << num << ", N2: " << num_inside3 << ", max: " << maxdist3 << ", av: " << avdist3 << std::endl;

      shared_ptr<HedgeSurface> hedge(new HedgeSurface(surf3, this));
      for (size_t kh=0; kh<associated_sf_.size(); ++kh)
	prevsfs.push_back(associated_sf_[kh]);
      setHedge(hedge.get());
      hedgesfs.push_back(hedge);
	  
    }
  if (!primary_.get() ||
      (num_inside3 >= num_in_primary_ && avdist3 < avdist_primary_))
    setPrimarySf(surf3, maxdist3, avdist3, num_inside3);

  return found;
}


//===========================================================================
shared_ptr<Plane> RevEngRegion::computePlane(vector<RevEngPoint*>& points,
					     const Point& norm_dir)
//===========================================================================
{
  Point normal1 = normalcone_.centre();
  Point normal(0.0, 0.0, 0.0);
  Point pos(0.0, 0.0, 0.0);
  double wgt = 1.0/(double)(points.size());
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      Point curr = points[ki]->getMongeNormal();
      Vector3D xyz = points[ki]->getPoint();
      normal += wgt*curr;
      pos += wgt*Point(xyz[0], xyz[1], xyz[2]);
    }
  
  impl_ = shared_ptr<ImplicitApprox>(new ImplicitApprox());
  impl_->approx(points, 1);
  Point pos3, normal3;
  impl_->projectPoint(pos, normal, pos3, normal3);
  if (normal3*norm_dir < 0.0)
    normal3 *= -1.0;
  
  shared_ptr<Plane> surf(new Plane(pos3, normal3));
  Point low = bbox_.low();
  Point high = bbox_.high();
  //double len = low.dist(high);
  //surf->setParameterBounds(-len, -len, len, len);
  
  bool plane_project = false;
  if (plane_project)
    {
      vector<pair<vector<RevEngPoint*>::iterator,
		  vector<RevEngPoint*>::iterator> > group;
      group.push_back(std::make_pair(points.begin(), points.end()));
      Point vec1, vec2;
      surf->getSpanningVectors(vec1, vec2);
      vector<Point> projected1, projected2;
      double maxdp1, avdp1, maxdp2, avdp2;
      RevEngUtils::projectToPlane(group, vec1, pos3, projected1, maxdp1, avdp1);
      RevEngUtils::projectToPlane(group, vec2, pos3, projected2, maxdp2, avdp2);

      std::ofstream ofp3("plane_project.g2");
      ofp3 << "400 1 0 4 255 0 0 255" << std::endl;
      ofp3 << projected1.size() << std::endl;
      for (size_t kr=0; kr<projected1.size(); ++kr)
	ofp3 << projected1[kr] << std::endl;
      ofp3 << "400 1 0 4 255 0 0 255" << std::endl;
      ofp3 << projected2.size() << std::endl;
      for (size_t kr=0; kr<projected2.size(); ++kr)
	ofp3 << projected2[kr] << std::endl;
    }
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
      Point curr = group_points_[ki]->getMongeNormal();
      avnorm += wgt*curr;
    }

  int nmb_in = 0;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Point curr = group_points_[ki]->getMongeNormal();
      double ang = curr.angle(avnorm);
      if (ang <= angtol)
	++nmb_in;
    }

  double frac = (double)nmb_in/(double)group_points_.size();
  return (frac >= inlim);
}

//===========================================================================
vector<RevEngRegion*> RevEngRegion::fetchAdjacentPlanar()
//===========================================================================
{
  vector<RevEngRegion*> adjacent_planar;
  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
    {
      bool found = false;
      if ((*it)->hasSurface())
	{
	  shared_ptr<ParamSurface> sf = (*it)->getSurface(0)->surface();
	  if (sf->instanceType() == Class_Plane)
	    adjacent_planar.push_back(*it);
	}
    }

  return adjacent_planar;
}

//===========================================================================
Point RevEngRegion::directionFromAdjacent(double angtol)
//===========================================================================
{
  Point dir;
  vector<std::pair<Point,int> > all_dir;
  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
    {
      if ((*it)->hasSurface())
	{
	  shared_ptr<ParamSurface> sf = (*it)->getSurface(0)->surface();
	  if (sf->instanceType() == Class_Plane || sf->instanceType() == Class_Cylinder)
	    {
	      shared_ptr<ElementarySurface> elem =
		dynamic_pointer_cast<ElementarySurface,ParamSurface>(sf);
	      Point curr_dir = elem->direction();
	      size_t kj;
	      for (kj=0; kj<all_dir.size(); ++kj)
		{
		  double ang = curr_dir.angle(all_dir[kj].first);
		  if (M_PI-ang < ang)
		    {
		      curr_dir *= -1.0;
		      ang = M_PI - ang;
		    }
		  if (ang < angtol)
		    {
		      int num = (*it)->numPoints();
		      double fac = 1.0/(double)(all_dir[kj].second+num);
		      curr_dir = fac*(all_dir[kj].second*all_dir[kj].first + num*curr_dir);
		      break;
		    }
		}
	      if (kj == all_dir.size())
		all_dir.push_back(std::make_pair(curr_dir, (*it)->numPoints()));
	    }
	}
    }

  for (size_t ki=0; ki<all_dir.size(); ++ki)
    for (size_t kj=ki+1; kj<all_dir.size(); ++kj)
      if (all_dir[kj].second > all_dir[ki].second)
	std::swap(all_dir[ki], all_dir[kj]);

  if (all_dir.size() > 0)
    dir = all_dir[0].first;

  return dir;
}

//===========================================================================
void RevEngRegion::getAdjacentElemInfo(vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> >& adj_elem,
				       vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> >& adj_elem_base,
				       vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> >& adj_primary)
//===========================================================================
{
  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
    {
      bool found = false;
      if ((*it)->hasSurface())
	{
	  shared_ptr<ParamSurface> sf = (*it)->getSurface(0)->surface();
	  shared_ptr<ElementarySurface> elem =
	    dynamic_pointer_cast<ElementarySurface,ParamSurface>(sf);
	  if (elem.get())
	    {
	      adj_elem.push_back(std::make_pair(elem,*it));
	      found = true;
	    }
	}
      
      if ((*it)->hasBaseSf() && (!found))
	{
	  shared_ptr<ParamSurface> sf = (*it)->getBase();
	  shared_ptr<ElementarySurface> elem =
	    dynamic_pointer_cast<ElementarySurface,ParamSurface>(sf);
	  if (elem.get())
	    {
	      adj_elem_base.push_back(std::make_pair(elem, *it));
	      found = true;
	    }
	}
      
      if ((*it)->hasPrimary() && (!found))
	{
	  shared_ptr<ParamSurface> sf = (*it)->getPrimary();
	  shared_ptr<ElementarySurface> elem =
	    dynamic_pointer_cast<ElementarySurface,ParamSurface>(sf);
	  if (elem.get())
	    adj_primary.push_back(std::make_pair(elem, *it));
	}
      
    }
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

      Point norm = group_points_[kr]->getMongeNormal();
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
  //double vark1 = ((double)group_points_.size()*k1_1)/fabs(d1) - fabs(k1_2);
  //double vark2 = ((double)group_points_.size()*k2_1)/fabs(d2) - fabs(k2_2);

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

  //double frac1 = (double)nmb1/(double)group_points_.size();
  //double frac2 = (double)nmb2/(double)group_points_.size();

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
      Point curr = group_points_[ki]->getMongeNormal();
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
      vec[ki] = group_points_[ki]->getMongeNormal();
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
  //int num = (int)vec.size();
  vector<double> distance;
  RevEngUtils::distToSurf(vec, surf, tol, maxdist, avdist, num_inside,
			  distance);

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
bool RevEngRegion::feasiblePlane(double zero_H, double zero_K) const
//===========================================================================
{
  double angtol = 0.2;
  double in_lim = 0.9;
  double in_lim2 = 0.5;
  double fac = 5.0;
  if (basesf_ && basesf_->instanceType() == Class_Plane)
    return true;

  if (normalcone_.angle() < angtol && (!normalcone_.greaterThanPi()))
    return true;

  int nmb_in = 0;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      Point normal = group_points_[kr]->getMongeNormal();
      if (avnorm_.angle(normal) <= angtol)
	nmb_in++;
    }
  double in_frac = (double)nmb_in/(double)group_points_.size();
  if (in_frac > in_lim)
    return true;

  if (in_frac < in_lim2)
    return false;
  
  if (std::max(MAH_, MAK_) > fac*std::max(zero_H, zero_K))
    return false;

  return true;
}

//===========================================================================
bool RevEngRegion::feasibleCylinder(double zero_H, double zero_K) const
//===========================================================================
{
  double fac1 = 5.0;
  double fac2 = 0.5;
  if (basesf_ && basesf_->instanceType() == Class_Cylinder)
    return true;
  if (MAK_ > fac1*zero_H)
    return false;
  if (MAK_ > fac2*MAH_)
    return false;
  // if (fabs(avK_) < fac2*MAK_)
  //   return false;

  return true;
}


//===========================================================================
bool RevEngRegion::extractCylinder(double tol, int min_pt, double angtol,
				   double mean_edge_len,
				   int prefer_elementary_,
				   vector<shared_ptr<HedgeSurface> >& hedgesfs,
				   vector<HedgeSurface*>& prevsfs,
				   vector<vector<RevEngPoint*> >& out_groups,
				   bool& repeat)
//===========================================================================
{
  // std::ofstream of("curr_region.g2");
  // writeRegionInfo(of);
  
  // std::ofstream of2("curr_normals.g2");
  // writeUnitSphereInfo(of2);

  bool found = false;
  int min_nmb = 20;
  //double eps = 1.0e-6;
  if ((int)group_points_.size() < min_nmb)
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

  vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> > adj_elem, adj_elem_base, adj_primary;
  getAdjacentElemInfo(adj_elem, adj_elem_base, adj_primary);

  vector<vector<RevEngPoint*> > configs;
  shared_ptr<Cylinder> cyl = computeCylinder(group_points_, tol, configs);
  std::ofstream ofs("cyl.g2");
  shared_ptr<Cylinder> cyl2 = cyl;
  double diag = bbox_.low().dist(bbox_.high());
  cyl2->setParamBoundsV(-0.5*diag,0.5*diag);
  cyl2->writeStandardHeader(ofs);
  cyl2->write(ofs);
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
  vector<pair<double, double> > dist_ang;
  vector<double> parvals;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			  cyl, tol, maxd, avd, num2, inpt, outpt, 
			  parvals, dist_ang, angtol);

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

  // Point axis_in, Cx_in, Cy_in, pos_in;
  // double rad_in;
  // Point low = bbox_.low();
  // Point high = bbox_.high();
  // vector<pair<vector<RevEngPoint*>::iterator,
  // 	      vector<RevEngPoint*>::iterator> > group_in;
  // group_in.push_back(std::make_pair(inpt.begin(), inpt.end()));
  // RevEngUtils::computeAxis(group_in, axis_in, Cx_in, Cy_in);
  // RevEngUtils::computeCylPosRadius(group_in, low, high, axis_in, Cx_in,
  // 				   Cy_in, pos_in, rad_in);
  // if (false) //inpt.size() > group_points_.size()/2 && (int)inpt.size() > min_nmb)
  //   {
  //     vector<vector<RevEngPoint*> > configs2;
  //     shared_ptr<Cylinder> cyl_in = computeCylinder(inpt, tol, configs2);
  //     vector<RevEngPoint*> inpt_in, outpt_in; //, inpt2, outpt2;
  //     vector<pair<double, double> > dist_ang_in;
  //     vector<double> parvals_in;
  //     double maxd_in, avd_in;
  //     int num2_in;
  //     RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
  // 			      cyl_in, tol, maxd_in, avd_in, num2_in,
  // 			      inpt_in, outpt_in, parvals_in, dist_ang_in, angtol);
  
  //     std::ofstream ofd2("in_out_cyl2.g2");
  //     ofd2 << "400 1 0 4 155 50 50 255" << std::endl;
  //     ofd2 << inpt_in.size() << std::endl;
  //     for (size_t kr=0; kr<inpt_in.size(); ++kr)
  // 	ofd2 << inpt_in[kr]->getPoint() << std::endl;
  //     ofd2 << "400 1 0 4 50 155 50 255" << std::endl;
  //     ofd2 << outpt_in.size() << std::endl;
  //     for (size_t kr=0; kr<outpt_in.size(); ++kr)
  // 	ofd2 << outpt_in[kr]->getPoint() << std::endl;

  //     if (num2_in > num2 || avd_in < avd)
  // 	{
  // 	  std::swap(cyl, cyl_in);
  // 	  std::swap(num2, num2_in);
  // 	  std::swap(avd, avd_in);
  // 	  std::swap(maxd, maxd_in);
  // 	  std::swap(parvals, parvals_in);
  // 	  std::swap(dist_ang, dist_ang_in);
  // 	  // std::cout << "Cylinder swap" << std::endl;
  // 	  // std::cout << group_points_.size() << " " << num2_in;
  // 	  // std::cout << " " << num2 << " " << avd_in << " " << avd;
  // 	  // std::cout << " " << maxd_in << " " << maxd << std::endl;
  // 	}
  //   }
  if (configs.size() > 1 && (!hasSurface()))
    {
      // Split point group and try again
      int keep_ix = -1;
      int keep_nmb = 0;
      for (size_t ki=0; ki<configs.size(); ++ki)
	if ((int)configs[ki].size() > keep_nmb)
	  {
	    keep_nmb = (int)configs[ki].size();
	    keep_ix = (int)ki;
	  }
      
      for (size_t ki=0; ki<configs.size(); ++ki)
	{
	  if ((int)ki == keep_ix)
	    continue;

	  extractSpesPoints(configs[ki], out_groups);
	}

      // Check that the remaing point cloud is connected
      splitRegion(out_groups);

      if (hasSurface())
	{
	  for (size_t ki=0; ki<associated_sf_.size(); ++ki)
	    prevsfs.push_back(associated_sf_[ki]);
	  clearSurface();
	}
      
      updateInfo();
      repeat = true;
    }
  else
    {
      int num = (int)group_points_.size();
      if (num2 > min_pt && num2 > num/2)
	{
	  // Check for deviant points at the boundary
	  vector<RevEngPoint*> dist_points;
	  identifyDistPoints(dist_ang, tol, maxd, avd, dist_points);
	  extractSpesPoints(dist_points, out_groups, true);
	  if (out_groups.size() > 0)
	    {
	      vector<vector<RevEngPoint*> > configs2;
	      shared_ptr<Cylinder> cyl_in = computeCylinder(group_points_, tol, configs2);
	      vector<RevEngPoint*> inpt_in, outpt_in; 
	      vector<pair<double, double> > dist_ang_in;
	      vector<double> parvals_in;
	      double maxd_in, avd_in;
	      int num2_in;
	      RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
				      cyl_in, tol, maxd_in, avd_in, num2_in,
				      inpt_in, outpt_in, parvals_in, dist_ang_in, angtol);

	      
	      std::ofstream ofd2("in_out_cylinder2.g2");
	      ofd2 << "400 1 0 4 155 50 50 255" << std::endl;
	      ofd2 << inpt_in.size() << std::endl;
	      for (size_t kr=0; kr<inpt_in.size(); ++kr)
		ofd2 << inpt_in[kr]->getPoint() << std::endl;
	      ofd2 << "400 1 0 4 50 155 50 255" << std::endl;
	      ofd2 << outpt_in.size() << std::endl;
	      for (size_t kr=0; kr<outpt_in.size(); ++kr)
		ofd2 << outpt_in[kr]->getPoint() << std::endl;

	      std::swap(cyl, cyl_in);
	      std::swap(num2, num2_in);
	      std::swap(avd, avd_in);
	      std::swap(maxd, maxd_in);
	      std::swap(parvals, parvals_in);
	      std::swap(dist_ang, dist_ang_in);
	    }
	}
	  
      double maxd_init, avd_init;
      int num_init;
      getAccuracy(maxd_init, avd_init, num_init);
      if (accuracyOK(min_pt, tol, num2, avd))
	{
	  bool OK = true;
	  double acc_fac = 1.5;
	  if (associated_sf_.size() > 0)
	    {
	      int sfcode;
	      int sftype = associated_sf_[0]->instanceType(sfcode);
	      double ang = (sftype == Class_Plane) ?
		normalcone_.angle() : M_PI;
	      double ang_lim = 0.1*M_PI;
	      
	      // Check with current approximating surface
	      if (prefer_elementary_ == ALWAYS_ELEM ||
		  prefer_elementary_ == PREFER_ELEM)
		{
		  if (!(ang > ang_lim && (num2 < num_init ||
					  (avd < avd_init &&
					   num2 < acc_fac*num_init))))
		    OK = false;
		}
	      else
		{
		  if (!(num2 < num_init ||
			(avd < avd_init && num2 < acc_fac*num_init)))
		    OK = true;
		}
	    }

	  if (OK)
	    {
	      found = true;
	      for (size_t kh=0; kh<group_points_.size(); ++kh)
		{
		  group_points_[kh]->setPar(Vector2D(parvals[2*kh],parvals[2*kh+1]));
		  group_points_[kh]->setSurfaceDist(dist_ang[kh].first, dist_ang[kh].second);
		}
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
	      setHedge(hedge.get());
	      hedgesfs.push_back(hedge);
	      // for (int ka=0; ka<divcyl[0]->nmbEntities(); ++ka)
	      // 	{
	      // 	  shared_ptr<ParamSurface> cyl2 = divcyl[0]->getSurface(ka);
	      // 	  shared_ptr<HedgeSurface> hedge(new HedgeSurface(cyl2, this));
	      // 	  for (size_t kh=0; kh<associated_sf_.size(); ++kh)
	      // 	    prevsfs.push_back(associated_sf_[kh]);
	      // 	  associated_sf_.push_back(hedge.get());
	      // 	  hedgesfs.push_back(hedge);
	  
	    }
	}
    }
  if (!primary_.get() ||
      (num2 >= num_in_primary_ && avd < avdist_primary_))
    setPrimarySf(cyl, maxd, avd, num2);
  int stop_break0 = 1;
  return found;
}

//===========================================================================
shared_ptr<Cylinder>
RevEngRegion::computeCylinder(vector<RevEngPoint*>& points, double tol,
			      vector<vector<RevEngPoint*> >& configs)
//===========================================================================
{
  // Cylinder orientation by covariance matrix of normal vectors
  Point axis;
  Point Cx;
  Point Cy;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(points.begin(), points.end()));
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

  Point first = points[0]->getMongeNormal();
  vector<Point> rest;
  rest.reserve(points.size());
  Vector3D xyz = points[0]->getPoint();
  Point mid(xyz[0], xyz[1], xyz[2]);
  double wgt = 1.0/(double)points.size();
  mid *= wgt;
  for (size_t ki=1; ki<points.size(); ++ki)
    {
      rest.push_back(points[ki]->getMongeNormal());
      xyz = points[ki]->getPoint();
      mid += wgt*Point(xyz[0], xyz[1], xyz[2]);
    }

  double lambda[3];
  double eigen[3][3];
  RevEngUtils::principalAnalysis(first, rest, lambda, eigen);

  // rest.push_back(first);
  // ImplicitApprox impl;
  // impl.approxPoints(rest, 1);

  // Point axisnorm, projpos;
  // impl.projectPoint(pnt, axis, projpos, axisnorm);
  // axisnorm.normalize_checked();
  // if (axis*axisnorm < 0.0)
  //   axisnorm *= -1;
  // Point Cynorm = Cy - (Cy*axisnorm)*axisnorm;
  // Cynorm.normalize();
  // Point Cxnorm = axisnorm.cross(Cynorm);
  // Cxnorm.normalize();
  // double radi;
  // Point pnti;
  // RevEngUtils::computeCylPosRadius(group, low, high,
  // 				   axisnorm, Cxnorm, Cynorm,
  // 				   pnti, radi);
  
  vector<Point> rotated;
  RevEngUtils::rotateToPlane(group, Cx, axis, pnt, rotated);
  // vector<Point> rotatedi;
  // RevEngUtils::rotateToPlane(group, Cxnorm, axisnorm, pnti, rotatedi);
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
  // of3 << "400 1 0 4 0 255 0 255" << std::endl;
  // of3 << rotatedi.size() << std::endl;
  // for (size_t kr=0; kr<rotatedi.size(); ++kr)
  //   of3 << rotatedi[kr] << std::endl;
  // of3 << "410 1 0 4 0 255 0 255" << std::endl;
  // of3 << "1" << std::endl;
  // of3 << pnti-0.5*len*axisnorm << " " << pnti+0.5*len*axisnorm << std::endl;

  shared_ptr<Line> line(new Line(pnt, axis));
  line->setParameterInterval(-len, len);
  shared_ptr<SplineCurve> line2;
  curveApprox(rotated, line, 2, 2, line2);
  vector<Point> der(2);
  line2->point(der, 0.0, 1);
  
  shared_ptr<Circle> circ(new Circle(rad, pnt, axis, Cy));
  //shared_ptr<Circle> circi(new Circle(radi, pnti, axisnorm, Cxnorm));
  vector<Point> projected;
  double maxdp, avdp;
  RevEngUtils::projectToPlane(group, axis, pnt, projected, maxdp, avdp);
  // vector<Point> projectedi;
  // double maxdpi, avdpi;
  // RevEngUtils::projectToPlane(group, axisnorm, pnti, projectedi, maxdpi, avdpi);
  ofp3 << "400 1 0 4 255 0 0 255" << std::endl;
  ofp3 << projected.size() << std::endl;
  for (size_t kr=0; kr<projected.size(); ++kr)
    ofp3 << projected[kr] << std::endl;
  // ofp3 << "400 1 0 4 0 255 0 255" << std::endl;
  // ofp3 << projectedi.size() << std::endl;
  // for (size_t kr=0; kr<projectedi.size(); ++kr)
  //   ofp3 << projectedi[kr] << std::endl;
  circ->writeStandardHeader(ofp3);
  circ->write(ofp3);
  // circi->writeStandardHeader(ofp3);
  // circi->write(ofp3);
  shared_ptr<SplineCurve> spl;
  Point xpos;
  vector<double> param;
  curveApprox(projected, tol, circ, param, spl, xpos);
  spl->writeStandardHeader(ofp3);
  spl->write(ofp3);
  double maxdc, avdc, maxds, avds;
  int num_inc, num_ins;
  RevEngUtils::distToCurve(projected, circ, tol, maxdc, avdc, num_inc);
  RevEngUtils::distToCurve(projected, spl, tol, maxds, avds, num_ins);

  if (xpos.dimension() == pnt.dimension())
    {
      Point vec = xpos - pnt;
      if (vec.length() > tol)
	{
	  vec.normalize();
	  Cy = vec;
	}
    }
  shared_ptr<Cylinder> cyl(new Cylinder(rad, pnt, axis, Cy));
  //cyl->setParamBoundsV(-len, len);

  // if (maxds < maxdc && num_ins > num_inc && 
  //     num_inc > 0.25*points.size() && rad < 2.0*len)
  if (maxds < maxdc && num_ins > num_inc && 
      num_ins > (int)points.size()/4 && rad < 2.0*len)
    {
      // Investigate point configuration
      //configSplit(points, param, cyl, spl, tol, configs);
      configSplit(points, param, cyl, spl, maxds, configs);
      std::ofstream ofconf("conf_groups.g2");
      for (size_t ki=0; ki<configs.size(); ++ki)
	{
	  ofconf << "400 1 0 0" << std::endl;
	  ofconf << configs[ki].size() << std::endl;
	  for (size_t kr=0; kr<configs[ki].size(); ++kr)
	    ofconf << configs[ki][kr]->getPoint() << std::endl;
	}
      int stop_out = 1;
    }
  
  if (configs.size() <= 1 && points.size() == group_points_.size() &&
      (num_ins >= num_inc && avds <= avdc) &&
      (num_ins > (int)group_points_.size()/2 && avds < tol) &&
      (!sweep_.get() || sweep_->type_ != 1 ||
       (sweep_->num_in_ < num_ins && sweep_->avdist_ > avds)))
    {
      // Possible linear sweep
      double len = bbox_.low().dist(bbox_.high());
      Point pt1 = pnt - len*axis; //startpt - len*axis; //pnt - len*axis;
      Point pt2 = pnt + len*axis; //startpt + len*axis; //pnt + len*axis;
      sweep_ = shared_ptr<SweepData>(new SweepData(1, spl, pt1, pt2, maxds, avds, num_ins));
    }


  return cyl;
}

//===========================================================================
bool RevEngRegion::extractLinearSweep(double tol, int min_pt, double angtol,
				      int prefer_elementary_,
				      std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
				      std::vector<HedgeSurface*>& prevsfs)
//===========================================================================
{
  // std::ofstream of("curr_region.g2");
  // writeRegionInfo(of);
  
  // std::ofstream of2("curr_normals.g2");
  // writeUnitSphereInfo(of2);

  if (!sweep_.get())
    return false;

  if (sweep_->type_ != 1)
    return false;
  
  bool found = false;
  shared_ptr<SplineSurface> surf;

  vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> > adj_elem, adj_elem_base, adj_primary;
  getAdjacentElemInfo(adj_elem, adj_elem_base, adj_primary);
  
  SweepSurfaceCreator sweep;
  shared_ptr<SplineCurve> along(new SplineCurve(sweep_->location_, sweep_->added_info_));
  Point mid = 0.5*(sweep_->location_ + sweep_->added_info_);
  surf = shared_ptr<SplineSurface>(sweep.linearSweptSurface(*along,
							  *sweep_->profile_, mid));
  if (!surf.get())
    return false;
  std::ofstream of3("sweep.g2");
  along->writeStandardHeader(of3);
  along->write(of3);
  sweep_->profile_->writeStandardHeader(of3);
  sweep_->profile_->write(of3);
  surf->writeStandardHeader(of3);
  surf->write(of3);
  of3 << "400 1 0 4 255 0 0 255" << std::endl;
  of3 << "1" << std::endl;
  of3 << sweep_->location_ << std::endl;

  // Check accuracy
  double maxd, avd; 
  int num2; 
  vector<RevEngPoint*> inpt, outpt;
  vector<pair<double, double> > dist_ang;
  vector<double> parvals;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
  			  surf, tol, maxd, avd, num2, inpt, outpt,
			  parvals, dist_ang, angtol);
  std::ofstream ofd("in_out_surf.g2");
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
  if (accuracyOK(min_pt, tol, num2, avd))
    {
     bool OK = true;
      if (associated_sf_.size() > 0)
	{
	  // Check with current approximating surface
	  double acc_fac1 = 1.25;
	  double acc_fac2 = 0.75;
	  if (prefer_elementary_ == ALWAYS_ELEM)
	    OK = false;
	  else if (prefer_elementary_ == PREFER_ELEM &&
		   ((double)num2 < acc_fac1*num_init ||
		    avd > acc_fac2*avd_init))
	    OK = false;
	  else if (prefer_elementary_ == BEST_ACCURACY &&
		   (num2 < num_init || avd > avd_init))
	    OK = false;
	}

      if (OK)
	{
	  found = true;
	  for (size_t kh=0; kh<group_points_.size(); ++kh)
	    {
	      group_points_[kh]->setPar(Vector2D(parvals[2*kh],parvals[2*kh+1]));
	      group_points_[kh]->setSurfaceDist(dist_ang[kh].first, dist_ang[kh].second);
	    }
	  setAccuracy(maxd, avd, num2);
	  std::cout << "Linear swept surface. N1: " << num << ", N2: " << num2 << ", max: " << maxd << ", av: " << avd << std::endl;

	  shared_ptr<HedgeSurface> hedge(new HedgeSurface(surf, this));
	  hedge->setLinearSweepInfo(sweep_->profile_, sweep_->location_, sweep_->added_info_);
	  for (size_t kh=0; kh<associated_sf_.size(); ++kh)
	    prevsfs.push_back(associated_sf_[kh]);
	  setHedge(hedge.get());
	  hedgesfs.push_back(hedge);
	}
    }

   return found;
 }

//===========================================================================
shared_ptr<SplineSurface>
RevEngRegion::computeLinearSwept(double tol, shared_ptr<SplineCurve>& profile,
				 Point& pt1, Point& pt2)
//===========================================================================
{
  // Axis by covariance matrix of normal vectors
  Point axis;
  Point Cx;
  Point Cy;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(group_points_.begin(), group_points_.end()));
  RevEngUtils::computeAxis(group, axis, Cx, Cy);

  // Point on axis
  Point low = bbox_.low();
  Point high = bbox_.high();
  double len = low.dist(high);
  double rad;
  Point pnt;
  RevEngUtils::computeCylPosRadius(group, low, high,
				   axis, Cx, Cy, pnt, rad);

  // Project the points onto the defined plane 
  vector<Point> projected;
  double maxdp, avdp;
  RevEngUtils::projectToPlane(group, axis, pnt, projected, maxdp, avdp);

  // Approximate the projected point cloud with a spline curve
  shared_ptr<Circle> circ(new Circle(rad, pnt, axis, Cx));
  Point xpos;
  vector<double> param;
  curveApprox(projected, tol, circ, param, profile, xpos);

  pt1 = pnt - len*axis;
  pt2 = pnt + len*axis;
  Point mid = 0.5*(pt1 + pt2);
  
  SweepSurfaceCreator sweep;
  shared_ptr<SplineCurve> along(new SplineCurve(pt1, pt2));
  shared_ptr<SplineSurface> swept_surf(sweep.linearSweptSurface(*along, *profile, mid));
				       
  return swept_surf;
}

//===========================================================================
bool RevEngRegion::extractSphere(double tol, int min_pt, double angtol,
				 double mean_edge_len,
				 int prefer_elementary_,
				 vector<shared_ptr<HedgeSurface> >& hedgesfs,
				 vector<HedgeSurface*>& prevsfs,
				 vector<vector<RevEngPoint*> >& out_groups)
//===========================================================================
{
  // std::ofstream of("curr_region.g2");
  // writeRegionInfo(of);
  
  // std::ofstream of2("curr_normals.g2");
  // writeUnitSphereInfo(of2);

  bool found = false;
  int min_nmb = 20;
  //double eps = 1.0e-6;
  if ((int)group_points_.size() < min_nmb)
    return false;
  
  vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> > adj_elem, adj_elem_base, adj_primary;
  getAdjacentElemInfo(adj_elem, adj_elem_base, adj_primary);
  
  shared_ptr<Sphere> sphere = computeSphere(group_points_);
  if (!sphere.get())
    return false;
  std::ofstream ofs("sph.g2");
  sphere->writeStandardHeader(ofs);
  sphere->write(ofs);
  
  // Check accuracy
  double maxd, avd; 
  int num2; 
  vector<RevEngPoint*> inpt, outpt;
  vector<pair<double, double> > dist_ang;
  vector<double> parvals;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
  			  sphere, tol, maxd, avd, num2, inpt, outpt,
			  parvals, dist_ang, angtol);
  std::ofstream ofd("in_out_sphere.g2");
  ofd << "400 1 0 4 155 50 50 255" << std::endl;
  ofd << inpt.size() << std::endl;
  for (size_t kr=0; kr<inpt.size(); ++kr)
    ofd << inpt[kr]->getPoint() << std::endl;
  ofd << "400 1 0 4 50 155 50 255" << std::endl;
  ofd << outpt.size() << std::endl;
  for (size_t kr=0; kr<outpt.size(); ++kr)
    ofd << outpt[kr]->getPoint() << std::endl;

  if (inpt.size() > group_points_.size()/2 && (int)inpt.size() > min_nmb)
    {
      shared_ptr<Sphere> sphere_in = computeSphere(inpt);
      vector<RevEngPoint*> inpt_in, outpt_in; //, inpt2, outpt2;
      vector<pair<double, double> > dist_ang_in;
      vector<double> parvals_in;
      double maxd_in, avd_in;
      int num2_in;
      RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			      sphere_in, tol, maxd_in, avd_in, num2_in,
			      inpt_in, outpt_in, parvals_in, dist_ang_in, angtol);
  
      std::ofstream ofd2("in_out_sphere2.g2");
      ofd2 << "400 1 0 4 155 50 50 255" << std::endl;
      ofd2 << inpt_in.size() << std::endl;
      for (size_t kr=0; kr<inpt_in.size(); ++kr)
	ofd2 << inpt_in[kr]->getPoint() << std::endl;
      ofd2 << "400 1 0 4 50 155 50 255" << std::endl;
      ofd2 << outpt_in.size() << std::endl;
      for (size_t kr=0; kr<outpt_in.size(); ++kr)
	ofd2 << outpt_in[kr]->getPoint() << std::endl;

      if (num2_in > num2 || avd_in < avd)
	{
	  std::swap(sphere, sphere_in);
	  std::swap(num2, num2_in);
	  std::swap(avd, avd_in);
	  std::swap(maxd, maxd_in);
	  std::swap(parvals, parvals_in);
	  std::swap(dist_ang, dist_ang_in);
	  // std::cout << "Sphere swap" << std::endl;
	  // std::cout << group_points_.size() << " " << num2_in;
	  // std::cout << " " << num2 << " " << avd_in << " " << avd;
	  // std::cout << " " << maxd_in << " " << maxd << std::endl;
	}
    }

  int num = (int)group_points_.size();
  double maxd_init, avd_init;
  int num_init;
  getAccuracy(maxd_init, avd_init, num_init);
  if (accuracyOK(min_pt, tol, num2, avd))
    {
     bool OK = true;
      if (associated_sf_.size() > 0)
	{
	  // Check with current approximating surface
	  double acc_fac1 = 1.25;
	  double acc_fac2 = 0.75;
	  if (prefer_elementary_ == ALWAYS_ELEM)
	    OK = false;
	  else if (prefer_elementary_ == PREFER_ELEM &&
		   ((double)num2 < acc_fac1*num_init ||
		    avd > acc_fac2*avd_init))
	    OK = false;
	  else if (prefer_elementary_ == BEST_ACCURACY &&
		   (num2 < num_init || avd > avd_init))
	    OK = false;
	}

      if (OK)
	{
	  found = true;
	  for (size_t kh=0; kh<group_points_.size(); ++kh)
	    {
	      group_points_[kh]->setPar(Vector2D(parvals[2*kh],parvals[2*kh+1]));
	      group_points_[kh]->setSurfaceDist(dist_ang[kh].first, dist_ang[kh].second);
	    }
	  setAccuracy(maxd, avd, num2);
	  std::cout << "Sphere. N1: " << num << ", N2: " << num2 << ", max: " << maxd << ", av: " << avd << std::endl;
	  shared_ptr<HedgeSurface> hedge(new HedgeSurface(sphere, this));
	  for (size_t kh=0; kh<associated_sf_.size(); ++kh)
	    prevsfs.push_back(associated_sf_[kh]);
	  setHedge(hedge.get());
	  hedgesfs.push_back(hedge);
	}
    }
  if (!primary_.get() ||
      (num2 >= num_in_primary_ && avd < avdist_primary_))
    setPrimarySf(sphere, maxd, avd, num2);
  int stop_break0 = 1;
  return found;
}

//===========================================================================
shared_ptr<Sphere> RevEngRegion::computeSphere(vector<RevEngPoint*>& points)
//===========================================================================
{
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(points.begin(), points.end()));
  Point centre;
  double radius;
  try {
    RevEngUtils::computeSphereProp(group, centre, radius);
  }
  catch (...)
    {
      shared_ptr<Sphere> dummy;
      return dummy;
    }

  Point y_axis = normalcone_.centre();
  double lambda[3];
  Point eigen1, eigen2, eigen3;
  getPCA(lambda, eigen1, eigen2, eigen3);
  Point z_axis = eigen2;
  Point x_axis = y_axis.cross(z_axis);
  
  shared_ptr<Sphere> sph(new Sphere(radius, centre, z_axis, x_axis));
  //cyl->setParamBoundsV(-len, len);
  return sph;
}

//===========================================================================
bool RevEngRegion::extractCone(double tol, int min_pt, double angtol,
			       double mean_edge_len,
			       int prefer_elementary_,
			       vector<shared_ptr<HedgeSurface> >& hedgesfs,
			       vector<HedgeSurface*>& prevsfs,
			       vector<vector<RevEngPoint*> >& out_groups)
//===========================================================================
{
  // std::ofstream of("curr_region.g2");
  // writeRegionInfo(of);
  
  // std::ofstream of2("curr_normals.g2");
  // writeUnitSphereInfo(of2);

  bool found = false;
  int min_nmb = 20;
  //double eps = 1.0e-6;
  double minang = 0.05;
  if ((int)group_points_.size() < min_nmb)
    return false;

  vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> > adj_elem, adj_elem_base, adj_primary;
  getAdjacentElemInfo(adj_elem, adj_elem_base, adj_primary);
  
  Point apex;
  shared_ptr<Cone> cone = computeCone(group_points_, apex);
  std::ofstream ofs("conesf.g2");
  cone->writeStandardHeader(ofs);
  cone->write(ofs);
  
  // Check accuracy
  double maxd, avd; 
  int num2; 
  vector<RevEngPoint*> inpt, outpt;
  vector<pair<double, double> > dist_ang;
  vector<double> parvals;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
  			  cone, tol, maxd, avd, num2, inpt, outpt,
			  parvals, dist_ang, angtol);
  std::ofstream ofd("in_out_cone.g2");
  ofd << "400 1 0 4 155 50 50 255" << std::endl;
  ofd << inpt.size() << std::endl;
  for (size_t kr=0; kr<inpt.size(); ++kr)
    ofd << inpt[kr]->getPoint() << std::endl;
  ofd << "400 1 0 4 50 155 50 255" << std::endl;
  ofd << outpt.size() << std::endl;
  for (size_t kr=0; kr<outpt.size(); ++kr)
    ofd << outpt[kr]->getPoint() << std::endl;

  if (inpt.size() > group_points_.size()/2 && (int)inpt.size() > min_nmb)
    {
      Point apex2;
      shared_ptr<Cone> cone_in = computeCone(inpt, apex2);
      vector<RevEngPoint*> inpt_in, outpt_in; //, inpt2, outpt2;
      vector<pair<double, double> > dist_ang_in;
      vector<double> parvals_in;
      double maxd_in, avd_in;
      int num2_in;
      RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			      cone_in, tol, maxd_in, avd_in, num2_in,
			      inpt_in, outpt_in, parvals_in, dist_ang_in, angtol);
  
      std::ofstream ofd2("in_out_cone2.g2");
      ofd2 << "400 1 0 4 155 50 50 255" << std::endl;
      ofd2 << inpt_in.size() << std::endl;
      for (size_t kr=0; kr<inpt_in.size(); ++kr)
	ofd2 << inpt_in[kr]->getPoint() << std::endl;
      ofd2 << "400 1 0 4 50 155 50 255" << std::endl;
      ofd2 << outpt_in.size() << std::endl;
      for (size_t kr=0; kr<outpt_in.size(); ++kr)
	ofd2 << outpt_in[kr]->getPoint() << std::endl;

      if (num2_in > num2 || avd_in < avd)
	{
	  std::swap(cone, cone_in);
	  std::swap(num2, num2_in);
	  std::swap(avd, avd_in);
	  std::swap(maxd, maxd_in);
	  std::swap(parvals, parvals_in);
	  std::swap(dist_ang, dist_ang_in);
	  std::swap(apex, apex2);
	  // std::cout << "Cone swap" << std::endl;
	  // std::cout << group_points_.size() << " " << num2_in;
	  // std::cout << " " << num2 << " " << avd_in << " " << avd;
	  // std::cout << " " << maxd_in << " " << maxd << std::endl;
	}
    }

  int num = (int)group_points_.size();
  double maxd_init, avd_init;
  int num_init;
  getAccuracy(maxd_init, avd_init, num_init);
  double cone_ang = cone->getConeAngle();
  if (0.5*M_PI-fabs(cone_ang) > minang && (!bbox_.containsPoint(apex, tol)) &&
      accuracyOK(min_pt, tol, num2, avd))
    {
     bool OK = true;
      if (associated_sf_.size() > 0)
	{
	  // Check with current approximating surface
	  double acc_fac1 = 1.25;
	  double acc_fac2 = 0.75;
	  if (prefer_elementary_ == ALWAYS_ELEM)
	    OK = false;
	  else if (prefer_elementary_ == PREFER_ELEM &&
		   ((double)num2 < acc_fac1*num_init ||
		    avd > acc_fac2*avd_init))
	    OK = false;
	  else if (prefer_elementary_ == BEST_ACCURACY &&
		   (num2 < num_init || avd > avd_init))
	    OK = false;
	}

      if (OK)
	{
	  found = true;
	  for (size_t kh=0; kh<group_points_.size(); ++kh)
	    {
	      group_points_[kh]->setPar(Vector2D(parvals[2*kh],parvals[2*kh+1]));
	      group_points_[kh]->setSurfaceDist(dist_ang[kh].first, dist_ang[kh].second);
	    }
	  setAccuracy(maxd, avd, num2);
	  std::cout << "Cone. N1: " << num << ", N2: " << num2 << ", max: " << maxd << ", av: " << avd << std::endl;
	  shared_ptr<HedgeSurface> hedge(new HedgeSurface(cone, this));
	  for (size_t kh=0; kh<associated_sf_.size(); ++kh)
	    prevsfs.push_back(associated_sf_[kh]);
	  setHedge(hedge.get());
	  hedgesfs.push_back(hedge);
	}
    }
  if (!primary_.get() ||
      (num2 >= num_in_primary_ && avd < avdist_primary_))
    setPrimarySf(cone, maxd, avd, num2);
  int stop_break0 = 1;
  return found;
}


//===========================================================================
shared_ptr<Cone> RevEngRegion::computeCone(vector<RevEngPoint*>& points, Point& apex)
//===========================================================================
{
  Point low = bbox_.low();
  Point high = bbox_.high();
  double len = low.dist(high);

  Point mid(0.0, 0.0, 0.0);
  double wgt = 1.0/(double)points.size();
  for (size_t kr=0; kr<points.size(); ++kr)
    {
      Vector3D xyz = points[kr]->getPoint();
      Point pnt(xyz[0], xyz[1], xyz[2]);
      mid += wgt*pnt;
    }
  Point axis;
  Point Cx;
  Point Cy;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(points.begin(), points.end()));
  RevEngUtils::coneAxis(group, axis, Cx, Cy);

  //Point apex;
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
  //cone->setParamBoundsV(-0.2*len, 0.2*len);
  
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
bool RevEngRegion::contextTorus(Point mainaxis[3],
				double tol, int min_pt, double angtol,
				double mean_edge_len,
				int prefer_elementary_,
				vector<shared_ptr<HedgeSurface> >& hedgesfs,
				vector<HedgeSurface*>& prevsfs,
				vector<vector<RevEngPoint*> >& out_groups)
//===========================================================================
{
  bool found = false;
  int min_nmb = 20;
  if ((int)group_points_.size() < min_nmb)
    return false;
  
  // std::ofstream of("curr_region.g2");
  // writeRegionInfo(of);

  vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> > adj_elem, adj_elem_base, adj_primary;
  getAdjacentElemInfo(adj_elem, adj_elem_base, adj_primary);

  double angtol2 = 0.1;
  Point pos, axis, Cx;
  vector<size_t> adj_ix;
  double R1, R2;
  double cyl_dom[4];
  bool outer;
  int plane_ix, cyl_ix;
  bool foundaxis = analyseTorusContext(adj_elem, tol, angtol2, adj_ix, plane_ix, cyl_ix,
				       pos, axis, Cx, R1, R2, cyl_dom, outer);
  if (!foundaxis)
    return false;

  // bool update = false;
  // if (update)
  //   {
  // // Move points to adjacent if appropriate. First sort adjacent according to number of points
  // for (size_t ki=0; ki<adj_ix.size(); ++ki)
  //   for (size_t kj=ki+1; kj<adj_ix.size(); ++kj)
  //     if (adj_elem[adj_ix[kj]].second->numPoints() > adj_elem[adj_ix[ki]].second->numPoints())
  // 	std::swap(adj_ix[kj], adj_ix[ki]);

  // for (size_t ki=0; ki<adj_ix.size(); ++ki)
  //   {
  //     RevEngRegion* next_reg = adj_elem[adj_ix[ki]].second;

  //     // Get seed points
  //     vector<RevEngPoint*> adj_pts = extractNextToAdjacent(next_reg);

  //     // Grow adjacent
  //     next_reg->growFromNeighbour(adj_pts, tol, this);
  //   }
  // updateInfo();
  
  // std::ofstream of("updated_region.g2");
  // writeRegionInfo(of);
  //   }

  double eps = 1.0e-2;
  Point pos2 = pos - R2*axis;
  Point Cy = axis.cross(Cx);
  shared_ptr<Torus> tor(new Torus(R1, R2, pos2, axis, Cx));
  
#ifdef DEBUG_TORUSCONTEXT
  // Rotate point cloud
  Point low = bbox_.low();
  Point high = bbox_.high();
  double len = low.dist(high);

  vector<Point> rotated;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(group_points_.begin(), group_points_.end()));
  RevEngUtils::rotateToPlane(group, Cx, axis, pos, rotated);
  std::ofstream of3("rotated_pts_tor.g2");
  of3 << "400 1 0 4 255 0 0 255" << std::endl;
  of3 << rotated.size() << std::endl;
  for (size_t kr=0; kr<rotated.size(); ++kr)
    of3 << rotated[kr] << std::endl;
  of3 << "410 1 0 4 0 0 255 255" << std::endl;
  of3 << "1" << std::endl;
  of3 << pos-0.5*len*axis << " " << pos+0.5*len*axis << std::endl;
  of3 << "410 1 0 4 0 255 0 255" << std::endl;
  of3 << "1" << std::endl;
  of3 << pos-0.5*len*Cx << " " << pos+0.5*len*Cx << std::endl;

  // Point cpos;
  // double crad;
  // RevEngUtils::computeCircPosRadius(rotated, Cy, Cx, axis, cpos, crad);
  // shared_ptr<Circle> circ(new Circle(crad, cpos, Cy, Cx));
  // circ->writeStandardHeader(of3);
  // circ->write(of3);

  Point cpos2 = pos + R1*Cx - R2*axis;
  shared_ptr<Circle> circ2(new Circle(R2, cpos2, Cy, Cx));
  circ2->writeStandardHeader(of3);
  circ2->write(of3);

  std::ofstream oft("torus_context.g2");
  tor->writeStandardHeader(oft);
  tor->write(oft);
#endif

  // Move points as appropriate
  // First check accuracy with respect to torus
  vector<RevEngPoint*> inpt, outpt;
  vector<pair<double, double> > dist_ang;
  double maxd, avd;
  int num_in;
  vector<double> parvals;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			  tor, tol, maxd, avd, num_in,
			  inpt, outpt, parvals, dist_ang,
			  angtol);
  
  // Extract planar and cylindrical points from torus
  vector<RevEngPoint*> planar;
  vector<RevEngPoint*> cyl_pts;
  vector<RevEngPoint*> remaining;
  shared_ptr<ElementarySurface> cyl = adj_elem[cyl_ix].first;
  shared_ptr<ElementarySurface> plane = adj_elem[plane_ix].first;
  Point axis2 = cyl->direction();
  double pihalf = 0.5*M_PI;
  double tor_dom[4];   // Torus domain
  tor_dom[0] = tor_dom[2] = std::numeric_limits<double>::max();
  tor_dom[1] = tor_dom[3] = std::numeric_limits<double>::lowest();
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Vector3D xyz = group_points_[ki]->getPoint();
      Point ptpos(xyz[0], xyz[1], xyz[2]);
      Point norm = group_points_[ki]->getMongeNormal();
      double dist = pos.dist(ptpos);
      double ang1 = axis.angle(norm);
      double ang2 = axis2.angle(norm);
      if (((dist < R1 && outer) || (dist > R1 && (!outer))) && ang1 <= angtol2)
	planar.push_back(group_points_[ki]);
      else if (fabs(pihalf-ang2) <= dist_ang[ki].second || dist_ang[ki].first > tol)
	{
	  double upar, vpar, tdist;
	  Point close;
	  cyl->closestPoint(ptpos, upar, vpar, close, tdist, eps);
	  if (upar >= cyl_dom[0] && upar <= cyl_dom[1] &&
	      vpar >= cyl_dom[2] && vpar <= cyl_dom[3] && tdist <= dist_ang[ki].first)
	    cyl_pts.push_back(group_points_[ki]);
	  else
	    {
	      remaining.push_back(group_points_[ki]);
	      tor_dom[0] = std::min(tor_dom[0], parvals[2*ki]);
	      tor_dom[1] = std::max(tor_dom[1], parvals[2*ki]);
	      tor_dom[2] = std::min(tor_dom[2], parvals[2*ki+1]);
	      tor_dom[3] = std::max(tor_dom[3], parvals[2*ki+1]);
	    }
	}
      else
	{
	  remaining.push_back(group_points_[ki]);
	  tor_dom[0] = std::min(tor_dom[0], parvals[2*ki]);
	  tor_dom[1] = std::max(tor_dom[1], parvals[2*ki]);
	  tor_dom[2] = std::min(tor_dom[2], parvals[2*ki+1]);
	  tor_dom[3] = std::max(tor_dom[3], parvals[2*ki+1]);
	}
    }

  // Plane to torus
  vector<RevEngPoint*>  plan2tor;
  for (auto it=adj_elem[plane_ix].second->pointsBegin();
       it!=adj_elem[plane_ix].second->pointsEnd(); ++it)
    {
      Vector3D xyz = (*it)->getPoint();
      Point ptpos(xyz[0], xyz[1], xyz[2]);
      Point norm = (*it)->getMongeNormal();
      double dist = pos.dist(ptpos);
      double ang = axis.angle(norm);
      ang = std::min(ang, M_PI-ang);
      if ((dist > R1 && outer) || (dist < R1 && (!outer)))
	{
	  double upar, vpar, tdist;
	  Point close;
	  tor->closestPoint(ptpos, upar, vpar, close, tdist, eps);
	  Point tnorm;
	  tor->normal(tnorm, upar, vpar);
	  double ang2 = tnorm.angle(norm);
	  if (tdist < (*it)->getSurfaceDist() && ang2 <= angtol2 &&
	      upar >= tor_dom[0] && upar <= tor_dom[1])
	    plan2tor.push_back(*it);
	}
    }
  
  // Cylinder to torus
  vector<RevEngPoint*>  cyl2tor;
  for (auto it=adj_elem[cyl_ix].second->pointsBegin();
       it!=adj_elem[cyl_ix].second->pointsEnd(); ++it)
    {
      Vector3D xyz = (*it)->getPoint();
      Vector2D uv = (*it)->getPar();
      Point ptpos(xyz[0], xyz[1], xyz[2]);
      Point norm = (*it)->getMongeNormal();
      if (uv[0] < cyl_dom[0] || uv[0] > cyl_dom[1] ||
	  uv[1] < cyl_dom[2] || uv[1] > cyl_dom[3])
	{
	  double upar, vpar, tdist;
	  Point close;
	  tor->closestPoint(ptpos, upar, vpar, close, tdist, eps);
	  Point tnorm;
	  tor->normal(tnorm, upar, vpar);
	  double ang2 = tnorm.angle(norm);
	  if (tdist < (*it)->getSurfaceDist() && ang2 <= angtol2)
	    cyl2tor.push_back(*it);
	}
    }
  

  
#ifdef DEBUG_TORUSCONTEXT
  if (planar.size() > 0)
    {
      std::ofstream ofp("planar_move.g2");
      ofp << "400 1 0 4 255 0 0 255" << std::endl;
      ofp << planar.size() << std::endl;
      for (size_t kr=0; kr<planar.size(); ++kr)
	ofp << planar[kr]->getPoint() << std::endl;
    }

  if (plan2tor.size() > 0)
    {
      std::ofstream oftp("plan_tor_move.g2");
      oftp << "400 1 0 4 0 255 0 255" << std::endl;
      oftp << plan2tor.size() << std::endl;
      for (size_t kr=0; kr<plan2tor.size(); ++kr)
	oftp << plan2tor[kr]->getPoint() << std::endl;
    }
  
  if (cyl_pts.size() > 0)
    {
      std::ofstream ofc("cyl_move.g2");
      ofc << "400 1 0 4 255 0 0 255" << std::endl;
      ofc << cyl_pts.size() << std::endl;
      for (size_t kr=0; kr<cyl_pts.size(); ++kr)
	ofc << cyl_pts[kr]->getPoint() << std::endl;
    }

  if (cyl2tor.size() > 0)
    {
      std::ofstream oftc("cyl_tor_move.g2");
      oftc << "400 1 0 4 0 255 0 255" << std::endl;
      oftc << cyl2tor.size() << std::endl;
      for (size_t kr=0; kr<cyl2tor.size(); ++kr)
	oftc << cyl2tor[kr]->getPoint() << std::endl;
    }
#endif

  // Recompute accuracy
  vector<RevEngPoint*> inpt2, outpt2;
  vector<pair<double, double> > dist_ang2;
  double maxd2, avd2;
  int num_in2;
  vector<double> parvals2;
  RevEngUtils::distToSurf(remaining.begin(), remaining.end(),
			  tor, tol, maxd2, avd2, num_in2,
			  inpt2, outpt2, parvals2, dist_ang2,
			  angtol);
  std::ofstream ofd2("in_out_tor_context.g2");
  ofd2 << "400 1 0 4 155 50 50 255" << std::endl;
  ofd2 << inpt2.size() << std::endl;
  for (size_t kr=0; kr<inpt2.size(); ++kr)
    ofd2 << inpt2[kr]->getPoint() << std::endl;
  ofd2 << "400 1 0 4 50 155 50 255" << std::endl;
  ofd2 << outpt2.size() << std::endl;
  for (size_t kr=0; kr<outpt2.size(); ++kr)
    ofd2 << outpt2[kr]->getPoint() << std::endl;

  if (num_in2 > min_pt && num_in2 > (int)remaining.size()/2)
    {
      // Check for deviant points at the boundary
      shared_ptr<RevEngRegion> reg(new RevEngRegion(classification_type_,
						    edge_class_type_, remaining));
      
      vector<RevEngPoint*> dist_points;
      reg->identifyDistPoints(dist_ang2, tol, maxd2, avd2, dist_points);
      reg->extractSpesPoints(dist_points, out_groups, true);

      vector<RevEngPoint*> remaining2;
      remaining2 = reg->getPoints();
      for (size_t ki=0; ki<remaining2.size(); ++ki)
	remaining2[ki]->setRegion(this);

      std::swap(remaining, remaining2);
    }

  if (plan2tor.size() > 0)
    remaining.insert(remaining.end(), plan2tor.begin(), plan2tor.end());
  if (cyl2tor.size() > 0)
    remaining. insert(remaining.end(), cyl2tor.begin(), cyl2tor.end());
  
  // Recompute accuracy
  vector<RevEngPoint*> inpt3, outpt3;
  vector<pair<double, double> > dist_ang3;
  double maxd3, avd3;
  int num_in3;
  vector<double> parvals3;
  RevEngUtils::distToSurf(remaining.begin(), remaining.end(),
			  tor, tol, maxd3, avd3, num_in3,
			  inpt3, outpt3, parvals3, dist_ang3,
			  angtol);
  std::ofstream ofd3("in_out_tor_context3.g2");
  ofd3 << "400 1 0 4 155 50 50 255" << std::endl;
  ofd3 << inpt3.size() << std::endl;
  for (size_t kr=0; kr<inpt3.size(); ++kr)
    ofd3 << inpt3[kr]->getPoint() << std::endl;
  ofd3 << "400 1 0 4 50 155 50 255" << std::endl;
  ofd3 << outpt3.size() << std::endl;
  for (size_t kr=0; kr<outpt3.size(); ++kr)
    ofd3 << outpt3[kr]->getPoint() << std::endl;

  // Move points
  std::swap(group_points_, remaining);
  bool OKsurf = accuracyOK(min_pt, tol, num_in3, avd3);
  if (OKsurf)
    {
      for (size_t ki=0; ki<group_points_.size(); ++ki)
	{
	  group_points_[ki]->setPar(Vector2D(parvals3[2*ki],parvals3[2*ki+1]));
	  group_points_[ki]->setSurfaceDist(dist_ang3[ki].first, dist_ang3[ki].second);
	}
      setAccuracy(maxd3, avd3, num_in3);
      shared_ptr<HedgeSurface> hedge(new HedgeSurface(tor, this));
      setHedge(hedge.get());
      hedgesfs.push_back(hedge);
    }
  updateInfo();

  if (planar.size() > 0)
    {
      // Parameterize
      vector<RevEngPoint*> inpt4, outpt4;
      vector<pair<double, double> > dist_ang4;
      double maxd4, avd4;
      int num_in4;
      vector<double> parvals4;
      RevEngUtils::distToSurf(planar.begin(), planar.end(),
			      plane, tol, maxd4, avd4, num_in4,
			      inpt4, outpt4, parvals4, dist_ang4,
			      angtol);
       for (size_t ki=0; ki<planar.size(); ++ki)
	{
	  planar[ki]->setPar(Vector2D(parvals4[2*ki],parvals4[2*ki+1]));
	  planar[ki]->setSurfaceDist(dist_ang4[ki].first, dist_ang4[ki].second);
	  adj_elem[plane_ix].second->addPoint(planar[ki]);
	}
       adj_elem[plane_ix].second->updateInfo();
    }
  
  if (cyl_pts.size() > 0)
    {
      // Parameterize
      vector<RevEngPoint*> inpt4, outpt4;
      vector<pair<double, double> > dist_ang4;
      double maxd4, avd4;
      int num_in4;
      vector<double> parvals4;
      RevEngUtils::distToSurf(cyl_pts.begin(), cyl_pts.end(),
			      plane, tol, maxd4, avd4, num_in4,
			      inpt4, outpt4, parvals4, dist_ang4,
			      angtol);
       for (size_t ki=0; ki<cyl_pts.size(); ++ki)
	{
	  cyl_pts[ki]->setPar(Vector2D(parvals4[2*ki],parvals4[2*ki+1]));
	  cyl_pts[ki]->setSurfaceDist(dist_ang4[ki].first, dist_ang4[ki].second);
	  adj_elem[cyl_ix].second->addPoint(cyl_pts[ki]);
	}
       adj_elem[cyl_ix].second->updateInfo();
    }
  
  return OKsurf;
 }

//===========================================================================
void RevEngRegion::growFromNeighbour(vector<RevEngPoint*>& seed, double tol,
				     RevEngRegion *neighbour)
//===========================================================================
{
  if (!hasSurface())
    return;   // No growing is possible

  shared_ptr<ParamSurface> surf = getSurface(0)->surface();
  shared_ptr<ElementarySurface> elem =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf);
  Point axis;
  if (elem.get())
    axis = elem->direction();
  double axis_ang = 0.0;
  if (surf->instanceType() == Class_Cylinder)
    axis_ang = 0.5*M_PI;   // To be extended as appropriate

  std::ofstream ofs1("seed1.g2");
  ofs1 << "400 1 0 4 0 200 55 255" << std::endl;
  ofs1 << seed.size() << std::endl;
  for (size_t ki=0; ki<seed.size(); ++ki)
    {
      Vector3D xyz = seed[ki]->getPoint();
      ofs1 << xyz << std::endl;
    }
  
  // Check initial points
  double eps = 1.0e-6;
  double tol2 = std::max(tol, 0.75*maxdist_);
  double tol3 = 0.5*tol;
  double angtol = 0.1;
  vector<RevEngPoint*> next_pts;
  double upar, vpar, dist, ang;
  Point close;
  for (size_t ki=0; ki<seed.size(); ++ki)
    {
      seed[ki]->setVisited();
      Vector3D xyz = seed[ki]->getPoint();
      Point norm = seed[ki]->getMongeNormal();
      surf->closestPoint(Point(xyz[0],xyz[1],xyz[2]), upar, vpar, close, dist, eps);
      if (!elem.get())
	surf->normal(axis, upar, vpar);
      ang = axis.angle(norm);
      if (dist <= tol3 || (dist <= tol2 && fabs(ang-axis_ang) <= angtol))
	next_pts.push_back(seed[ki]);
    }

  std::ofstream ofs2("seed2.g2");
  ofs2 << "400 1 0 4 0 200 55 255" << std::endl;
  ofs2 << next_pts.size() << std::endl;
  for (size_t ki=0; ki<next_pts.size(); ++ki)
    {
      Vector3D xyz = next_pts[ki]->getPoint();
      ofs2 << xyz << std::endl;
    }
  
  // Grow
  for (size_t ki=0; ki<next_pts.size(); ++ki)
    {
      vector<ftSamplePoint*> next2 = next_pts[ki]->getNeighbours();
      for (size_t kj=0; kj<next2.size(); ++kj)
	{
	  RevEngPoint *curr = dynamic_cast<RevEngPoint*>(next2[kj]);
	  if (curr->visited())
	    continue;
	  RevEngRegion *adj_reg = curr->region();
	  if (adj_reg != neighbour)
	    continue;
	  curr->setVisited();
	  Vector3D xyz = curr->getPoint();
	  Point norm = curr->getMongeNormal();
	  surf->closestPoint(Point(xyz[0],xyz[1],xyz[2]), upar, vpar, close, dist, eps);
	  if (!elem.get())
	    surf->normal(axis, upar, vpar);
	  ang = axis.angle(norm);
	  if (dist <= tol3 || (dist <= tol2 && fabs(ang-axis_ang) <= angtol))
	    next_pts.push_back(curr);
	}
    }
  for (int ka=0; ka<neighbour->numPoints(); ++ka)
    neighbour->getPoint(ka)->unsetVisited();

  if (surf->instanceType() == Class_Cylinder)
    {
      double dom[4];
      Vector2D uv = group_points_[0]->getPar();
      dom[0] = dom[1] = uv[0];
      dom[2] = dom[3] = uv[1];
      double dist, ang, avang;
      double fac = 1.0/(double)group_points_.size();
      group_points_[0]->getSurfaceDist(dist, ang);
      avang = fac*ang;
      for (size_t ki=1; ki<group_points_.size(); ++ki)
	{
	  Vector2D uv = group_points_[ki]->getPar();
	  group_points_[ki]->getSurfaceDist(dist, ang);
	  dom[0] = std::min(dom[0], uv[0]);
	  dom[1] = std::max(dom[1], uv[0]);
	  dom[2] = std::min(dom[2], uv[1]);
	  dom[3] = std::max(dom[3], uv[1]);
	  avang += fac*ang;
	}
      
      vector<RevEngPoint*> adjpts;
      vector<double> par_and_dist;
      double avd, ava;
      int nn;
      getAdjInsideDist(surf, dom, tol, neighbour, avd, ava, nn, adjpts, par_and_dist);

      std::ofstream of2("in_cyl_pts.g2");
      of2 << "400 1 0 4 75 75 75 255" << std::endl;
      of2 << adjpts.size() << std::endl;
      for (size_t kh=0; kh<adjpts.size(); ++kh)
	of2 << adjpts[kh]->getPoint() << std::endl;
      int stop_cyl = 1;
    }
  
  // Move points
  for (size_t ki=0; ki<next_pts.size(); ++ki)
    {
      addPoint(next_pts[ki]);
      neighbour->removePoint(next_pts[ki]);
    }

  std::ofstream of("updated_region_adj.g2");
  writeRegionInfo(of);

  // Update this region
  if (next_pts.size() > 0)
    checkReplaceSurf(tol, true);
  
}

//===========================================================================
vector<RevEngPoint*> RevEngRegion::extractNextToAdjacent(RevEngRegion* reg)
//===========================================================================
{
  vector<RevEngPoint*> result;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      vector<ftSamplePoint*> next = group_points_[ki]->getNeighbours();
      for (size_t kj=0; kj<next.size(); ++kj)
	{
	  RevEngPoint *curr = dynamic_cast<RevEngPoint*>(next[kj]);
	  RevEngRegion *adj_reg = curr->region();
	  if (adj_reg == reg)
	    {
	      result.push_back(group_points_[ki]);
	      break;
	    }
	}
    }
  return result;
}

//===========================================================================
bool
RevEngRegion::analyseTorusContext(vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> >& adj,
				  double tol, double angtol, vector<size_t>& adjacent_ix,
				  int &plane_ix, int& cyl_ix,
				  Point& pos, Point& axis, Point& Cx, double& R1,
				  double& R2, double dom[4], bool& outer)
//===========================================================================
{
  for (size_t ki=0; ki<adj.size(); ++ki)
    for (size_t kj=ki+1; kj<adj.size(); ++kj)
      if (adj[kj].second->numPoints() > adj[ki].second->numPoints())
	std::swap(adj[ki], adj[kj]);
  
  bool plane=false, cyl=false;
  plane_ix=-1;
  cyl_ix=-1;
  vector<Point> dir;
  double radius;
  for (size_t ki=0; ki<adj.size(); ++ki)
    {
      shared_ptr<ElementarySurface> elem = adj[ki].first;
      if (!(elem->instanceType() == Class_Plane || elem->instanceType() == Class_Cylinder))
	continue;

      Point curr_dir = elem->direction();
      dir.push_back(curr_dir);
      adjacent_ix.push_back(ki);
      if (elem->instanceType() == Class_Plane)
	{
	  plane = true;
	  plane_ix = (int)ki;
	}
      else
	{
	  cyl = true;
	  pos = elem->location();
	  Cx = elem->direction2();
	  radius = elem->radius(0.0, 0.0);
	  cyl_ix = (int)ki;
	}
	  
      for (size_t kj=ki+1; kj<adj.size(); ++kj)
	{
	  if (!(adj[kj].first->instanceType() == Class_Plane ||
		adj[kj].first->instanceType() == Class_Cylinder))
	    continue;
	  curr_dir = adj[kj].first->direction();
	  double ang = dir[0].angle(curr_dir);
	  ang = std::min(ang, M_PI-ang);
	  if (ang < angtol)
	    {
	      dir.push_back(curr_dir);
	      adjacent_ix.push_back(kj);
	      if (adj[kj].first->instanceType() == Class_Plane)
		{
		  plane = true;
		  if (plane_ix < 0)
		    plane_ix = (int)kj;
		}
	      else
		{
		  if (!cyl)
		    {
		      pos = adj[kj].first->location();
		      Cx = adj[kj].first->direction2();
		      radius = adj[kj].first->radius(0.0, 0.0);
		      cyl_ix = (int)kj;
		    }
		  cyl = true;
		}
	    }
	}

      if (plane && cyl)
	break;
      else
	{
	  plane = false;
	  cyl = false;
	  dir.clear();
	  adjacent_ix.clear();
	  plane_ix = cyl_ix = -1;
	}
    }

  if (plane && cyl)
    {
      // Compute axis
      for (size_t ki=1; ki<dir.size(); ++ki)
	if (dir[0]*dir[ki] < 0.0)
	  dir[ki] *= -1;

      // axis = Point(0.0, 0.0, 0.0);
      // double frac = 1.0/(double)dir.size();
      // for (size_t ki=0; ki<dir.size(); ++ki)
      // 	axis += frac*dir[ki];
      outer = true;
      axis = adj[plane_ix].first->direction();
      Point axis2 = adj[cyl_ix].first->direction();

#ifdef DEBUG_TORUSCONTEXT
      std::ofstream of("adj_elem.g2");
      for (size_t ki=0; ki<adjacent_ix.size(); ++ki)
	{
	  adj[adjacent_ix[ki]].first->writeStandardHeader(of);
	  adj[adjacent_ix[ki]].first->write(of);
	  adj[adjacent_ix[ki]].second->writeRegionInfo(of);
	}
#endif
      
#ifdef DEBUG_TORUSCONTEXT
      // Check accuracy of alternative plane
      Point loc = adj[cyl_ix].first->location();
      Point mid(0.0, 0.0, 0.0);
      double wgt = 1.0/(double)(adj[plane_ix].second->numPoints());
      for (auto it=adj[plane_ix].second->pointsBegin();
	   it!=adj[plane_ix].second->pointsEnd(); ++it)
	{
	  Vector3D pos0 = (*it)->getPoint();
	  Point xpos(pos0[0], pos0[1], pos0[2]);
	  Point vec = xpos - loc;
	  Point pos2 = loc + (vec*axis2)*axis2;
	  mid += wgt*pos2;
	}
      shared_ptr<Plane> plane(new Plane(mid, axis2));
      vector<RevEngPoint*> inptp, outptp;
      vector<pair<double, double> > dist_angp;
      double maxdp, avdp;
      int num_inp;
      vector<double> parvalsp;
      RevEngUtils::distToSurf(adj[plane_ix].second->pointsBegin(),
			      adj[plane_ix].second->pointsEnd(),
			      plane, tol, maxdp, avdp, num_inp,
			      inptp, outptp, parvalsp, dist_angp, -1.0);

      std::ofstream ofd1("in_out_alt_plane.g2");
      ofd1 << "400 1 0 4 155 50 50 255" << std::endl;
      ofd1 << inptp.size() << std::endl;
      for (size_t kr=0; kr<inptp.size(); ++kr)
	ofd1 << inptp[kr]->getPoint() << std::endl;
      ofd1 << "400 1 0 4 50 155 50 255" << std::endl;
      ofd1 << outptp.size() << std::endl;
      for (size_t kr=0; kr<outptp.size(); ++kr)
	ofd1 << outptp[kr]->getPoint() << std::endl;
      
      // Check accuracy of alternative cylinder
      Point Cx = adj[cyl_ix].first->direction2();
      Point Cy = axis.cross(Cx);
      Cx = Cy.cross(axis);  // Project to plane
      
      vector<pair<vector<RevEngPoint*>::iterator,
		  vector<RevEngPoint*>::iterator> > points;
      points.push_back(std::make_pair(adj[cyl_ix].second->pointsBegin(),
				      adj[cyl_ix].second->pointsEnd()));
      BoundingBox bb = adj[cyl_ix].second->getBbox();
      Point xpos;
      double rad;
      Point low = bb.low();
      Point high = bb.high();
      RevEngUtils::computeCylPosRadius(points, low, high, axis, Cx, Cy,
				       xpos, rad);
      shared_ptr<Cylinder> cyl(new Cylinder(rad, xpos, axis, Cy));
      vector<RevEngPoint*> inpt, outpt;
      vector<pair<double, double> > dist_ang;
      double maxd, avd;
      int num_in;
      vector<double> parvals;
      RevEngUtils::distToSurf(points[0].first, points[0].second,
			      cyl, tol, maxd, avd, num_in,
			      inpt, outpt, parvals, dist_ang, -1.0);

      std::ofstream ofd2("in_out_alt_cyl.g2");
      ofd2 << "400 1 0 4 155 50 50 255" << std::endl;
      ofd2 << inpt.size() << std::endl;
      for (size_t kr=0; kr<inpt.size(); ++kr)
	ofd2 << inpt[kr]->getPoint() << std::endl;
      ofd2 << "400 1 0 4 50 155 50 255" << std::endl;
      ofd2 << outpt.size() << std::endl;
      for (size_t kr=0; kr<outpt.size(); ++kr)
	ofd2 << outpt[kr]->getPoint() << std::endl;
#endif
      
      // Bound cylinder
      int num_pt = adj[cyl_ix].second->numPoints();
      double dist, ang;
      double angtol = 0.1;
      int ka;
      for (ka=0; ka<num_pt; ++ka)
	{
	  Vector2D uv = adj[cyl_ix].second->getPoint(ka)->getPar();
	  adj[cyl_ix].second->getPoint(ka)->getSurfaceDist(dist, ang);
	  if (ang <= angtol)
	    {
	      dom[0] = dom[1] = uv[0];
	      dom[2] = dom[3] = uv[1];
	      break;
	    }
	}
      for (; ka<num_pt; ++ka)
	{
	  Vector2D uv = adj[cyl_ix].second->getPoint(ka)->getPar();
	  adj[cyl_ix].second->getPoint(ka)->getSurfaceDist(dist, ang);
	  if (ang > angtol)
	    continue;
	  dom[0] = std::min(dom[0], uv[0]);
	  dom[1] = std::max(dom[1], uv[0]);
	  dom[2] = std::min(dom[2], uv[1]);
	  dom[3] = std::max(dom[3], uv[1]);
	}
      if (dom[1] - dom[2] > 2*M_PI)
	{
	  dom[0] = 0.0;
	  dom[1] = 2*M_PI;
	}
      
#ifdef DEBUG_TORUSCONTEXT
      shared_ptr<ElementarySurface> elem2(adj[cyl_ix].first->clone());
      elem2->setParameterBounds(dom[0], dom[2], dom[1], dom[3]);
      std::ofstream ofcyl("bd_cyl.g2");
      elem2->writeStandardHeader(ofcyl);
      elem2->write(ofcyl);
#endif

      // Compute distance from cylinder to plane
      Point pos1 = adj[cyl_ix].first->point(0.5*(dom[0]+dom[1]), dom[2]);
      Point pos2 = adj[cyl_ix].first->point(0.5*(dom[0]+dom[1]), dom[3]);
      double upar1, upar2, upar3, vpar1, vpar2, vpar3, dist1, dist2, dist3;
      Point close1, close2, close3;
      double eps = 1.0e-6;
      adj[plane_ix].first->closestPoint(pos1, upar1, vpar1, close1, dist1, eps);
      adj[plane_ix].first->closestPoint(pos2, upar2, vpar2, close2, dist2, eps);
      Point pos3 = (dist1 < dist2) ? pos + dom[2]*axis2 : pos + dom[3]*axis2;
      adj[plane_ix].first->closestPoint(pos3, upar3, vpar3, close3, dist3, eps);

      R2 = std::min(dist1, dist2);
      R1 = radius - R2;
      pos = close3;

      // Check if the torus in inwards or outwards
      double cdist;
      RevEngPoint *closest_planar = adj[plane_ix].second->closestPoint(pos, cdist);
      if (cdist > 0.99*radius)
	{
	  outer = false;
	  R1 = radius + R2;
	}

#ifdef DEBUG_TORUSCONTEXT
      // Bound plane
      double dom2[4];
      int num_pt2 = adj[plane_ix].second->numPoints();
      Vector2D uv2 = adj[plane_ix].second->getPoint(0)->getPar();
      dom2[0] = dom2[1] = uv2[0];
      dom2[2] = dom2[3] = uv2[1];

      for (int ka=0; ka<num_pt2; ++ka)
	{
	  Vector2D uv2 = adj[plane_ix].second->getPoint(ka)->getPar();
	  dom2[0] = std::min(dom2[0], uv2[0]);
	  dom2[1] = std::max(dom2[1], uv2[0]);
	  dom2[2] = std::min(dom2[2], uv2[1]);
	  dom2[3] = std::max(dom2[3], uv2[1]);
	}
      shared_ptr<ElementarySurface> elem3(adj[plane_ix].first->clone());
      elem3->setParameterBounds(dom2[0], dom2[2], dom2[1], dom2[3]);
      std::ofstream ofpl("bd_plane.g2");
      elem3->writeStandardHeader(ofpl);
      elem3->write(ofpl);
#endif
      
      return true;
    }
  
  return false;  // Requested configuration not found   
}

//===========================================================================
RevEngPoint* RevEngRegion::closestPoint(const Point& pos, double& dist)
//===========================================================================
{
    dist = std::numeric_limits<double>::max();
    int ix = -1;
    for (size_t ki=0; ki<group_points_.size(); ++ki)
      {
	Vector3D xyz = group_points_[ki]->getPoint();
	Point pt(xyz[0], xyz[1], xyz[2]);
	double curr_dist = pos.dist(pt);
	if (curr_dist < dist)
	  {
	    dist = curr_dist;
	    ix = (int)ki;
	  }
      }

    if (ix >= 0)
      return group_points_[ix];
    else
      return 0;
}

//===========================================================================
bool RevEngRegion::extractTorus(double tol, int min_pt, double angtol,
				double mean_edge_len,
				int prefer_elementary_,
				vector<shared_ptr<HedgeSurface> >& hedgesfs,
				vector<HedgeSurface*>& prevsfs,
				vector<vector<RevEngPoint*> >& out_groups)
//===========================================================================
{
  bool found = false;
  int min_nmb = 20;
  if ((int)group_points_.size() < min_nmb)
    return false;
  
  // std::ofstream of("curr_region.g2");
  // writeRegionInfo(of);

  vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> > adj_elem, adj_elem_base, adj_primary;
  getAdjacentElemInfo(adj_elem, adj_elem_base, adj_primary);
  
  shared_ptr<Torus> tor2;
  shared_ptr<Torus> tor1 = computeTorus(group_points_, tol, tor2);
  if (!tor1.get())
    return false;
  
  // Check accuracy
  double maxd1, avd1, maxd2, avd2;
  int num1, num2;
  //shared_ptr<ParamSurface> surf = cyl;
  vector<RevEngPoint*> in1, out1,in2, out2;
  vector<pair<double, double> > dist_ang1;
  vector<pair<double, double> > dist_ang2;
  vector<double> parvals1;
  vector<double> parvals2;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			  tor1, tol, maxd1, avd1, num1, in1, out1,
			  parvals1, dist_ang1, angtol);
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
  			  tor2, tol, maxd2, avd2, num2, in2, out2,
			  parvals2, dist_ang2, angtol);

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

  size_t insize = std::max(in1.size(), in2.size());
  if (insize > group_points_.size()/2 && (int)insize > min_nmb)
    {
      shared_ptr<Torus> tor_in2;
      shared_ptr<Torus> tor_in1 =
	computeTorus(in1.size() > in2.size() ? in1 : in2, tol, tor_in2);
if (tor_in1.get())
{
      vector<RevEngPoint*> inpt_in1, outpt_in1, inpt_in2, outpt_in2; 
      vector<pair<double, double> > dist_ang_in1, dist_ang_in2;
      double maxd_in1, avd_in1, maxd_in2, avd_in2;
      int num1_in, num2_in;
      vector<double> parvals_in1;
      vector<double> parvals_in2;
      RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			      tor_in1, tol, maxd_in1, avd_in1, num1_in,
			      inpt_in1, outpt_in1, parvals_in1, dist_ang_in1,
			      angtol);
  
      RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			      tor_in2, tol, maxd_in2, avd_in2, num2_in,
			      inpt_in2, outpt_in2, parvals_in2,
			      dist_ang_in2, angtol);
  
      std::ofstream ofd2("in_out_tor2.g2");
      ofd2 << "400 1 0 4 155 50 50 255" << std::endl;
      ofd2 << inpt_in1.size() << std::endl;
      for (size_t kr=0; kr<inpt_in1.size(); ++kr)
	ofd2 << inpt_in1[kr]->getPoint() << std::endl;
      ofd2 << "400 1 0 4 50 155 50 255" << std::endl;
      ofd2 << outpt_in1.size() << std::endl;
      for (size_t kr=0; kr<outpt_in1.size(); ++kr)
	ofd2 << outpt_in1[kr]->getPoint() << std::endl;
      ofd2 << "400 1 0 4 155 50 50 255" << std::endl;
      ofd2 << inpt_in2.size() << std::endl;
      for (size_t kr=0; kr<inpt_in2.size(); ++kr)
	ofd2 << inpt_in2[kr]->getPoint() << std::endl;
      ofd2 << "400 1 0 4 50 155 50 255" << std::endl;
      ofd2 << outpt_in2.size() << std::endl;
      for (size_t kr=0; kr<outpt_in2.size(); ++kr)
	ofd2 << outpt_in2[kr]->getPoint() << std::endl;

      if (num1_in > num1 || avd_in1 < avd1)
	{
	  std::swap(tor1, tor_in1);
	  std::swap(num1, num1_in);
	  std::swap(avd1, avd_in1);
	  std::swap(maxd1, maxd_in1);
	  std::swap(parvals1, parvals_in1);
	  std::swap(dist_ang1, dist_ang_in1);
	  // std::cout << "Torus swap 1" << std::endl;
	  // std::cout << group_points_.size() << " " << num1_in << " ";
	  // std::cout << num1 << " " << avd_in1 << " " << avd1;
	  // std::cout << " " << maxd_in1 << " " << maxd1 << std::endl;
	}
      if (num2_in > num2 || avd_in2 < avd2)
	{
	  std::swap(tor2, tor_in2);
	  std::swap(num2, num2_in);
	  std::swap(avd2, avd_in2);
	  std::swap(maxd2, maxd_in2);
	  std::swap(parvals2, parvals_in2);
	  std::swap(dist_ang2, dist_ang_in2);
	  // std::cout << "Torus swap 2" << std::endl;
	  // std::cout << group_points_.size() << " " << num2_in << " ";
	  // std::cout << num2 << " " << avd_in2 << " " << avd2;
	  // std::cout << " " << maxd_in2 << " " << maxd2 << std::endl;
	}
    }
    }
  
  int num = (int)group_points_.size();
  double maxd_init, avd_init;
  int num_init;
  getAccuracy(maxd_init, avd_init, num_init);
  if (accuracyOK(min_pt, tol, num1, avd1) && num1 > num2) 
    {
      bool OK = true;
      if (associated_sf_.size() > 0)
	{
	  // Check with current approximating surface
	  double acc_fac1 = 1.25;
	  double acc_fac2 = 0.75;
	  if (prefer_elementary_ == ALWAYS_ELEM)
	    OK = false;
	  else if (prefer_elementary_ == PREFER_ELEM &&
		   ((double)num1 < acc_fac1*num_init ||
		    avd1 > acc_fac2*avd_init))
	    OK = false;
	  else if (prefer_elementary_ == BEST_ACCURACY &&
		   (num1 < num_init || avd1 > avd_init))
	    OK = false;
	}

      if (OK)
	{
	  found = true;
	  for (size_t kh=0; kh<group_points_.size(); ++kh)
	    {
	      group_points_[kh]->setPar(Vector2D(parvals1[2*kh],parvals1[2*kh+1]));
	      group_points_[kh]->setSurfaceDist(dist_ang1[kh].first, dist_ang1[kh].second);
	    }
	  setAccuracy(maxd1, avd1, num1);      
      
	  std::cout << "Torus 1. N1: " << num << ", N2: " << num1 << ", max: " << maxd1 << ", av: " << avd1 << std::endl;
	  // associated_sf_.clear();
	  // for (int ka=0; ka<divmod[0]->nmbEntities(); ++ka)
	  // 	{
	  //shared_ptr<ParamSurface> tor1_2 = divmod[0]->getSurface(ka);
	  shared_ptr<HedgeSurface> hedge(new HedgeSurface(tor1, this));
	  for (size_t kh=0; kh<associated_sf_.size(); ++kh)
	    prevsfs.push_back(associated_sf_[kh]);
	  setHedge(hedge.get());
	  hedgesfs.push_back(hedge);

	}
	// }
	}
  else if (accuracyOK(min_pt, tol, num2, avd2) && num2 >= num1) 
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
	  for (size_t kh=0; kh<group_points_.size(); ++kh)
	    {
	      group_points_[kh]->setPar(Vector2D(parvals2[2*kh],parvals2[2*kh+1]));
	      group_points_[kh]->setSurfaceDist(dist_ang2[kh].first, dist_ang2[kh].second);
	    }
	  setAccuracy(maxd2, avd2, num2);
      
	  std::cout << "Torus 2. N1: " << num << ", N2: " << num2 << ", max: " << maxd2 << ", av: " << avd2 << std::endl;
	  // associated_sf_.clear();
	  // for (int ka=0; ka<divmod[0]->nmbEntities(); ++ka)
	  // 	{
	  // 	  shared_ptr<ParamSurface> tor2_2 = divmod[0]->getSurface(ka);
	  shared_ptr<HedgeSurface> hedge(new HedgeSurface(tor2, this));
	  for (size_t kh=0; kh<associated_sf_.size(); ++kh)
	    prevsfs.push_back(associated_sf_[kh]);
	  setHedge(hedge.get());
	  hedgesfs.push_back(hedge);
	  // }
	}
    }
  if (!primary_.get() ||
      (num1 >= num_in_primary_ && avd1 < avdist_primary_))
    setPrimarySf(tor1, maxd1, avd1, num1);
  if (!primary_.get() ||
      (num1 >= num_in_primary_ && avd1 < avdist_primary_))
    setPrimarySf(tor2, maxd2, avd2, num2);
  int stop_break = 1;
  return found;
}
 
//===========================================================================
  shared_ptr<Torus> RevEngRegion::computeTorus(vector<RevEngPoint*>& points,
					       double tol,
					       shared_ptr<Torus>& torus2)
//===========================================================================
{
  std::ofstream of("torus_compute.g2");
  
  // Compute mean curvature and initial point in plane
  double k2mean = 0.0;
  double wgt = 1.0/(double)points.size();
  for (size_t kr=0; kr<points.size(); ++kr)
    {
      double kmax = points[kr]->maxPrincipalCurvature();
      k2mean += wgt*kmax;
    }
  double rd = 1.0/k2mean;
  
  vector<Point> centr(points.size());
  Point mid(0.0, 0.0, 0.0);
  for (size_t kr=0; kr<points.size(); ++kr)
    {
      double kmax = points[kr]->maxPrincipalCurvature();
      k2mean += wgt*kmax;

      Vector3D xyz = points[kr]->getPoint();
      Point xyz2(xyz[0], xyz[1], xyz[2]);
      Point norm = points[kr]->getMongeNormal();
      centr[kr] = xyz2 + rd*norm;
      mid += wgt*centr[kr];
    }
  
  of << "400 1 0 4 155 100 0 255" << std::endl;
  of << points.size() << std::endl;
  for (size_t kr=0; kr<points.size(); ++kr)
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


  // std::ofstream ofimpl("impl_plane.g2");
  // impl->visualize(points, ofimpl);
  
  Vector3D xyz = points[0]->getPoint();
  Point xyz2(xyz[0], xyz[1], xyz[2]);
  Point Cx = centr[0] - xyz2;
  Cx -= (Cx*normal)*normal;
  Cx.normalize();
  Point Cy = Cx.cross(normal);
  
  double rad;
  Point pnt;
  try {
    RevEngUtils::computeCircPosRadius(centr, normal, Cx, Cy, pnt, rad);
  }
  catch (...)
    {
      shared_ptr<Torus> dummy;
      return dummy;
    }
  pnt -= ((pnt - pos)*normal)*normal;
  shared_ptr<Circle> circ(new Circle(rad, pnt, normal, Cx));
  circ->writeStandardHeader(of);
  circ->write(of);

  vector<Point> rotated;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(points.begin(), points.end()));
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
  Point xpos;
  vector<double> param;
  try {
    curveApprox(rotated, tol, circ2, param, spl, xpos);
  }
  catch (...)
    {
    }
  if (spl.get())
    {
      spl->writeStandardHeader(of3);
      spl->write(of3);
    }
  
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
  //std::cout << "Torus small radius: " << fabs(rd) << ", " << crad << std::endl;

  torus2 = tor2;
  return tor1;
}

//===========================================================================
bool RevEngRegion::tryOtherSurf(int prefer_elementary_, bool replace)
//===========================================================================
{
  if (associated_sf_.size() > 0 && (!replace))
    return false;
  
  if (prefer_elementary_ == BEST_ACCURACY ||
      associated_sf_.size() == 0)
    return true;

  if (prefer_elementary_ == ALWAYS_ELEM)
    {
      int sfcode;
      ClassType type = associated_sf_[0]->instanceType(sfcode);
      double fac = 5.0;
      if (type == Class_Plane && MAH_ > fac*MAK_)
	return true;
      else
	return false;
    }
  
  // Check if the current accuracy is sufficient
  int num = (int)group_points_.size();
  double maxd_init, avd_init;
  int num_init;
  getAccuracy(maxd_init, avd_init, num_init);
  if ((double)num_init > 0.75*num)
    return false;
  else
    return true;
  
}


//===========================================================================
bool RevEngRegion::extractFreeform(double tol, int min_pt, double angtol,
				   double mean_edge_len,
				   int prefer_elementary_,
				   vector<shared_ptr<HedgeSurface> >& hedgesfs,
				   vector<HedgeSurface*>& prevsfs,
				   vector<vector<RevEngPoint*> >& out_groups)
//===========================================================================
{
  // std::ofstream of("curr_region.g2");
  // writeRegionInfo(of);
  
  // std::ofstream of2("curr_normals.g2");
  // writeUnitSphereInfo(of2);

  bool found = false;
  int min_nmb = 20;
  //double eps = 1.0e-6;
  if ((int)group_points_.size() < min_nmb)
    return false;

  vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> > adj_elem, adj_elem_base, adj_primary;
  getAdjacentElemInfo(adj_elem, adj_elem_base, adj_primary);
  
  shared_ptr<SplineSurface> spl = computeFreeform(tol);
  if (!spl.get())
    return false;
  std::ofstream ofs("spl.g2");
  spl->writeStandardHeader(ofs);
  spl->write(ofs);
  
  // Check accuracy
  double maxd, avd; 
  int num2; 
  vector<RevEngPoint*> inpt, outpt;
  vector<pair<double, double> > dist_ang;
  vector<double> parvals;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
  			  spl, tol, maxd, avd, num2, inpt, outpt,
			  parvals, dist_ang, angtol);
  std::ofstream ofd("in_out_spl.g2");
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
  if (associated_sf_.size() > 0 ||  (num2 > min_pt && num2 > num/2))
    {
      extractOutPoints(dist_ang, tol, out_groups);
      if (out_groups.size() > 0)
	{
	  // Some points has been removed from the group. Redo surface
	  // generation
	  spl = computeFreeform(tol);
	  if (!spl.get())
	    return false;
	  std::ofstream ofs2("spl2.g2");
	  spl->writeStandardHeader(ofs2);
	  spl->write(ofs2);

	  num = (int)group_points_.size();
	  dist_ang.clear();
	  parvals.clear();
	  inpt.clear();
	  outpt.clear();
	  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
				  spl, tol, maxd, avd, num2, inpt, outpt,
				  parvals, dist_ang, angtol);
	  std::ofstream ofd2("in_out_spl2.g2");
	  ofd2 << "400 1 0 4 155 50 50 255" << std::endl;
	  ofd2 << inpt.size() << std::endl;
	  for (size_t kr=0; kr<inpt.size(); ++kr)
	    ofd2 << inpt[kr]->getPoint() << std::endl;
	  ofd2 << "400 1 0 4 50 155 50 255" << std::endl;
	  ofd2 << outpt.size() << std::endl;
	  for (size_t kr=0; kr<outpt.size(); ++kr)
	    ofd2 << outpt[kr]->getPoint() << std::endl;
	}
      int stop_break_out = 1;
    }

  if (accuracyOK(min_pt, tol, num2, avd))
    {
     bool OK = true;
      if (associated_sf_.size() > 0)
	{
	  // Check with current approximating surface
	  double acc_fac1 = 1.25;
	  double acc_fac2 = 0.75;
	  if (prefer_elementary_ == ALWAYS_ELEM)
	    OK = false;
	  else if (prefer_elementary_ == PREFER_ELEM &&
		   ((double)num2 < acc_fac1*num_init ||
		    avd > acc_fac2*avd_init))
	    OK = false;
	  else if (prefer_elementary_ == BEST_ACCURACY &&
		   (num2 < num_init || avd > avd_init))
	    OK = false;
	}

      if (OK)
	{
	  found = true;
	  for (size_t kh=0; kh<group_points_.size(); ++kh)
	    {
	      group_points_[kh]->setPar(Vector2D(parvals[2*kh],parvals[2*kh+1]));
	      group_points_[kh]->setSurfaceDist(dist_ang[kh].first, dist_ang[kh].second);
	    }
	  setAccuracy(maxd, avd, num2);
	  std::cout << "Spline. N1: " << num << ", N2: " << num2 << ", max: " << maxd << ", av: " << avd << std::endl;
	  shared_ptr<HedgeSurface> hedge(new HedgeSurface(spl, this));
	  for (size_t kh=0; kh<associated_sf_.size(); ++kh)
	    prevsfs.push_back(associated_sf_[kh]);
	  setHedge(hedge.get());
	  hedgesfs.push_back(hedge);
	}
    }
  int stop_break0 = 1;
  return found;
}


//===========================================================================
 void RevEngRegion::extractOutPoints(vector<pair<double, double> >& dist_ang,
				     double tol,
				     vector<vector<RevEngPoint*> >& out_groups)
//===========================================================================
 {
   vector<RevEngPoint*> move;
   for (size_t kr=0; kr<dist_ang.size(); ++kr)
     {
       if (dist_ang[kr].first > tol &&
	   (associated_sf_.size() == 0 ||
	    group_points_[kr]->getSurfaceDist() > tol))
	 move.push_back(group_points_[kr]);
     }
   extractSpesPoints(move, out_groups, true);
   if (out_groups.size() > 0)
     updateInfo();
 }

//===========================================================================
 void RevEngRegion::identifyAngPoints(vector<pair<double, double> >& dist_ang,
				      double tol,
				      vector<RevEngPoint*>& ang_points)
//===========================================================================
 {
   for (size_t kr=0; kr<dist_ang.size(); ++kr)
     {
       if (dist_ang[kr].second > tol)
	 ang_points.push_back(group_points_[kr]);
     }
 }

//===========================================================================
 void RevEngRegion::identifyDistPoints(vector<pair<double, double> >& dist_ang,
				       double tol, double maxd, double avd,
				       vector<RevEngPoint*>& dist_points)
//===========================================================================
 {
   double tol2 = (avd > 0.9*tol) ? 3.0*tol : 2.0*tol;
   for (size_t kr=0; kr<dist_ang.size(); ++kr)
     {
       if (dist_ang[kr].first > tol2)
	 dist_points.push_back(group_points_[kr]);
     }
 }

//===========================================================================
 void RevEngRegion::extractSpesPoints(vector<RevEngPoint*>& move,
				      vector<vector<RevEngPoint*> >& out_groups,
				      bool outer)
//===========================================================================
 {
   // Note: derived information in this group is not updated!!
   shared_ptr<RevEngRegion> dummy_reg(new RevEngRegion(edge_class_type_));
      
   for (size_t kr=0; kr<move.size(); ++kr)
     move[kr]->setRegion(dummy_reg.get());  // Preliminary

   vector<vector<RevEngPoint*> > sep_move;
   for (size_t kr=0; kr<move.size(); ++kr)
     {
       if (move[kr]->visited())
	 continue;
       vector<RevEngPoint*> curr_move;
       move[kr]->fetchConnected(dummy_reg.get(), (int)move.size(), curr_move);
       sep_move.push_back(curr_move);
     }
      
   for (size_t kr=0; kr<move.size(); ++kr)
     move[kr]->unsetVisited();
      
   std::ofstream m1("move_group.g2");
   std::ofstream m2("not_move_group.g2");
   for (size_t kr=0; kr<sep_move.size(); ++kr)
     {
       size_t kh;
       for (kh=0; kh<sep_move[kr].size(); ++kh)
	 {
	   vector<ftSamplePoint*> next = sep_move[kr][kh]->getNeighbours();
	   size_t kh1;
	   for (kh1=0; kh1<next.size(); ++kh1)
	     {
	       RevEngPoint *pt = dynamic_cast<RevEngPoint*>(next[kh1]);
		  
	       RevEngRegion *reg = pt->region();
	       if (reg != this && reg != dummy_reg.get())
		 break;
	     }
	   if (kh1 < next.size())
	     break;
	 }
       if (kh < sep_move[kr].size() || (!outer))
	 {
	   out_groups.push_back(sep_move[kr]);
	   m1 << "400 1 0 0" << std::endl;
	   m1 << sep_move[kr].size() << std::endl;
	   for (size_t kh1=0; kh1<sep_move[kr].size(); ++kh1)
	     m1 << sep_move[kr][kh1]->getPoint() << std::endl;
	 }
       else
	 {
	   m2 << "400 1 0 0" << std::endl;
	   m2 << sep_move[kr].size() << std::endl;
	   for (size_t kh1=0; kh1<sep_move[kr].size(); ++kh1)
	     m2 << sep_move[kr][kh1]->getPoint() << std::endl;
	 }
     }
   for (size_t kr=0; kr<move.size(); ++kr)
     move[kr]->setRegion(this);

   // Extract group of out points from current group
   for (size_t ki=0; ki<out_groups.size(); ++ki)
     for (size_t kj=0; kj<out_groups[ki].size(); ++kj)
       {
	 out_groups[ki][kj]->unsetRegion();
	 out_groups[ki][kj]->addMove();
	 vector<RevEngPoint*>::iterator it = std::find(group_points_.begin(),
						 group_points_.end(),
						 out_groups[ki][kj]);
	 if (it != group_points_.end())
	   {
	     std::swap(*it, group_points_[group_points_.size()-1]);
	     group_points_.pop_back();
	     //group_points_.erase(it);
	   }
       }

   int stop_break = 1;
 }
 
//===========================================================================
shared_ptr<SplineSurface> RevEngRegion::updateFreeform(double tol)
//===========================================================================
{
  shared_ptr<SplineSurface> dummy;
  if (associated_sf_.size() == 0)
    return dummy;

  shared_ptr<ParamSurface> surf0 = associated_sf_[0]->surface();
  shared_ptr<SplineSurface> surf =
    dynamic_pointer_cast<SplineSurface,ParamSurface>(surf0);
  if (!surf.get())
    return dummy;

  vector<double> data, param;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Vector3D xyz = group_points_[ki]->getPoint();
      Vector2D uv = group_points_[ki]->getPar();
      data.insert(data.end(), xyz.begin(), xyz.end());
      param.insert(param.end(), uv.begin(), uv.end());
    }

  int dim = surf->dimension();
  ApproxSurf approx(surf, data, param, dim, tol);
  //ApproxSurf approx(surf3, data, param, dim, tol, 0, false, true, 0, true);
  approx.setMBA(true);
  approx.setFixBoundary(false);
  int max_iter = 1;
  double maxd, avd;
  int num_out;
  shared_ptr<SplineSurface> surf2;
  try {
    surf2 = approx.getApproxSurf(maxd, avd, num_out, max_iter);
  }
  catch (...)
  {
    std::cout << "Surface update failed" << std::endl;
  }
  
 return surf2;
}
  

//===========================================================================
shared_ptr<SplineSurface> RevEngRegion::computeFreeform(double tol)
//===========================================================================
{
// Parameterize
  vector<double> data, param;
  int inner1=0, inner2=0;
      
  // Compute PCA axes
  double lambda[3];
  Point eigen1, eigen2, eigen3;
  getPCA(lambda, eigen1, eigen2, eigen3);

  //bool usebasesf = false;
  bool done = false;
  bool close1 = false, close2 = false;
  if (associated_sf_.size() > 0)
    {
      done = parameterizeOnSurf(associated_sf_[0]->surface(), data,
				param, inner1, inner2, close1, close2);
    }
  if (!done && primary_.get() && avdist_primary_ <= 5.0*tol)
    {
      done = parameterizeOnSurf(primary_, data, param, inner1, inner2, close1, close2);
    }
  if (!done)
    {
     // Parameterize on plane
      RevEngUtils::parameterizeWithPlane(group_points_, bbox_, eigen1,
					 eigen2, data, param);
    }
  std::ofstream ofpar("parpoints.g2");
  int nmbpar = (int)param.size()/2;
  ofpar << "400 1 0 0" << std::endl;
  ofpar << nmbpar << std::endl;
  for (int ka=0; ka<nmbpar; ++ka)
    ofpar << param[2*ka] << " " << param[2*ka+1] << "  0.0" << std::endl;

  vector<double> param2;
  double umin, umax, vmin, vmax;
  bool repar = reparameterize(param, param2, umin, umax, vmin, vmax);
  //std::cout << "repar: " << repar << std::endl;
  if (repar)
    {
      vector<double> p_edg1, p_edg2;
      for (size_t ki=0; ki<group_points_.size(); ++ki)
	{
	  vector<ftSamplePoint*> next = group_points_[ki]->getNeighbours();
	  for (size_t kj=0; kj<next.size(); ++kj)
	    {
	      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(next[kj]);
	      if (pt->region() != this)
		continue;
	      std::vector<RevEngPoint*>::iterator it = std::find(group_points_.begin(),
								 group_points_.end(), pt);
	      if (it != group_points_.end())
		{
		  int ix = (int)(it - group_points_.begin());
		  p_edg1.push_back(param[2*ki]);
		  p_edg1.push_back(param[2*ki+1]);
		  p_edg2.push_back(param2[2*ki]);
		  p_edg2.push_back(param2[2*ki+1]);
		  p_edg1.push_back(param[2*ix]);
		  p_edg1.push_back(param[2*ix+1]);
		  p_edg2.push_back(param2[2*ix]);
		  p_edg2.push_back(param2[2*ix+1]);
		}
	    }
	}

      std::ofstream ofe1("par_edgs1.g2");
      std::ofstream ofe2("par_edgs2.g2");
      ofe1 << "410 1 0 4 255 0 0 255" << std::endl;
      ofe2 << "410 1 0 4 255 0 0 255" << std::endl;
      ofe1 << p_edg1.size()/4 << std::endl;
      ofe2 << p_edg2.size()/4 << std::endl;
      for (size_t ki=0; ki<p_edg1.size(); ki+=4)
	{
	  ofe1 << p_edg1[ki] << " " << p_edg1[ki+1] << " 0.0 " << p_edg1[ki+2];
	  ofe1 << " " << p_edg1[ki+3] << " 0.0" << std::endl;
	  ofe2 << p_edg2[ki] << " " << p_edg2[ki+1] << " 0.0 " << p_edg2[ki+2];
	  ofe2 << " " << p_edg2[ki+3] << " 0.0" << std::endl;
	}

      double plenfac = 5.0;
      if (umax - umin > plenfac*(vmax - vmin))
	{
	  if (inner1 <= inner2)
	    inner1 = inner2 + 1;
	}
      else if (vmax - vmin > plenfac*(umax - umin))
	{
	  if (inner2 <= inner1)
	    inner2 = inner1 + 1;
	}
      std::swap(param, param2);
    }

  // Extend with extra points in corners far from the point cloud
  size_t nmb_prev_extend = param.size();
  if ((!close1) && (!close2))
    {
      try {
	extendInCorner(data, param, umin, umax, vmin, vmax);
      }
      catch (...)
	{
	  std::cout << "Corner extend failed" << std::endl;
	}
    }
  
  // Approximate
  double maxd, avd;
  int num_out;
  int max_iter = (associated_sf_.size() > 0) ? 3 : 6; //2;
  int dim = bbox_.low().dimension();
  int order = 4;
  shared_ptr<SplineSurface> surf;
  double del = 0.01;
  vector<double> parvals;
  try {
    surf = RevEngUtils::surfApprox(data, dim, param, order, order, 
				   order+inner1,  order+inner2,
				   close1, close2, max_iter, 
				   tol, maxd, avd, num_out, parvals, del);
  }
  catch (...)
    {
       std::cout << "Surface approximation failed" << std::endl;
   }

  if (surf.get())
    {
      std::ofstream of("spline_approx.g2");
      surf->writeStandardHeader(of);
      surf->write(of);
    }
  else
    int no_surf = 1;

  return surf;
}
											 
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

//===========================================================================
bool RevEngRegion::reparameterize(vector<double>& param, vector<double>& param2,
				  double& umin, double& umax, double& vmin, double& vmax)
//===========================================================================
{
  // Define raster
  int nmb_div = 15;
  vector<vector<int> > raster;
  vector<double> param_copy(param.begin(), param.end());  // Raster definition
  // changes sequence of parameter points
  defineRaster(param_copy, nmb_div, raster, umin, umax, vmin, vmax);
  int div2 = (int)raster.size();
  if (div2 == 0)
    return false;
  int div1 = (int)raster[0].size();
  double udel = (umax - umin)/(double)(div1);
  double vdel = (vmax - vmin)/(double)(div2);

  // Count number of empty cells
  int nmb_zero = 0;
  for (int kb=0; kb<div2; ++kb)
    for (int ka=0; ka<div1; ++ka)
      if (raster[kb][ka] == 0)
	++nmb_zero;

  if (nmb_zero < (int)(0.5*div1*div2))
    return false;

  // Identify points on edges
  double eps = 1.0e-10;
  vector<double> edge_pts;
  for (size_t kr=0; kr<param.size(); kr+=2)
    {
      if (param[kr]-umin < eps || umax-param[kr] < eps)
	{
	  edge_pts.push_back((param[kr]-umin < eps) ? 0.0 : (double)div1);
	  edge_pts.push_back(0.5*(int)(2*(param[kr+1]-vmin)/vdel));
	}
      if (param[kr+1]-vmin < eps || vmax-param[kr+1] < eps)
	{
	  edge_pts.push_back(0.5*(int)(2*(param[kr]-umin)/udel));
	  edge_pts.push_back((param[kr+1]-vmin < eps) ? 0.0 : (double)div2);
	}
    }

  // Find candidates for endpoints of a mid curve through the paramaeter points
  vector<Point> cand_par;
  vector<double> wd;
  int i1, i2, i3, i4, ix, kc;
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
	    size_t kr;
	    for (kr=0; kr<edge_pts.size(); kr+=2)
	      {
		if (((kc == 0 && edge_pts[kr+1] == 0) ||
		     (kc == 1 && edge_pts[kr+1] == div2)) &&
		    (edge_pts[kr] >= i1 && edge_pts[kr] <= i2))
		  break;
	      }
	    
	    if (kr < edge_pts.size())
	      {
		Point currpar(2);
		currpar[0] = edge_pts[kr]; //0.5*(i1+i2);
		currpar[1] = (kc == 0) ? 0 : (double)div2;
		cand_par.push_back(currpar);
		wd.push_back(0.5*udel*(i2-i1));
	      }
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
	    size_t kr;
	    for (kr=0; kr<edge_pts.size(); kr+=2)
	      {
		if (((kc == 0 && edge_pts[kr] == 0) ||
		     (kc == 1 && edge_pts[kr] == div1)) &&
		    (edge_pts[kr+1] >= i1 && edge_pts[kr+1] <= i2))
		  break;
	      }
	    
	    if (kr < edge_pts.size())
	      {
		Point currpar(2);
		currpar[0] = (kc == 0) ? 0 : (double)div1;
		currpar[1] = edge_pts[kr+1]; //0.5*(i1+i2);
		cand_par.push_back(currpar);
		wd.push_back(0.5*vdel*(i2-i1));
	      }
	  }
      }

  if (cand_par.size() < 2)
    return false;
  
  // Decide on endpoints
  int ix1=-1, ix2=-1;
  double maxd = std::numeric_limits<double>::lowest();
  for (size_t kr=0; kr<cand_par.size(); ++kr)
    for (size_t kh=kr+1; kh<cand_par.size(); ++kh)
      {
	double dd = cand_par[kr].dist(cand_par[kh]);
	dd -= (wd[kr] + wd[kh]);
	if (dd > maxd)
	  {
	    ix1 = (int)kr;
	    ix2 = (int)kh;
	    maxd = dd;
	  }
      }

  // Compute points along a sceleton in the parameter point cloud
  double curr1[2], curr2[2];
  double wd1, wd2;
  vector<double> ptpar;
  vector<double> width;
  ptpar.push_back(cand_par[ix1][0]);
  ptpar.push_back(cand_par[ix1][1]);
  width.push_back(wd[ix1]);

  Point prev;
  int sgn1=0, sgn2=0;
  bool loop = false;
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

      curr1[0] = curr2[0] = ptpar[kc];
      curr1[1] = curr2[1] = ptpar[kc+1];
      if (sgn1 == 0 || curr1[0] == prev[0])
	sgn1 = (cand_par[ix2][0] > curr1[0]) ? 1 : -1;
      if (sgn2 == 0 || curr2[1] == prev[1])
	sgn2 = (cand_par[ix2][1] > curr2[1]) ? 1 : -1;
      double fac1 = (curr1[0] == 0.0) ? 0.5 : 1.0;
      double fac2 = (curr2[1] == 0.0) ? 0.5 : 1.0;
      curr1[0] += fac1*sgn1;
      curr2[1] += fac2*sgn2;
      prev = Point(ptpar[kc], ptpar[kc+1]);
      if (prev.dist(cand_par[ix2]) <= 1.0)
	break;
      curr1[0] = std::min(curr1[0], div1-0.5);
      curr1[1] = std::min(curr1[1], div2-0.5);
      curr2[0] = std::min(curr2[0], div1-0.5);
      curr2[1] = std::min(curr2[1], div2-0.5);
      curr1[0] = std::max(curr1[0], 0.0);
      curr1[1] = std::max(curr1[1], 0.0);
      curr2[0] = std::max(curr2[0], 0.0);
      curr2[1] = std::max(curr2[1], 0.0);

      // Extent with constant first parameter direction
      getParExtent(curr1, 0, raster, i1, i2);
      curr1[1] = 0.5*(double)(i1+i2+1);
      wd1 = 0.5*(i1 - i2);

      // Constant second parameter direction
      getParExtent(curr2, 1, raster, i3, i4);
      curr2[0] = 0.5*(double)(i3+i4+1);
      wd2 = 0.5*(i3 - i4);

      // Select continuation
      Point lastpar(ptpar[ptpar.size()-2], ptpar[ptpar.size()-1]);
      Point c1(curr1[0], curr1[1]);
      Point c2(curr2[0], curr2[1]);
      if (i1-i2 < i3-i4 || c2.dist(lastpar) < eps) 
	{
	  ptpar.push_back(curr1[0]);
	  ptpar.push_back(curr1[1]);
	  width.push_back(udel*wd1);
	  lastpar = c1;
	}
      else
	{
	  ptpar.push_back(curr2[0]);
	  ptpar.push_back(curr2[1]);
	  width.push_back(vdel*wd2);
	  lastpar = c2;
	}
	
      // Check for loops
      for (size_t kh=2; kh<ptpar.size()-2; kh+=2)
	{
	  Point currpar(ptpar[kh], ptpar[kh+1]);
	  if (currpar.dist(lastpar) < eps)
	    {
	      loop = true;
	      break;
	    }
	}
      if (loop)
	break;
      
      int stop_break0 = 1;      
    }
  if (loop)
    return false;
  
  Point last(ptpar[ptpar.size()-2], ptpar[ptpar.size()-1]);
  if (last.dist(cand_par[ix2]) > 0.0)
    {
      ptpar.push_back(cand_par[ix2][0]);
      ptpar.push_back(cand_par[ix2][1]);
      width.push_back(wd[ix2]);
    }

  if (ptpar.size() <= 2)
    return false;
  
  // Translate to real coordinates
  vector<double> ptpar2(ptpar.size());
  for (size_t kr=0; kr<ptpar.size(); kr+=2)
    {
      ptpar2[kr] = umin+udel*ptpar[kr];
      ptpar2[kr+1] = vmin+vdel*ptpar[kr+1];
    }

  // Extend
  double npar1[2], npar2[2];
  for (int ka=0; ka<2; ++ka)
    {
      npar1[ka] = ptpar2[ka] - 0.1*(ptpar2[2+ka]-ptpar2[ka]);
      npar2[ka] = ptpar2[ptpar2.size()-2+ka] + 0.1*(ptpar2[ptpar2.size()-2+ka]-
						   ptpar2[ptpar2.size()-4+ka]);
    }
  ptpar2.insert(ptpar2.begin(), &npar1[0], &npar1[0]+2);
  ptpar2.insert(ptpar2.end(), &npar2[0], &npar2[0]+2);
  width.insert(width.begin(), width[0]);
  width.push_back(width[width.size()-1]);
  

  std::ofstream ofp("parpt.g2");
  ofp << "400 1 0 4 200 55 0 255" << std::endl;
  ofp << ptpar2.size()/2 << std::endl;
  for (size_t kr=0; kr<ptpar2.size(); kr+=2)
    ofp << ptpar2[kr] << " " << ptpar2[kr+1] << " 0.0" << std::endl;

  // Define "mid" curve
  vector<double> par2(ptpar2.size()/2);
  par2[0] = 0.0;
  for (size_t kr=2; kr<ptpar2.size(); kr+=2)
    {
      double dd = Utils::distance_squared(&ptpar2[kr-2], &ptpar2[kr], &ptpar2[kr]);
      par2[kr/2] = par2[(kr-2)/2] + sqrt(dd);
    }

  int in = 4, ik = 4;
  int dim = 2;
  double tol = 0.1;
  ApproxCurve approx(ptpar2, par2, dim, tol, in, ik);
  
  double maxdist, avdist;
  shared_ptr<SplineCurve> midcv = approx.getApproxCurve(maxdist, avdist, 1);

  ApproxCurve approx2(width, par2, 1, tol, in, ik);
  
  double maxdistw, avdistw;
  shared_ptr<SplineCurve> widthcv = approx2.getApproxCurve(maxdistw, avdistw, 1);

  std::ofstream ofmid("midcv.g2");
  midcv->writeStandardHeader(ofmid);
  midcv->write(ofmid);

  double avw = 0.0;
  for (size_t kr=0; kr<width.size(); ++kr)
    avw += width[kr];
  avw /= (double)width.size();

  param2.resize(param.size());
  double tmin = midcv->startparam();
  double tmax = midcv->endparam();
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
      double tmp = fabs(der[1]*vec);
      // if (tmp > 1.0e-6)
      // 	std::cout << "inner: " << tmp << std::endl;
      if (vec2*der[1] < 0.0)
	dist *= -1.0;
      Point wdt = widthcv->ParamCurve::point(tpar);
      //dist *= (avw/wdt[0]);
      param2[kr] = tpar;
      param2[kr+1] = dist;
      int stop_break0 = 1;
    }
  
  std::ofstream ofpar2("parpoints_repar.g2");
  int nmbpar = (int)param2.size()/2;
  ofpar2 << "400 1 0 0" << std::endl;
  ofpar2 << nmbpar << std::endl;
  for (int ka=0; ka<nmbpar; ++ka)
    ofpar2 << param2[2*ka] << " " << param2[2*ka+1] << "  0.0" << std::endl;
  
  vector<vector<int> > raster2;
  double umin2, umax2, vmin2, vmax2;
  vector<double> param2_copy(param2.begin(), param2.end());  // Raster definition
  // changes sequence of parameter points
  defineRaster(param2_copy, nmb_div, raster2, umin2, umax2, vmin2, vmax2);

  // Count number of empty cells
  int nmb_zero2 = 0;
  for (size_t kr=0; kr<raster2.size(); ++kr)
    for (size_t kh=0; kh<raster2[kr].size(); ++kh)
      if (raster2[kr][kh] == 0)
	++nmb_zero2;

  if (nmb_zero2 >= nmb_zero)
    return false;

  umin = umin2;
  umax = umax2;
  vmin = vmin2;
  vmax = vmax2;
  return true;
}

//===========================================================================
void RevEngRegion::getParExtent(double curr[2], int pdir, vector<vector<int> >& raster,
				int& i1, int& i2)
//===========================================================================
{
  int div = (pdir == 0) ? (int)raster.size() : (int)raster[0].size();
  int j1 = (int)curr[pdir];
  int pdir2 = 1 - pdir;
  
  for (i1=(int)curr[pdir2]; i1<div; ++i1)
    {
      int r1 = (pdir == 0) ? raster[i1][j1] : raster[j1][i1];
      if (r1 > 0)
	break;
    }
  if (i1 == div)
    i1=(int)curr[pdir2];
  for (; i1<div; ++i1)
    {
      int r1 = (pdir == 0) ? raster[i1][j1] : raster[j1][i1];
      if (r1 == 0)
	break;
    }
  for (i2=(int)curr[pdir2]; i2>=0; --i2)
    {
      int r1 = (pdir == 0) ? raster[i2][j1] : raster[j1][i2];
      if (r1 > 0)
	break;
    }
  if (i2 < 0)
    i2=(int)curr[pdir2];
  for (; i2>=0; --i2)
    {
      int r1 = (pdir == 0) ? raster[i2][j1] : raster[j1][i2];
      if (r1 == 0)
      break;
    }
  }


//===========================================================================
void RevEngRegion::defineRaster(vector<double>& param, int nmb_div,
				vector<vector<int> >& raster, double& umin,
				double& umax, double& vmin, double& vmax)
//===========================================================================
{
  umin = std::numeric_limits<double>::max();
  umax = std::numeric_limits<double>::lowest();
  vmin = std::numeric_limits<double>::max();
  vmax = std::numeric_limits<double>::lowest();
  for (size_t kr=0; kr<param.size(); kr+=2)
    {
      umin = std::min(umin, param[kr]);
      umax = std::max(umax, param[kr]);
      vmin = std::min(vmin, param[kr+1]);
      vmax = std::max(vmax, param[kr+1]);
    }
  
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
  
  raster.resize(div2);
  for (size_t kr=0; kr<raster.size(); ++kr)
    {
      raster[kr].resize(div1, 0);
    }

  int nmb_pts = (int)param.size()/2;
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
}

//===========================================================================
void RevEngRegion::extendInCorner(vector<double>& data, vector<double>& param,
				  double umin, double umax, double vmin, double vmax)
//===========================================================================
{
  Point corner[4];
  corner[0] = Point(umin, vmin);
  corner[1] = Point(umax, vmin);
  corner[2] = Point(umin, vmax);
  corner[3] = Point(umax, vmax);
  double cdist[4];
  cdist[0] = cdist[1] = cdist[2] = cdist[3] = std::numeric_limits<double>::max();
  for (size_t kr=0; kr<param.size(); kr+=2)
    {
      Point currpar(param[kr], param[kr+1]);
      for (int ka=0; ka<4; ++ka)
	{
	  double dd = corner[ka].dist(currpar);
	  cdist[ka] = std::min(cdist[ka], dd);
	}
    }

  int div = 15;
  double minc = std::min((umax-umin)/(double)div, (vmax-vmin)/(double)div);
  vector<double> corner_data;
  vector<double> corner_par;
  int min_nmb = 20;
  for (int ka=0; ka<4; ++ka)
    {
      if (cdist[ka] > minc)
      {
	// Estimate corner point
	// Collect points in the vicinity of the corner
	vector<double> data2;
	vector<double> param2;
	double umin2 = std::numeric_limits<double>::max();
	double umax2 = std::numeric_limits<double>::lowest();
	double vmin2 = std::numeric_limits<double>::max();
	double vmax2 = std::numeric_limits<double>::lowest();
	double lim = 1.1*cdist[ka];  //2.0
	for (size_t kr=0; kr<param.size()/2; ++kr)
	  {
	    Point currpar(param[2*kr], param[2*kr+1]);
	    if (currpar.dist(corner[ka]) < lim)
	      {
		data2.insert(data2.end(), data.begin()+3*kr, data.begin()+3*(kr+1));
		param2.insert(param2.end(), param.begin()+2*kr, param.begin()+2*(kr+1));
		umin2 = std::min(umin2, param[2*kr]);
		umax2 = std::max(umax2, param[2*kr]);
		vmin2 = std::min(vmin2, param[2*kr+1]);
		vmax2 = std::max(vmax2, param[2*kr+1]);
	      }
	  }

	while ((int)data2.size()/3 < min_nmb)
	  {
	    data2.clear();
	    param2.clear();
	    lim += 0.5*cdist[ka];
	    for (size_t kr=0; kr<param.size()/2; ++kr)
	      {
		Point currpar(param[2*kr], param[2*kr+1]);
		if (currpar.dist(corner[ka]) < lim)
		  {
		    data2.insert(data2.end(), data.begin()+3*kr, data.begin()+3*(kr+1));
		    param2.insert(param2.end(), param.begin()+2*kr, param.begin()+2*(kr+1));
		    umin2 = std::min(umin2, param[2*kr]);
		    umax2 = std::max(umax2, param[2*kr]);
		    vmin2 = std::min(vmin2, param[2*kr+1]);
		    vmax2 = std::max(vmax2, param[2*kr+1]);
		  }
	      }
	  }

	// Approximate sub cloud with a planar surface
	umin2 = std::min(umin2, corner[ka][0]);
	umax2 = std::max(umax2, corner[ka][0]);
	vmin2 = std::min(vmin2, corner[ka][1]);
	vmax2 = std::max(vmax2, corner[ka][1]);
	umin2 -= 0.1*(umax2 - umin2);
	umax2 += 0.1*(umax2 - umin2);
	vmin2 -= 0.1*(vmax2 - vmin2);
	vmax2 += 0.1*(vmax2 - vmin2);

	shared_ptr<SplineSurface> plane =
	  RevEngUtils::surfApprox(data2, 3, param2, 2, 2, 2, 2, umin, umax,
				  vmin, vmax);
	Point pos = plane->ParamSurface::point(corner[ka][0], corner[ka][1]);
	corner_data.insert(corner_data.end(), pos.begin(), pos.end());
	corner_par.insert(corner_par.end(), corner[ka].begin(), corner[ka].end());
      }
    }

  if (corner_data.size() > 0)
    {
      data.insert(data.end(), corner_data.begin(), corner_data.end());
      param.insert(param.end(), corner_par.begin(), corner_par.end());

      std::ofstream of("added_corners.g2");
      for (size_t ki=0; ki<corner_data.size()/3; ++ki)
	{
	  of << "400 1 0 4 255 0 0 255" << std::endl;
	  of << "1" << std::endl;
	  for (int ka=0; ka<3; ++ka)
	    of << corner_data[3*ki+ka] << " ";
	  of << std::endl;
	}
    }
}

//===========================================================================
void RevEngRegion::implicitizeSplit()
//===========================================================================
{
  double lambda[2];
  Point eigen1, eigen2, eigen3;
  getPCA(lambda, eigen1, eigen2, eigen3);
  vector<Point> pos_and_der(3*group_points_.size());
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      Vector3D xyz = group_points_[kr]->getPoint();
      pos_and_der[3*kr] = Point(xyz[0], xyz[1], xyz[2]);
      Point norm = group_points_[kr]->getMongeNormal();
      Point vec1 = eigen1 - (eigen1*norm)*norm;
      vec1.normalize();
      pos_and_der[3*kr+1] = vec1;
      Point vec2 = norm.cross(vec1);
      vec2.normalize();
      pos_and_der[3*kr+2] = vec2;
    }

  int degree = 2;
  vector<double> coefs;
  ImplicitApprox approx;
  approx.polynomialSurf(pos_and_der, degree, coefs);

  double maxfield, avfield, maxdist, avdist, maxang, avang;
  int ndiv;
  approx.polynomialSurfAccuracy(pos_and_der, degree, coefs,
				maxfield, avfield, maxdist, avdist,
				ndiv, maxang, avang);
  int stop_break = 1;
  
 }

//===========================================================================
bool RevEngRegion::parameterizeOnSurf(shared_ptr<ParamSurface> surf, 
				      vector<double>& data, vector<double>& param,
				      int& inner1, int& inner2, bool& close1, bool& close2)
//===========================================================================
{
  ClassType classtype = surf->instanceType();
  if (classtype == Class_Cone || classtype == Class_Torus)
    {
      // Check angle. An almost plane cone is not appropriate for
      // parametrization
      shared_ptr<ParamSurface> sf = surf;
      shared_ptr<BoundedSurface> bdsf =
	dynamic_pointer_cast<BoundedSurface,ParamSurface>(sf);
      if (bdsf.get())
	sf = bdsf->underlyingSurface();
      if (classtype == Class_Cone)
	{
	  shared_ptr<Cone> cone =
	    dynamic_pointer_cast<Cone,ParamSurface>(sf);
	  double angle = cone->getConeAngle();
	  if (angle > 0.3*M_PI)
	    return false;
	}
      else  if (classtype == Class_Torus)
	{
	  shared_ptr<Torus> torus =
	    dynamic_pointer_cast<Torus,ParamSurface>(sf);
	  double rad1 = torus->getMajorRadius();
	  double rad2 = torus->getMinorRadius();
	  if (rad2 > 0.9*rad1)
	    return false;
	}
    }
  
  bool OK = RevEngUtils::parameterizeOnPrimary(group_points_, surf, data, param,
					       inner1, inner2, close1, close2);
  return OK;
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
      //std::cout << "Exception in SVD" << std::endl;
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
RevEngRegion::peelOffRegions(int min_point, double tol,
			     vector<shared_ptr<HedgeSurface> >& hedgesfs,
			     vector<vector<RevEngPoint*> >& out_groups,
			     vector<RevEngPoint*>& single_pts)
//===========================================================================
{
  std::ofstream of1("region_to_peel.g2");
  writeRegionInfo(of1);

  // Check if parts of the region can be represented by a plane
  double angtol0 = 0.2; //0.1*M_PI;
  double angtol = 0.1;
  double dfrac = 0.75;

  // Count fraction of normals closer to the centre than the tolerance
  double in_frac;
  if (normalcone_.angle() <= angtol0)
    in_frac = 1.0;
  else
    {
      int nmb_in = 0;
      for (size_t kr=0; kr<group_points_.size(); ++kr)
	{
	  Point normal = group_points_[kr]->getMongeNormal();
	  if (avnorm_.angle(normal) <= angtol0)
	    nmb_in++;
	}
      in_frac = (double)nmb_in/(double)group_points_.size();
    }
  frac_norm_in_ = in_frac;

  bool apply_plane = false, apply_cyl = false;
  
  double in_lim = 0.9;
  double maxdist1, avdist1, maxdist2, avdist2;
  int num_inside1, num_inside2;
  shared_ptr<Plane> plane = computePlane(group_points_, avnorm_);
  vector<RevEngPoint*> inpt1, outpt1, inpt2, outpt2;
  vector<pair<double, double> > dist_ang1, dist_ang2;
  vector<double> parvals1, parvals2;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			  plane, tol, maxdist1, avdist1, num_inside1, inpt1,
			  outpt1, parvals1, dist_ang1, angtol);

  vector<vector<RevEngPoint*> > configs;
  shared_ptr<Cylinder> cyl = computeCylinder(group_points_, tol, configs);
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			  cyl, tol, maxdist2, avdist2, num_inside2, inpt2,
			  outpt2, parvals2, dist_ang2, angtol);

  vector<RevEngPoint*> ang_points1;
  vector<RevEngPoint*> ang_points2;
  double peel_frac = 0.2;
  if ((in_frac >= in_lim && (!(MAH_ > 10.0*MAK_ && MAH_ > 0.1))) ||
      accuracyOK(min_point, tol, num_inside1, avdist1))
    {
      std::ofstream ofd("in_out_plane.g2");
      ofd << "400 1 0 4 155 50 50 255" << std::endl;
      ofd << inpt1.size() << std::endl;
      for (size_t kr=0; kr<inpt1.size(); ++kr)
	ofd << inpt1[kr]->getPoint() << std::endl;
      ofd << "400 1 0 4 50 155 50 255" << std::endl;
      ofd << outpt1.size() << std::endl;
      for (size_t kr=0; kr<outpt1.size(); ++kr)
	ofd << outpt1[kr]->getPoint() << std::endl;
      
      if (num_inside1 > group_points_.size()/3)
	  identifyAngPoints(dist_ang1, angtol, ang_points1);

      if (ang_points1.size() > 0)
	{
	  std::ofstream ofa("ang_points1.g2");
	  ofa << "400 1 0 4 50 50 155 255" << std::endl;
	  ofa << ang_points1.size() << std::endl;
	  for (size_t kr=0; kr<ang_points1.size(); ++kr)
	    ofa << ang_points1[kr]->getPoint() << std::endl;
	}
      apply_plane = true;
   }

  if (MAH_ > 3*MAK_ || accuracyOK(min_point, tol, num_inside2, avdist2))
    {
      std::ofstream ofd("in_out_cyl.g2");
      ofd << "400 1 0 4 155 50 50 255" << std::endl;
      ofd << inpt2.size() << std::endl;
      for (size_t kr=0; kr<inpt2.size(); ++kr)
	ofd << inpt2[kr]->getPoint() << std::endl;
      ofd << "400 1 0 4 50 155 50 255" << std::endl;
      ofd << outpt2.size() << std::endl;
      for (size_t kr=0; kr<outpt2.size(); ++kr)
	ofd << outpt2[kr]->getPoint() << std::endl;
      
      if (num_inside2 > group_points_.size()/3)
	{
	  Point axis = cyl->direction();
	  double pihalf = 0.5*M_PI;
	  for (size_t ki=0; ki<group_points_.size(); ++ki)
	    {
	      if (dist_ang2[ki].second > angtol && dist_ang2[ki].first > dfrac*tol)
		ang_points2.push_back(group_points_[ki]);
	    }
	}

      if (ang_points2.size() > 0)
	{
	  std::ofstream ofa("ang_points2.g2");
	  ofa << "400 1 0 4 50 50 155 255" << std::endl;
	  ofa << ang_points2.size() << std::endl;
	  for (size_t kr=0; kr<ang_points2.size(); ++kr)
	    ofa << ang_points2[kr]->getPoint() << std::endl;
	}
      apply_cyl = true;
   }

  if (ang_points1.size() > 0 && ang_points1.size() < peel_frac*group_points_.size() &&
      (!(ang_points2.size() > 0 && ang_points2.size() < ang_points1.size())))
    {
      extractSpesPoints(ang_points1, out_groups, true);
      if (out_groups.size() > 0)
	{
	  updateInfo();
	  inpt1.clear();
	  outpt1.clear();
	  parvals1.clear();
	  dist_ang1.clear();
	  plane = computePlane(group_points_, avnorm_);
	  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
				  plane, tol, maxdist1, avdist1, num_inside1, inpt1,
				  outpt1, parvals1, dist_ang1, angtol);
	}
    }
  else if (ang_points2.size() > 0 && ang_points2.size() < peel_frac*group_points_.size())
    {
      extractSpesPoints(ang_points2, out_groups, true);
      if (out_groups.size() > 0)
	{
	  updateInfo();
	  configs.clear();
	  inpt2.clear();
	  outpt2.clear();
	  parvals2.clear();
	  dist_ang2.clear();
	  cyl = computeCylinder(group_points_, tol, configs);
	  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
				  cyl, tol, maxdist2, avdist2, num_inside2, inpt2,
				  outpt2, parvals2, dist_ang2, angtol);
	}
    }

  if (apply_plane && apply_cyl)
    {
      if ((num_inside2 > num_inside1 && avdist2 < avdist1) ||
	  (avdist2 < avdist1 && num_inside2 > (int)(0.9*num_inside1)))
	apply_plane = false;
      else
	apply_cyl = false;
    }

  if (apply_plane && num_inside1 > group_points_.size()/3 || avdist1 < tol)
    {
      setBaseSf(plane, maxdist1, avdist1, num_inside1);
      
      if (accuracyOK(min_point, tol, num_inside1, avdist1))
	{
	  for (size_t kh=0; kh<group_points_.size(); ++kh)
	    {
	      group_points_[kh]->setPar(Vector2D(parvals1[2*kh],parvals1[2*kh+1]));
	      group_points_[kh]->setSurfaceDist(dist_ang1[kh].first, dist_ang1[kh].second);
	    }
	  setAccuracy(maxdist1, avdist1, num_inside1);
	  shared_ptr<HedgeSurface> hedge(new HedgeSurface(plane, this));
	  setHedge(hedge.get());
	  hedgesfs.push_back(hedge);
	}
    }
  else if (apply_cyl && num_inside2 > group_points_.size()/3 || avdist2 < tol)
    {
    
      if (configs.size() > 1)
	{
	  // Split point group 
	  int keep_ix = -1;
	  int keep_nmb = 0;
	  for (size_t ki=0; ki<configs.size(); ++ki)
	    if ((int)configs[ki].size() > keep_nmb)
	      {
		keep_nmb = (int)configs[ki].size();
		keep_ix = (int)ki;
	      }
      
	  for (size_t ki=0; ki<configs.size(); ++ki)
	    {
	      if ((int)ki == keep_ix)
		continue;

	      extractSpesPoints(configs[ki], out_groups);
	    }

	  // Check that the remaing point cloud is connected
	  splitRegion(out_groups);
	  updateInfo();
	}
      else
	{
	  setBaseSf(cyl, maxdist2, avdist2, num_inside2);

	  if (accuracyOK(min_point, tol, num_inside2, avdist2))
	    {
	      for (size_t kh=0; kh<group_points_.size(); ++kh)
		{
		  group_points_[kh]->setPar(Vector2D(parvals2[2*kh],parvals2[2*kh+1]));
		  group_points_[kh]->setSurfaceDist(dist_ang2[kh].first,
						    dist_ang2[kh].second);
		}
	      setAccuracy(maxdist2, avdist2, num_inside2);
	      shared_ptr<HedgeSurface> hedge(new HedgeSurface(cyl, this));
	      setHedge(hedge.get());
	      hedgesfs.push_back(hedge);
	    }
  	}
    }
      

   for (size_t kr=0; kr<out_groups.size(); )
    {
      if (out_groups[kr].size() == 1)
  	{
  	  out_groups[kr][0]->unsetRegion();
  	  single_pts.push_back(out_groups[kr][0]);
	  out_groups.erase(out_groups.begin()+kr);
  	}
      else
	++kr;
    }
  
  int stop_break = 1;
}


  
//===========================================================================
void
RevEngRegion::splitComposedRegions(int classtype,
				   vector<shared_ptr<RevEngRegion> >& added_groups,
				   vector<RevEngPoint*>& single_pts)
//===========================================================================
{
  std::ofstream of1("regions_to_split.g2");
  writeRegionInfo(of1);
  
  vector<RevEngPoint*> deviant;
  vector<vector<RevEngPoint*> > connected;
  splitFromSurfaceNormals(deviant, connected);
  for (size_t kj=0; kj<connected.size(); ++kj)
    {
      if (connected[kj].size() == 1)
	{
	  connected[kj][0]->unsetRegion();
	  single_pts.push_back(connected[kj][0]);
	}
      else
	{
	  shared_ptr<RevEngRegion> reg(new RevEngRegion(classtype,
							edge_class_type_,
							connected[kj]));
	  added_groups.push_back(reg);
	}
    }
      
  if (deviant.size() > 0)
    {
      for (size_t kh=0; kh<deviant.size(); ++kh)
	deviant[kh]->setGaussRad(1.0);
      if (deviant.size() == 1)
	single_pts.push_back(deviant[0]);
      else
	{
	  shared_ptr<RevEngRegion> reg2(new RevEngRegion(classtype,
							 edge_class_type_,
							 deviant));
	  vector<vector<RevEngPoint*> > connected2;
	  reg2->splitRegion(connected2);
	  if (reg2->numPoints() == 1)
	    {
	      reg2->group_points_[0]->unsetRegion();
	      single_pts.push_back(reg2->group_points_[0]);
	    }
	  else
	    added_groups.push_back(reg2);
	  for (size_t kj=0; kj<connected2.size(); ++kj)
	    {
	      if (connected2[kj].size() == 1)
		{
		  connected2[kj][0]->unsetRegion();
		  single_pts.push_back(connected2[kj][0]);
		}
	      else
		{
		  shared_ptr<RevEngRegion> reg3(new RevEngRegion(classtype,
								 edge_class_type_,
								 connected2[kj]));
		  added_groups.push_back(reg3);
		}
	    }
	}
      // vector<RevEngPoint*> deviant2;
      // reg2->splitFromSurfaceNormals(deviant2, connected2);
      // added_groups.push_back(reg2);
      // for (size_t kj=0; kj<connected2.size(); ++kj)
      // 	{
      // 	  shared_ptr<RevEngRegion> reg3(new RevEngRegion(classtype,
      // 							 connected2[kj]));
      // 	  added_groups.push_back(reg3);
      // 	}
      // if (deviant2.size() > 0)
      // 	{
      // 	  shared_ptr<RevEngRegion> reg4(new RevEngRegion(classtype,
      // 							 deviant2));
      // 	  //vector<vector<RevEngPoint*> > deviant4;
      // 	  //reg4->splitRegion(deviant4);
      // 	  added_groups.push_back(reg4);
      // 	  // for (size_t kh=0; kh<deviant4.size(); ++kh)
      // 	  //   {
      // 	  //     shared_ptr<RevEngRegion> reg5(new RevEngRegion(classtype,
      // 	  // 						     deviant4[kh]));
      // 	  //     added_groups.push_back(reg5);
      // 	  //   }
      // 	}
    }

  std::ofstream of2("split_regions.g2");
  writeRegionInfo(of2);
  for (size_t ki=0; ki<added_groups.size(); ++ki)
    added_groups[ki]->writeRegionInfo(of2);

  int stop_break = 1;
}

//===========================================================================
void
RevEngRegion::splitWithShapeIndex(vector<shared_ptr<RevEngRegion> >& updated_regions)
//===========================================================================
{
  // Collect new regions based on shape index
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      if (group_points_[ki]->region() != this)
	continue;
      shared_ptr<RevEngRegion> curr(new RevEngRegion(CLASSIFICATION_SHAPEINDEX,
						     edge_class_type_));
      updated_regions.push_back(curr);
      group_points_[ki]->setRegion(curr.get());
      curr->collect(group_points_[ki], this);
    }
  
  if (updated_regions.size() <= 1)
    {
      // Reset current region pointer
      for (size_t ki=0; ki<group_points_.size(); ++ki)
	group_points_[ki]->setRegion(this);
    }
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
      //RevEngPoint *pt = dynamic_cast<RevEngPoint*>(group_points_[kr]);
      double Grad = group_points_[kr]->getGaussRad();
      if (Grad < radlim)
	{
	  group_points_[kr]->unsetRegion();
	  group_points_[kr]->addMove();
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
      of3 << "410 1 0 4 0 55 200 255" << std::endl;
      of3 << group_points_.size() << std::endl;
      for (size_t kr=0; kr<group_points_.size(); ++kr)
	{
	  Vector3D xyz = group_points_[kr]->getPoint();
	  Point xyz2(xyz[0], xyz[1], xyz[2]);
	  Point norm = group_points_[kr]->getMongeNormal();
	  of3 << xyz2 << " " << xyz2 + norm << std::endl;
	}
    }
  if (smallrad.size() > 0)
    {
      of3 << "400 1 0 4 0 155 100 255" << std::endl;
      of3 << smallrad.size() << std::endl;
      for (size_t kr=0; kr<smallrad.size(); ++kr)
	of3 << smallrad[kr]->getPoint() << std::endl;
      of3 << "410 1 0 4 100 55 100 255" << std::endl;
      of3 << smallrad.size() << std::endl;
      for (size_t kr=0; kr<smallrad.size(); ++kr)
	{
	  Vector3D xyz = smallrad[kr]->getPoint();
	  Point xyz2(xyz[0], xyz[1], xyz[2]);
	  Point norm = smallrad[kr]->getMongeNormal();
	  of3 << xyz2 << " " << xyz2 + norm << std::endl;
	}
    }

  // Split remaining region points into disjunct groups
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

  if (connected.size() <= 1)
    return;
  
  group_points_.clear();
  group_points_ = connected[0];
  
  // Update bounding box and principal curvature summary
  updateInfo();

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
  double fac = 1.0/(double)group_points_.size();
  for  (size_t kj=0; kj<group_points_.size(); ++kj)
    {
      double k1 = group_points_[kj]->minPrincipalCurvature();
      double k2 = group_points_[kj]->maxPrincipalCurvature();
      double H = group_points_[kj]->meanCurvature();
      double K = group_points_[kj]->GaussCurvature();
      mink1_ = std::min(mink1_, fabs(k1));
      maxk1_ = std::max(maxk1_, fabs(k1));
      mink2_ = std::min(mink2_, fabs(k2));
      maxk2_ = std::max(maxk2_, fabs(k2));
      avH_ += fac*H;
      avK_ += fac*K;
      MAH_ += fac*fabs(H);
      MAK_ += fac*fabs(K);
      Vector3D point = group_points_[kj]->getPoint();
      Point point2(point[0], point[1], point[2]);
      bbox_.addUnionWith(point2);
    }

  normalcone_ = DirectionCone(group_points_[0]->getMongeNormal());
  avnorm_ = Point(0.0, 0.0, 0.0);
  for  (size_t kj=0; kj<group_points_.size(); ++kj)
    {
      Point norm = group_points_[kj]->getMongeNormal();
      normalcone_.addUnionWith(norm);
      avnorm_ += fac*norm;
    }

  double anglim = 0.2;
  if (normalcone_.angle() <= anglim)
    frac_norm_in_ = 1.0;
  else
    {
      int nmb_in = 0;
      for (size_t kr=0; kr<group_points_.size(); ++kr)
	{
	  Point normal = group_points_[kr]->getMongeNormal();
	  if (avnorm_.angle(normal) <= anglim)
	    nmb_in++;
	}
      frac_norm_in_ = (double)nmb_in/(double)group_points_.size();
    }
  
  if (hasSurface() && group_points_.size() > 0)
    {
      Vector2D par = group_points_[0]->getPar();
      domain_[0] = domain_[1] = par[0];
      domain_[2] = domain_[3] = par[1];
      for  (size_t kj=1; kj<group_points_.size(); ++kj)
	{
	  Vector2D par = group_points_[kj]->getPar();
	  domain_[0] = std::min(domain_[0], par[0]);
	  domain_[1] = std::max(domain_[1], par[0]);
	  domain_[2] = std::min(domain_[2], par[1]);
	  domain_[3] = std::max(domain_[3], par[1]);
 	}
    }
}

//===========================================================================
void RevEngRegion::addPoint(RevEngPoint* point)
//===========================================================================
{
  int nmb = (int)group_points_.size();
  group_points_.push_back(point);
  point->setRegion(this);
  Vector3D point2 = point->getPoint();
  Point point3(point2[0], point2[1], point2[2]);
  bbox_.addUnionWith(point3);
  Point norm = point->getMongeNormal();
  normalcone_.addUnionWith(norm);
  double k1 = point->minPrincipalCurvature();
  double k2 = point->maxPrincipalCurvature();
  mink1_ = std::min(mink1_, fabs(k1));
  maxk1_ = std::max(maxk1_, fabs(k1));
  mink2_ = std::min(mink2_, fabs(k2));
  maxk2_ = std::max(maxk2_, fabs(k2));
  avnorm_ = (nmb*avnorm_ + norm)/(double)(nmb+1);
  double H = point->meanCurvature();
  double K = point->GaussCurvature();
  avH_ = (nmb*avH_ + H)/(double)(nmb+1);
  avK_ = (nmb*avK_ + K)/(double)(nmb+1);
  MAH_ = (nmb*MAH_ + fabs(H))/(double)(nmb+1);
  MAK_ = (nmb*MAK_ + fabs(K))/(double)(nmb+1);

  if (hasSurface())
    {
      Vector2D par = point->getPar();
      domain_[0] = std::min(domain_[0], par[0]);
      domain_[1] = std::max(domain_[1], par[0]);
      domain_[2] = std::min(domain_[2], par[1]);
      domain_[3] = std::max(domain_[3], par[1]);
    }
 
}

//===========================================================================
void RevEngRegion::computeDomain()
//===========================================================================
{
  if (hasSurface() && group_points_.size() > 0)
    {
      Vector2D par = group_points_[0]->getPar();
      domain_[0] = domain_[1] = par[0];
      domain_[2] = domain_[3] = par[1];
      for  (size_t kj=1; kj<group_points_.size(); ++kj)
	{
	  Vector2D par = group_points_[kj]->getPar();
	  domain_[0] = std::min(domain_[0], par[0]);
	  domain_[1] = std::max(domain_[1], par[0]);
	  domain_[2] = std::min(domain_[2], par[1]);
	  domain_[3] = std::max(domain_[3], par[1]);
 	}
    }
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
  //double del1 = (p2-p1)/(double)(n1);
  //double del2 = (p4-p3)/(double)(n2);

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

struct integrate_info
{
  RevEngRegion *adjacent;
  int nmb_pt_adj;
  double maxd, avd, maxd_adj, avd_adj;
  int nmb_in, nmb_in_adj;
  double min_ang, max_ang, min_dist, max_dist;
  bool outlier;

  integrate_info()
  {
    adjacent = 0;
    nmb_pt_adj = nmb_in = nmb_in_adj = 0;
    maxd = avd = maxd_adj = avd_adj = 0.0;
    min_ang = max_ang = min_dist = max_dist = -1.0;
    outlier = false;
  }

  void setAdjacent(RevEngRegion* adj)
  {
    adjacent = adj;
  }

  void setOutlier()
  {
    outlier = true;
  }

  void setNmbAdj(int nmb)
  {
    nmb_pt_adj = nmb;
  }

  void setAccuracy(double maxd1, double avd1, int nmb_in1, double maxd2,
		   double avd2, int nmb_in2)
  {
    maxd = maxd1;
    avd = avd1;
    nmb_in = nmb_in1;
    maxd_adj = maxd2;
    avd_adj = avd2;
    nmb_in_adj = nmb_in2;
  }

  void setMinMax(double min_a, double max_a, double min_d, double max_d)
  {
    min_ang = min_a;
    max_ang = max_a;
    min_dist = min_d;
    max_dist = max_d;
  }
};
  
//===========================================================================
bool RevEngRegion::integrateInAdjacent(double mean_edge_len, int min_next,
				       int max_next, double tol, double angtol,
				       int max_nmb_outlier,
				       RevEngRegion* taboo)
//===========================================================================
{
  // if (adjacent_regions_.size() != 1)
  //   return false;   // To be removed

  if ((int)group_points_.size() > max_next/2)
    return false;

#ifdef DEBUG_INTEGRATE
  std::ofstream of("curr_integrate.g2");
  of << "400 1 0 4 0 255 0 255" << std::endl;
  of << group_points_.size() << std::endl;
  for (size_t kh=0; kh<group_points_.size(); ++kh)
    of << group_points_[kh]->getPoint() << std::endl;
#endif 
  size_t kj=0;
  double local_len = group_points_[0]->getMeanEdgLen(10.0*mean_edge_len);
  double radius = 3.0*(local_len + mean_edge_len);
  radius = std::min(radius, 20.0*mean_edge_len);
  size_t adjsize = adjacent_regions_.size();
  vector<integrate_info> info(adjsize);
  
  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it, ++kj)
    info[kj].setAdjacent(*it);
  
  if (group_points_.size() <= max_nmb_outlier)
    {
      // Simplified test
      vector<RevEngRegion*> adj_reg;
      vector<pair<double,double> > adj_info;
      double lentol = 2.0*mean_edge_len;
      for (size_t ki=0; ki<group_points_.size(); ++ki)
	{
	  vector<RevEngRegion*> pt_adj_reg;
	  vector<RevEngPoint*> pt_adj_pt;
	  group_points_[ki]->getAdjInfo(mean_edge_len, pt_adj_reg, pt_adj_pt);
	  Point monge1 = group_points_[ki]->getMongeNormal();
	  for (size_t kj=0; kj<pt_adj_pt.size(); ++kj)
	    {
	      if (pt_adj_reg[kj] == this)
		continue;
	      double len = group_points_[ki]->pntDist(pt_adj_pt[kj]);
	      if (len > lentol)
		continue;
	      Point monge2 = pt_adj_pt[kj]->getMongeNormal();
	      if (monge1*monge2 < 0.0 || monge1.angle(monge2) > angtol)
		continue;
	      adj_reg.push_back(pt_adj_reg[kj]);
	      adj_info.push_back(std::make_pair(len, monge1.angle(monge2)));
	    }
	}

      for (size_t kj=0; kj<info.size(); ++kj)
	{
	  RevEngRegion *reg2 = info[kj].adjacent;
	  for (size_t ki=0; ki<adj_reg.size(); ++ki)
	    {
	      if (adj_reg[ki] == reg2)
		{
		  if (info[kj].max_dist < 0.0)
		    {
		      info[kj].min_dist = info[kj].max_dist = adj_info[ki].first;
		      info[kj].min_ang = info[kj].max_ang = adj_info[ki].second;
		    }
		  else
		    {
		      info[kj].min_dist = std::min(info[kj].min_dist,adj_info[ki].first);
		      info[kj].max_dist = std::max(info[kj].max_dist,adj_info[ki].first);
		      info[kj].min_ang = std::min(info[kj].min_ang,adj_info[ki].second);
		      info[kj].max_ang = std::max(info[kj].max_ang,adj_info[ki].second);
		    }
		}
	    }
	}
    }

  kj = 0;
  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it, ++kj)
    {
      double maxd1, maxd2, avd1, avd2;
      maxd1 = maxd2 = avd1 = avd2 = std::numeric_limits<double>::max();
      int nmb_in1 = 0, nmb_in2 = 0;
      bool local_approx = (info[kj].max_dist < 0.0);
      bool outlier = false;
      int nmb_pt_adj;
      bool computed = computeIntegrateInfo(group_points_, *it, tol, angtol, radius,
					   local_approx, min_next, max_next, max_nmb_outlier, 
					   outlier, nmb_pt_adj, maxd2, avd2, nmb_in2, maxd1, 
					   avd1, nmb_in1);
      info[kj].setNmbAdj(nmb_pt_adj);
      if (outlier)
	info[kj].setOutlier();
      if (!computed)
	continue;
      // shared_ptr<ParamSurface> surf;
      // if ((*it)->hasSurface())
      // 	{
      // 	  surf = (*it)->getSurface(0)->surface();
      // 	  (*it)->getAccuracy(maxd1, avd1, nmb_in1);
      // 	}
      // else if ((*it)->hasBaseSf())
      // 	(*it)->getBase(surf, maxd1, avd1, nmb_in1);

      // if (surf.get())
      // 	{
      // 	  // Fetch accuracy of current surface
      // 	  info[kj].setNmbAdj((*it)->numPoints());
	  
      // 	  // Check accuracy of new points
      // 	  vector<RevEngPoint*> in2, out2;
      // 	  vector<pair<double,double> > dist_ang2;
      // 	  vector<double> parvals2;
      // 	  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
      // 				  surf, tol, maxd2, avd2, nmb_in2, in2, out2,
      // 				  parvals2, dist_ang2, angtol);
      // 	}
      // else if (info[kj].max_dist < 0.0)
      // 	{
      // 	  // Fetch nearby points
      // 	  vector<RevEngPoint*> nearpts;
      // 	  if ((*it)->numPoints() < min_next)
      // 	    nearpts.insert(nearpts.end(), (*it)->pointsBegin(), (*it)->pointsEnd());
      // 	  else
      // 	    {
      // 	      for (size_t ki=0; ki<group_points_.size(); ++ki)
      // 		{
      // 		  nearpts.clear();
      // 		  group_points_[ki]->fetchClosePoints2(radius, min_next,
      // 						       max_next-(int)group_points_.size(),
      // 						       nearpts, *it);
      // 		  if (nearpts.size() > max_nmb_outlier)
      // 		    break;
      // 		}
      // 	    }
      // 	  size_t nearnmb = nearpts.size();
      // 	  info[kj].setNmbAdj((int)nearnmb);
      // 	  if (((int)(nearnmb+group_points_.size()) <= max_nmb_outlier) ||
      // 	      ((int)nearnmb <= max_nmb_outlier && (*it)->numPoints() > min_next))
      // 	    {
      // 	      if ((int)group_points_.size() <= max_nmb_outlier)
      // 		info[kj].setOutlier();
      // 	      continue;
      // 	    }

      // 	  nearpts.insert(nearpts.end(), group_points_.begin(), group_points_.end());
      // 	  BoundingBox bbox = bbox_;
      // 	  for (size_t ki=0; ki<nearnmb; ++ki)
      // 	    {
      // 	      Vector3D xyz = nearpts[ki]->getPoint();
      // 	      Point xyz2(xyz[0], xyz[1], xyz[2]);
      // 	      bbox.addUnionWith(xyz2);
      // 	    }
      // 	  surf = surfApprox(nearpts, bbox);

      // 	  // Check accuracy
      // 	  vector<RevEngPoint*> in1, in2, out1, out2;
      // 	  vector<pair<double,double> > dist_ang1, dist_ang2;
      // 	  vector<double> parvals1, parvals2;
      // 	  RevEngUtils::distToSurf(nearpts.begin(), nearpts.begin()+nearnmb, surf,
      // 				  tol, maxd1, avd1, nmb_in1, in1, out1,
      // 				  parvals1, dist_ang1, angtol);
      // 	  RevEngUtils::distToSurf(nearpts.begin()+nearnmb, nearpts.end(), surf,
      // 				  tol, maxd2, avd2, nmb_in2, in2, out2,
      // 				  parvals2, dist_ang2, angtol);
      // 	}
#ifdef DEBUG_INTEGRATE
      int num = (*it)->numPoints();
      of << "400 1 0 4 255 0 0 255" << std::endl;
      of << num << std::endl;
      for (int ka=0; ka<num; ++ka)
      	of << (*it)->getPoint(ka)->getPoint() << std::endl;

      // if (surf.get())
      // 	{
      // 	  surf->writeStandardHeader(of);
      // 	  surf->write(of);
      // 	}
#endif
      info[kj].setAccuracy(maxd2, avd2, nmb_in2, maxd1, avd1, nmb_in1);
    }

  // Select adjacent
  bool outlier = true;
  int ix = -1;
  double score = 0.0;
  int num = (int)group_points_.size();
  for (size_t kj=0; kj<info.size(); ++kj)
    {
      if (!info[kj].outlier)
	outlier = false;
      if (taboo && info[kj].adjacent == taboo /*&& info.size()>1*/)
	continue;
      if (prev_region_ && info[kj].adjacent == prev_region_ && info.size()>1)
	continue;
      if (info[kj].avd > tol || (info[kj].nmb_in < num/2 && info[kj].maxd > tol))
	continue;
      double frac1 = (double)info[kj].nmb_in/(double)num;
      double frac2 = (double)info[kj].nmb_in_adj/(double)info[kj].nmb_pt_adj;
      if (frac1 < 0.9*frac2 && info[kj].maxd > tol)
	continue;
      double avH = avH_*info[kj].adjacent->avH_;
      double avK = avK_*info[kj].adjacent->avK_;
      double curr_score = (tol/std::max(1.0e-6, info[kj].avd)) + /* frac2/frac1 +*/
	(info[kj].adjacent->hasSurface()) + (avH > 0.0) + (avK > 0.0);
      if (curr_score > score)
	{
	  ix = (int)kj;
	  score = curr_score;
	}
    }

  if (ix < 0 && (!outlier))
    {
      double div = std::numeric_limits<double>::max();
      for (size_t kj=0; kj<info.size(); ++kj)
	{
	  if (info[kj].max_dist < 0)
	    continue;
	  if (taboo && info[kj].adjacent == taboo /*&& info.size()>1*/)
	    continue;
	  if (prev_region_ && info[kj].adjacent == prev_region_ && info.size()>1)
	    continue;
	  double curr_div = info[kj].min_dist + info[kj].min_ang;
	  if (curr_div < div)
	    {
	      div = curr_div;
	      ix = (int)kj;
	    }
	}
    }
  
  if (outlier)
    {
     for (size_t ki=0; ki<group_points_.size(); ++ki)
       {
	group_points_[ki]->unsetRegion();
	group_points_[ki]->setOutlier();
	group_points_[ki]->addMove();
       }
     for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
       (*it)->removeAdjacentRegion(this);
     //updateInfo();
     return true;
    }
  else if (ix >= 0)
    {
      if (info[ix].adjacent->hasSurface())
	{
	  // Set parameter value, distance and angle difference
	  shared_ptr<ParamSurface> surf = info[ix].adjacent->getSurface(0)->surface();
	  double maxd2, avd2;
	  int nmb_in2;
	  vector<RevEngPoint*> in2, out2;
	  vector<pair<double,double> > dist_ang2;
	  vector<double> parvals2;
	  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
				  surf, tol, maxd2, avd2, nmb_in2, in2, out2,
				  parvals2, dist_ang2, angtol);
	  for (size_t ki=0; ki<group_points_.size(); ++ki)
	    {
	      group_points_[ki]->setPar(Vector2D(parvals2[2*ki],parvals2[2*ki+1]));
	      group_points_[ki]->setSurfaceDist(dist_ang2[ki].first, dist_ang2[ki].second);
	    }
	}
     for (size_t ki=0; ki<group_points_.size(); ++ki)
       group_points_[ki]->addMove();
      vector<pair<double, double> > dummy;
      info[ix].adjacent->addRegion(this, dummy);
      for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
	(*it)->removeAdjacentRegion(this);
      updateInfo();
      return true;
      }
  return false;
}

//===========================================================================
bool RevEngRegion::computeIntegrateInfo(vector<RevEngPoint*>& points, RevEngRegion *adj_reg,
					double tol, double angtol, double radius,
					bool local_approx, int min_next, int max_next,
					int max_nmb_outlier, bool& outlier, int& nmb_pt_adj,
					double& maxdist, double& avdist, int& nmb_in,
					double& maxdist_adj, double& avdist_adj, int& nmb_in_adj)
//===========================================================================
{
  outlier = false;
  nmb_pt_adj = 0;
  
  shared_ptr<ParamSurface> surf;
  if (adj_reg->hasSurface())
    {
      surf = adj_reg->getSurface(0)->surface();
      adj_reg->getAccuracy(maxdist_adj, avdist_adj, nmb_in_adj);
    }
  else if (adj_reg->hasBaseSf())
    adj_reg->getBase(surf, maxdist_adj, avdist_adj, nmb_in_adj);

  if (surf.get())
    {
      // Fetch accuracy of current surface
      nmb_pt_adj = adj_reg->numPoints();
	  
      // Check accuracy of new points
      vector<RevEngPoint*> in, out;
      vector<pair<double,double> > dist_ang;
      vector<double> parvals;
      RevEngUtils::distToSurf(points.begin(), points.end(),
			      surf, tol, maxdist, avdist, nmb_in, in, out,
			      parvals, dist_ang, angtol);
    }
  else if (local_approx)
    {
      // Fetch nearby points
      vector<RevEngPoint*> nearpts;
      if (adj_reg->numPoints() < min_next)
	nearpts.insert(nearpts.end(), adj_reg->pointsBegin(), adj_reg->pointsEnd());
      else
	{
	  for (size_t ki=0; ki<points.size(); ++ki)
	    {
	      nearpts.clear();
	      points[ki]->fetchClosePoints2(radius, min_next,
					    max_next-(int)points.size(),
					    nearpts, adj_reg);
	      if (nearpts.size() > max_nmb_outlier)
		break;
	    }
	}
      size_t nearnmb = nearpts.size();
      nmb_pt_adj = (int)nearnmb;
      if (((int)(nearnmb+points.size()) <= max_nmb_outlier) ||
	  ((int)nearnmb <= max_nmb_outlier && adj_reg->numPoints() > min_next))
	{
	  if ((int)points.size() <= max_nmb_outlier)
	    outlier = true;
	  return false;
	}

      nearpts.insert(nearpts.end(), points.begin(), points.end());
      BoundingBox bbox = bbox_;
      for (size_t ki=0; ki<nearnmb; ++ki)
	{
	  Vector3D xyz = nearpts[ki]->getPoint();
	  Point xyz2(xyz[0], xyz[1], xyz[2]);
	  bbox.addUnionWith(xyz2);
	}
      surf = surfApprox(nearpts, bbox);

      // Check accuracy
      vector<RevEngPoint*> in1, in2, out1, out2;
      vector<pair<double,double> > dist_ang1, dist_ang2;
      vector<double> parvals1, parvals2;
      RevEngUtils::distToSurf(nearpts.begin(), nearpts.begin()+nearnmb, surf,
			      tol, maxdist_adj, avdist_adj, nmb_in_adj, in1, out1,
			      parvals1, dist_ang1, angtol);
      RevEngUtils::distToSurf(nearpts.begin()+nearnmb, nearpts.end(), surf,
			      tol, maxdist, avdist, nmb_in, in2, out2,
			      parvals2, dist_ang2, angtol);
    }

  return true;
}

//===========================================================================
bool RevEngRegion::adjustWithCylinder(double tol, double angtol, int min_point_region,
				      vector<vector<RevEngPoint*> >& out_groups,
				      vector<RevEngRegion*>& grown_regions,
				      vector<HedgeSurface*>& adj_surfs)
//===========================================================================
{
  if (!hasSurface())
    return false;

  HedgeSurface *hedge = getSurface(0);
  int code;
  ClassType type = hedge->instanceType(code);
  if (type != Class_Cylinder &&
      (!(type == Class_SplineSurface && code == LINEARSWEPT_SURF)))
    return false;
  shared_ptr<ParamSurface> surf = hedge->surface();

  std::ofstream of("cylinder_adjust.g2");
  writeRegionInfo(of);
  
  // Make bounding parameter domain
  int ixd = (type == Class_Cylinder) ? 0 : 2;
  //int ixd2 = (ixd == 0) ? 2 : 0;
  int ixp = ixd/2;
  double dom[4];
  Vector2D uv = group_points_[0]->getPar();
  dom[0] = dom[1] = uv[0];
  dom[2] = dom[3] = uv[1];
  double dist, ang, avang;
  double fac = 1.0/(double)group_points_.size();
  group_points_[0]->getSurfaceDist(dist, ang);
  avang = fac*ang;
  for (size_t ki=1; ki<group_points_.size(); ++ki)
    {
      Vector2D uv = group_points_[ki]->getPar();
      group_points_[ki]->getSurfaceDist(dist, ang);
      dom[0] = std::min(dom[0], uv[0]);
      dom[1] = std::max(dom[1], uv[0]);
      dom[2] = std::min(dom[2], uv[1]);
      dom[3] = std::max(dom[3], uv[1]);
      avang += fac*ang;
    }

  std::ofstream ofp("projected_pts.g2");
  Point axis, pnt;
  shared_ptr<ParamCurve> crv;
  if (type == Class_Cylinder)
    {
      shared_ptr<Cylinder> cyl = dynamic_pointer_cast<Cylinder,ParamSurface>(surf);
      axis = cyl->getAxis();
      pnt = cyl->getLocation();
      shared_ptr<Circle> circ(new Circle(cyl->radius(0,0), pnt, axis, cyl->direction2()));
      crv = circ;
    }
  else
    {
      vector<Point> der(3);
      double upar = 0.5*(dom[0]+dom[1]);
      double vpar = 0.5*(dom[2]+dom[3]);
      surf->point(der, upar, vpar, 1);
      axis = der[1];
      axis.normalize();
      pnt = der[0];
      vector<shared_ptr<ParamCurve> > cvs = surf->constParamCurves(vpar, true);
      crv = cvs[0];
    }
  
  vector<Point> projected;
  double maxdp, avdp;
  RevEngUtils::projectToPlane(group_points_, axis, pnt, projected, maxdp, avdp);
  ofp << "400 1 0 4 255 0 0 255" << std::endl;
  ofp << projected.size() << std::endl;
  for (size_t kr=0; kr<projected.size(); ++kr)
    ofp << projected[kr] << std::endl;
  crv->writeStandardHeader(ofp);
  crv->write(ofp);
  
  // Reduce domain from the start
  vector<RevEngPoint*> outer;
  double del = 0.02;
  double del2 = del*(dom[ixd+1] - dom[ixd]);
  double par;
  int knmb = 10;
  int ka;
  double dfac = 2.0;
  //double afac = 2.0;
  double pfac = 0.5;
  int part = (int)(del*(double)group_points_.size());
  for (ka=0, par=dom[ixd]+del2; ka<knmb; ++ka, par+=del2)
    {
      vector<RevEngPoint*> curr_pts;
      double avd = 0.0, ava = 0.0;
      int nn = 0;
      for (size_t ki=0; ki<group_points_.size(); ++ki)
	{
	  Vector2D uv = group_points_[ki]->getPar();
	  if (uv[ixp] >= par-del2 && uv[ixp]<par)
	    {
	      curr_pts.push_back(group_points_[ki]);
	      group_points_[ki]->getSurfaceDist(dist, ang);
	      avd += dist;
	      ava += ang;
	      ++nn;
	    }
	}
      avd /= (double)nn;
      ava /= (double)nn;

      // Check with adjacent regions
      std::ofstream ofo("outer_cand.g2");
      ofo << "400 1 0 4 0 255 0 255" << std::endl;
      ofo << curr_pts.size() << std::endl;
     for (size_t ki=0; ki<curr_pts.size(); ++ki)
       ofo << curr_pts[ki]->getPoint() << std::endl;
     
      vector<RevEngRegion*> next_reg;
      vector<int> nmb_next;
      for (size_t ki=0; ki<curr_pts.size(); ++ki)
	{
	  vector<RevEngRegion*> adjr;
	  curr_pts[ki]->adjacentRegions(adjr);
	  for (size_t kj=0; kj<adjr.size(); ++kj)
	    {
	      size_t kr=0;
	      for (kr=0; kr<next_reg.size(); ++kr)
		if (next_reg[kr] == adjr[kj])
		  break;
	      if (kr == next_reg.size())
		{
		  next_reg.push_back(adjr[kj]);
		  nmb_next.push_back(1);
		}
	      else
		nmb_next[kr]++;
	    }
	}

      
      if (avd > dfac*avdist_ && (ava > dfac*avang || nn < pfac*part))
	{
	  outer.insert(outer.end(), curr_pts.begin(), curr_pts.end());
	  dom[ixd] = par;
	}
      else
	break;
      int stop_break = 1;
    }

  // Reduce domain from the end
  for (ka=0, par=dom[ixd+1]-del2; ka<knmb; ++ka, par-=del2)
    {
      vector<RevEngPoint*> curr_pts;
       double avd = 0.0, ava = 0.0;
      int nn = 0;
      for (size_t ki=0; ki<group_points_.size(); ++ki)
	{
	  Vector2D uv = group_points_[ki]->getPar();
	  if (uv[ixp] > par && uv[ixp]<=par+del2)
	    {
	      curr_pts.push_back(group_points_[ki]);
	      group_points_[ki]->getSurfaceDist(dist, ang);
	      avd += dist;
	      ava += ang;
	      ++nn;
	    }
	}
      avd /= (double)nn;
      ava /= (double)nn;
      // Check with adjacent regions
      std::ofstream ofo("outer_cand.g2");
      ofo << "400 1 0 4 0 255 0 255" << std::endl;
      ofo << curr_pts.size() << std::endl;
     for (size_t ki=0; ki<curr_pts.size(); ++ki)
       ofo << curr_pts[ki]->getPoint() << std::endl;
     
      vector<RevEngRegion*> next_reg;
      vector<int> nmb_next;
      for (size_t ki=0; ki<curr_pts.size(); ++ki)
	{
	  vector<RevEngRegion*> adjr;
	  curr_pts[ki]->adjacentRegions(adjr);
	  for (size_t kj=0; kj<adjr.size(); ++kj)
	    {
	      size_t kr=0;
	      for (kr=0; kr<next_reg.size(); ++kr)
		if (next_reg[kr] == adjr[kj])
		  break;
	      if (kr == next_reg.size())
		{
		  next_reg.push_back(adjr[kj]);
		  nmb_next.push_back(1);
		}
	      else
		nmb_next[kr]++;
	    }
	}

      
      if (avd > dfac*avdist_ && (ava > dfac*avang || nn < pfac*part))
	{
	  outer.insert(outer.end(), curr_pts.begin(), curr_pts.end());
	  dom[ixd+1] = par;
	}
      else
	break;
      int stop_break = 1;
    }

  // Integrate points from adjacent regions
  vector<vector<RevEngPoint*> > adjpts;
  vector<RevEngRegion*> adjreg;
  vector<vector<double> > par_and_dist;
  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
    {
      std::ofstream of2("adj_group.g2");
      int num = (*it)->numPoints();
      of2 << "400 1 0 4 255 0 0 255" << std::endl;
      of2 << num << std::endl;
      for (int ka=0; ka<num; ++ka)
      	of2 << (*it)->getPoint(ka)->getPoint() << std::endl;

      vector<Point> projected2;
      double maxdp2, avdp2;
      vector<RevEngPoint*> curr_pts = (*it)->getPoints();
      RevEngUtils::projectToPlane(curr_pts, axis, pnt, projected2, maxdp2, avdp2);
      std::ofstream ofp2("projected_pts_adj.g2");
      ofp2 << "400 1 0 4 100 155 0 255" << std::endl;
      ofp2 << projected2.size() << std::endl;
      for (size_t kr=0; kr<projected2.size(); ++kr)
	ofp2 << projected2[kr] << std::endl;

      if ((*it)->hasSurface())
	{
	  shared_ptr<ParamSurface> surf2((*it)->getSurface(0)->surface()->clone());
	  double upar2, vpar2, dist;
	  Point close;
	  surf2->closestPoint(pnt, upar2, vpar2, close, dist, tol);
	  if (!surf2->isBounded())
	    {
	      BoundingBox bb = getBbox();
	      double len = bb.low().dist(bb.high());
	      shared_ptr<Cylinder> elem1 =
		dynamic_pointer_cast<Cylinder,ParamSurface>(surf2);
	      shared_ptr<Plane> elem2 =
		dynamic_pointer_cast<Plane,ParamSurface>(surf2);
	      if (elem1.get())
		elem1->setParamBoundsV(-len, len);
	      else if (elem2.get())
		elem2->setParameterBounds(-len, -len, len, len);
	    }

	  if (surf2->isBounded())
	    {
	      vector<shared_ptr<ParamCurve> > cvs2_1 = surf2->constParamCurves(upar2, false);
	      vector<shared_ptr<ParamCurve> > cvs2_2 = surf2->constParamCurves(vpar2, true);
	      cvs2_1[0]->writeStandardHeader(ofp2);
	      cvs2_1[0]->write(ofp2);
	      cvs2_2[0]->writeStandardHeader(ofp2);
	      cvs2_2[0]->write(ofp2);
	    }
	}

      vector<RevEngPoint*> curr_adjpts;
       vector<double> curr_par_and_dist;
       double avd, ava;
       int nn;
       getAdjInsideDist(surf, dom, tol, *it, avd, ava, nn, curr_adjpts, curr_par_and_dist);

      if (avd < dfac*avdist_ /*&& (ava < dfac*avang || nn == num)*/)
	{
	  adjpts.push_back(curr_adjpts);
	  adjreg.push_back(*it);
	  par_and_dist.push_back(curr_par_and_dist);
	}
     }

  std::ofstream of2("move2pts.g2");
  for (size_t ki=0; ki<adjpts.size(); ++ki)
    {
      of2 << "400 1 0 4 75 75 75 255" << std::endl;
      of2 << adjpts[ki].size() << std::endl;
      for (size_t kh=0; kh<adjpts[ki].size(); ++kh)
	of2 << adjpts[ki][kh]->getPoint() << std::endl;
    }

  // Adjust point groups
  if (outer.size() > 0)
    {
      extractSpesPoints(outer, out_groups);
      splitRegion(out_groups);
      }

  for (size_t ki=0; ki<adjpts.size(); ++ki)
    {
      for (size_t kh=0; kh<adjpts[ki].size(); ++kh)
	{
	  adjpts[ki][kh]->setMoved();
	  adjreg[ki]->removePoint(adjpts[ki][kh]);
	  adjpts[ki][kh]->setRegion(this);
	  adjpts[ki][kh]->setPar(Vector2D(par_and_dist[ki][4*kh],
					  par_and_dist[ki][4*kh+1]));
	  adjpts[ki][kh]->setSurfaceDist(par_and_dist[ki][4*kh+2], par_and_dist[ki][4*kh+3]);
	  addPoint(adjpts[ki][kh]);
	}

      if (adjreg[ki]->numPoints() == 0)
	{
	  for (auto it=adjreg[ki]->adjacent_regions_.begin();
	       it!=adjreg[ki]->adjacent_regions_.end(); ++it)
	    {
	      if (*it != this)
		{
		  (*it)->addAdjacentRegion(this);
		  (*it)->removeAdjacentRegion(adjreg[ki]);
		}
	    }
	  grown_regions.push_back(adjreg[ki]);
	  int num_sf = adjreg[ki]->numSurface();
	  for (int kb=0; kb<num_sf; ++kb)
	    adj_surfs.push_back(adjreg[ki]->getSurface(kb));
	      
	  removeAdjacentRegion(adjreg[ki]);
	}
      else
	{
	  if (adjreg[ki]->hasSurface())
	    {
	      if (adjreg[ki]->numPoints() >= min_point_region)
		adjreg[ki]->checkReplaceSurf(tol);
	      else
		{
		  int num_sf = adjreg[ki]->numSurface();
		  for (int kb=0; kb<num_sf; ++kb)
		    adj_surfs.push_back(adjreg[ki]->getSurface(kb));
		  adjreg[ki]->clearSurface();
		}
	    }
	  else
	    adjreg[ki]->updateInfo();
	}
    }
  
  //updateInfo();
  checkReplaceSurf(tol);

   
  return (outer.size() > 0 || adjpts.size() > 0);
}

//===========================================================================
void RevEngRegion::getAdjInsideDist(shared_ptr<ParamSurface> surf, double dom[4],
				    double tol, RevEngRegion* reg,
				    double& avd, double& ava, int& nn,
				    vector<RevEngPoint*>& adjpts,
				    vector<double>& par_and_dist)
//===========================================================================
{
  double maxdd, avdd;
  int num_inside;
  vector<pair<double, double> > dist_ang;
  vector<double> parvals;
  vector<RevEngPoint*> in, out;
  RevEngUtils::distToSurf(reg->pointsBegin(),
			  reg->pointsEnd(), surf, tol, maxdd, avdd, 
			  num_inside, in, out, parvals, dist_ang);
  avd = 0.0;
  ava = 0.0;
  nn = 0;
  for (int kh=0; kh<(int)dist_ang.size(); ++kh)
    {
      if (parvals[2*kh] >= dom[0] && parvals[2*kh] <= dom[1] &&
	  parvals[2*kh+1] >= dom[2] && parvals[2*kh+1] <= dom[3])
	{
	  adjpts.push_back(reg->getPoint(kh));
	  par_and_dist.push_back(parvals[2*kh]);
	  par_and_dist.push_back(parvals[2*kh+1]);
	  par_and_dist.push_back(dist_ang[kh].first);
	  par_and_dist.push_back(dist_ang[kh].second);
	  avd += dist_ang[kh].first;
	  ava += dist_ang[kh].second;
	  ++nn;
	}
    }
  if (nn > 0)
    {
      avd /= (double)nn;
      ava /= (double)nn;
    }
}

//===========================================================================
void RevEngRegion::addRegion(RevEngRegion* reg,
			     vector<pair<double, double> >& dist_ang,
			     double maxd, double avd, int num_inside)
//===========================================================================
{
  int num = reg->numPoints();
  bbox_.addUnionWith(reg->boundingBox());
  normalcone_.addUnionWith(reg->getNormalCone());
  if (num_inside >= 0)
    {
      maxdist_ = std::max(maxdist_, maxd);
      double div = (double)((int)group_points_.size() + num);
      avdist_ = ((double)(group_points_.size())*avdist_ + num*avd)/div;
      num_inside_ += num_inside;
    }

  double mink1, maxk1, mink2, maxk2;
  reg->getPrincipalCurvatureInfo(mink1, maxk1, mink2, maxk2);
  mink1_ = std::min(mink1_, mink1);
  maxk1_ = std::max(maxk1_, maxk1);
  mink2_ = std::min(mink2_, mink2);
  maxk2_ = std::max(maxk2_, maxk2);

  int num_all = numPoints() + reg->numPoints();
  double fac1 = (double)numPoints()/(double)num_all;
  double fac2 = (double)reg->numPoints()/(double)num_all;
  double avH, avK, MAH, MAK;
  reg->getAvCurvatureInfo(avH, avK, MAH, MAK);
  Point avnorm = reg->getMeanNormal();
  avH_ = fac1*avH_ + fac2*avH;
  avK_ = fac1*avK_ + fac2*avK;
  MAH_ = fac1*MAH_ + fac2*MAH;
  MAK_ = fac1*MAK_ + fac2*MAK;
  avnorm_ = fac1*avnorm_ + fac2*avnorm;

  size_t kr=0;
  for (auto it=reg->pointsBegin(); it != reg->pointsEnd(); ++it, ++kr)
    {
      if (dist_ang.size() > 0)
	(*it)->setSurfaceDist(dist_ang[kr].first, dist_ang[kr].second);
      (*it)->setRegion(this);
    }
  group_points_.insert(group_points_.end(), reg->pointsBegin(),
		       reg->pointsEnd());

  if (hasSurface())
    computeDomain();
  
}


//===========================================================================
void RevEngRegion::growWithSurf(int max_nmb, double tol,
				vector<RevEngRegion*>& grown_regions,
				vector<HedgeSurface*>& adj_surfs)
//===========================================================================
{
  //double eps = 1.0e-6;
  if (associated_sf_.size() == 0)
    return;  // No surface with which to check growt

  int sfcode;
  ClassType classtype = associated_sf_[0]->instanceType(sfcode);
  shared_ptr<ParamSurface> surf = associated_sf_[0]->surface();
  shared_ptr<BoundedSurface> bdsurf =
    dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf);
  if (bdsurf.get())
    surf = bdsurf->underlyingSurface();

  shared_ptr<ParamSurface> primary;
  if (primary_.get() && avdist_primary_ < tol &&
      num_in_primary_ > (int)group_points_.size()/2)
    primary = primary_;

  int write_extend = 1;
  bool changed = false;
  vector<RevEngRegion*> adj_reg;
  adj_reg.insert(adj_reg.end(), adjacent_regions_.begin(), adjacent_regions_.end());
  for (size_t ki=0; ki<adj_reg.size(); ++ki)
    for (size_t kj=ki+1; kj<adj_reg.size(); ++kj)
      if (adj_reg[kj]->numPoints() > adj_reg[ki]->numPoints())
	std::swap(adj_reg[ki], adj_reg[kj]);

  for (size_t ki=0; ki<adj_reg.size(); ++ki)
    {
      if (adj_reg[ki]->prev_region_ && adj_reg[ki]->prev_region_ == this)
	{
	  continue;
	}

      int num = adj_reg[ki]->numPoints();
      
      std::ofstream of("curr_extend.g2");
       if (write_extend)
	{
	  of << "400 1 0 4 0 255 0 255" << std::endl;
	  of << group_points_.size() << std::endl;
	  for (size_t kh=0; kh<group_points_.size(); ++kh)
	    of << group_points_[kh]->getPoint() << std::endl;
      
	  of << "400 1 0 4 255 0 0 255" << std::endl;
	  of << num << std::endl;
	  for (int ka=0; ka<num; ++ka)
	    of << adj_reg[ki]->getPoint(ka)->getPoint() << std::endl;
	  surf->writeStandardHeader(of);
	  surf->write(of);
	  if (primary)
	    {
	      primary->writeStandardHeader(of);
	      primary->write(of);
	    }
	}
      
       if (adj_reg[ki]->hasSurface())
	 {
	  double anglim = 0.1;
	  double score;
	  bool compatible = associated_sf_[0]->isCompatible(adj_reg[ki]->getSurface(0),
							    anglim, tol, classtype,
							    score);
	  int sfcode2;
	  ClassType classtype2 = adj_reg[ki]->getSurface(0)->instanceType(sfcode);
	  bool check_pair = false;
	  //if (check_pair)
	  if (compatible || (classtype == classtype2 && adj_reg[ki]->numAdjacentRegions() == 1))
	    {
	      BoundingBox bb = bbox_;
	      bb.addUnionWith(adj_reg[ki]->getBbox());
	      vector<pair<vector<RevEngPoint*>::iterator,vector<RevEngPoint*>::iterator> > points(2);
	      points[0] = make_pair(group_points_.begin(), group_points_.end());
	      points[1] = make_pair(adj_reg[ki]->pointsBegin(), adj_reg[ki]->pointsEnd());
	      vector<int> nmbpts(2);
	      nmbpts[0] = numPoints();
	      nmbpts[1] = adj_reg[ki]->numPoints();
	      shared_ptr<ParamSurface> merged;
	      if (classtype == Class_Cylinder)
		merged = RevEngUtils::doMergeCylinders(points, bb, nmbpts, false);
	      else if (classtype == Class_Plane)
		merged = RevEngUtils::doMergePlanes(points, bb, nmbpts, false);
	      if (merged.get())
		{
		  double maxd1, maxd2, avd1, avd2;
		  int num_in1, num_in2;
		  vector<RevEngPoint*> in1, out1, in2, out2;
		  vector<pair<double, double> > dist_ang1, dist_ang2;
		  vector<double> parvals1, parvals2;
		  RevEngUtils::distToSurf(points[0].first, points[0].second, merged, tol,
					  maxd1, avd1, num_in1, in1, out1, parvals1,
					  dist_ang1);
		  RevEngUtils::distToSurf(points[1].first, points[1].second, merged, tol,
					  maxd2, avd2, num_in2, in2, out2, parvals2,
					  dist_ang2);
		  if (write_extend)
		    {
		      merged->writeStandardHeader(of);
		      merged->write(of);
		    }

		  if (num_in1 > nmbpts[0]/2 && avd1 <= tol && num_in2 >= nmbpts[1]/2 &&
		      avd2 <= tol)
		    {
		      // Merge regions and replace surface
		      // Update current
		      for (size_t kr=0; kr<group_points_.size(); ++kr)
			{
			  group_points_[kr]->setPar(Vector2D(parvals1[2*kr],parvals1[2*kr+1]));
			  group_points_[kr]->setSurfaceDist(dist_ang1[kr].first, dist_ang1[kr].second);
			}
		      setAccuracy(maxd1, avd1, num_in1);
		      associated_sf_[0]->replaceSurf(merged);

		      // Include adjacent region in present
		      vector<RevEngRegion*> added_adjacent;
		      includeAdjacentRegion(adj_reg[ki], maxd2, avd2, num_in2, parvals2,
					    dist_ang2, added_adjacent);
		      for (size_t kr=0; kr<added_adjacent.size(); ++kr)
			{
			  if (std::find(adj_reg.begin(), adj_reg.end(), added_adjacent[kr]) != adj_reg.end())
			    {
			      size_t kj;
			      for (kj=ki+1; kj<adj_reg.size(); ++kj)
				if (adj_reg[kj]->numPoints() < added_adjacent[kr]->numPoints())
				  {
				    adj_reg.insert(adj_reg.begin()+kj, added_adjacent[kr]);
				    break;
				  }
			      if (kj == adj_reg.size())
				adj_reg.push_back(added_adjacent[kr]);
			    }
			}
		      grown_regions.push_back(adj_reg[ki]);
		      
		      int num_sf = adj_reg[ki]->numSurface();
		      for (int kb=0; kb<num_sf; ++kb)
			adj_surfs.push_back(adj_reg[ki]->getSurface(kb));
		      
		      removeAdjacentRegion(adj_reg[ki]);
		    }
		  int stop_compare0 = 1;
		  continue;
		}
	    }
	  int stop_compare = 1;
	}
      

	      
      // Check accuracy
      double maxd, avd;
      int num_inside;
      double maxd_init, avd_init;
      int num_inside_init;
      adj_reg[ki]->getAccuracy(maxd_init, avd_init, num_inside_init);

      shared_ptr<ParamSurface> curr_sf;
      int ka;
      bool curr_changed = false;
      for (curr_sf=primary, ka=0; ka<2; curr_sf=surf, ++ka)
	{
	  if (!curr_sf.get())
	    continue;
	  vector<RevEngPoint*> in, out;
	  vector<pair<double, double> > dist_ang;
	  vector<double> parvals;
	  RevEngUtils::distToSurf(adj_reg[ki]->pointsBegin(),
				  adj_reg[ki]->pointsEnd(), curr_sf, tol,
				  maxd, avd, num_inside, in, out, parvals, dist_ang);
	  vector<Vector3D> better1, worse1;
	  for (size_t kh=0; kh<dist_ang.size(); ++kh)
	    {
	      double ptdist, ptang;
	      adj_reg[ki]->group_points_[kh]->getSurfaceDist(ptdist, ptang);
	      if (dist_ang[kh].first <= ptdist && dist_ang[kh].second < ptang)
		better1.push_back(adj_reg[ki]->group_points_[kh]->getPoint());
	      else
		worse1.push_back(adj_reg[ki]->group_points_[kh]->getPoint());
	    }
	  if (write_extend)
	    {
	  std::ofstream ofc1("better_worse1.g2");
	  ofc1 << "400 1 0 4 155 50 50 255" << std::endl;
	  ofc1 << better1.size() << std::endl;
	  for (size_t kr=0; kr<better1.size(); ++kr)
	    ofc1 << better1[kr] << std::endl;
	  ofc1 << "400 1 0 4 50 155 50 255" << std::endl;
	  ofc1 << worse1.size() << std::endl;
	  for (size_t kr=0; kr<worse1.size(); ++kr)
	    ofc1 << worse1[kr] << std::endl;
      
      
	  std::ofstream ofd("in_out_extend.g2");
	  ofd << "400 1 0 4 155 50 50 255" << std::endl;
	  ofd << in.size() << std::endl;
	  for (size_t kr=0; kr<in.size(); ++kr)
	    ofd << in[kr]->getPoint() << std::endl;
	  ofd << "400 1 0 4 50 155 50 255" << std::endl;
	  ofd << out.size() << std::endl;
	  for (size_t kr=0; kr<out.size(); ++kr)
	    ofd << out[kr]->getPoint() << std::endl;
	    }
      
	  if ((adj_reg[ki]->hasSurface() && num_inside > std::max(3*num_inside_init/4, 2*num/3) &&
	       avd < tol /*avd_init*/) ||
	      ((!adj_reg[ki]->hasSurface()) && num_inside > 2*num/3 && avd < tol))
	    {
	      // Criteria must be updated
	      changed = true;
	      curr_changed = true;
	      
	      // Include adjacent region in present
	      vector<RevEngRegion*> added_adjacent;
	      includeAdjacentRegion(adj_reg[ki], maxd, avd, num_inside, parvals, dist_ang,
				    added_adjacent);
	      for (size_t kr=0; kr<added_adjacent.size(); ++kr)
		{
		  if (std::find(adj_reg.begin(), adj_reg.end(), added_adjacent[kr]) != adj_reg.end())
		    {
		      size_t kj;
		      for (kj=ki+1; kj<adj_reg.size(); ++kj)
			if (adj_reg[kj]->numPoints() < added_adjacent[kr]->numPoints())
			  {
			    adj_reg.insert(adj_reg.begin()+kj, added_adjacent[kr]);
			    break;
			  }
		      if (kj == adj_reg.size())
			adj_reg.push_back(added_adjacent[kr]);
		    }
		}
	      grown_regions.push_back(adj_reg[ki]);

	      int num_sf = adj_reg[ki]->numSurface();
	      for (int kb=0; kb<num_sf; ++kb)
		adj_surfs.push_back(adj_reg[ki]->getSurface(kb));
	      
	      removeAdjacentRegion(adj_reg[ki]);

	    }
	  if (curr_changed)
	    break;
	}
    }
  

  // Check if the surface should be updated
  if (changed)
    checkReplaceSurf(tol);
}

//===========================================================================
void RevEngRegion::mergeAdjacentSimilar(double tol,
					vector<RevEngRegion*>& grown_regions,
					vector<HedgeSurface*>& adj_surfs)
//===========================================================================
{
  if (associated_sf_.size() == 0)
    return;  // No surface with which to check growt

  int sfcode;
  ClassType classtype = associated_sf_[0]->instanceType(sfcode);
  HedgeSurface *hedge = associated_sf_[0];
  vector<RevEngRegion*> adj_reg;
  vector<double> score;
  double anglim = 0.1;
  double frac2 = 0.75;
  double frac3 = 2.0;
  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
    {
      if ((*it)->hasSurface())
	{
	  double curr_score;
	  bool compatible = associated_sf_[0]->isCompatible((*it)->getSurface(0),
							    anglim, tol,
							    classtype, curr_score);
	  if (compatible)
	    {
	      adj_reg.push_back(*it);
	      score.push_back(curr_score);
	    }
	}
    }

  if (adj_reg.size() == 0)
    return; // Nothing with which to merge

  for (size_t ki=0; ki<adj_reg.size(); ++ki)
    for (size_t kj=ki+1; kj<adj_reg.size(); ++kj)
      if (score[kj] > score[ki])
	{
	  std::swap(score[ki], score[kj]);
	  std::swap(adj_reg[ki], adj_reg[kj]);
	}

  shared_ptr<ParamSurface> surf;
  int ka;
  double maxdist=0.0, avdist=0.0;
  int num_in = 0;
  vector<double> maxd(adj_reg.size()+1, 0.0);
  vector<double> avd(adj_reg.size()+1, 0.0);
  vector<int> ninside(adj_reg.size()+1, 0);
  vector<vector<double> > parvals(adj_reg.size()+1);
  vector<vector<pair<double,double> > > dist_ang(adj_reg.size()+1);
  for (ka=adj_reg.size(); ka>=1; --ka)
    {
      // Create surface from combined point set
      vector<pair<vector<RevEngPoint*>::iterator,
		  vector<RevEngPoint*>::iterator> > points;
      BoundingBox bbox(3);
      vector<int> nmbpts;

      nmbpts.push_back(numPoints());
      points.push_back(std::make_pair(pointsBegin(), pointsEnd()));
      bbox = boundingBox();
      
      for (int kb=0; kb<ka; ++kb)
	{
	  nmbpts.push_back(adj_reg[kb]->numPoints());
	  points.push_back(std::make_pair(adj_reg[kb]->pointsBegin(),
					  adj_reg[kb]->pointsEnd()));
	  bbox.addUnionWith(adj_reg[kb]->boundingBox());
	}
      
      int num_all = 0;
      for (int kb=0; kb<=ka; ++kb)
	num_all += nmbpts[kb];
      double frac = 1.0/(double)num_all;
      if (classtype == Class_Plane)
	{
	  surf = RevEngUtils::doMergePlanes(points, bbox, nmbpts);
	}
      else if (classtype == Class_Cylinder)
	{
	  surf = RevEngUtils::doMergeCylinders(points, bbox, nmbpts);
	}
      else if (classtype == Class_Sphere)
	{
	  Point normal = frac*numPoints()*avnorm_;
	  
	  for (int kb=1; kb<=ka; ++kb)
	    normal += frac*adj_reg[kb-1]->numPoints()*adj_reg[kb-1]->getMeanNormal();
	  normal.normalize_checked();
	  surf = RevEngUtils::doMergeSpheres(points, bbox, nmbpts, normal);
	}
      else if (classtype == Class_Torus)
	{
	  surf = RevEngUtils::doMergeTorus(points, bbox, nmbpts);
	}
      if (!surf.get())
	continue;

      // Check accuracy
       std::ofstream of("in_out_adj.g2");
       maxdist = avdist = 0.0;
       num_in = 0;
      for (int kb=0; kb<=ka; ++kb)
	{
	  maxd[kb] = avd[kb] = 0.0;
	  ninside[kb] = 0;
	  parvals[kb].clear();
	  dist_ang[kb].clear();
	  vector<RevEngPoint*> in, out;
	  RevEngUtils::distToSurf(points[kb].first, points[kb].second, surf, tol,
				  maxd[kb], avd[kb], ninside[kb], in, out, parvals[kb],
				  dist_ang[kb]);
	  maxdist = std::max(maxdist, maxd[kb]);
	  avdist += frac*nmbpts[kb]*avd[kb];
	  num_in += ninside[kb];

	  of << "400 1 0 4 255 0 0 255" << std::endl;
	  of << in.size() << std::endl;
	  for (size_t kj=0; kj<in.size(); ++kj)
	    of << in[kj]->getPoint() << std::endl;
	  of << "400 1 0 4 0 255 0 255" << std::endl;
	  of << out.size() << std::endl;
	  for (size_t kj=0; kj<out.size(); ++kj)
	    of << out[kj]->getPoint() << std::endl;
	}

      int num = numPoints();
      double init_max = frac*num*maxdist_;
      double init_av = frac*num*avdist_;
      double init_in = frac*num_inside_;
      for (int kb=0; kb<ka; ++kb)
	{
	  int num2 = adj_reg[kb]->numPoints();
	  double max2, av2;
	  int num_in2;
	  adj_reg[kb]->getAccuracy(max2, av2, num_in2);
	  init_max += frac*num2*max2;
	  init_av += frac*num2*av2;
	  init_in += frac*num_in2;
	}
      if (num_in > num_all/2 && avdist < tol &&
	  frac*num_in > frac2*init_in && avdist < frac3*init_av)
	break;

      // Swap adjacent regions to skip the least accurate region
      if (ka > 1 && avd[1] > avd[ka])  // The test should be made more accurate
	{
	  std::swap(adj_reg[0], adj_reg[ka-1]);
	  std::swap(score[0], score[ka-1]);
	}
    }

  if (!surf.get())
    return;
  if (ka >= 1)
    {
      setAccuracy(maxd[0], avd[0], ninside[0]);
      for (size_t ki=0; ki<group_points_.size(); ++ki)
	{
	  group_points_[ki]->setPar(Vector2D(parvals[0][2*ki],parvals[0][2*ki+1]));
	  group_points_[ki]->setSurfaceDist(dist_ang[0][ki].first, dist_ang[0][ki].second);
	}
      
      for (int kb=0; kb<ka; ++kb)
	{
	  vector<RevEngRegion*> added_adjacent;
	  includeAdjacentRegion(adj_reg[kb], maxd[kb+1], avd[kb+1], ninside[kb+1],
				parvals[kb+1], dist_ang[kb+1], added_adjacent);
	  grown_regions.push_back(adj_reg[kb]);
	  int num_sf = adj_reg[kb]->numSurface();
	  for (int kc=0; kc<num_sf; ++kc)
	    adj_surfs.push_back(adj_reg[kb]->getSurface(kc));
	  removeAdjacentRegion(adj_reg[kb]);
	}
      updateInfo();
      associated_sf_[0]->replaceSurf(surf);
      computeDomain();
    }
  
}

//===========================================================================
void RevEngRegion::includeAdjacentRegion(RevEngRegion* reg, double maxd, double avd,
			   int num_inside, vector<double>& parvals,
			   vector<pair<double, double> >& dist_ang,
			   vector<RevEngRegion*>& added_adjacent)
//===========================================================================
{
  // First update parameter values
  size_t kr=0;
  for (auto it1=reg->pointsBegin(); it1!=reg->pointsEnd(); ++it1, kr+=2)
    {
      (*it1)->addMove();
      (*it1)->setPar(Vector2D(parvals[kr],parvals[kr+1]));
    }
  addRegion(reg, dist_ang, maxd, avd, num_inside);


  // Update adjacent regions
  for (auto it2=reg->adjacent_regions_.begin(); it2!=reg->adjacent_regions_.end(); ++it2)
    {
      if (*it2 != this)
	{
	  size_t nmb_adj = adjacent_regions_.size();
	  addAdjacentRegion(*it2);
	  if (adjacent_regions_.size() > nmb_adj)
	    {
	      added_adjacent.push_back(*it2);
	      (*it2)->addAdjacentRegion(this);
	    }
	  (*it2)->removeAdjacentRegion(reg);
	}
    }
}

//===========================================================================
void RevEngRegion::adjustBoundaries(double mean_edge_len, double tol, double angtol)
//===========================================================================
{
  std::ofstream of1("curr_regions_adjust.g2");
  writeRegionInfo(of1);

  int min_nmb = 100;
  size_t min_bd = 10;
  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
    {
      if (!(*it)->hasSurface())
	continue;
      if ((*it)->numPoints() < min_nmb)
	continue;

      std::ofstream of2("adj_regions_adjust.g2");
      (*it)->writeRegionInfo(of2);

      // Get boundary points
      vector<RevEngPoint*> adj_pts1 = extractNextToAdjacent(*it);
      vector<RevEngPoint*> adj_pts2 = (*it)->extractNextToAdjacent(this);
      if (adj_pts1.size() < min_bd || adj_pts2.size() < min_bd)
	continue;

      of1 << "400 1 0 4 0 255 0 255" << std::endl;
      of1 << adj_pts1.size() << std::endl;
      for (size_t ki=0; ki<adj_pts1.size(); ++ki)
	of1 << adj_pts1[ki]->getPoint() << std::endl;
      
      of2 << "400 1 0 4 0 255 0 255" << std::endl;
      of2 << adj_pts2.size() << std::endl;
      for (size_t ki=0; ki<adj_pts2.size(); ++ki)
	of2 << adj_pts2[ki]->getPoint() << std::endl;

      vector<RevEngPoint*> end_pts1, end_pts2;
      for (size_t ki=0; ki<adj_pts1.size(); ++ki)
	{
	  int num_adj_reg = adj_pts1[ki]->numAdjacentRegions();
	  if (num_adj_reg > 2)
	    end_pts1.push_back(adj_pts1[ki]);
	}
      for (size_t ki=0; ki<adj_pts2.size(); ++ki)
	{
	  int num_adj_reg = adj_pts2[ki]->numAdjacentRegions();
	  if (num_adj_reg > 2)
	    end_pts2.push_back(adj_pts2[ki]);
	}
      
      of1 << "400 1 0 4 255 0 0 255" << std::endl;
      of1 << end_pts1.size() << std::endl;
      for (size_t ki=0; ki<end_pts1.size(); ++ki)
	of1 << end_pts1[ki]->getPoint() << std::endl;
      
      of2 << "400 1 0 4 255 0 0 255" << std::endl;
      of2 << end_pts2.size() << std::endl;
      for (size_t ki=0; ki<end_pts2.size(); ++ki)
	of2 << end_pts2[ki]->getPoint() << std::endl;

      vector<RevEngPoint*> sorted1 = sortPtsSeq(mean_edge_len, adj_pts1, end_pts1);
      vector<RevEngPoint*> sorted2 = sortPtsSeq(mean_edge_len, adj_pts2, end_pts2);

      int degree = 3;
      double tol2 = 10.0*tol;
      int maxiter = 8;
      shared_ptr<SplineCurve> cv1 = RevEngUtils::createCurve(sorted1, degree,
							     tol2, maxiter);
      if (cv1.get())
	{
	  cv1->writeStandardHeader(of1);
	  cv1->write(of1);
	}
      
      shared_ptr<SplineCurve> cv2 = RevEngUtils::createCurve(sorted2, degree,
							     tol2, maxiter);
      if (cv2.get())
	{
	  cv2->writeStandardHeader(of2);
	  cv2->write(of2);
	}

      shared_ptr<SplineCurve> mid = RevEngUtils::midCurve(cv1, cv2);
      if (mid.get())
	{
	  mid->writeStandardHeader(of1);
	  mid->write(of1);
	  mid->writeStandardHeader(of2);
	  mid->write(of2);
	}

      vector<RevEngPoint*> out1 = extractBdOutPoints(mid, adj_pts1, tol);
      vector<RevEngPoint*> out2 = (*it)->extractBdOutPoints(mid, adj_pts2, tol);
      if (out1.size() > 0)
	{
	  std::ofstream of3("bd_out_pt1.g2");
	  of3 << "400 1 0 4 255 0 0 255" << std::endl;
	  of3 << out1.size() << std::endl;
	  for (size_t ki=0; ki<out1.size(); ++ki)
	    of3 << out1[ki]->getPoint() << std::endl;
	}
      if (out2.size() > 0)
	{
	  std::ofstream of4("bd_out_pt2.g2");
	  of4 << "400 1 0 4 255 0 0 255" << std::endl;
	  of4 << out2.size() << std::endl;
	  for (size_t ki=0; ki<out2.size(); ++ki)
	    of4 << out2[ki]->getPoint() << std::endl;
	}
       
      int stop_break = 1;
    }
}

//===========================================================================
vector<RevEngPoint*>  RevEngRegion::extractBdOutPoints(shared_ptr<SplineCurve>& crv,
						       vector<RevEngPoint*>& seq_pts,
						       double tol)
//===========================================================================
{
  vector<RevEngPoint*> left, right, upper, lower;
  vector<RevEngPoint*> points;
  points.insert(points.end(), seq_pts.begin(), seq_pts.end());
  int seq_size = seq_pts.size();

  for (size_t ki=0; ki<points.size(); ++ki)
    points[ki]->setVisited();

  double maxdist = 0.0;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      if (ki == seq_size)
	maxdist *= 2.0;
      Vector3D xyz = points[ki]->getPoint();
      Point pos(xyz[0], xyz[1], xyz[2]);
      double tpar, dist;
      Point close;
      crv->ParamCurve::closestPoint(pos, tpar, close, dist);
      if (dist > maxdist)
	{
	  if (ki < seq_size)
	    maxdist = dist;
	  else
	    continue;
	}
      
      if (tpar <= crv->startparam())
	lower.push_back(points[ki]);
      else if (tpar >= crv->endparam())
	upper.push_back(points[ki]);
      else
	{
	  vector<Point> der(2);
	  crv->point(der, tpar, 1);
	  Point norm = points[ki]->getMongeNormal();
	  Point vec = pos - close;
	  Point vec2 = vec.cross(der[1]);
	  if (vec2*norm >= 0)
	    left.push_back(points[ki]);
	  else
	    right.push_back(points[ki]);
	}
      
      vector<ftSamplePoint*> next = points[ki]->getNeighbours();
      for (size_t kj=0; kj<next.size(); ++kj)
	{
	  RevEngPoint *curr = dynamic_cast<RevEngPoint*>(next[kj]);
	  RevEngRegion *adj_reg = curr->region();
	  if (curr->visited())
	    continue;
	  if (adj_reg != this)
	    continue;
	  curr->setVisited();
	  points.push_back(curr);
	}
    }

  for (size_t ki=0; ki<points.size(); ++ki)
    points[ki]->unsetVisited();
  
  return (left.size() < right.size()) ? left : right;
}

//===========================================================================
vector<RevEngPoint*>  RevEngRegion::sortPtsSeq(double mean_edge_len,
					       vector<RevEngPoint*>& seq_pts,
					       vector<RevEngPoint*>& sub_pts)
//===========================================================================
{
  vector<RevEngPoint*> sub_pts2;
  if (sub_pts.size() < 2)
    {
      // Check for maximum distance between input points
      vector<RevEngPoint*> dummy;
      if (seq_pts.size() < 2)
	return dummy;

      double fac = 10.0;
      size_t ix1 = 0, ix2 = 1;
      ix2 = std::min(ix2, seq_pts.size()-1);
      double maxlen = seq_pts[ix1]->pntDist(seq_pts[ix2]);
      for (size_t ki=0; ki<seq_pts.size(); ++ki)
	for (size_t kj=ki+1; kj<seq_pts.size(); ++kj)
	  {
	    double len = seq_pts[ki]->pntDist(seq_pts[kj]);
	    if (len > maxlen)
	      {
		maxlen = len;
		ix1 = ki;
		ix2 = kj;
	      }
	  }

      if (sub_pts.size() == 1)
	{
	  sub_pts2.push_back(sub_pts[0]);
	  double len1 = sub_pts[0]->pntDist(seq_pts[ix1]);
	  double len2 = sub_pts[0]->pntDist(seq_pts[ix2]);
	  if (len1 >= len2)
	    sub_pts2.push_back(seq_pts[ix1]);
	  else
	    sub_pts2.push_back(seq_pts[ix2]);
	}
      else
	{
	  sub_pts2.push_back(seq_pts[ix1]);
	  if (maxlen > fac*mean_edge_len)
	    sub_pts2.push_back(seq_pts[ix2]);
	}
    }
  else
    sub_pts2.insert(sub_pts2.end(), sub_pts.begin(), sub_pts.end());
  vector<RevEngPoint*> seq_pts2(seq_pts.begin(), seq_pts.end());
  vector<vector<RevEngPoint*> > all_seq;
  
  while (sub_pts2.size() > 0)
    {
      size_t ix1 = 0, ix2 = 1;
      ix2 = std::min(ix2, sub_pts2.size()-1);
      double maxlen = sub_pts2[ix1]->pntDist(sub_pts2[ix2]);
      for (size_t ki=0; ki<sub_pts2.size(); ++ki)
	for (size_t kj=ki+1; kj<sub_pts2.size(); ++kj)
	  {
	    double len = sub_pts2[ki]->pntDist(sub_pts2[kj]);
	    if (len > maxlen)
	      {
		maxlen = len;
		ix1 = ki;
		ix2 = kj;
	      }
	  }

      vector<RevEngPoint*> sortseq;
      sortseq.push_back(sub_pts2[ix1]);
      bool more_pts = true;
      RevEngPoint *curr = sub_pts2[ix1];
      vector<RevEngPoint*> prev_pts;
      while (more_pts)
	{
	  vector<RevEngPoint*> next_pts;
	  vector<ftSamplePoint*> next = curr->getNeighbours();
	  for (size_t ki=0; ki<next.size(); ++ki)
	    {
	      size_t kj;
	      for (kj=0; kj<next_pts.size(); ++kj)
		if (next_pts[kj] == next[ki])
		  break;
	      if (kj == next_pts.size())
		{
		  size_t kh;
		  for (kh=0; kh<sortseq.size(); ++kh)
		    if (sortseq[kh] == next[ki])
		      break;
		  if (kh == sortseq.size())
		    {
		      size_t kr;
		      for (kr=0; kr<seq_pts2.size(); ++kr)
			if (seq_pts2[kr] == next[ki])
			  break;
		      if (kr < seq_pts2.size())
			next_pts.push_back(dynamic_cast<RevEngPoint*>(next[ki]));
		    }
		}
	    }
	  for (size_t ki=0; ki<next_pts.size(); ++ki)
	    for (size_t kj=ki+1; kj<next_pts.size(); )
	      {
		if (next_pts[ki] == next_pts[kj])
		  next_pts.erase(next_pts.begin()+kj);
		else
		  ++kj;
	      }

	  if (next_pts.size() == 0)
	    next_pts = prev_pts;
	  
	  if (next_pts.size() == 0)
	    break;
	  if (next_pts.size() == 1)
	    {
	      prev_pts.clear();
	      sortseq.push_back(next_pts[0]);
	    }
	  else
	    {
	      double minlen = curr->pntDist(next_pts[0]);
	      size_t ix3 = 0;
	      for (size_t kr=1; kr<next_pts.size(); ++kr)
		{
		  double len = curr->pntDist(next_pts[kr]);
		  if (len < minlen)
		    {
		      minlen = len;
		      ix3 = kr;
		    }
		}
	      prev_pts.clear();
	      sortseq.push_back(next_pts[ix3]);
	      for (size_t kr=0; kr<next_pts.size(); ++kr)
		if (kr != ix3)
		  prev_pts.push_back(next_pts[kr]);
	    }
	  curr = sortseq[sortseq.size()-1];
	}

      std::ofstream of("sorted_seq.g2");
      for (size_t ki=0; ki<sortseq.size(); ++ki)
	{
	  of << "400 1 0 4 0 0 0 255" << std::endl;
	  of << "1" << std::endl;
	  of << sortseq[ki]->getPoint() << std::endl;
	}
      int stop_break = 1;

      for (size_t ki=0; ki<sortseq.size(); ++ki)
	{
	  for (size_t kj=0; kj<sub_pts2.size(); ++kj)
	    {
	      if (sortseq[ki] == sub_pts2[kj])
		{
		  sub_pts2.erase(sub_pts2.begin()+kj);
		  break;
		}
	      for (size_t kj=0; kj<seq_pts2.size(); ++kj)
		{
		  if (sortseq[ki] == seq_pts2[kj])
		    {
		      seq_pts2.erase(seq_pts2.begin()+kj);
		      break;
		    }
		}
	    }
	}
      all_seq.push_back(sortseq);
    }

  // Join sub sequences
  vector<RevEngPoint*> sorted;
  if (all_seq.size() > 0)
    {
      size_t all_size = all_seq.size();
      while (all_size > 1)
	{
	  double t1min = std::numeric_limits<double>::max();
	  double t2min = std::numeric_limits<double>::max();
	  int x1 = -1, x2 = -1;
	  bool turn1 = false, turn2 = false;
	  for (size_t kj=1; kj<all_size; ++kj)
	    {
	      double l1 = all_seq[0][0]->pntDist(all_seq[kj][0]);
	      double l2 = all_seq[0][0]->pntDist(all_seq[kj][all_seq[kj].size()-1]);
	      double l3 = all_seq[0][all_seq[0].size()-1]->pntDist(all_seq[kj][0]);
	      double l4 = all_seq[0][all_seq[0].size()-1]->pntDist(all_seq[kj][all_seq[kj].size()-1]);
	      if (std::min(l1, l2) < t1min)
		{
		  t1min = std::min(l1, l2);
		  x1 = (int)kj;
		  turn1 = (l1 < l2);
		}
	    
	      if (std::min(l3, l4) < t2min)
		{
		  t2min = std::min(l3, l4);
		  x2 = (int)kj;
		  turn2 = (l4 < l3);
		}
	    }
	  if (t1min < t2min)
	    {
	      if (turn1)
		{
		  size_t last = all_seq[x1].size()-1;
		  for (size_t kr=0; kr<all_seq[x1].size()/2; ++kr)
		    std::swap(all_seq[x1][kr], all_seq[x1][last-kr]);
		}
	      all_seq[0].insert(all_seq[0].begin(), all_seq[x1].begin(),
				all_seq[x1].end());
	      if (x1 < all_size-1)
		std::swap(all_seq[x1],all_seq[all_size-1]);
	    }
	  else
	    {
	      if (turn2)
		{
		  size_t last = all_seq[x2].size()-1;
		  for (size_t kr=0; kr<all_seq[x2].size()/2; ++kr)
		    std::swap(all_seq[x2][kr], all_seq[x2][last-kr]);
		}
	      all_seq[0].insert(all_seq[0].end(), all_seq[x2].begin(),
				all_seq[x2].end());
	      if (x2 < all_size-1)
		std::swap(all_seq[x2],all_seq[all_size-1]);
	    }
	  all_size--;
	}
    }

  bool include_missing = false;
  if (include_missing)
    {
      // Include missing input points
      for (size_t ki=0; ki<seq_pts2.size(); ++ki)
	{
	  double minlen = seq_pts2[ki]->pntDist(all_seq[0][0]);
	  size_t min_ix = 0;
	  for (size_t kj=1; kj<all_seq[0].size(); ++kj)
	    {
	      double len = seq_pts2[ki]->pntDist(all_seq[0][kj]);
	      if (len < minlen)
		{
		  min_ix = kj;
		  minlen = len;
		}
	    }

	  double len1 = (min_ix == 0) ? 0.0 : seq_pts2[ki]->pntDist(all_seq[0][min_ix-1]);
	  double len2 = (min_ix == all_seq[0].size()-1) ? 0.0 :
	    seq_pts2[ki]->pntDist(all_seq[0][min_ix+1]);
	  size_t ix = min_ix + (len2 > len1);
	  all_seq[0].insert(all_seq[0].begin()+ix, seq_pts2[ki]);
    }
    }
  return all_seq[0];
}

//===========================================================================
void RevEngRegion::adjustWithSurf(double tol, double angtol)
//===========================================================================
{
  //double eps = 1.0e-6;
  if (associated_sf_.size() == 0)
    return;  // No surface with which to check growt

  //int sfcode;
  //ClassType classtype = associated_sf_[0]->instanceType(sfcode);
  shared_ptr<ParamSurface> surf = associated_sf_[0]->surface();
  shared_ptr<BoundedSurface> bdsurf =
    dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf);
  if (bdsurf.get())
    surf = bdsurf->underlyingSurface();

  std::ofstream res1("residuals_source.txt");
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    res1 << group_points_[ki]->getPoint() << " " << group_points_[ki]->getSurfaceDist() << std::endl;
  
  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
    {
      // Check criteria for growt
      int num = (*it)->numPoints();
      // if (num > group_points_.size())
      // 	continue;
      
      std::ofstream res2("residuals_adj.txt");
      for (size_t kh=0; kh<(*it)->group_points_.size(); ++kh)
	res2 << (*it)->group_points_[kh]->getPoint() << " " << (*it)->group_points_[kh]->getSurfaceDist() << std::endl;
  
       std::ofstream of("curr_extend.g2");
      of << "400 1 0 4 0 255 0 255" << std::endl;
      of << group_points_.size() << std::endl;
      for (size_t kh=0; kh<group_points_.size(); ++kh)
      	of << group_points_[kh]->getPoint() << std::endl;
      
      of << "400 1 0 4 255 0 0 255" << std::endl;
      of << num << std::endl;
      for (int ka=0; ka<num; ++ka)
      	of << (*it)->getPoint(ka)->getPoint() << std::endl;
      
      surf->writeStandardHeader(of);
      surf->write(of);
	      
      vector<RevEngPoint*> in, out;
      vector<pair<double, double> > dist_ang;
      vector<double> parvals;
      double maxd, avd;
      int num_inside;
      RevEngUtils::distToSurf((*it)->pointsBegin(),
			      (*it)->pointsEnd(), surf, tol,
			      maxd, avd, num_inside, in, out,
			      parvals, dist_ang);
      vector<RevEngPoint*> better1, worse1;
      std::ofstream res3("residuals_adj2.txt");
      for (size_t kh=0; kh<(*it)->group_points_.size(); ++kh)
	res3 << (*it)->group_points_[kh]->getPoint() << " " << dist_ang[kh].first << std::endl;
  
      for (size_t kh=0; kh<dist_ang.size(); ++kh)
	{
	  double ptdist = tol, ptang = angtol;
	  if ((*it)->hasSurface())
	    (*it)->group_points_[kh]->getSurfaceDist(ptdist, ptang);
	  if (dist_ang[kh].first <= ptdist && dist_ang[kh].second < angtol) //ptang)
	  //if (dist_ang[kh].first <= tol)
	    better1.push_back((*it)->group_points_[kh]);
	  else
	    worse1.push_back((*it)->group_points_[kh]);
	}
      std::ofstream ofc1("better_worse1.g2");
      ofc1 << "400 1 0 4 155 50 50 255" << std::endl;
      ofc1 << better1.size() << std::endl;
      for (size_t kr=0; kr<better1.size(); ++kr)
	ofc1 << better1[kr]->getPoint() << std::endl;
      ofc1 << "400 1 0 4 50 155 50 255" << std::endl;
      ofc1 << worse1.size() << std::endl;
      for (size_t kr=0; kr<worse1.size(); ++kr)
	ofc1 << worse1[kr]->getPoint() << std::endl;
      
      std::ofstream ofd("in_out_extend.g2");
      ofd << "400 1 0 4 155 50 50 255" << std::endl;
      ofd << in.size() << std::endl;
      for (size_t kr=0; kr<in.size(); ++kr)
	ofd << in[kr]->getPoint() << std::endl;
      ofd << "400 1 0 4 50 155 50 255" << std::endl;
      ofd << out.size() << std::endl;
      for (size_t kr=0; kr<out.size(); ++kr)
	ofd << out[kr]->getPoint() << std::endl;

      vector<RevEngPoint*> better2, worse2;
      if ((*it)->hasSurface())
	{
	  shared_ptr<ParamSurface> surf2 = (*it)->getSurface(0)->surface();
	  double maxd2, avd2;
	  int num_inside2;
	  vector<RevEngPoint*> in2, out2;
	  vector<pair<double, double> > dist_ang2;
	  vector<double> parvals2;
	  RevEngUtils::distToSurf(pointsBegin(),
				  pointsEnd(), surf2, tol,
				  maxd2, avd2, num_inside2, in2, out2,
				  parvals2, dist_ang2);
	  for (size_t kh=0; kh<dist_ang2.size(); ++kh)
	    {
	      double ptdist, ptang;
	      group_points_[kh]->getSurfaceDist(ptdist, ptang);
	      if (dist_ang2[kh].first <= ptdist && dist_ang2[kh].second < angtol) //ptang)
		{
		  better2.push_back(group_points_[kh]);
		}
	      else
		worse2.push_back(group_points_[kh]);
	}
	  std::ofstream ofc2("better_worse2.g2");
	  ofc2 << "400 1 0 4 155 50 50 255" << std::endl;
	  ofc2 << better2.size() << std::endl;
	  for (size_t kr=0; kr<better2.size(); ++kr)
	    ofc2 << better2[kr]->getPoint() << std::endl;
	  ofc2 << "400 1 0 4 50 155 50 255" << std::endl;
	  ofc2 << worse2.size() << std::endl;
	  for (size_t kr=0; kr<worse2.size(); ++kr)
	    ofc2 << worse2[kr]->getPoint() << std::endl;
	}

      // Move points
      bool move2 = false;
      size_t b1size = 2*better1.size();
      while (better1.size() < b1size)
	{
	  b1size = better1.size();
	  for (size_t ki=0; ki<better1.size(); )
	    {
	      if (better1[ki]->isNeighbour(this))
		{
		  move2 = true;
		  (*it)->removePoint(better1[ki]);
		  better1[ki]->setRegion(this);
		  better1[ki]->addMove();
		  group_points_.push_back(better1[ki]);
		  better1.erase(better1.begin()+ki);
		}
	      else
		++ki;
	    }
	}

      if (false)
	{
      size_t b2size = 2*better2.size();
      while (better2.size() < b2size)
	{
	  b2size = better2.size();
	  for (size_t ki=0; ki<better2.size(); )
	    {
	      if (better2[ki]->isNeighbour(*it))
		{
		  move2 = true;
		  removePoint(better2[ki]);
		  better2[ki]->setRegion(*it);
		  better2[ki]->addMove();
		  (*it)->addPoint(better2[ki]);
		  better2.erase(better2.begin()+ki);
		}
	      else
		++ki;
	    }
	}
	}
     
      if ((*it)->hasSurface() && (*it)->numPoints() > 5 && move2)
	{
	  // Check if the surface should be updated
	  (*it)->checkReplaceSurf(tol);
	}
    }
  

  // Check if the surface should be updated
  checkReplaceSurf(tol);
  std::ofstream res4("residuals_res.txt");
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    res4 << group_points_[ki]->getPoint() << " " << group_points_[ki]->getSurfaceDist() << std::endl;
  int stop_break = 1;
  
}
  
//===========================================================================
void RevEngRegion::checkReplaceSurf(double tol, bool always)
//===========================================================================
{
  int sfcode;
  ClassType classtype[2];
  classtype[0] = Class_Unknown;
  classtype[1] = associated_sf_[0]->instanceType(sfcode);
  shared_ptr<ParamSurface> primary;
  if (primary_.get() && avdist_primary_ < tol &&
      num_in_primary_ > (int)group_points_.size()/2)
    {
      primary = primary_;
      classtype[0] = primary->instanceType();
    }

  shared_ptr<SplineCurve> profile;
  Point pt1, pt2;
  for (int ka=0; ka<2; ++ka)
    {
      if (classtype[ka] == Class_Unknown)
	continue;
  
      shared_ptr<ParamSurface> updated, updated2;
      if (classtype[ka] == Class_Plane)
	{
	  updated = computePlane(group_points_, avnorm_);
	}
      else if (classtype[ka] == Class_Cylinder)
	{
	  vector<vector<RevEngPoint*> > configs;
	  updated = computeCylinder(group_points_, tol, configs);
	}
      else if (classtype[ka] == Class_Sphere)
	{
	  updated = computeSphere(group_points_);
	}
      else if (classtype[ka] == Class_Cone)
	{
	  Point apex;
	  updated = computeCone(group_points_, apex);
	}
      else if (classtype[ka] == Class_Torus)
	{
	  shared_ptr<Torus> torus2;
	  updated = computeTorus(group_points_, tol, torus2);
	  updated2 = torus2;
	}
      else if (classtype[ka] == Class_SplineSurface && sfcode == 1)
	{
	  updated = computeLinearSwept(tol, profile, pt1, pt2);
	}
      else if (classtype[ka] == Class_SplineSurface)
	{
	  updated = updateFreeform(tol);
	  if (!updated.get())
	    updated = computeFreeform(tol);
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
	  vector<pair<double, double> > dist_ang;
	  vector<double> parvals;
	  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
				  updated, tol, maxd, avd, num_inside,
				  inpt, outpt, parvals, dist_ang);

	  if (updated2.get())
	    {
	      double maxd2, avd2;
	      int num_inside2;
	      vector<RevEngPoint*> inpt2, outpt2;
	      vector<pair<double, double> > dist_ang2;
	      vector<double> parvals2;
	      RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
				      updated2, tol, maxd2, avd2, num_inside2,
				      inpt2, outpt2, parvals2, dist_ang2);
	      if ((num_inside2 > num_inside || (num_inside2 == num_inside && avd2 < avd))
		  && (((num_inside2 > num_inside_ || (num_inside2 == num_inside_ && avd2 < avdist_))
		       && avd2 < tol) || always))
		{
		  replacesurf = updated2;
		  for (size_t kh=0; kh<group_points_.size(); ++kh)
		    {
		      group_points_[kh]->setPar(Vector2D(parvals2[2*kh],parvals2[2*kh+1]));
		      group_points_[kh]->setSurfaceDist(dist_ang2[kh].first, dist_ang2[kh].second);
		    }
		  setAccuracy(maxd2, avd2, num_inside2);
		}
	      if (ka == 0 && num_inside2 >= num_in_primary_ && avd2 < avdist_primary_)
		{
		  setPrimarySf(updated2, maxd2, avd2, num_inside2);
		}
	      int stop_break2 = 1;
	    }
	  if ((!replacesurf.get()) &&
	      (((num_inside > num_inside_ || (num_inside == num_inside_ && avd < avdist_)) &&
		avd < tol) || always))
	    {
	      replacesurf = updated;
	      for (size_t kh=0; kh<group_points_.size(); ++kh)
		{
		  group_points_[kh]->setPar(Vector2D(parvals[2*kh],parvals[2*kh+1]));
		  group_points_[kh]->setSurfaceDist(dist_ang[kh].first, dist_ang[kh].second);
		}
	      setAccuracy(maxd, avd, num_inside);
	    }
	  if (ka == 0 && num_inside >= num_in_primary_ && avd < avdist_primary_)
	    {
	      setPrimarySf(updated, maxd, avd, num_inside);
	    }
	  int stop_break = 1;
	}

      if (replacesurf.get())
	{
	  associated_sf_[0]->replaceSurf(replacesurf);
	  if (sfcode == 1 && profile.get())
	    {
	      // A linear swept surface
	      associated_sf_[0]->setLinearSweepInfo(profile, pt1, pt2);
	    }
	  computeDomain();
	}
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


int compare_t(const void* el1, const void* el2)
{
  if (((double*)el1)[0] < ((double*)el2)[0])
    return -1;
  else if (((double*)el1)[0] > ((double*)el2)[0])
    return 1;
  else
    return 0;
}

//===========================================================================
void  RevEngRegion::curveApprox(vector<Point>& points, double tol,
				shared_ptr<Circle> circle, vector<double>& parval,
				shared_ptr<SplineCurve>& curve, Point& xpos)
//===========================================================================
{
  double eps = 0.001;
  vector<double> pts;
  vector<double> param;
  double tmin = circle->startparam();
  double tmax = circle->endparam();
  double tdel = tmax - tmin;
  double tmin2 = tmax;
  double tmax2 = tmin;
  vector<double> tmppts;
  tmppts.reserve(4*points.size());
  parval.resize(points.size());
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      double tpar, dist;
      Point close;
      circle->closestPoint(points[ki], tmin, tmax, tpar, close, dist);
      // if (ki > 0 && ((tpar-tmin < tmin2-tpar && tmax-tmax2 < tmin2-tmin) ||
      // 		     (tmax-tpar < tpar-tmax2 && tmin2-tmin < tmax-tmax2)))
      if (ki > 0 && (tpar<tmin2 || tpar>tmax2) &&
	  std::min(fabs(tmin2-tpar+tdel),fabs(tpar+tdel-tmax2)) <
	  std::min(fabs(tpar-tmax2),fabs(tpar-tmin2)))
	//(fabs(tmin2-tpar+tdel) < fabs(tpar-tmax2) || fabs(tpar+tdel-tmax2) < fabs(tpar-tmax2)))
	{
	  if (tpar-tmin < tmax-tpar)
	    tpar += tdel;
	  else
	    tpar -= tdel;
	}

      parval[ki] = tpar;
      tmppts.push_back(tpar);
      tmppts.insert(tmppts.end(), points[ki].begin(), points[ki].end());
      // pts.insert(pts.end(), points[ki].begin(), points[ki].end());
      // param.push_back(tpar);
      tmin2 = std::min(tmin2, tpar);
      tmax2 = std::max(tmax2, tpar);
    }

  qsort(&tmppts[0], points.size(), 4*sizeof(double), compare_t);
  pts.resize(3*points.size());
  param.resize(points.size());
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      param[ki] = tmppts[4*ki];
      for (size_t ka=0; ka<3; ++ka)
	pts[3*ki+ka] = tmppts[4*ki+ka+1];
    }

  if (tmax2 - tmin2 < 2*M_PI-eps && (tmin2 < -eps || tmax2 > 2*M_PI+eps))
    {
      double tpar = (tmax2 > 2*M_PI+eps) ? 0.5*(tmax2 + tmin2 - 2*M_PI)
		      : 0.5*(tmin2 + tmax2 + 2*M_PI);
      xpos = circle->ParamCurve::point(tpar);
    }
    
  int inner = (int)(2.0*(tmax2 - tmin2)/M_PI);
  int ik = 4;
  int in = ik + inner;
  // double tdel = (tmax2 - tmin2)/(double)(in - ik + 1);
  // double et[12];
  // for (int ka=0; ka<ik; ++ka)
  //   {
  //     et[ka] = tmin2;
  //     et[in+ka] = tmax2;
  //   }
  // for (int ka=ik; ka<in; ++ka)
  //   et[ka] = tmin2 + (ka-ik+1)*tdel;

  double smoothwgt = 1.0e-9; //0.001;
  ApproxCurve approx(pts, param, 3, tol, in, ik);
  approx.setSmooth(smoothwgt);
  int maxiter = 4; //3;
  double maxdist, avdist;
  curve = approx.getApproxCurve(maxdist, avdist, maxiter);
  // vector<double> ecoef(3*in, 0.0);
  // shared_ptr<SplineCurve> cv(new SplineCurve(in, ik, et, &ecoef[0], 3));

  // SmoothCurve smooth(3);
  // vector<int> cfn(in, 0);
  // vector<double> wgts(param.size(), 1.0);
  // smooth.attach(cv, &cfn[0]);

  // double wgt1 = 0.0, wgt2 = 0.1, wgt3 = 0.1;
  // double approxwgt = 1.0 - wgt1 - wgt2 - wgt3;
  // smooth.setOptim(wgt1, wgt2, wgt3);
  // smooth.setLeastSquares(pts, param, wgts, approxwgt);

  // shared_ptr<SplineCurve> curve0;
  // smooth.equationSolve(curve0);
  std::ofstream of("points_and_tangents.txt");
  of << points.size() << std::endl;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      of << pts[3*ki] << " " << pts[3*ki+1] << " ";
      vector<Point> der(2);
      curve->point(der, param[ki], 1);
      of << der[1][0] << " " << der[1][1] << std::endl;
    }
  int stop_break = 1;
}

//===========================================================================
void  RevEngRegion::curveApprox(vector<Point>& points,
				shared_ptr<ParamCurve> cvin,
				int ik, int in, 
				shared_ptr<SplineCurve>& curve)
//===========================================================================
{
  vector<double> pts;
  vector<double> param;
  double tmin = cvin->startparam();
  double tmax = cvin->endparam();
  double tmin2 = tmax;
  double tmax2 = tmin;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      double tpar, dist;
      Point close;
      cvin->closestPoint(points[ki], tmin, tmax, tpar, close, dist);
      pts.insert(pts.end(), points[ki].begin(), points[ki].end());
      param.push_back(tpar);
      tmin2 = std::min(tmin2, tpar);
      tmax2 = std::max(tmax2, tpar);
    }

  double tdel = (tmax2 - tmin2)/(double)(in - ik + 1);
  vector<double> et(ik+in);
  for (int ka=0; ka<ik; ++ka)
    {
      et[ka] = tmin2;
      et[in+ka] = tmax2;
    }
  for (int ka=ik; ka<in; ++ka)
    et[ka] = tmin2 + (ka-ik+1)*tdel;

  vector<double> ecoef(3*in, 0.0);
  shared_ptr<SplineCurve> cv(new SplineCurve(in, ik, &et[0], &ecoef[0], 3));

  SmoothCurve smooth(3);
  vector<int> cfn(in, 0);
  vector<double> wgts(param.size(), 1.0);
  smooth.attach(cv, &cfn[0]);

  smooth.setOptim(0.0, 0.001, 0.001);
  smooth.setLeastSquares(pts, param, wgts, 0.998);

  smooth.equationSolve(curve);
  int stop_break = 1;
}

//===========================================================================
void RevEngRegion::configSplit(vector<RevEngPoint*>& points,
			       vector<double>& param,
			       shared_ptr<Cylinder> cyl,
			       shared_ptr<SplineCurve> spl, double tol,
			       vector<vector<RevEngPoint*> >& configs)
//===========================================================================
{
  // Check cylinder axis
  Point axis = cyl->direction();
  double angtol = 0.2;
  double inlim = 0.8;
  int nmb_in = 0;
  double mpi2 = 0.5*M_PI;
  vector<Vector3D> low, high;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      Point normal = group_points_[kr]->getMongeNormal();
      double ang = axis.angle(normal);
      if (ang < mpi2)
	low.push_back(group_points_[kr]->getPoint());
      else
	high.push_back(group_points_[kr]->getPoint());
      if (fabs(ang-mpi2) <= angtol)
	nmb_in++;
    }
  double in_frac = (double)nmb_in/(double)group_points_.size();
  if (in_frac < inlim)
    return;

  std::ofstream of1("low_axis.g2");
  of1 << "400 1 0 4 100 155 0 255" << std::endl;
  of1 << low.size() << std::endl;
  for (size_t kr=0; kr<low.size(); ++kr)
    of1 << low[kr] << std::endl;
  
  std::ofstream of2("high_axis.g2");
  of2 << "400 1 0 4 0 155 100 255" << std::endl;
  of2 << high.size() << std::endl;
  for (size_t kr=0; kr<high.size(); ++kr)
    of2 << high[kr] << std::endl;
  
  // Compute all intersections between the cylinder and the spline curve
  double eps = std::min(1.0e-6, 0.1*tol);
  vector<double> intpar;
  vector<pair<double,double> > int_cvs;
  intersectCurveCylinder(spl.get(), cyl->getLocation(), cyl->getAxis(),
			 cyl->getRadius(), eps, intpar, int_cvs);
  for (size_t ki=0; ki<int_cvs.size(); ++ki)
    {
      intpar.push_back(int_cvs[ki].first);
      intpar.push_back(int_cvs[ki].second);
    }
  if (intpar.size() == 0)
    return;

  // Define parameter intervals of different configurations
  std::sort(intpar.begin(), intpar.end());
  vector<double> delpar;

  // Check startpoint, endpoint and points in the middle of intervals
  double upar, vpar, dist;
  Point close;
  Point pos = spl->ParamCurve::point(spl->startparam());
  cyl->closestPoint(pos, upar, vpar, close, dist, eps);
  delpar.push_back(spl->startparam()-tol);
  if (dist > tol)
    {
      delpar.push_back(intpar[0]);
    }
  for (size_t ki=1; ki<intpar.size(); ++ki)
    {
      double tpar = 0.5*(intpar[ki-1] + intpar[ki]);
      pos = spl->ParamCurve::point(tpar);
      cyl->closestPoint(pos, upar, vpar, close, dist, eps);
      if (dist > tol)
	delpar.push_back(intpar[ki]);
    }
  pos = spl->ParamCurve::point(spl->endparam());
  cyl->closestPoint(pos, upar, vpar, close, dist, eps);
  if (dist > tol)
    {
      if (delpar.size() > 0 &&
	  intpar[intpar.size()-1] > delpar[delpar.size()-1])
	delpar.push_back(intpar[intpar.size()-1]);
    }
  delpar.push_back(spl->endparam());

  if (delpar.size() == 0)
    return;
  
  // Divide point set according to configuration
  configs.resize(delpar.size()-1);
  for (size_t kj=0; kj<points.size(); ++kj)
    {
      for (size_t ki=1; ki<delpar.size(); ++ki)
	if (param[kj] > delpar[ki-1] && param[kj] <= delpar[ki])
	  {
	    configs[ki-1].push_back(points[kj]);
	    break;
	  }
    }
  int stop_break = 1;
}

//===========================================================================
bool  RevEngRegion::hasEdgeBetween(RevEngRegion* adj)
//===========================================================================
{
  if (adj->numPoints() < (int)group_points_.size())
    return adj->hasEdgeBetween(this);

  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      vector<ftSamplePoint*> next = group_points_[ki]->getNeighbours();
      for (size_t kj=0; kj<next.size(); ++kj)
	{
	  RevEngPoint *pt = dynamic_cast<RevEngPoint*>(next[kj]);
	  if (!pt->isEdge(edge_class_type_))
	    continue;

	  vector<ftSamplePoint*> next2 = pt->getNeighbours();
	  for (size_t kr=0; kr<next2.size(); ++kr)
	    {
	      RevEngPoint *pt2 = dynamic_cast<RevEngPoint*>(next2[kr]);
	      if (pt2->region() == adj)
		return true;
	    }
	}
    }
  return false;
}

//===========================================================================
void RevEngRegion::store(std::ostream& os) const
//===========================================================================
{
  os << group_points_.size() << std::endl;
for (size_t ki=0; ki<group_points_.size(); ++ki)
  os << group_points_[ki]->getIndex() << " ";
 os << std::endl;
 os << classification_type_ << " " << frac_norm_in_ << std::endl;
 os << maxdist_ << " " << avdist_ << " " << num_inside_ << std::endl;
 int base = basesf_.get() ? 1 : 0;
 os << base << std::endl;
 if (base)
   {
     basesf_->writeStandardHeader(os);
     basesf_->write(os);
     os << maxdist_base_ << " " << avdist_base_ << " " << num_in_base_ << std::endl;
   }

 int primary = primary_.get() ? 1 : 0;
 os << primary << std::endl;
 if (primary)
   {
     primary_->writeStandardHeader(os);
     primary_->write(os);
     os << maxdist_primary_ << " " << avdist_primary_ << " " << num_in_primary_ << std::endl;
   }
 os << associated_sf_.size() << std::endl;
 for (size_t ki=0; ki<associated_sf_.size(); ++ki)
     os << associated_sf_[ki]->getId() << " ";
}

//===========================================================================
void RevEngRegion::read(std::istream& is,
			shared_ptr<ftPointSet>& tri_sf,
			vector<int>& associated_sf_id)
//===========================================================================
{
  GoTools::init();
  int num_points;
  is >> num_points;
  group_points_.resize(num_points);
  int ix;
  for (int ki=0; ki<num_points; ++ki)
    {
      is >> ix;
      RevEngPoint* pt = dynamic_cast<RevEngPoint*>((*tri_sf)[ix]);
      pt->setRegion(this);
      group_points_[ki] = pt;
      }
  is >> classification_type_ >> frac_norm_in_;
  is >> maxdist_ >> avdist_ >> num_inside_;
  int base;
  is >> base;

  if (base)
    {
      ObjectHeader header;
      header.read(is);
      shared_ptr<GeomObject> obj(Factory::createObject(header.classType()));
      obj->read(is);
      basesf_ = dynamic_pointer_cast<ParamSurface,GeomObject>(obj);
      is >> maxdist_base_ >> avdist_base_ >> num_in_base_;
    }

  int primary;
  is >> primary;
  if (primary)
    {
      ObjectHeader header;
      header.read(is);
      shared_ptr<GeomObject> obj(Factory::createObject(header.classType()));
      obj->read(is);
      primary_ = dynamic_pointer_cast<ParamSurface,GeomObject>(obj);
      is >> maxdist_primary_ >> avdist_primary_ >> num_in_primary_;
    }

  int num_sf;
  is >> num_sf;
  for (int ki=0; ki<num_sf; ++ki)
    {
      int sf_id;
      is >> sf_id;
      associated_sf_id.push_back(sf_id);
    }
  if (num_sf > 0)
    computeDomain();

  // Bounding box and principal curvature summary
  maxk2_ = std::numeric_limits<double>::lowest();
  mink2_ = std::numeric_limits<double>::max();
  maxk1_ = std::numeric_limits<double>::lowest();
  mink1_ = std::numeric_limits<double>::max();
  MAH_ = MAK_ = avH_ = avK_ = 0.0;
  bbox_ = BoundingBox(3);
  double fac = 1.0/(double)group_points_.size();
  for  (size_t kj=0; kj<group_points_.size(); ++kj)
    {
      double k1 = group_points_[kj]->minPrincipalCurvature();
      double k2 = group_points_[kj]->maxPrincipalCurvature();
      double H = group_points_[kj]->meanCurvature();
      double K = group_points_[kj]->GaussCurvature();
      mink1_ = std::min(mink1_, fabs(k1));
      maxk1_ = std::max(maxk1_, fabs(k1));
      mink2_ = std::min(mink2_, fabs(k2));
      maxk2_ = std::max(maxk2_, fabs(k2));
      avH_ += fac*H;
      avK_ += fac*K;
      MAH_ += fac*fabs(H);
      MAK_ += fac*fabs(K);
      Vector3D point = group_points_[kj]->getPoint();
      Point point2(point[0], point[1], point[2]);
      bbox_.addUnionWith(point2);
    }
  
  normalcone_ = DirectionCone(group_points_[0]->getMongeNormal());
  avnorm_ = Point(0.0, 0.0, 0.0);
  for  (size_t kj=1; kj<group_points_.size(); ++kj)
    {
      Point norm = group_points_[kj]->getMongeNormal();
      normalcone_.addUnionWith(norm);
      avnorm_ += fac*norm;
    }
  
}

void RevEngRegion::writeSubTriangulation(std::ostream& of)
{
  of << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      vector<ftSamplePoint*> next = group_points_[kr]->getNeighbours();
      size_t nmb_next = next.size();
      for (int ka=(int)(next.size()-1); ka>=0; --ka)
	{
	  RevEngPoint *pr = dynamic_cast<RevEngPoint*>(next[ka]);
	  if (pr->region() != this)
	    next.erase(next.begin()+ka);
	}

      int bd = (next.size() < nmb_next) ? 1 : 0;
      of << kr << " " << group_points_[kr]->getPoint() << " " << bd << std::endl;
      of << next.size() << " ";
      for (size_t kh=0; kh<next.size(); ++kh)
	{
	  vector<RevEngPoint*>::iterator it = std::find(group_points_.begin(),
							group_points_.end(), next[kh]);
	  if (it == group_points_.end())
	    std::cout << "writeSubTriangulation: missing connection" << std::endl;
	  else
	    {
	      size_t ix = it - group_points_.begin();
	      of << ix << " ";
	    }
	}
      of << std::endl;
    }
}

void RevEngRegion::writeSurface(std::ostream& of)
{
  if (associated_sf_.size() == 0)
    return;
  shared_ptr<ParamSurface> surf = associated_sf_[0]->surface();
  shared_ptr<ElementarySurface> elemsf =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf);
  if (elemsf.get())
    {
      // double umin = std::numeric_limits<double>::max();
      // double umax = std::numeric_limits<double>::lowest();
      // double vmin = std::numeric_limits<double>::max();
      // double vmax = std::numeric_limits<double>::lowest();
      // for (size_t ki=0; ki<group_points_.size(); ++ki)
      // 	{
      // 	  Vector2D par = group_points_[ki]->getPar();
      // 	  umin = std::min(umin, par[0]);
      // 	  umax = std::max(umax, par[0]);
      // 	  vmin = std::min(vmin, par[1]);
      // 	  vmax = std::max(vmax, par[1]);
      // 	}
      double umin = domain_[0];
      double umax = domain_[1];
      double vmin = domain_[2];
      double vmax = domain_[3];
      shared_ptr<ElementarySurface> elemsf2(elemsf->clone());
      if (elemsf2->instanceType() != Class_Plane && umax-umin > 2*M_PI)
	{
	  umin = 0;
	  umax = 2*M_PI;
	}
      if (elemsf2->instanceType() != Class_Plane && elemsf2->instanceType() != Class_Sphere &&
	  (umin < -2.0*M_PI || umax > 2.0*M_PI))
	{
	  umin = 0;
	  umax = 2*M_PI;
	}
      if (elemsf2->instanceType() == Class_Sphere && (umin < 0.0 || umax > 2.0*M_PI))
	{
	  umin = 0;
	  umax = 2*M_PI;
	}
      
     if (elemsf2->instanceType() == Class_Sphere &&
	 (vmax-vmin > M_PI || vmin < -0.5*M_PI || vmax > 0.5*M_PI))
       {
	 vmin = -0.5*M_PI;
	 vmax = 0.5*M_PI;
       }
      if (elemsf2->instanceType() == Class_Torus &&
	  (vmax - vmin > 2*M_PI || vmin < -2.0*M_PI || vmax > 2.0*M_PI))
       {
	 vmin = 0;
	 vmax = 2*M_PI;
       }
	 
      elemsf2->setParameterBounds(umin, vmin, umax, vmax);
      elemsf2->writeStandardHeader(of);
      elemsf2->write(of);
    }
  else
    {
      surf->writeStandardHeader(of);
      surf->write(of);
    }
}

void RevEngRegion::writeRegionInfo(std::ostream& of)
{
  double len = bbox_.low().dist(bbox_.high());
  double ll = 0.05*len;
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
      of << xyz2 << " " << xyz2 + ll*norm << std::endl;
    }
  
  of << "410 1 0 4 0 100  155 255" << std::endl;
  of << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      Vector3D xyz = group_points_[kr]->getPoint();
      Point xyz2(xyz[0], xyz[1], xyz[2]);
      Point vec = group_points_[kr]->minCurvatureVec();
      of << xyz2 << " " << xyz2 + ll*vec << std::endl;
    }

  double lambda[2];
  Point eigen1, eigen2, eigen3;
  getPCA(lambda, eigen1, eigen2, eigen3);
  std::ofstream of2("points_two_tangents.txt");
  of2 << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      Vector3D xyz = group_points_[kr]->getPoint();
      Point norm = group_points_[kr]->getMongeNormal();
      Point vec1 = eigen1 - (eigen1*norm)*norm;
      vec1.normalize();
      Point vec2 = norm.cross(vec1);
      vec2.normalize();
      of2 << xyz << " " << vec1 << " " << vec2 << std::endl;
    }
}


void RevEngRegion::writeUnitSphereInfo(std::ostream& of)
{
  of << "400 1 0 4 100  0 155 255" << std::endl;
  of << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(group_points_[kr]);
      Point norm = pt->getMongeNormal();
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
