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

#ifndef _REVENGREGION_H
#define _REVENGREGION_H

#include "GoTools/compositemodel/RevEngPoint.h"
#include "GoTools/compositemodel/ImplicitApprox.h"
#include "GoTools/utils/BoundingBox.h"
//#include "GoTools/compositemodel/HedgeSurface.h"
#include <set>


namespace Go
{
  class HedgeSurface;
  //class RevEngPoint;
  class Circle;
  class SplineCurve;
  class Plane;
  class Cylinder;
  class Sphere;
  class Cone;
  class Torus;
  class SplneSurface;
  
  enum
    {
     CLASSIFICATION_UNDEF, CLASSIFICATION_CURVATURE, CLASSIFICATION_SHAPEINDEX, CLASSIFICATION_POINTASSOCIATION
    };

  struct SweepData
  {
    int type_;  // Linear = 1, rotational = 2
    shared_ptr<SplineCurve> profile_;
    Point location_;
    Point added_info_;
    double maxdist_;
    double avdist_;
    int num_in_;

    SweepData(int type, shared_ptr<SplineCurve> profile, Point location, 
	      Point info2, double maxdist, double avdist, int num_in)
    {
      type_ = type;
      profile_ = profile;
      location_ = location;
      added_info_ = info2;
      maxdist_ = maxdist;
      avdist_ = avdist;
      num_in_ = num_in;
    }
  };

  
  class RevEngRegion
  {
  public:
    // Constructor
    RevEngRegion();

    RevEngRegion(int classification_type);

    RevEngRegion(int classification_type,
		 std::vector<RevEngPoint*> & points);

    ~RevEngRegion();
    
     int getClassification()
    {
      return classification_type_;
    }

    bool isCompatible(ClassType classtype, int sfcode);
    
    // Extend group
    void addPoint(RevEngPoint* point);

    // Remove. NB! Leaves overview information invalid.
    void removePoint(RevEngPoint* point);

    void addRegion(RevEngRegion* reg,
		   std::vector<std::pair<double, double> >& dist_ang,
		   double maxd=0.0, double avd=0.0, int num_inside=-1);
    
    // Update overview information
    void updateInfo();

    // Extend region with adjacent points having the same classification
    void collect(RevEngPoint *pt, RevEngRegion* prev=0);

    int numPoints()
    {
      return (int)group_points_.size();
    }

    RevEngPoint* getPoint(int ix)
    {
      if (ix < 0 || ix >= (int)group_points_.size())
	return 0;
      else
	return group_points_[ix];
    }

    std::vector<RevEngPoint*>::iterator pointsBegin()
    {
      return group_points_.begin();
    }
    
    std::vector<RevEngPoint*>::iterator pointsEnd()
    {
      return group_points_.end();
    }
    
    std::vector<RevEngPoint*> getPoints()
    {
      return group_points_;
    }

    const BoundingBox& getBbox()
    {
      return bbox_;
    }

    BoundingBox getParameterBox();
    
    RevEngPoint* seedPoint();
    
    // void growLocal(RevEngPoint* seed, double tol, double radius, int min_close,
    // 		   std::vector<RevEngPoint*>& out);

    void growWithSurf(int max_nmb, double tol,
		      std::vector<RevEngRegion*>& grown_regions,
		      std::vector<HedgeSurface*>& adj_surfs);

    void
    splitComposedRegions(int classtype,
			 std::vector<shared_ptr<RevEngRegion> >& added_groups);
    
    void
    splitWithShapeIndex(std::vector<shared_ptr<RevEngRegion> >& updated_regions);
    
    void
    splitFromSurfaceNormals(std::vector<RevEngPoint*>& smallrad,
			    std::vector<std::vector<RevEngPoint*> >& separate_group);
    void splitRegion(std::vector<std::vector<RevEngPoint*> >& separate_groups);

    void updateRegion(double approx_tol, double anglim,
		      std::vector<RevEngRegion*>& adapted_regions,
		      std::vector<shared_ptr<RevEngRegion> >& outdiv_regions);

     void joinRegions(double approx_tol, double anglim,
		     std::vector<RevEngRegion*>& adapted_regions);

    void extractOutPoints(std::vector<std::pair<double, double> >& dist_ang,
			  double tol,
			  std::vector<std::vector<RevEngPoint*> >& out_groups);
    
   bool cylindertype()
    {
      if (classification_type_ == CLASSIFICATION_CURVATURE)
	return (group_points_[0]->C1_surf() == C1_RIDGE ||
		group_points_[0]->C1_surf() == C1_VALLEY);
      else if (classification_type_ == CLASSIFICATION_SHAPEINDEX)
	return (group_points_[0]->SI_surf() == SI_RUT ||
		group_points_[0]->SI_surf() == SI_RID);
      else
	return false;
    }
    
   bool planartype()
    {
      if (classification_type_ == CLASSIFICATION_CURVATURE)
	return (group_points_[0]->C1_surf() == C1_FLAT);
      else if (classification_type_ == CLASSIFICATION_SHAPEINDEX)
	return (group_points_[0]->SI_surf() == SI_PLANE);
      else
	return false;
    }
    
    bool extractPlane(double tol, int min_pt, double angtol,
		      std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
		      std::vector<HedgeSurface*>& prevsfs,
		      std::vector<std::vector<RevEngPoint*> >& out_groups,
		      std::ostream& fileout);

    bool extractCylinder(double tol, int min_pt, double angtol,
			 double mean_edge_len,
			 std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			 std::vector<HedgeSurface*>& prevsfs,
			 std::vector<std::vector<RevEngPoint*> >& out_groups,
			 std::ostream& fileout);

    bool extractSphere(double tol, int min_pt, double angtol,
		       double mean_edge_len,
		       std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
		       std::vector<HedgeSurface*>& prevsfs,
		       std::vector<std::vector<RevEngPoint*> >& out_groups,
		       std::ostream& fileout);

    bool extractLinearSweep(double tol, int min_pt, double angtol,
			    std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			    std::vector<HedgeSurface*>& prevsfs,
			    std::ostream& fileout);

    bool extractCone(double tol, int min_pt, double angtol,
		     double mean_edge_len,
		     std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
		     std::vector<HedgeSurface*>& prevsfs,
		     std::vector<std::vector<RevEngPoint*> >& out_groups,
		     std::ostream& fileout);

    bool extractTorus(double tol, int min_pt, double angtol,
		      double mean_edge_len,
		      std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
		      std::vector<HedgeSurface*>& prevsfs,
		      std::vector<std::vector<RevEngPoint*> >& out_groups,
		      std::ostream& fileout);

    bool extractFreeform(double tol, int min_pt, double angtol,
			 double mean_edge_len,
			 std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			 std::vector<HedgeSurface*>& prevsfs,
			 std::vector<std::vector<RevEngPoint*> >& out_groups,
			 std::ostream& fileout);

    void implicitizeSplit();
    
    void setHedge(HedgeSurface* surface)
    {
      associated_sf_.clear();
      associated_sf_.push_back(surface);
    }

    bool hasSurface()
    {
      return (associated_sf_.size() > 0);
    }

    int numSurface()
    {
      return (int)associated_sf_.size();
    }

    HedgeSurface* getSurface(int ix)
    {
      return (ix < 0 || ix>=(int)associated_sf_.size()) ? 0 : associated_sf_[ix];
    }

    const BoundingBox& boundingBox()
    {
      return bbox_;
    }

    const DirectionCone& getNormalCone()
    {
      return normalcone_;
    }

    void getPrincipalCurvatureInfo(double& mink1, double& maxk1, double& mink2, double& maxk2)
    {
      mink1 = mink1_;
      maxk1 = maxk1_;
      mink2 = mink2_;
      maxk2 = maxk2_;
    }
    
    void setAccuracy(double maxdist, double avdist, int num_inside)
    {
      maxdist_ = maxdist;
      avdist_ = avdist;
      num_inside_ = num_inside;
    }

    void getAccuracy(double& maxdist, double& avdist, int& num_inside)
    {
      maxdist = maxdist_;
      avdist = avdist_;
      num_inside = num_inside_;
    }

    double getMaxSfDist()
    {
      return maxdist_;
    }

    void setVisited(bool visited)
    {
      visited_ = visited;
    }

    bool visited()
    {
      return visited_;
    }

    bool possiblePlane(double angtol, double inlim);
    bool possibleCylinder(double angtol, double inlim);
    bool possibleCone(double angtol, double inlim);
    bool possibleTorus(double angtol, double inlim);

    bool hasSweepInfo()
    {
      return (sweep_.get() != 0);
    }

    int sweepType()
    {
      return (sweep_.get() ? sweep_->type_ : 0);
    }

    bool hasDivideInfo()
    {
      return false;
    }

    void setRegionAdjacency();

    bool integrateInAdjacent(double mean_edge_len, int min_next,
			     int max_next, double tol, double angtol,
			     int max_nmb_outlier);

    void addAdjacentRegion(RevEngRegion* adj_reg)
    {
      //adj_reg->addAdjacentRegion(this);
      adjacent_regions_.insert(adj_reg);
    }
    
    void removeAdjacentRegion(RevEngRegion* adj_reg)
    {
      //adj_reg->removeAdjacentRegion(this);
      adjacent_regions_.erase(adj_reg);
    }
    
    void clearRegionAdjacency()
    {
      adjacent_regions_.clear();
    }
    
    bool hasAdjacentRegion(RevEngRegion* adj_reg)
    {
      return (adjacent_regions_.find(adj_reg) != adjacent_regions_.end());
    }
    
    void adjustWithSurf(double tol, double angtol);

    bool hasBaseSf()
    {
      return basesf_.get();
    }

    void setBaseSf(shared_ptr<ParamSurface> base, double maxd, double avd,
		   int num_in)
    {
      basesf_ = base;
      maxdist_base_ = maxd;
      avdist_base_ = avd;
      num_in_base_ = num_in;
    }

    void getBase(shared_ptr<ParamSurface>& base, double& maxd, double& avd,
		 int& num_in)
    {
      base = basesf_;
      maxd = maxdist_base_;
      avd = avdist_base_;
      num_in = num_in_base_;
    }
    
    void getBaseDist(double& maxd, double& avd, int& num_in)
    {
      maxd = maxdist_base_;
      avd = avdist_base_;
      num_in = num_in_base_;
    }

    bool hasEdgeBetween(RevEngRegion* adj);

    bool hasPrimary()
    {
      return (primary_.get() != 0);
    }

    shared_ptr<ParamSurface> getPrimary()
    {
      return primary_;
    }

    void getPrimaryInfo(double& maxdist_primary,
			double& avdist_primary, int& num_in_primary)
    {
      maxdist_primary = maxdist_primary_;
      avdist_primary = avdist_primary_;
      num_in_primary = num_in_primary_;
    }

    
    void writeRegionInfo(std::ostream& of);
    void writeUnitSphereInfo(std::ostream& of);
    void writeSubTriangulation(std::ostream& of);
    void writeSurface(std::ostream& of);

    void store(std::ostream& os) const;
    void read(std::istream& is, shared_ptr<ftPointSet>& tri_sf);
	      
  private:
    std::vector<RevEngPoint*> group_points_;   // Points belonging to classified segment
    int classification_type_;
    std::vector<HedgeSurface*> associated_sf_;  // Can be two due to split along
    // seam of closed surface (should be fixed)
    shared_ptr<ImplicitApprox> impl_;
    double mink1_, maxk1_, mink2_, maxk2_;
    BoundingBox bbox_;
    DirectionCone normalcone_;
    double maxdist_, avdist_;
    int num_inside_;
    std::set<RevEngRegion*> adjacent_regions_;
    shared_ptr<ParamSurface> basesf_;
    double maxdist_base_, avdist_base_;
    int num_in_base_;
    shared_ptr<ParamSurface> primary_;
    double maxdist_primary_, avdist_primary_;
    int num_in_primary_;

    shared_ptr<SweepData> sweep_;
    bool visited_;

    const Point& pluckerAxis();
    void extendWithGaussRad();
    void extendWithGaussRad2();
    void analyseNormals(double tol, Point& normal, Point& centre, double& radius);
    void analysePlaneProperties(Point avnorm, double angtol,
				std::vector<RevEngPoint*>& in,
				std::vector<RevEngPoint*> out);
    void analyseCylinderProperties(Point avvec, double angtol,
				   std::vector<RevEngPoint*>& in,
				   std::vector<RevEngPoint*> out);
    void  curveApprox(std::vector<Point>& points, double tol,
		      shared_ptr<Circle> circle,
		      shared_ptr<SplineCurve>& curve, Point& xpos);
    void  curveApprox(std::vector<Point>& points,
		      shared_ptr<ParamCurve> cvin,
		      int ik, int in, 
		      shared_ptr<SplineCurve>& curve);
    shared_ptr<Plane> computePlane(std::vector<RevEngPoint*>& points);
    shared_ptr<Cylinder> computeCylinder(std::vector<RevEngPoint*>& points,
					 double tol);
    shared_ptr<Sphere> computeSphere(std::vector<RevEngPoint*>& points);
    shared_ptr<Cone> computeCone(std::vector<RevEngPoint*>& points, Point& apex);
    shared_ptr<Torus> computeTorus(std::vector<RevEngPoint*>& points,
				   double tol, shared_ptr<Torus>& torus2);
    shared_ptr<SplineSurface> computeLinearSwept(double tol, shared_ptr<SplineCurve>& profile,
						 Point& pt1, Point& pt2);
    shared_ptr<SplineSurface> computeFreeform(double tol);
    shared_ptr<SplineSurface> updateFreeform(double tol);
    void getPCA(double lambda[3], Point& eigen1, Point& eigen2, Point& eigen3);
    shared_ptr<SplineSurface> surfApprox(vector<RevEngPoint*>& points,
					 const BoundingBox& bbox);
    void splitCylinderRad(const Point& pos, const Point& axis,
			  const Point& Cx, const Point& Cy,
			  int nmb_split, std::vector<Point>& centr,
			  std::vector<double>& rad);
    void approximationAccuracy(std::vector<RevEngPoint*>& points,
			       shared_ptr<ParamSurface> surf,
			       double tol, double angtol,
			       double& maxd, double& avd,
			       std::vector<RevEngPoint*>& in,
			       std::vector<RevEngPoint*>& out);
    void checkReplaceSurf(double tol);
    bool parameterizeOnSurf(shared_ptr<ParamSurface> surf,
			    std::vector<double>& data,
			    std::vector<double>& param,
			    int& inner1, int& inner2, bool& close1, bool& close2);
    void setPrimarySf(shared_ptr<ParamSurface> surf, double maxd, double avd,
		      int num_in)
    {
      primary_ = surf;
      maxdist_primary_ = maxd;
      avdist_primary_ = avd;
      num_in_primary_ = num_in;
    }

    bool reparameterize(std::vector<double>& param, std::vector<double>& param2,
			double& umin, double& umax, double& vmin, double& vmax);

    void getParExtent(double curr[2], int pdir, std::vector<std::vector<int> >& raster,
		      int& i1, int& i2);

    void defineRaster(std::vector<double>& param, int nmb_div,
		      std::vector<std::vector<int> >& raster, double& umin,
		      double& umax, double& vmin, double& vmax);

    void extendInCorner(std::vector<double>& data, std::vector<double>& param,
			double umin, double umax, double vmin, double vmax);
  };
}

#endif
