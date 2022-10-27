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


namespace Go
{
  class HedgeSurface;
  //class RevEngPoint;

  enum
    {
     CLASSIFICATION_UNDEF, CLASSIFICATION_CURVATURE, CLASSIFICATION_SHAPEINDEX
    };
  
  class RevEngRegion
  {
  public:
    // Constructor
    RevEngRegion();

    RevEngRegion(int classification_type);

    RevEngRegion(int classification_type,
		 std::vector<RevEngPoint*> & points);

     int getClassification()
    {
      return classification_type_;
    }

    bool isCompatible(ClassType classtype, int sfcode);
    
    // Extend group
    void addPoint(RevEngPoint* point);

    // Extend region with adjacent points having the same classification
    void collect(RevEngPoint *pt);

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
    
    RevEngPoint* seedPoint();
    
    void growLocal(RevEngPoint* seed, double tol, double radius, int min_close,
		   std::vector<RevEngPoint*>& out);

    void growWithSurf(int max_nmb, double tol,
		      std::vector<RevEngRegion*>& grown_regions,
		      std::vector<HedgeSurface*>& adj_surfs);

    void
    splitFromSurfaceNormals(std::vector<RevEngPoint*>& smallrad,
			    std::vector<std::vector<RevEngPoint*> >& separate_group);
    void splitRegion(std::vector<std::vector<RevEngPoint*> >& separate_groups);

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
    
    bool extractPlane(double tol, int min_pt,
		      std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
		      std::vector<HedgeSurface*>& prevsfs,
		      std::ostream& fileout);

    bool extractCylinder(double tol, int min_pt, double mean_edge_len,
			 std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			 std::vector<HedgeSurface*>& prevsfs,
			 std::ostream& fileout);

    bool extractTorus(double tol, int min_pt, double mean_edge_len,
		      std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
		      std::vector<HedgeSurface*>& prevsfs,
		      std::ostream& fileout);

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

    void setVisited(bool visited)
    {
      visited_ = visited;
    }

    bool visited()
    {
      return visited_;
    }

  private:
    vector<RevEngPoint*> group_points_;   // Points belonging to classified segment
    int classification_type_;
    std::vector<HedgeSurface*> associated_sf_;  // Can be two due to split along
    // seam of closed surface (should be fixed)
    shared_ptr<ImplicitApprox> impl_;
    double mink1_, maxk1_, mink2_, maxk2_;
    BoundingBox bbox_;
    DirectionCone normalcone_;
    double maxdist_, avdist_;
    int num_inside_;

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
  };
}

#endif
