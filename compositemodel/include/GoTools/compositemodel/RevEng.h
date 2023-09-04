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

#ifndef _REVENG_H
#define _REVENG_H

#include "GoTools/compositemodel/ftPointSet.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/RevEngRegion.h"
#include "GoTools/utils/Point.h"
#include "GoTools/utils/BoundingBox.h"

namespace Go
{
  class RevEngPoint;
  //class RevEngRegion;
  class HedgeSurface;

  // Elementary surface types to recognize (omitting currently
  // ellipsoid, elliptic cylinder, ...)
  enum
  {
   PLANE, CYLINDER, SPHERE, CONE, TORUS
  };

  struct SurfaceProperties
  {
    int sfix_;
    ClassType type_;
    ClassType prev_type_;
    Point dir_, loc_;
    double rad1_, rad2_;
    int num_points_;

    SurfaceProperties(int ix, ClassType type, int num, Point& dir, Point& loc,
		      ClassType prev_type=Class_Unknown, double rad1=-1,
		      double rad2=-1)
    {
      sfix_ = ix;
      type_ = type;
      num_points_ = num;
      dir_ = dir;
      loc_ = loc;
      prev_type_ = prev_type;
      rad1_ = rad1;
      rad2_ = rad2;
    }
  };
  
  class RevEng
  {
  public:
    RevEng();
    
    RevEng(shared_ptr<ftPointSet> tri_sf, double mean_edge_len);

     ~RevEng();

    // Should this class have an option to run all operations in one
     // sequence without being started from outside?
      
      void enhancePoints();
      void enhancePoints2();

    void edgeClassification();
      void classifyPoints();
      
      void growRegions();

      void recognizeElementary();
      
      // void recognizePlanes();
      
      // void recognizeCylinders();
      
       void mergePlanes(size_t first, size_t last);

    // Could be private
    void mergeCylinders(size_t first, size_t last);

      void trimSurfaces();

    shared_ptr<SurfaceModel> createModel();

    void storeClassified(std::ostream& os) const;
    void readClassified(std::istream& is);
    void storeGrownRegions(std::ostream& os) const;
    void readGrownRegions(std::istream& is);
    void curvatureFilter();

    double getInitApproxTol();
    void setApproxTolerance();
    void setApproxTol(double eps)
    {
      approx_tol_= eps;
    }
    
    double getApproxTol()
    {
      return approx_tol_;
    }

    void setEdgeClassificationParams();
    int getEdgeClassificationType()
    {
      return edge_class_type_;
    }
    
    void setEdgeClassificationType(int edge_class_type)
    {
      edge_class_type_ = edge_class_type;
    }
    
    double getCfac()
    {
      return cfac_;
    }
    void setCfac(double cfac)
    {
      cfac_ = cfac;
    }

    double getPCAlim();
    void setPCAlim(double pca_lim)
    {
      pca_lim_ = pca_lim;
    }
    
    double getCnesslim();
    void setCnesslim(double cness_lim)
    {
      cness_lim_ = cness_lim;
    }
    
    double getRPfac()
    {
      return rpfac_;
    }
    void setRPfac(double rpfac)
    {
      rpfac_ = rpfac;
    }
    
    void setClassificationParams();
    int getClassificationType()
    {
      return classification_type_;
    }
    void setClassificationType(int classification_type)
    {
      classification_type_ = classification_type;
    }

    double getMeanCurvatureZero()
    {
      return zero_H_;
    }

    void setMeanCurvatureZero(double zero_H)
    {
      zero_H_ = zero_H;
    }

    
    double getGaussCurvatureZero()
    {
      return zero_K_;
    }

    void setGaussCurvatureZero(double zero_K)
    {
      zero_K_ = zero_K;
    }

    double getShapeIndexZero()
    {
      return zero_si_;
    }

    void setShapeIndexZero(double zero_si)
    {
      zero_si_ = zero_si;
    }

    int getElementaryPreferLevel()
    {
      return prefer_elementary_;
    }

    void setElementaryPreferLevel(int preferlevel)
    {
      prefer_elementary_ = preferlevel;
    }
    
    

  private:
    shared_ptr<ftPointSet> tri_sf_;
    double mean_edge_len_;
    std::vector<shared_ptr<RevEngRegion> > regions_;
    std::vector<shared_ptr<HedgeSurface> > surfaces_;  // I think the 
    // surfaces must be collected here to have a stable storage
    // The surfaces can be freeform as well as primary. The collection
    // will be build gradually. The number of surfaces will increase and
    // decrease based on recognition, merging and splitting by trimming
    shared_ptr<SurfaceModel> sfmodel_;
    BoundingBox bbox_;
    int min_next_;  // Minimum number of neighbouring points
    int max_next_;  // Estimate for maximum number of neighbouring points
    double rfac_;   // Factor for radius in which to search for neighbouring points
    int edge_class_type_ = CNESS_EDGE;
    int classification_type_ = CLASSIFICATION_CURVATURE;
    double cfac_;   // Edge points from curvature is given by
    // cfac_ times the average length of triangulation edges in a vertex
    double pca_lim_; // Limit for edge classification from surface variation
    double cness_lim_; // Limit for edge classification from curvedness
    double norm_ang_lim_; // Limit for when the cone angle corresponding
    // to triangle normals indicate an edge
    double norm_plane_lim_;  // Limit for when the cone angle corresponding
    // to triangle normals indicate a plane
    double zero_H_;  // When mean curvature is considered zero
    double zero_K_;  // When Gauss curvature is considered zero
    double zero_si_; // When shape index is considered zero
    int min_point_region_;
    double approx_tol_;  // Approximation tolerance in region growing
    double int_tol_;  // Intersection tolerance
    double anglim_;
    int max_nmb_outlier_;
    int rpix_;
    double rpfac_, ffac_, sfac_;

    int prefer_elementary_; // 0 = always, 1 = preferred, 2 = best accuracy

    Point dirvec_[3];
    
    void initParameters();
    void growSurface(size_t& ix);
    void mergeSurfaces();
    void mergeSplineSurfaces();
    shared_ptr<HedgeSurface> doMerge(std::vector<size_t>& cand_ix,
				     std::vector<double>& cand_score,
				     ClassType type);
    shared_ptr<ParamSurface> approxMergeSet(std::vector<size_t>& cand_ix,
					    std::vector<size_t>& select_ix,
					    ClassType type);
    
    shared_ptr<ParamSurface> doMergePlanes(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
					   std::vector<RevEngPoint*>::iterator> > points,
					   const BoundingBox& bbox,
					   std::vector<int>& nmbpts,
					   bool set_bound = true);
    shared_ptr<ParamSurface> doMergeCylinders(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
					      std::vector<RevEngPoint*>::iterator> > points,
					      const BoundingBox& bbox,
					      std::vector<int>& nmbpts,
					      bool set_bound = true);
    shared_ptr<ParamSurface> doMergeTorus(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
					  std::vector<RevEngPoint*>::iterator> > points,
					  const BoundingBox& bbox,
					  std::vector<int>& nmbpts);

    void adaptToMainAxis();
    void cylinderFit(std::vector<int>& sf_ix, Point normal);
    
    void collectAxis(std::vector<SurfaceProperties>& sfprop);

    shared_ptr<HedgeSurface> doMergeSpline(shared_ptr<HedgeSurface>& surf1,
					   DirectionCone& cone1,
					   shared_ptr<HedgeSurface>& surf2,
					   DirectionCone& cone2);

    double computeCylRadius(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
				std::vector<RevEngPoint*>::iterator> >& points,
			    Point mid, Point vec1, Point vec2);

    void computeMonge(RevEngPoint* pt, std::vector<RevEngPoint*>& points,
		      Point& vec1, Point& vec2, double radius);

    void setRp(RevEngPoint* first, double rp[2]);
    double computeRp(RevEngPoint* first, std::vector<std::vector<int> >& tri);
    int setSmallRegionNumber();
    
    void storeParams(std::ostream& os) const;
    void readParams(std::istream& is);
    
    void writeRegionStage(std::ostream& of, std::ostream& ofs) const;
  };

} // namespace Go

#endif // _REVENG_H
