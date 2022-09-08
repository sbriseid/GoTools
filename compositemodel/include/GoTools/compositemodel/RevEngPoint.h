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

#ifndef _REVENGPOINT_H
#define _REVENGPOINT_H

#include "GoTools/compositemodel/ftPointSet.h"
//#include "GoTools/compositemodel/RevEngRegion.h"
#include "GoTools/utils/Array.h"
#include "GoTools/utils/Point.h"
#include "GoTools/utils/DirectionCone.h"

namespace Go
{
  class RevEngRegion;
  
  // Enumerations representing classification results
  enum
  {
   PCA_EDGE_UNDEF, PCA_NOT_EDGE, PCA_CLOSE_EDGE, PCA_EDGE
  };

  enum
  {
   C1_EDGE_UNDEF, C1_NOT_EDGE, C1_CLOSE_EDGE, C1_EDGE
  };

  enum
  {
   C2_EDGE_UNDEF, C2_NOT_EDGE, C2_CLOSE_EDGE, C2_EDGE
  };

  enum
  {
   PCA_UNDEF, PCA_PLANAR, PCA_LINEAR, PCA_OTHER
  };

  enum
  {
   C1_UNDEF, C1_PEAK, C1_RIDGE, C1_SRIDGE, C1_NONE, C1_FLAT, C1_MINSURF, C1_PIT, C1_VALLEY, C1_SVALLEY 
  };
    
  enum
  {
   SI_UNDEF, SI_PLANE, SI_SCUP, SI_TRO, SI_RUT, SI_SRUT, SI_SAD, SI_SRID, SI_RID, SI_DOM, SI_SCAP
   // Plane, spherical cup, trough, rut, saddle rut, saddle, saddle ridge,
   // ridge, dome, spherical cap
  };

  /** Subclass of ftSamplePoint; Enhanced point for use in reverse engineering
   */
  class RevEngPoint : public ftSamplePoint
  {

  public:
    /// Constructor
    RevEngPoint();
    
   RevEngPoint(Vector3D xyz, int bnd);

    virtual ~RevEngPoint();

    void computeTriangNormal();

    void addCovarianceEigen(Point& eigen1, double lambda1, Point& eigen2,
			    double lambda2, Point& eigen3, double lambda3);

    void addMongeInfo(Point& norm, Point& mincvec, double minc, Point& maxcvec,
		      double maxc, double currdist, double avdist, double eps);

    const Point& getTriangNormal()
    {
      return normalcone_.centre();
    }

    double getTriangAngle()
    {
      return normalcone_.greaterThanPi() ? 2.0*M_PI : normalcone_.angle();
    }

    double getMeanEdgLen();

    const Point& getMongeNormal()
    {
      return Mongenormal_;
    }

    Point fetchClosePoints(double radius, int min_nmb,
			   std::vector<Point>& nearpts);


    void fetchClosePoints2(double radius, int min_nmb,
			   std::vector<RevEngPoint*>& nearpts);


    void setVisited()
    {
      visited_ = 1;
    }

    
    void unsetVisited()
    {
      visited_ = 0;
    }

    bool visited()
    {
      return (visited_ > 0);
    }

    double getSurfaceVariation()
    {
      return sfvariation_;
    }

    double getShapeIndex()
    {
      return shapeindex_;
    }

    double maxPrincipalCurvature()
    {
      return kmax_;
    }

    double minPrincipalCurvature()
    {
      return kmin_;
    }

    double GaussCurvature()
    {
      return gausscurv_;
    }

    double meanCurvature()
      {
	return meancurv_;
      }

    double GaussCurvature0()
    {
      return gausscurv0_;
    }

    double meanCurvature0()
      {
	return meancurv0_;
      }

    void setMeanCurvature(double mean)
    {
      meancurv_ = mean;
    }

    void setGaussCurvature(double Gauss)
    {
      gausscurv_ = Gauss;
    }

    void updateCurvature()
    {
      meancurv0_ = meancurv_;
      gausscurv0_ = gausscurv_;
    }

    double getCurvedness()
    {
      return curvedness_;
    }

    void setClassification(int ctype, int c1_edge, int si_type,
			   int c2_edge, int pca_edge)
    {
      edge_[0] = pca_edge;
      edge_[1] = c1_edge;
      edge_[2] = c2_edge;
      surf_[0] = PCA_UNDEF;
      surf_[1] = ctype;
      surf_[2] = si_type;
    }

    bool isEdge()
    {
      return (edge_[0] == PCA_EDGE || edge_[1] == C1_EDGE ||
	      edge_[2] == C2_EDGE);
    }

    bool closeEdge()
    {
      return (edge_[0] >= PCA_CLOSE_EDGE || edge_[1] >= C1_CLOSE_EDGE ||
	      edge_[2] >= C2_CLOSE_EDGE);
    }

    bool notEdge()
    {
      return (edge_[0] <= PCA_NOT_EDGE && edge_[1] <= C1_NOT_EDGE &&
	      edge_[2] <= C2_NOT_EDGE);
    }

    bool isolatedEdge();

    void setEdgeUndef()
    {
      edge_[0] = PCA_EDGE_UNDEF;
      edge_[1] = C1_EDGE_UNDEF;
      edge_[2] = C2_EDGE_UNDEF;
    }

    int surfaceClassification(int classification_type);
    
    void adjustWithTriangNorm(double anglim);

    bool hasRegion()
    {
      return (region_ != 0);
    }

    RevEngRegion* region()
    {
      return region_;
    }

    void setRegion(RevEngRegion* region)
    {
      region_ = region;
    }

    void store(std::ostream& os) const;
    void read(std::istream& is, double eps, vector<int>& next_ix);
    
  private:
    double avedglen_;
    
    // How to mark that the properties are not computed? Can check dimension
    // of eigenvectors and principal curvature vectors. Compute derived info
    // once the primary one is in place
    Point eigen1_, eigen2_, eigen3_;  // Eigenvectors of covariance matrix
    double lambda1_, lambda2_, lambda3_; // Eigenvalues of covariance matrix
    // lambda1_>= lambda2_ >= lambda3_
    Point Mongenormal_, kvecmin_, kvecmax_;  // Normal and principal curvature vectors
    double kmin_, kmax_;  // Principal curvatures from Monge's patch
    double ptdist_, avdist_;
    DirectionCone normalcone_;  // Span of surface normals computed from triangulation

    double meancurv0_, meancurv_;
    double gausscurv0_, gausscurv_;
    double sfvariation_;  // Surface variation computed from eigenvalues
    double curvedness_;   // Curvedness computed from principal curvatures
    double shapeindex_;   // Computed from principal curvatures

    // Parameters corresponding to classification results
    int edge_[3];   // Results of edge classification.
    // Sequence: Surface variation (PCA), curvature (C1), curveness (C2)
    int surf_[3];   // Results of surface classification
    // Sequence: PCA, curvature (C1) shape index (SI)
    
    // Group (segment) of classified points
    RevEngRegion* region_;

    mutable int visited_;

    void getNearby(Vector3D xyz, double radius,
		   std::vector<RevEngPoint*>& near);
  };
}

#endif
