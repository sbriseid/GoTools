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

#include "GoTools/compositemodel/RevEng.h"
#include "GoTools/compositemodel/RevEngPoint.h"
#include "GoTools/compositemodel/RevEngRegion.h"
#include "GoTools/compositemodel/RevEngUtils.h"
#include "GoTools/compositemodel/HedgeSurface.h"
#include "GoTools/utils/DirectionCone.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/BoundedSurface.h"
#include <vector>
#include <fstream>
#include <iostream> // @@ debug

using namespace Go;
using std::vector;
using std::pair;
using std::istream;
using std::ostream;


#define MAX_COLORS 12
int colors[MAX_COLORS][3] = {
  {255, 0, 0},
  {0, 255, 0},
  {0, 0, 255},
  {255, 255, 0},
  {255, 0, 255},
  {0, 255, 255},
  {128, 255, 0},
  {255, 128, 0},
  {128, 0, 255},
  {255, 0, 128},
  {0, 128, 255},
  {0, 255, 128},
};

//===========================================================================
RevEng::RevEng(shared_ptr<ftPointSet> tri_sf, double mean_edge_len)
  : tri_sf_(tri_sf), mean_edge_len_(mean_edge_len)
//===========================================================================
{
  // Set default parameters
  initParameters();
}


//===========================================================================
RevEng::RevEng()
//===========================================================================
{
  // Empty infrastructure for reading stage
  initParameters();
}


//===========================================================================
RevEng::~RevEng()
//===========================================================================
{

}


//===========================================================================
void RevEng::enhancePoints()
//===========================================================================
{
  int writepoints = 0;
  int nmbpt = tri_sf_->size();
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]); 

      // Compute surface normal from triangulation
      pt->computeTriangNormal();

      double avlen = pt->getMeanEdgLen();

      //Fetch nearby points
      vector<Point> nearpts;
      double local_len = pt->getMeanEdgLen();
      double radius = 0.5*rfac_*(local_len + mean_edge_len_);
      Point curr = pt->fetchClosePoints(radius, min_next_, nearpts);

      std::ofstream of("nearpts.g2");
      if (writepoints)
	{
	  of << "400 1 0 4 255 0 0 255" << std::endl;
	  of << "1" << std::endl;
	  of << curr << std::endl << std::endl;
	  of << "400 1 0 4 0 255 0 255" << std::endl;
	  of << nearpts.size() << std::endl;
	  for (size_t kr=0; kr<nearpts.size(); ++kr)
	    of << nearpts[kr] << std::endl;
	}
      
      // Compute eigenvectors and values of covariance matrix
      double lambda[3];
      double eigenvec[3][3];
      RevEngUtils::principalAnalysis(curr, nearpts, lambda, eigenvec);
      Point eigen1(eigenvec[0][0], eigenvec[0][1], eigenvec[0][2]);
      Point eigen2(eigenvec[1][0], eigenvec[1][1], eigenvec[1][2]);
      Point eigen3(eigenvec[2][0], eigenvec[2][1], eigenvec[2][2]);
      if (eigen3*pt->getTriangNormal() <  0.0)
	{
	  eigen2 *= -1;
	  eigen3 *= -1;
	}
      pt->addCovarianceEigen(eigen1, lambda[0], eigen2, lambda[1],
			     eigen3, lambda[2]);
      
      if (writepoints)
	{
	  for (int ki=0; ki<3; ++ki)
	    {
	      Point vec(eigenvec[ki][0], eigenvec[ki][1], eigenvec[ki][2]);
	      of << "410 1 0 4 0 0 0 255" << std::endl;
	      of << "1" << std::endl;
	      of << curr << " " << curr+0.1*vec << std::endl;
	    }
	}
      // Compute normal and curvature using Monge patch
      Point normal, mincvec, maxcvec;
      double minc, maxc;
      double currdist, avdist;
      // RevEngUtils::computeMonge(curr, nearpts, eigen1, eigen2, normal, mincvec, minc,
      // 				maxcvec, maxc, currdist, avdist);
      RevEngUtils::computeMonge(curr, nearpts, eigen1, eigen3, normal, mincvec, minc,
				maxcvec, maxc, currdist, avdist);

      // Orient vectors with respect to triangulation normal
      // The normal vectors should be OK. Curvature vectors are not necessarily
      // consistent with regard to orientation

      pt->addMongeInfo(normal, mincvec, minc, maxcvec, maxc, currdist, avdist,
		       zero_si_);
      int stop_break = 1;
    }

  std::cout << "Start curvature filter" << std::endl;
  
  curvatureFilter();

  std::cout << "Finish curvature filter" << std::endl;
  
  vector<Vector3D> triangcorners;
  vector<Vector3D> triangplane;
  vector<Vector3D> curvaturecorners;
  std::ofstream of("triangnorm.g2");
  std::ofstream ofM("Mongenorm.g2");
  of << "410 1 0 4 0 0 255 255" << std::endl;
  of << nmbpt << std::endl;
  ofM << "410 1 0 4 255 0 0 255" << std::endl;
  ofM << nmbpt << std::endl;
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]); 
      Vector3D xyz = pt->getPoint();
      Point xyz2(xyz[0], xyz[1], xyz[2]);
      Point norm = pt->getTriangNormal();
      double ang = pt->getTriangAngle();
      of << xyz2 << " " << xyz2+0.1*norm << std::endl;
      if (pt->getTriangAngle() > norm_ang_lim_)
	triangcorners.push_back(xyz);
      Point Mnorm = pt->getMongeNormal();
      ofM << xyz2 << " " << xyz2+Mnorm << std::endl;

      if (ang <= norm_plane_lim_)
	triangplane.push_back(xyz);
      double avlen = pt->getMeanEdgLen();
      double maxpc = std::max(fabs(pt->maxPrincipalCurvature()),
			      fabs(pt->minPrincipalCurvature()));
      double crvrad = 1.0/maxpc; //fabs(maxpc);
      if (crvrad < cfac_*avlen)
	curvaturecorners.push_back(xyz);
    }

  std::ofstream of2("triangcorners.g2");
  of2 << "400 1 0 4 0 0 0 255" << std::endl;
  of2 << triangcorners.size() << std::endl;
  for (size_t kj=0; kj<triangcorners.size(); ++kj)
    of2 << triangcorners[kj] << std::endl;

  std::ofstream of3("curvaturecorners.g2");
  of3 << "400 1 0 4 10 10 10 255" << std::endl;
  of3 << curvaturecorners.size() << std::endl;
  for (size_t kj=0; kj<curvaturecorners.size(); ++kj)
    of3 << curvaturecorners[kj] << std::endl;

  std::ofstream of4("triangplane.g2");
  of4 << "400 1 0 4 200 0 200 255" << std::endl;
  of4 << triangplane.size() << std::endl;
  for (size_t kj=0; kj<triangplane.size(); ++kj)
    of4 << triangplane[kj] << std::endl;
}

//===========================================================================
void RevEng::curvatureFilter()
//===========================================================================
{
  int writepoints = 0;
  int nmbpt = tri_sf_->size();
  double radius_fac = 0.7;
  vector<vector<RevEngPoint*> > nearpts(nmbpt);
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      
      // Fetch nearby points
      double local_len = pt->getMeanEdgLen();
      double radius = 0.5*radius_fac*rfac_*(local_len + mean_edge_len_);
      pt->fetchClosePoints2(radius, min_next_, nearpts[ki]);
      vector<double> H0(nearpts[ki].size()+1), K0(nearpts[ki].size()+1);
      H0[0] = pt->meanCurvature0();
      K0[0] = pt->GaussCurvature0();
      for (size_t kj=0; kj<nearpts[ki].size(); ++kj)
	{
	  H0[kj+1] = nearpts[ki][kj]->meanCurvature0();
	  K0[kj+1] = nearpts[ki][kj]->GaussCurvature0();
	}
      std::sort(H0.begin(), H0.end());
      std::sort(K0.begin(), K0.end());
      pt->setMeanCurvature(0.5*H0[H0.size()/2] + H0[(H0.size()+1)/2]);
      pt->setGaussCurvature(0.5*K0[K0.size()/2] + K0[(K0.size()+1)/2]);
    }
  
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      pt->updateCurvature();
    }

  int nmbsmooth = 2;
  for (int ka=0; ka<nmbsmooth; ++ka)
    {
      for (int ki=0; ki<nmbpt; ++ki)
	{
	  RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      
	  // Fetch nearby points
	  double Hmean = pt->meanCurvature0();
	  double Kmean = pt->GaussCurvature0();
	  for (size_t kj=0; kj<nearpts[ki].size(); ++kj)
	    {
	      Hmean += nearpts[ki][kj]->meanCurvature0();
	      Kmean += nearpts[ki][kj]->GaussCurvature0();
	    }
	  Hmean /= (double)(nearpts[ki].size()+1);
	  Kmean /= (double)(nearpts[ki].size()+1);
	  pt->setMeanCurvature(Hmean);
	  pt->setGaussCurvature(Kmean);
	}
      for (int ki=0; ki<nmbpt; ++ki)
	{
	  RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
	  pt->updateCurvature();
	}
   }
}

//===========================================================================
void RevEng::setClassificationParams()
//===========================================================================
{
  int nmbpts = tri_sf_->size();
  
  vector<double> sfvar(nmbpts);
  vector<double> curvedness(nmbpts);
  double varmean = 0.0, curvedmean = 0.0;
  double varh = 0.0, curvh = 0.0;
  double hfac = 6.0/(nmbpts*7.0);
  for (int ki=0; ki<nmbpts; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]); 
      sfvar[ki] = pt->getSurfaceVariation();
      curvedness[ki] = pt->getCurvedness();
      varmean += sfvar[ki]/(double)nmbpts;
      curvedmean += curvedness[ki]/(double)nmbpts;
      varh += hfac*sfvar[ki];
      curvh += hfac*curvedness[ki];
    }
  pca_lim_ = varh;
  cness_lim_ = curvh;

  // // Surface variation
  // std::sort(sfvar.begin(), sfvar.end());
  // double varmed = 0.5*(sfvar[nmbpts/2]+sfvar[(nmbpts+1)/2]);
  // double varQ1 = sfvar[nmbpts/4];
  // double varQ2 = sfvar[3*nmbpts/4];

  // double stepmean = (sfvar[sfvar.size() - 1] - sfvar[0])/(double)(nmbpts-1);
  // double stepfac = 20.0; //50.0; //100.0;
  // double steplim1 = stepfac*stepmean;
  // vector<int> largestep;
  // for (int ki=0; ki<nmbpts-1; ++ki)
  //   {
  //     double diff = sfvar[ki+1] - sfvar[ki];
  //     if (diff > steplim1)
  // 	largestep.push_back(ki);
  //   }

  // // Curvedness
  // std::sort(curvedness.begin(), curvedness.end());
  // double stepmean2 = (curvedness[curvedness.size() - 1] - curvedness[0])/(double)(nmbpts-1);
  // double steplim2 = stepfac*stepmean2;
  // double curvedmed = 0.5*(curvedness[nmbpts/2] + curvedness[nmbpts/2+1]);
  // double curvedQ1 = curvedness[nmbpts/4];
  // double curvedQ2 = curvedness[3*nmbpts/4];
  // vector<int> largecurved;
  // for (int ki=0; ki<nmbpts-1; ++ki)
  //   {
  //     double diff = curvedness[ki+1] - curvedness[ki];
  //     if (diff > steplim2)
  // 	largecurved.push_back(ki);
  //   }

  // Testing
  // std::ofstream of1("sfvar_1.g2");
  // of1 << "400 1 0 4 155 0 100 255" << std::endl;
  // of1 << nmbpts-largestep[0]-1 << std::endl;
  // std::ofstream of2("curvedness_1.g2");
  // of2 << "400 1 0 4 0 155 100 255" << std::endl;
  // of2 << nmbpts-largecurved[0]-1 << std::endl;
  vector<Vector3D> pts_v, pts_c;
  for (int ka=0; ka<nmbpts; ++ka)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ka]);
      double var = pt->getSurfaceVariation();
      double curved = pt->getCurvedness();
      Vector3D xyz = pt->getPoint();
      // if (var > sfvar[largestep[0]])
      // 	of1 << xyz << std::endl;
      // if (curved > curvedness[largecurved[0]])
      // 	of2 << xyz << std::endl;
      if (var > varh)
	pts_v.push_back(xyz);
      if (curved > curvh)
	pts_c.push_back(xyz);
    }
  std::ofstream of3("sfvar.g2");
  of3 << "400 1 0 4 155 0 100 255" << std::endl;
  of3 << pts_v.size() << std::endl;
  for (size_t kr=0; kr<pts_v.size(); ++kr)
    of3 << pts_v[kr] << std::endl;
  std::ofstream of5("curvedness.g2");
  of5 << "400 1 0 4 0 155 100 255" << std::endl;
  of5 << pts_c.size() << std::endl;
  for (size_t kr=0; kr<pts_c.size(); ++kr)
    of5 << pts_c[kr] << std::endl;

int stop_break = 1;
}

//===========================================================================
void RevEng::classifyPoints()
//===========================================================================
{
  // Fetch relevant values for all points

  // Update limits
  setClassificationParams();

  vector<vector<Vector3D> > class_pts(9);
  vector<vector<Vector3D> > class_shape(10);
  int nmbpts = tri_sf_->size();
  for (int ki=0; ki<nmbpts; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      Vector3D xyz = pt->getPoint();

      // Curvature surface classification
      int ctype = C1_UNDEF;
      double gausscurv = pt->GaussCurvature();
      double meancurv = pt->meanCurvature();
      if (meancurv < -zero_H_)
	{
	  if (gausscurv > zero_K_)
	    {
	      ctype = C1_PEAK;
	      class_pts[0].push_back(xyz);
	    }
	  else if (gausscurv < -zero_K_)
	    {
	      ctype = C1_SRIDGE;
	      class_pts[2].push_back(xyz);
	    }
	  else
	    {
	      ctype = C1_RIDGE;
	      class_pts[1].push_back(xyz);
	    }
	}
      else if (meancurv > zero_H_)
	{
	  if (gausscurv > zero_K_)
	    {
	      ctype = C1_PIT;
	      class_pts[6].push_back(xyz);
	    }
	  else if (gausscurv < -zero_K_)
	    {
	      ctype = C1_SVALLEY;
	      class_pts[8].push_back(xyz);
	    }
	  else
	    {
	      ctype = C1_VALLEY;
	      class_pts[7].push_back(xyz);
	    }
	}
      else
	{
	  if (gausscurv > zero_K_)
	    {
	      ctype = C1_NONE;
	      class_pts[3].push_back(xyz);
	    }
	  else if (gausscurv < -zero_K_)
	    {
	      ctype = C1_MINSURF;
	      class_pts[5].push_back(xyz);
	    }
	  else
	    {
	      ctype = C1_FLAT;
	      class_pts[4].push_back(xyz);
	    }
	}

      // Curvature edge classification
      double avlen = pt->getMeanEdgLen();
      double maxpc = std::max(fabs(pt->maxPrincipalCurvature()),
			      fabs(pt->minPrincipalCurvature()));
      double crvrad = 1.0/maxpc; 
      int c1_edge = (crvrad < cfac_*avlen) ? C1_EDGE : C1_NOT_EDGE;

      // Surface classification with shape index
      int si_type = SI_UNDEF;
      double shapeindex = pt->getShapeIndex();
      if (shapeindex < -0.875)
	{
	  si_type = SI_SCUP;
	  class_shape[0].push_back(xyz);
	}
      else  if (shapeindex < -0.625)
	{
	  si_type = SI_TRO;
	  class_shape[1].push_back(xyz);
	}
      else  if (shapeindex < -0.375)
	{
	  si_type = SI_RUT;
	  class_shape[2].push_back(xyz);
	}
       else  if (shapeindex < -0.125)
	 {
	   si_type = SI_SRUT;
	   class_shape[3].push_back(xyz);
	 }
       else  if (shapeindex < 0.125)
	 {
	   si_type = SI_SAD;
	   class_shape[4].push_back(xyz);
	 }
       else  if (shapeindex < 0.375)
	 {
	   si_type = SI_SRID;
	   class_shape[5].push_back(xyz);
	 }
      else  if (shapeindex < 0.625)
	{
	  si_type = SI_RID;
	  class_shape[6].push_back(xyz);
	}
      else  if (shapeindex < 0.875)
	{
	  si_type = SI_DOM;
	  class_shape[7].push_back(xyz);
	}
      else  if (shapeindex <= 1.0)
	{
	  si_type = SI_SCAP;
	  class_shape[8].push_back(xyz);
	}
      else
	{
	  si_type = SI_PLANE;
	  class_shape[9].push_back(xyz);
	}

      // Edge classification with curvedness
      double curved = pt->getCurvedness();
      int c2_edge = (curved > cness_lim_) ? C2_EDGE : C2_NOT_EDGE;

      // Edge classification with surface variation
      double var = pt->getSurfaceVariation();
      int pca_edge = (var > pca_lim_) ? PCA_EDGE : PCA_NOT_EDGE;

      // Store classification in point
      pt->setClassification(ctype, c1_edge, si_type, c2_edge, pca_edge);
   }

  std::ofstream of("curvature_segments.g2");
  for (int ka=0; ka<3; ++ka)
    for (int kb=0; kb<3; ++kb)
      {
	of << "400 1 0 4 ";
	for (int kc=0; kc<3; ++kc)
	  of << colors[3*ka+kb][kc] << " ";
	of << "255" << std::endl;
	of << class_pts[3*ka+kb].size() << std::endl;
	for (size_t kr=0; kr<class_pts[3*ka+kb].size(); ++kr)
	  of << class_pts[3*ka+kb][kr] << std::endl;
      }


  // Shape index
  std::ofstream of2("shapeindex_segments.g2");
  for (int ka=0; ka<10; ++ka)
    {
      of2 << "400 1 0 4 ";
      for (int kc=0; kc<3; ++kc)
	of2 << colors[ka][kc] << " ";
      of2 << "255" << std::endl;
      of2 << class_shape[ka].size() << std::endl;
      for (size_t kr=0; kr<class_shape[ka].size(); ++kr)
	of2 << class_shape[ka][kr] << std::endl;
      }
  
  // Specify/clean edge classification
  bool closeedge = false;
  int nmbedge = 2;
  for (int ka=0; ka<nmbpts; ++ka)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ka]);
      if (pt->isEdge())
	{
	  if (pt->isolatedEdge(nmbedge, closeedge))
	    pt->setEdgeUndef();

	  else 
	    // If the angular difference between triangle normals is less
	    // then the limit, classify the point as CLOSE_EDGE.
	    pt->adjustWithTriangNorm(norm_ang_lim_);
	}
    }

  vector<Vector3D> edgepts;
  for (int ka=0; ka<nmbpts; ++ka)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ka]);
      if (pt->isEdge())
	edgepts.push_back(pt->getPoint());
    }

  if (edgepts.size() > 0)
    {
      std::ofstream of3("edgepts.g2");
      of3 << "400 1 0 4 255 0 0 255" << std::endl;
      of3 << edgepts.size() << std::endl;
      for (size_t kr=0; kr<edgepts.size();++kr)
	of3 << edgepts[kr] << std::endl;
    }
  int stop_break = 1;
}

struct
{
  bool operator()(shared_ptr<RevEngRegion> a, shared_ptr<RevEngRegion> b)
  {
    return (a->numPoints() > b->numPoints());
    //   return 1;
    // else if (a->numPoints() == b->numPoints())
    //   return 2;
    // else
    //   return 3;
  }
} sort_region;

//===========================================================================
void RevEng::growRegions(int classification_type)
//===========================================================================
{
  // First pass: Collect continous regions
  // Prepare approximation tolerance for second pass
  int nmbpts = tri_sf_->size();
  vector<double> pointdist;
  for (int ki=0; ki<nmbpts; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      if (!pt->isEdge())
	pointdist.push_back(pt->getPointDistance());
      if (pt->hasRegion())
	continue;
      if (pt->closeEdge())
	continue;

      shared_ptr<RevEngRegion> region(new RevEngRegion(classification_type));
      regions_.push_back(region);
      pt->setRegion(region.get());
      region->collect(pt);
    }
  std::sort(pointdist.begin(), pointdist.end());
  double dlim = 0.93; //0.75;
  int dix = (int)(dlim*pointdist.size());
  approx_tol_ = pointdist[dix];   //  This should probably be set earlier
  // to give the user a chance to update

  std::sort(regions_.begin(), regions_.end(), sort_region);
  if (regions_.size() > 0)
    {
      std::ofstream of("regions1.g2");
      for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  // BoundingBox bbox = regions_[kr]->boundingBox();
	  // if (bbox.low().dist(bbox.high()) < 0.1)
	  //   std::cout << "Small bounding box" << std::endl;
	  if (regions_[kr]->numPoints() < 5)
	    continue;
	  of << "400 1 0 0" << std::endl;
	  int nmb = regions_[kr]->numPoints();
	  of << nmb << std::endl;
	  for (int ki=0; ki<nmb; ++ki)
	    {
	      of << regions_[kr]->getPoint(ki)->getPoint() << std::endl;
	    }
	}
    }

  // Second pass:  Verify/update regions using the distribution of
  // normals associated to the points on the unit sphere
  size_t regsize = regions_.size();
  for (size_t kr=0; kr<regsize; ++kr)
    {
      if (regions_[kr]->numPoints() < min_point_region_)
	break;

      int classtype = regions_[kr]->getClassification();
      vector<RevEngPoint*> deviant;
      vector<vector<RevEngPoint*> > connected;
      regions_[kr]->splitFromSurfaceNormals(deviant, connected);
      for (size_t kj=0; kj<connected.size(); ++kj)
	{
	  shared_ptr<RevEngRegion> reg(new RevEngRegion(classtype,
							connected[kj]));
	  regions_.push_back(reg);
	}
      
      if (deviant.size() > 0)
	{
	  shared_ptr<RevEngRegion> reg2(new RevEngRegion(classtype,
							 deviant));
	  vector<vector<RevEngPoint*> > deviant2;
	  reg2->splitRegion(deviant2);
	  regions_.push_back(reg2);
	  for (size_t kh=0; kh<deviant2.size(); ++kh)
	    {
	      shared_ptr<RevEngRegion> reg3(new RevEngRegion(classtype,
							     deviant2[kh]));
	      regions_.push_back(reg3);
	    }
	}
    }
  
  std::sort(regions_.begin(), regions_.end(), sort_region);
  
  if (regions_.size() > 0)
    {
      std::ofstream of2("regions2.g2");
      for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  // BoundingBox bbox = regions_[kr]->boundingBox();
	  // if (bbox.low().dist(bbox.high()) < 0.1)
	  //   std::cout << "Small bounding box" << std::endl;
	  if (regions_[kr]->numPoints() < 5)
	    continue;
	  of2 << "400 1 0 0" << std::endl;
	  int nmb = regions_[kr]->numPoints();
	  of2 << nmb << std::endl;
	  for (int ki=0; ki<nmb; ++ki)
	    {
	      of2 << regions_[kr]->getPoint(ki)->getPoint() << std::endl;
	    }
	}
    }

  bool thirdpass = false;
  if (thirdpass)
    {
      // Third pass: Verify/update regions with surface approximation
      std::ofstream of2("seed_points.g2");
      for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  if (regions_[kr]->numPoints() < min_point_region_)
	    break;

	  // Fetch seed point
	  RevEngPoint *seed = regions_[kr]->seedPoint();
	  if (seed)
	    {
	      of2 << "400 1 0 4 255 0 0 255" << std::endl;
	      of2 << "1" << std::endl;
	      of2 << seed->getPoint() << std::endl;
	    }

	  vector<RevEngPoint*> outpts;
	  double local_len = seed->getMeanEdgLen();
	  double radius = 2.0*rfac_*local_len;
	  int mnext = std::max(min_next_, regions_[kr]->numPoints()/100);
	  regions_[kr]->growLocal(seed, approx_tol_, radius, mnext, outpts);
	  int stop_break = 1;
	}
  
      // Until all points are checked

      // Fetch a classified point (not edge or near edge). Which classification
      // to select? Can they be combined somehow

      // While not hit edge, or no new points of the same type is found
      // Fetch all ring 1 points
      // Collect the points that are of the same type
      // Continue from all collected points, stopping in points that are assigned
      // already (recusivity)

      // Depends strongly on correct classification

      // Will leave unassigned points
    }
}


//===========================================================================
void RevEng::recognizeElementary()
//===========================================================================
{
  // For each surface type in increasing complexity
  // Would it be better to try to recognize cone before sphere?
  for (int sftype=0; sftype<5; ++sftype)
    {
      switch (sftype)
      {
	case PLANE:
	  // Recognition
	  // Merge almost similar planes if the accuracy is not suffering
	  recognizePlanes();
	  break;
	  case CYLINDER:
	    recognizeCylinders();
	    break;
	    default:
	      // Do nothing
	      break;
      }
    }
  
}

//===========================================================================
void RevEng::recognizePlanes()
//===========================================================================
{
  std::ofstream planeout("plane.g2");
  double anglim = 0.1;
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      // Check type
      bool planar = regions_[ki]->planartype();
      DirectionCone normalcone = regions_[ki]->getNormalCone();
      if (!planar && (normalcone.angle() > anglim ||
		      normalcone.greaterThanPi()))
	continue;  // Not a probable plane

      // Try to fit the point cloud with a plane
      vector<shared_ptr<HedgeSurface> > plane_sfs;
      bool found = regions_[ki]->extractPlane(approx_tol_, min_point_region_,
					      plane_sfs, planeout);
      if (plane_sfs.size() > 0)
	surfaces_.insert(surfaces_.end(), plane_sfs.begin(), plane_sfs.end());
      // if not planar
      // continue;

      // apply method for recognition of plane
      // Should RevEngPoint be given as input or should the method take a
      // vector of Points as input to enforce independence on a triangulation?
      // The class HedgeSurface must lie in compositemodel due to the
      // connection to RevEngRegion
      // Should the computation take place in HedgeSurface? Or in RevEngRegion?
      // Or in this class?
      
      // some points may be disassembled

      // check if the number of points in the plane is large enough
      // if not, should the region be removed?
      int stop_break = 1;
    }
}

//===========================================================================
void RevEng::recognizeCylinders()
//===========================================================================
{
  std::ofstream cylout("cylinder.g2");
  size_t nmbsfs = surfaces_.size();
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      // Check type
      bool cyl = regions_[ki]->cylindertype();
      // if (!cyl)
      // 	continue;
      
      // So far, try to regognize cylinders
      vector<shared_ptr<HedgeSurface> > cyl_sfs;
      bool found = regions_[ki]->extractCylinder(approx_tol_, min_point_region_,
						 cyl_sfs, cylout);
      if (cyl_sfs.size() > 0)
	surfaces_.insert(surfaces_.end(), cyl_sfs.begin(), cyl_sfs.end());
      
      // if not planar
      // continue;

      // apply method for recognition of plane
      // Should RevEngPoint be given as input or should the method take a
      // vector of Points as input to enforce independence on a triangulation?
      // The class HedgeSurface must lie in compositemodel due to the
      // connection to RevEngRegion
      // Should the computation take place in HedgeSurface? Or in RevEngRegion?
      // Or in this class?
      
      // some points may be disassembled

      // check if the number of points in the plane is large enough
      // if not, should the region be removed?
    }

  // Now the surfaces between nmbsfs and surfaces_.size() are cylinder surfaces
  // and can be merged if they represent the same cylinder
  if (surfaces_.size() + nmbsfs > 1)
    mergeCylinders(nmbsfs, surfaces_.size());

  std::ofstream of("surfaces.g2");
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
      surf->writeStandardHeader(of);
      surf->write(of);
    }
}

//===========================================================================
void RevEng::mergePlanes()
//===========================================================================
{
  // vector<shared_ptr<HedgeSurface> > planes;

  // Collect all planes
  // while every plane is considered
  // {
  // Find the largest plane;
  // Copy region to current;
  // for (all planes not previously considered)
  //   {
  //     Check similarity (deviation between point and plane, plane normal cone);
  //     if (similar)
  // 	add region to current; remember region;
  //   }
  // recognize plane based on current;

  // for (all regions in current)
  //   {
  //     for (all points)
  // 	{
  // 	  check accuracy;
  // 	  remember points with good enough accuracy (or number of points);
  // 	}
  //     if (accuracy not good enough)
  // 	remove region from current;
  //   }

  // recognize plane based on current;
  // clean up (merge regions in current if close, attach several regions to new?updated HedgeSurfaces if not;
  // }

  // Should unassembled planar points be assigned to appropriate plane if
  // close engough?

  // Close enough and good enough must be defined
  // Points may be disassembled

  // A lot of the operations here are the same for all surface types.
  // Maybe an inheritance in HedgeSurface instead of a swith would
  // facilitate specialization on type keeping the main algorithm common?
}

//===========================================================================
void RevEng::mergeCylinders(size_t first, size_t last)
//===========================================================================
{
  // Sort according to number of associated points
  for (size_t ki=first; ki<last; ++ki)
    {
      int nmb1 = surfaces_[ki]->numPoints();
      for (size_t kj=ki+1; kj<last; ++kj)
	{
	  int nmb2 = surfaces_[kj]->numPoints();
	  if (nmb2 > nmb1)
	    std::swap(surfaces_[ki], surfaces_[kj]);
	}
    }

  // Find similar cylinders
  for (size_t ki=first; ki<last; ++ki)
    {
      vector<size_t> cand_ix;
      int code;
      ClassType type1 = surfaces_[ki]->instanceType(code);
      if (type1 != Class_Cylinder && type1 != Class_BoundedSurface)
	continue;

      cand_ix.push_back(ki); 
      shared_ptr<ParamSurface> surf1 = surfaces_[ki]->surface();
      shared_ptr<Cylinder> cyl1 =
	dynamic_pointer_cast<Cylinder,ParamSurface>(surf1);
      if (!cyl1.get())
	{
	  shared_ptr<BoundedSurface> bdsf1 =
	    dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf1);
	  if (bdsf1.get())
	    {
	      surf1 = bdsf1->underlyingSurface();
	      cyl1 = dynamic_pointer_cast<Cylinder,ParamSurface>(surf1);
	    }
	}
      if (!cyl1.get())
	continue; 

      Point axis1 = cyl1->getAxis();
      Point pos1 = cyl1->getLocation();
      double rad1 = cyl1->getRadius();
      double dlim = std::max(0.01*rad1, approx_tol_);
      for (size_t kj=ki+1; kj<last; ++kj)
	{
	  ClassType type2 = surfaces_[kj]->instanceType(code);
	  if (type2 != Class_Cylinder && type2 != Class_BoundedSurface)
	    continue;
	  
	  shared_ptr<ParamSurface> surf2 = surfaces_[kj]->surface();
	  shared_ptr<Cylinder> cyl2 =
	    dynamic_pointer_cast<Cylinder,ParamSurface>(surf2);
	  if (!cyl2.get())
	    {
	      shared_ptr<BoundedSurface> bdsf2 =
		dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf2);
	      if (bdsf2.get())
		{
		  surf2 = bdsf2->underlyingSurface();
		  cyl2 = dynamic_pointer_cast<Cylinder,ParamSurface>(surf2);
		}
	    }
 	  if (!cyl2.get())
	    continue;  // Should not happen
	  Point axis2 = cyl2->getAxis();
	  Point pos2 = cyl2->getLocation();
	  double rad2 = cyl2->getRadius();

	  double ang = axis1.angle(axis2);
	  ang = std::min(ang, M_PI-ang);
	  Point pos2_0 = pos1 + ((pos2-pos1)*axis1)*axis1;
	  double pdist = pos2.dist(pos2_0);

	  if (ang > anglim_ || pdist > 5.0*dlim ||
	      fabs(rad2-rad1) > dlim)
	    continue;
	  cand_ix.push_back(kj);
	}

      if (cand_ix.size() > 1)
	{
	  std::ofstream pre("pre_merge.g2");
	  for (size_t kr=0; kr<cand_ix.size(); ++kr)
	    {
	      shared_ptr<ParamSurface> surf = surfaces_[cand_ix[kr]]->surface();
	      surf->writeStandardHeader(pre);
	      surf->write(pre);
	    }
	  shared_ptr<HedgeSurface> merged_surf = doMergeCylinders(cand_ix);
	  if (merged_surf.get())
	    {
	      std::ofstream post("post_merge.g2");
	      for (size_t kr=0; kr<cand_ix.size(); ++kr)
		{
		  shared_ptr<ParamSurface> surf = surfaces_[cand_ix[kr]]->surface();
		  surf->writeStandardHeader(post);
		  surf->write(post);
		}
	      shared_ptr<ParamSurface> surfm = merged_surf->surface();
	      surfm->writeStandardHeader(post);
	      surfm->write(post);
		  
	      std::swap(surfaces_[cand_ix[0]], merged_surf);
	      if (cand_ix[0] > ki)
		std::swap(surfaces_[ki], surfaces_[cand_ix[0]]);
	      for (size_t kr=cand_ix.size()-1; kr>=1; --kr)
		{
		  surfaces_.erase(surfaces_.begin()+cand_ix[kr]);
		}
	      last -= (cand_ix.size()-1);
	    }
	}
      int stop_break = 1;
    }
}

//===========================================================================
shared_ptr<HedgeSurface> RevEng::doMergeCylinders(vector<size_t>& cand_ix)
//===========================================================================
{
  size_t candsize = cand_ix.size();
  shared_ptr<HedgeSurface> merged_surf;
  for (size_t kh=0; kh<cand_ix.size(); ++kh)
    {
      // Collect regions and point clouds
      vector<RevEngRegion*> regions;
      vector<pair<vector<RevEngPoint*>::iterator,
		  vector<RevEngPoint*>::iterator> > points;
      BoundingBox bbox(3);
      for (size_t ki=0; ki<cand_ix.size(); ++ki)
	{
	  HedgeSurface* surf = surfaces_[cand_ix[ki]].get();
	  vector<RevEngRegion*> reg = surf->getRegions();
	  regions.insert(regions.end(), reg.begin(), reg.end());
	  for (size_t kj=0; kj<reg.size(); ++kj)
	    {
	      points.push_back(std::make_pair(reg[kj]->pointsBegin(),
					      reg[kj]->pointsEnd()));
	      bbox.addUnionWith(reg[kj]->boundingBox());
	    }
	}
  

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
      shared_ptr<Cylinder> cyl(new Cylinder(rad, pnt, axis, Cy));
      double len = low.dist(high);
      cyl->setParamBoundsV(-len, len);

      // Check accuracy
      double dfac = 5.0;
      for (size_t ki=0; ki<regions.size(); ++ki)
	{
	  double maxd, avd;
	  int num2;
	  RevEngUtils::distToSurf(points[ki].first, points[ki].second,
				  cyl, approx_tol_, maxd, avd, num2);

	  double maxd_init, avd_init;
	  int num2_init;
	  regions[ki]->getAccuracy(maxd_init, avd_init, num2_init);
	  int num = regions[ki]->numPoints();

	  if (num2 < num/2 || maxd > dfac*maxd_init || avd > dfac*avd_init
	      || avd > approx_tol_)
	    {
	      size_t kj, kr;
	      for (kj=0; kj<cand_ix.size(); ++kj)
		{
		  HedgeSurface* surf = surfaces_[cand_ix[kj]].get();
		  vector<RevEngRegion*> reg = surf->getRegions();
		  for (kr=0; kr<reg.size(); ++kr)
		    if (reg[kr] == regions[ki])
		      break;
		  if (kr < reg.size())
		    break;
		}
	      if (kj < cand_ix.size())
		cand_ix.erase(cand_ix.begin() + kj);
	    }
	  int stop_break = 1;
	}
      if (cand_ix.size() <= 1)
	break;
      if (cand_ix.size() == candsize)
	{
	  merged_surf =
	    shared_ptr<HedgeSurface>(new HedgeSurface(cyl, regions));
	  for (size_t kj=0; kj<regions.size(); ++kj)
	    regions[kj]->setHedge(merged_surf.get());
	  break;
	}
      candsize = cand_ix.size();
    }

  return merged_surf;
}

//===========================================================================
void RevEng::trimPrimitives()
//===========================================================================
{
  // for (size_t ki=0; ki<surfaces_.size(); ++ki)
  //   {
  //     if (!surfaces_[ki]->hasSurface())
  // 	continue;
  //     shared_ptr<ParamSurface> sf1 = surfaces_[ki]->getSurface();
  //     for (size_t kj=0; kj<surfaces_.size(); ++kj)
  // 	{
  // 	  if (!surfaces_[kj]->hasSurface())
  // 	    continue;
  // 	  shared_ptr<ParamSurface> sf2 = surfaces_[kj]->getSurface();

  // 	  // Intersect sf1 and sf2
  // 	  // Remember intersection curves
  // 	}
  //   }
}

 //===========================================================================
void RevEng::initParameters()
//===========================================================================
{
  // Set default parameters
  min_next_ = 10;  // Minimum number of neighbouring points
  rfac_ = 3.0;  // Factor for radius in which to search for neighbouring points
  cfac_ = 8.0;  // Edge points from curvature is given by
  // cfac_ times the average length of triangulation edges in a vertex
  norm_plane_lim_= 0.005; // Limit for when the cone angle corresponding
  // to triangle normals indicate an edge
  zero_H_ = 0.0001;  // When mean curvature is considered zero
  zero_K_ = 0.0001;  // When Gauss curvature is considered zero
  zero_si_ = 0.001; // When shape index is considered zero
  norm_ang_lim_ = 0.025*M_PI; // Limit for when the cone angle corresponding
    // to triangle normals indicate an edge
  pca_lim_ = cness_lim_ = std::numeric_limits<double>::max();
  min_point_region_ = 100; //10;  // Should be updated with regard to the total
  // number of points
  approx_tol_ = 0.001;  // Very preliminary
  anglim_ = 0.01;
}

 //===========================================================================
void RevEng::storeClassified(ostream& os) const
//===========================================================================
{
  storeParams(os);
  int nmbpts = tri_sf_->size();
  os << nmbpts << std::endl;
  for (int ki=0; ki<nmbpts; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      pt->store(os);
    }
}

 //===========================================================================
void RevEng::readClassified(istream& is)
//===========================================================================
{
  readParams(is);
  int nmbpts;
  is >> nmbpts;
  tri_sf_ = shared_ptr<ftPointSet>(new ftPointSet());
  vector<vector<int> > next_ix(nmbpts);
  for (int ki=0; ki<nmbpts; ++ki)
    {
      shared_ptr<RevEngPoint> vertex(new RevEngPoint());
      vertex->read(is, zero_si_, next_ix[ki]);
      tri_sf_->addEntry(vertex);
    }

  // Add next information
  for (int ki=0; ki<nmbpts; ++ki)
    {
      ftSamplePoint* pt1 = (*tri_sf_)[ki];
      for (size_t kr=0; kr<next_ix[ki].size(); ++kr)
	{
	  int ix = next_ix[ki][kr];
	  ftSamplePoint* pt2 = (*tri_sf_)[ix];
	  pt1->addNeighbour(pt2);
	}
    }
  
}

 //===========================================================================
void RevEng::storeParams(ostream& os) const
//===========================================================================
{
  os << mean_edge_len_ << " " << min_next_ << " " << rfac_ << " " << cfac_;
  os << " " << pca_lim_ << " " << cness_lim_ << " " << norm_ang_lim_;
  os << " " << norm_plane_lim_ << " " << zero_H_ << " " << zero_K_;
  os << " " << zero_si_ << std::endl;
}

 //===========================================================================
void RevEng::readParams(istream& is)
//===========================================================================
{
  is >> mean_edge_len_ >> min_next_ >> rfac_ >> cfac_ >> pca_lim_ >> cness_lim_;
  is >> norm_ang_lim_ >> norm_plane_lim_ >> zero_H_ >> zero_K_ >> zero_si_;
}
