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
#include "GoTools/compositemodel/ImplicitApprox.h"
#include "GoTools/utils/DirectionCone.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/Plane.h"
#include "GoTools/geometry/Torus.h"
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
  std::ofstream of01("minc1.g2");
  std::ofstream of02("minc2.g2");
  std::ofstream of03("maxc1.g2");
  std::ofstream of04("maxc2.g2");
  of01 << "410 1 0 4 200 50 0 255" << std::endl;
  of01 << nmbpt << std::endl;
  of02 << "410 1 0 4 200 50 0 255" << std::endl;
  of02 << nmbpt << std::endl;
  of03 << "410 1 0 4 0 50 200 255" << std::endl;
  of03 << nmbpt << std::endl;
  of04 << "410 1 0 4 0 50 200 255" << std::endl;
  of04 << nmbpt << std::endl;
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]); 

      // Compute surface normal from triangulation
      pt->computeTriangNormal(100.0*mean_edge_len_);
      if (pt->isOutlier())
	continue;

      if (pt->getNmbNeighbour() == 0)
	{
	  pt->setOutlier();
	  continue;
	}

      double avlen = pt->getMeanEdgLen();

      //Fetch nearby points
      vector<Point> nearpts;
      double local_len = pt->getMeanEdgLen(10.0*mean_edge_len_);
      double radius = rfac_*(local_len + mean_edge_len_);
      radius = std::min(radius, 20.0*mean_edge_len_);
      //radius *= 1.5; // TEST 
      //double radius = 0.5*rfac_*(local_len + mean_edge_len_);
      Point curr = pt->fetchClosePoints(radius, min_next_, max_next_, nearpts);

      Point mincvec(0.0, 0.0, 0.0), maxcvec(0.0, 0.0, 0.0);
      Point mincvec2(0.0, 0.0, 0.0), maxcvec2(0.0, 0.0, 0.0);

      if (nearpts.size() >= 3)
	{

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
      Point tnorm = pt->getTriangNormal();
      if (tnorm.length() < 1.0e-10)
	{
	  int stop_norm = 1;
	}
      else if (eigen3*tnorm <  0.0)
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
      Point normal;//, mincvec, maxcvec;
      double minc, maxc;
      double currdist, avdist;
      // RevEngUtils::computeMonge(curr, nearpts, eigen1, eigen2, normal, mincvec, minc,
      // 				maxcvec, maxc, currdist, avdist);
      RevEngUtils::computeMonge(curr, nearpts, eigen1, eigen3, normal, mincvec, minc,
				maxcvec, maxc, currdist, avdist);
      // Orient vectors with respect to triangulation normal
      // The normal vectors should be OK. Curvature vectors are not necessarily
      // consistent with regard to orientation

      double minc2, maxc2;
      RevEngUtils::TaubinCurvature(curr, nearpts, eigen1, eigen3, mincvec2, minc2,
				   maxcvec2, maxc2);
      
      if (normal*tnorm < 0.0)
	normal *= -1;
      pt->addMongeInfo(normal, mincvec, minc, maxcvec, maxc, currdist, avdist,
		       zero_si_);
	}
      of01 << curr << " " << curr+mincvec << std::endl;
      of02 << curr << " " << curr+mincvec2 << std::endl;
      of03 << curr << " " << curr+maxcvec << std::endl;
      of04 << curr << " " << curr+maxcvec2 << std::endl;

      int stop_break = 1;
    }

  for (int ki=0; ki<nmbpt; ++ki)
    {
      // Check orientation of curvature
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]); 
      if (pt->isOutlier())
	continue;
      Point tnorm = pt->getTriangNormal();
      if (tnorm.length() < 1.0e-10)
	{
	  // Fetch triangulation normal in neighbour
	  vector<ftSamplePoint*> next = pt->getNeighbours();
	  double mindist = std::numeric_limits<double>::max();
	  int ix = -1;
	  for (size_t kr=0; kr<next.size(); ++kr)
	    {
	      double dist = pt->pntDist(next[kr]);
	      RevEngPoint *nextpt = dynamic_cast<RevEngPoint*>(next[kr]);
	      Point nextnorm = nextpt->getTriangNormal();
	      if (dist < mindist && nextnorm.length() > 1.0e-10)
		{
		  mindist = dist;
		  ix = (int)kr;
		}
	    }
	  if (ix >= 0)
	    {
	      RevEngPoint *nextpt = dynamic_cast<RevEngPoint*>(next[ix]);
	      tnorm = nextpt->getTriangNormal();
	      Point PCAnorm = pt->getPCANormal();
	      if (tnorm*PCAnorm < 0.0)
		pt->turnPCA();
	      Point Mongenorm = pt->getMongeNormal();
	      if (tnorm*Mongenorm < 0.0)
		pt->turnMongeNorm();
	    }
	  
	}
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
      if (pt->isOutlier())
	continue;
      Vector3D xyz = pt->getPoint();
      Point xyz2(xyz[0], xyz[1], xyz[2]);
      Point norm = pt->getTriangNormal();
      double ang = pt->getTriangAngle();
      double avlen = pt->getMeanEdgLen();
      of << xyz2 << " " << xyz2+avlen*norm << std::endl;
      if (pt->getTriangAngle() > norm_ang_lim_)
	triangcorners.push_back(xyz);
      Point Mnorm = pt->getMongeNormal();
      ofM << xyz2 << " " << xyz2+avlen*Mnorm << std::endl;

      if (ang <= norm_plane_lim_)
	triangplane.push_back(xyz);
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
      radius = std::min(radius, 20.0*mean_edge_len_);
      pt->fetchClosePoints2(radius, min_next_, max_next_, nearpts[ki]);
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

  int nmbsmooth = 5;
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
  std::cout << "Approx tol: " << approx_tol_ << std::endl;
  std::cout << "New tolerance: " << std::endl;
  std::cin >> approx_tol_;

  std::sort(regions_.begin(), regions_.end(), sort_region);
  if (regions_.size() > 0)
    {
      std::ofstream of("regions1.g2");
      vector<Vector3D> small;
      for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  // BoundingBox bbox = regions_[kr]->boundingBox();
	  // if (bbox.low().dist(bbox.high()) < 0.1)
	  //   std::cout << "Small bounding box" << std::endl;
	  int nmb = regions_[kr]->numPoints();
	  if (nmb < 50)
	    {
	      for (int ki=0; ki<nmb; ++ki)
		small.push_back(regions_[kr]->getPoint(ki)->getPoint());
	    }
	  else
	    {
	      of << "400 1 0 0" << std::endl;
	      of << nmb << std::endl;
	      for (int ki=0; ki<nmb; ++ki)
		{
		  of << regions_[kr]->getPoint(ki)->getPoint() << std::endl;
		}
	    }
	}
      std::ofstream ofs("small_regions.g2");
      ofs << "400 1 0 4 0 0 0 255" << std::endl;
      ofs << small.size() << std::endl;
      for (size_t kr=0; kr<small.size(); ++kr)
	ofs << small[kr] << std::endl;
    }

  bool applysecond = true; //false;
  if (applysecond)
    {
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
	  for (size_t kh=0; kh<deviant.size(); ++kh)
	    deviant[kh]->setGaussRad(1.0);
	  shared_ptr<RevEngRegion> reg2(new RevEngRegion(classtype,
							 deviant));
	  vector<RevEngPoint*> deviant2;
	  vector<vector<RevEngPoint*> > connected2;
	  reg2->splitFromSurfaceNormals(deviant2, connected2);
	  regions_.push_back(reg2);
	  for (size_t kj=0; kj<connected2.size(); ++kj)
	    {
	      shared_ptr<RevEngRegion> reg3(new RevEngRegion(classtype,
							    connected2[kj]));
	      regions_.push_back(reg3);
	    }
	  if (deviant2.size() > 0)
	    {
	      shared_ptr<RevEngRegion> reg4(new RevEngRegion(classtype,
							     deviant2));
	      vector<vector<RevEngPoint*> > deviant4;
	      reg4->splitRegion(deviant4);
	      regions_.push_back(reg4);
	      for (size_t kh=0; kh<deviant4.size(); ++kh)
		{
		  shared_ptr<RevEngRegion> reg5(new RevEngRegion(classtype,
								 deviant4[kh]));
		  regions_.push_back(reg5);
		}
	    }
	}
    }
  
  std::sort(regions_.begin(), regions_.end(), sort_region);
  
  if (regions_.size() > 0)
    {
      std::ofstream of2("regions2.g2");
       vector<Vector3D> small;
     for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  // BoundingBox bbox = regions_[kr]->boundingBox();
	  // if (bbox.low().dist(bbox.high()) < 0.1)
	  //   std::cout << "Small bounding box" << std::endl;
	  int nmb = regions_[kr]->numPoints();
	  if (nmb < 50)
	    {
	      for (int ki=0; ki<nmb; ++ki)
		small.push_back(regions_[kr]->getPoint(ki)->getPoint());
	    }
	  else
	    {
	      of2 << "400 1 0 0" << std::endl;
	      int nmb = regions_[kr]->numPoints();
	      of2 << nmb << std::endl;
	      for (int ki=0; ki<nmb; ++ki)
		{
		  of2 << regions_[kr]->getPoint(ki)->getPoint() << std::endl;
		}
	    }
	}
      std::ofstream ofs("small_regions2.g2");
      ofs << "400 1 0 4 0 0 0 255" << std::endl;
      ofs << small.size() << std::endl;
      for (size_t kr=0; kr<small.size(); ++kr)
	ofs << small[kr] << std::endl;
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

  // Sort regions according to number of points
  std::sort(regions_.begin(), regions_.end(), sort_region);

  // Identify adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

  // Simplify regions structure
  std::cout << "Number of regions pre integrate: " << regions_.size() << std::endl;
  for (int ka=(int)(regions_.size()-1); ka>=0; --ka)
    {
      bool to_be_removed = regions_[ka]->integrateInAdjacent(mean_edge_len_,
							     min_next_, max_next_,
							     approx_tol_, 0.5,
							     max_nmb_outlier_);
      if (to_be_removed)
	{
	  regions_.erase(regions_.begin()+ka);
	}
    }
  std::cout << "Number of regions post integrate: " << regions_.size() << std::endl;
  
  std::sort(regions_.begin(), regions_.end(), sort_region);
  
  if (regions_.size() > 0)
    {
      std::ofstream of4("regions2_1.g2");
       vector<Vector3D> small;
     for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  // BoundingBox bbox = regions_[kr]->boundingBox();
	  // if (bbox.low().dist(bbox.high()) < 0.1)
	  //   std::cout << "Small bounding box" << std::endl;
	  int nmb = regions_[kr]->numPoints();
	  if (nmb < 50)
	    {
	      for (int ki=0; ki<nmb; ++ki)
		small.push_back(regions_[kr]->getPoint(ki)->getPoint());
	    }
	  else
	    {
	      of4 << "400 1 0 0" << std::endl;
	      int nmb = regions_[kr]->numPoints();
	      of4 << nmb << std::endl;
	      for (int ki=0; ki<nmb; ++ki)
		{
		  of4 << regions_[kr]->getPoint(ki)->getPoint() << std::endl;
		}
	    }
	}
      std::ofstream ofs4("small_regions2_1.g2");
      ofs4 << "400 1 0 4 0 0 0 255" << std::endl;
      ofs4 << small.size() << std::endl;
      for (size_t kr=0; kr<small.size(); ++kr)
	ofs4 << small[kr] << std::endl;
     }

  std::cout << "Pre join. Number of regions: " << regions_.size() << std::endl;
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      vector<RevEngRegion*> adapted_regions;
      regions_[ki]->joinRegions(approx_tol_, 10.0*anglim_, adapted_regions);
      for (size_t kj=0; kj<adapted_regions.size(); ++kj)
	{
	  size_t kr=0;
	  for (kr=0; kr<regions_.size(); ++kr)
	    if (adapted_regions[kj] == regions_[kr].get())
	      break;
	  if (kr < regions_.size())
	    regions_.erase(regions_.begin()+kr);
	}
    }
  std::cout << "Post join. Number of regions: " << regions_.size() << std::endl;
  
  if (regions_.size() > 0)
    {
      std::ofstream of3_2("regions2_2.g2");
       vector<Vector3D> small;
     for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  // BoundingBox bbox = regions_[kr]->boundingBox();
	  // if (bbox.low().dist(bbox.high()) < 0.1)
	  //   std::cout << "Small bounding box" << std::endl;
	  int nmb = regions_[kr]->numPoints();
	  if (nmb < 50)
	    {
	      for (int ki=0; ki<nmb; ++ki)
		small.push_back(regions_[kr]->getPoint(ki)->getPoint());
	    }
	  else
	    {
	      of3_2 << "400 1 0 0" << std::endl;
	      int nmb = regions_[kr]->numPoints();
	      of3_2 << nmb << std::endl;
	      for (int ki=0; ki<nmb; ++ki)
		{
		  of3_2 << regions_[kr]->getPoint(ki)->getPoint() << std::endl;
		}
	    }
	}
      std::ofstream ofs3_2("small_regions2_2.g2");
      ofs3_2 << "400 1 0 4 0 0 0 255" << std::endl;
      ofs3_2 << small.size() << std::endl;
      for (size_t kr=0; kr<small.size(); ++kr)
	ofs3_2 << small[kr] << std::endl;
     }

  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }


bool updatereg = true; //false;
  if (updatereg)
    {
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      vector<RevEngRegion*> adapted_regions;
      vector<shared_ptr<RevEngRegion> > outdiv_regions;
      regions_[ki]->updateRegion(approx_tol_, 10.0*anglim_, adapted_regions, outdiv_regions);
      for (size_t kj=0; kj<adapted_regions.size(); ++kj)
	{
	  size_t kr=0;
	  for (kr=0; kr<regions_.size(); ++kr)
	    if (adapted_regions[kj] == regions_[kr].get())
	      break;
	  if (kr < regions_.size())
	    regions_.erase(regions_.begin()+kr);
	}
      if (outdiv_regions.size() > 0)
	regions_.insert(regions_.end(), outdiv_regions.begin(), outdiv_regions.end());
    }
	      
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

  if (regions_.size() > 0)
    {
      std::ofstream of3("regions2_3.g2");
       vector<Vector3D> small;
     for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  // BoundingBox bbox = regions_[kr]->boundingBox();
	  // if (bbox.low().dist(bbox.high()) < 0.1)
	  //   std::cout << "Small bounding box" << std::endl;
	  int nmb = regions_[kr]->numPoints();
	  if (nmb < 50)
	    {
	      for (int ki=0; ki<nmb; ++ki)
		small.push_back(regions_[kr]->getPoint(ki)->getPoint());
	    }
	  else
	    {
	      of3 << "400 1 0 0" << std::endl;
	      int nmb = regions_[kr]->numPoints();
	      of3 << nmb << std::endl;
	      for (int ki=0; ki<nmb; ++ki)
		{
		  of3 << regions_[kr]->getPoint(ki)->getPoint() << std::endl;
		}
	    }
	}
      std::ofstream ofs2("small_regions2_3.g2");
      ofs2 << "400 1 0 4 0 0 0 255" << std::endl;
      ofs2 << small.size() << std::endl;
      for (size_t kr=0; kr<small.size(); ++kr)
	ofs2 << small[kr] << std::endl;
     }
    }
  // Simplify regions structure
  std::cout << "Number of regions pre integrate: " << regions_.size() << std::endl;
  for (int ka=(int)(regions_.size()-1); ka>=0; --ka)
    {
      bool to_be_removed = regions_[ka]->integrateInAdjacent(mean_edge_len_,
							     min_next_, max_next_,
							     approx_tol_, 0.5,
							     max_nmb_outlier_);
      if (to_be_removed)
	{
	  regions_.erase(regions_.begin()+ka);
	}
    }
  std::cout << "Number of regions post integrate: " << regions_.size() << std::endl;
  
  std::sort(regions_.begin(), regions_.end(), sort_region);
  
  int stop_break2 = 1;
}


//===========================================================================
void RevEng::recognizeElementary()
//===========================================================================
{
  std::ofstream planeout("plane.g2");
  std::ofstream cylout("cylinder.g2");
  std::ofstream coneout("cone.g2");
  std::ofstream torout("torus.g2");
  std::ofstream splout("spline.g2");
  double frac = 0.75;   // Fraction of points with a certain property
  double angfac = 10.0;
  
  min_point_region_ = 200;  // Testing
  int min_point_in = 20;
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      if (regions_[ki]->numPoints() < min_point_region_)
	continue;
      
      std::ofstream of1("region.g2");
      regions_[ki]->writeRegionInfo(of1);
      std::ofstream of2("unitsphere.g2");
      regions_[ki]->writeUnitSphereInfo(of2);

      bool found1 = false, found2 = false, found3 = false, found4 = false, found5 = false;
      if (true) //regions_[ki]->possiblePlane(2.0*angfac*anglim_, frac))
      {
	vector<shared_ptr<HedgeSurface> > plane_sfs;
	vector<HedgeSurface*> prev_surfs;
	found1 = regions_[ki]->extractPlane(approx_tol_, min_point_in,
					    plane_sfs, prev_surfs, planeout);
      
	for (size_t kr=0; kr<prev_surfs.size(); ++kr)
	  {
	    size_t kj;
	    for (kj=0; kj<surfaces_.size(); ++kj)
	      if (surfaces_[kj].get() == prev_surfs[kr])
		break;
	    if (kj < surfaces_.size())
	      surfaces_.erase(surfaces_.begin()+kj);
	  }
	if (plane_sfs.size() > 0)
	  surfaces_.insert(surfaces_.end(), plane_sfs.begin(), plane_sfs.end());
      }

      
      if (true) //regions_[ki]->possibleCylinder(angfac*anglim_, frac))
	 {
	   vector<shared_ptr<HedgeSurface> > cyl_sfs;
	   vector<HedgeSurface*> prev_surfs;
	   found2 = regions_[ki]->extractCylinder(approx_tol_, min_point_in,
						  mean_edge_len_, cyl_sfs, prev_surfs,
						  cylout);

      
	   for (size_t kr=0; kr<prev_surfs.size(); ++kr)
	     {
	       size_t kj;
	       for (kj=0; kj<surfaces_.size(); ++kj)
		 if (surfaces_[kj].get() == prev_surfs[kr])
		   break;
	       if (kj < surfaces_.size())
		 surfaces_.erase(surfaces_.begin()+kj);
	     }
	   if (cyl_sfs.size() > 0)
	     surfaces_.insert(surfaces_.end(), cyl_sfs.begin(), cyl_sfs.end());
	 }

       
      if (regions_[ki]->possibleCone(angfac*anglim_, frac))
	 {
	   vector<shared_ptr<HedgeSurface> > cone_sfs;
	   vector<HedgeSurface*> prev_surfs;
	   found3 = regions_[ki]->extractCone(approx_tol_, min_point_in,
					      mean_edge_len_, cone_sfs, prev_surfs,
					      coneout);

      
	   for (size_t kr=0; kr<prev_surfs.size(); ++kr)
	     {
	       size_t kj;
	       for (kj=0; kj<surfaces_.size(); ++kj)
		 if (surfaces_[kj].get() == prev_surfs[kr])
		   break;
	       if (kj < surfaces_.size())
		 surfaces_.erase(surfaces_.begin()+kj);
	     }
	   if (cone_sfs.size() > 0)
	     surfaces_.insert(surfaces_.end(), cone_sfs.begin(), cone_sfs.end());
	 }

       
       if (regions_[ki]->possibleTorus(angfac*anglim_, frac))  // Not the final parameters
	{
	  vector<shared_ptr<HedgeSurface> > tor_sfs;
	  vector<HedgeSurface*> prev_surfs;
	  found4 = regions_[ki]->extractTorus(approx_tol_, min_point_in,
					      mean_edge_len_, tor_sfs, prev_surfs,
					      torout);

	  for (size_t kr=0; kr<prev_surfs.size(); ++kr)
	     {
	       size_t kj;
	       for (kj=0; kj<surfaces_.size(); ++kj)
		 if (surfaces_[kj].get() == prev_surfs[kr])
		   break;
	       if (kj < surfaces_.size())
		 surfaces_.erase(surfaces_.begin()+kj);
	     }
	   if (tor_sfs.size() > 0)
	     surfaces_.insert(surfaces_.end(), tor_sfs.begin(), tor_sfs.end());
	 }

       if (true)
	 {
	  vector<shared_ptr<HedgeSurface> > spl_sfs;
	  vector<HedgeSurface*> prev_surfs;
	  found5 = regions_[ki]->extractFreeform(approx_tol_, min_point_in,
						 mean_edge_len_, spl_sfs, prev_surfs,
						 splout);

	  for (size_t kr=0; kr<prev_surfs.size(); ++kr)
	     {
	       size_t kj;
	       for (kj=0; kj<surfaces_.size(); ++kj)
		 if (surfaces_[kj].get() == prev_surfs[kr])
		   break;
	       if (kj < surfaces_.size())
		 surfaces_.erase(surfaces_.begin()+kj);
	     }
	   if (spl_sfs.size() > 0)
	     surfaces_.insert(surfaces_.end(), spl_sfs.begin(), spl_sfs.end());
	 }
       
      if (regions_[ki]->hasDivideInfo())
	{
	}
      int stop_break = 1;
    }

  if (regions_.size() > 0)
    {
      std::ofstream of("regions3.g2");
      vector<Vector3D> small;
      for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  // BoundingBox bbox = regions_[kr]->boundingBox();
	  // if (bbox.low().dist(bbox.high()) < 0.1)
	  //   std::cout << "Small bounding box" << std::endl;
	  int nmb = regions_[kr]->numPoints();
	  if (nmb < 50)
	    {
	      for (int ki=0; ki<nmb; ++ki)
		small.push_back(regions_[kr]->getPoint(ki)->getPoint());
	    }
	  else
	    {
	      of << "400 1 0 0" << std::endl;
	      of << nmb << std::endl;
	      for (int ki=0; ki<nmb; ++ki)
		{
		  of << regions_[kr]->getPoint(ki)->getPoint() << std::endl;
		}
	      if (regions_[kr]->hasSurface())
		{
		  HedgeSurface *sf = regions_[kr]->getSurface(0);
		  shared_ptr<ParamSurface> sf2 = sf->surface();
		  sf2->writeStandardHeader(of);
		  sf2->write(of);
		}
	    }
	}
      std::ofstream ofs("small_regions3.g2");
      ofs << "400 1 0 4 0 0 0 255" << std::endl;
      ofs << small.size() << std::endl;
      for (size_t kr=0; kr<small.size(); ++kr)
	ofs << small[kr] << std::endl;
    }

  bool doGrow = true; //false;
  if (doGrow)
    {
      std::cout << "Number of regions, pre grow with surf: " << regions_.size() << std::endl;
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      if (regions_[ki]->hasSurface())
	growSurface(ki);
    }
      std::cout << "Number of regions, post grow with surf: " << regions_.size() << std::endl;

  if (regions_.size() > 0)
    {
      std::ofstream of2("regions4.g2");
      vector<Vector3D> small;
      for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  // BoundingBox bbox = regions_[kr]->boundingBox();
	  // if (bbox.low().dist(bbox.high()) < 0.1)
	  //   std::cout << "Small bounding box" << std::endl;
	  int nmb = regions_[kr]->numPoints();
	  if (nmb < 50)
	    {
	      for (int ki=0; ki<nmb; ++ki)
		small.push_back(regions_[kr]->getPoint(ki)->getPoint());
	    }
	  else
	    {
	      of2 << "400 1 0 0" << std::endl;
	      of2 << nmb << std::endl;
	      for (int ki=0; ki<nmb; ++ki)
		{
		  of2 << regions_[kr]->getPoint(ki)->getPoint() << std::endl;
		}
	      if (regions_[kr]->hasSurface())
		{
		  HedgeSurface *sf = regions_[kr]->getSurface(0);
		  shared_ptr<ParamSurface> sf2 = sf->surface();
		  sf2->writeStandardHeader(of2);
		  sf2->write(of2);
		}
	    }
	}
      std::ofstream ofs2("small_regions4.g2");
      ofs2 << "400 1 0 4 0 0 0 255" << std::endl;
      ofs2 << small.size() << std::endl;
      for (size_t kr=0; kr<small.size(); ++kr)
	ofs2 << small[kr] << std::endl;
    }
  
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      if (regions_[ki]->hasSurface())
	regions_[ki]->adjustWithSurf(approx_tol_, 10.0*anglim_);
    }
      std::cout << "Number of regions, post grow with surf: " << regions_.size() << std::endl;

  if (regions_.size() > 0)
    {
      std::ofstream of5("regions5.g2");
      vector<Vector3D> small;
      for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  // BoundingBox bbox = regions_[kr]->boundingBox();
	  // if (bbox.low().dist(bbox.high()) < 0.1)
	  //   std::cout << "Small bounding box" << std::endl;
	  int nmb = regions_[kr]->numPoints();
	  if (nmb < 50)
	    {
	      for (int ki=0; ki<nmb; ++ki)
		small.push_back(regions_[kr]->getPoint(ki)->getPoint());
	    }
	  else
	    {
	      of5 << "400 1 0 0" << std::endl;
	      of5 << nmb << std::endl;
	      for (int ki=0; ki<nmb; ++ki)
		{
		  of5 << regions_[kr]->getPoint(ki)->getPoint() << std::endl;
		}
	      if (regions_[kr]->hasSurface())
		{
		  HedgeSurface *sf = regions_[kr]->getSurface(0);
		  shared_ptr<ParamSurface> sf2 = sf->surface();
		  sf2->writeStandardHeader(of5);
		  sf2->write(of5);
		}
	    }
	}
      std::ofstream ofs5("small_regions5.g2");
      ofs5 << "400 1 0 4 0 0 0 255" << std::endl;
      ofs5 << small.size() << std::endl;
      for (size_t kr=0; kr<small.size(); ++kr)
	ofs5 << small[kr] << std::endl;
    }
      std::cout << "Number of regions, pre grow with surf: " << regions_.size() << std::endl;
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      if (regions_[ki]->hasSurface())
	growSurface(ki);
    }
      std::cout << "Number of regions, post grow with surf: " << regions_.size() << std::endl;

  if (regions_.size() > 0)
    {
      std::ofstream of2("regions6.g2");
      vector<Vector3D> small;
      for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  // BoundingBox bbox = regions_[kr]->boundingBox();
	  // if (bbox.low().dist(bbox.high()) < 0.1)
	  //   std::cout << "Small bounding box" << std::endl;
	  int nmb = regions_[kr]->numPoints();
	  if (nmb < 50)
	    {
	      for (int ki=0; ki<nmb; ++ki)
		small.push_back(regions_[kr]->getPoint(ki)->getPoint());
	    }
	  else
	    {
	      of2 << "400 1 0 0" << std::endl;
	      of2 << nmb << std::endl;
	      for (int ki=0; ki<nmb; ++ki)
		{
		  of2 << regions_[kr]->getPoint(ki)->getPoint() << std::endl;
		}
	      if (regions_[kr]->hasSurface())
		{
		  HedgeSurface *sf = regions_[kr]->getSurface(0);
		  shared_ptr<ParamSurface> sf2 = sf->surface();
		  sf2->writeStandardHeader(of2);
		  sf2->write(of2);
		}
	    }
	}
      std::ofstream ofs2("small_regions6.g2");
      ofs2 << "400 1 0 4 0 0 0 255" << std::endl;
      ofs2 << small.size() << std::endl;
      for (size_t kr=0; kr<small.size(); ++kr)
	ofs2 << small[kr] << std::endl;
    }
    }

  // Merge surface pieces representing the surface
  std::cout << "Pre merge: " << surfaces_.size() << " surfaces" << std::endl;
  mergeSurfaces();
  std::cout << "Post merge: " << surfaces_.size() << " surfaces" << std::endl;

  adaptToMainAxis();
  int stop_break = 1;
}

//===========================================================================
void RevEng::growSurface(size_t& ix)
//===========================================================================
{
  vector<RevEngRegion*> grown_regions;
  int min_nmb = 5*min_point_region_;  // Should be set from distribution of how many
	  // points the regions have
  vector<HedgeSurface*> adj_surfs;
  regions_[ix]->growWithSurf(min_nmb, approx_tol_, grown_regions, adj_surfs);
  if (grown_regions.size() > 0)
    {
      for (size_t kr=0; kr<grown_regions.size(); ++kr)
	{
	  size_t kj;
	  for (kj=0; kj<regions_.size(); )
	    {
	      if (kj == ix)
		{
		  ++kj;
		  continue;
		}

	      if (grown_regions[kr] == regions_[kj].get())
		{
		  regions_.erase(regions_.begin()+kj);
		  if (kj < ix)
		    --ix;
		}
	      else
		++kj;
	    }
	}
      for (size_t kr=0; kr<adj_surfs.size(); ++kr)
	{
	  size_t kj;
	  for (kj=0; kj<surfaces_.size(); ++kj)
	    if (surfaces_[kj].get() == adj_surfs[kr])
	      break;
	  if (kj < surfaces_.size())
	    surfaces_.erase(surfaces_.begin()+kj);
	}
    }

  for (size_t kj=0; kj<regions_.size(); ++kj)
    regions_[kj]->setVisited(false);
}

//===========================================================================
void RevEng::mergeSurfaces()
//===========================================================================
{
  // Sort according to number of associated points
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      int nmb1 = surfaces_[ki]->numPoints();
      for (size_t kj=ki+1; kj<surfaces_.size(); ++kj)
	{
	  int nmb2 = surfaces_[kj]->numPoints();
	  if (nmb2 > nmb1)
	    std::swap(surfaces_[ki], surfaces_[kj]);
	}
    }

  std::ofstream ofm("merged_sfs.g2");
  std::ofstream of("surfs0.g2");
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
      surf->writeStandardHeader(of);
      surf->write(of);

      int nreg = surfaces_[ki]->numRegions();
      for (int ka=0; ka<nreg; ++ka)
	{
	  RevEngRegion *reg =  surfaces_[ki]->getRegion(ka);
	  int nmb = reg->numPoints();
	  of << "400 1 0 4 255 0 0 255" << std::endl;
	  of << nmb << std::endl;
	  for (int kb=0; kb<nmb; ++kb)
	    {
	      of << reg->getPoint(kb)->getPoint() << std::endl;
	    }
	}
    }

  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      std::ofstream ofn("surfsn.g2");
      for (size_t kh=0; kh<surfaces_.size(); ++kh)
	{
	  shared_ptr<ParamSurface> surf = surfaces_[kh]->surface();
	  surf->writeStandardHeader(ofn);
	  surf->write(ofn);
	}
      
      // Identify possible merge candidates
      vector<size_t> cand_ix;
      cand_ix.push_back(ki);
      for (size_t kj=ki+1; kj<surfaces_.size(); ++kj)
	{
	  if (surfaces_[ki]->isCompatible(surfaces_[kj].get(), anglim_, approx_tol_))
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
	  shared_ptr<HedgeSurface> merged_surf = doMerge(cand_ix);
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

	      surfm->writeStandardHeader(ofm);
	      surfm->write(ofm);
	      int nreg = merged_surf->numRegions();
	      for (int ka=0; ka<nreg; ++ka)
		{
		  RevEngRegion *reg =  merged_surf->getRegion(ka);
		  int nmb = reg->numPoints();
		  ofm << "400 1 0 4 255 0 0 255" << std::endl;
		  ofm << nmb << std::endl;
		  for (int kb=0; kb<nmb; ++kb)
		    {
		      ofm << reg->getPoint(kb)->getPoint() << std::endl;
		    }
		}
	     
		  
	      std::swap(surfaces_[cand_ix[0]], merged_surf);
	      if (cand_ix[0] > ki)
		std::swap(surfaces_[ki], surfaces_[cand_ix[0]]);
	      for (size_t kr=cand_ix.size()-1; kr>=1; --kr)
		{
		  surfaces_.erase(surfaces_.begin()+cand_ix[kr]);
		}
	    }
	}
      int stop_break = 1;
    }
  std::ofstream of1("surfs1.g2");
  shared_ptr<std::ofstream> ofp(new std::ofstream("planes_2.g2"));
  shared_ptr<std::ofstream> ofc1(new std::ofstream("cylinders_2.g2"));
  shared_ptr<std::ofstream> ofc2(new std::ofstream("cones_2.g2"));
  shared_ptr<std::ofstream> oft(new std::ofstream("tori_2.g2"));
  shared_ptr<std::ofstream> ofsp(new std::ofstream("spline_2.g2"));
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      shared_ptr<std::ofstream> ofs;
      if (surfaces_[ki]->isPlane())
	ofs = ofp;
      else if (surfaces_[ki]->isCylinder())
	ofs = ofc1;
      else if (surfaces_[ki]->isCone())
	ofs = ofc2;
      else if (surfaces_[ki]->isTorus())
	ofs = oft;
      else if (surfaces_[ki]->isSpline())
	ofs = ofsp;
      shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
      surf->writeStandardHeader(of1);
      surf->write(of1);
      surf->writeStandardHeader(*ofs);
      surf->write(*ofs);
      
      int nreg = surfaces_[ki]->numRegions();
      int nmb = 0;
      for (int ka=0; ka<nreg; ++ka)
	{
	  RevEngRegion *reg =  surfaces_[ki]->getRegion(ka);
	  nmb += reg->numPoints();
	}
      of1 << "400 1 0 4 255 0 0 255" << std::endl;
      of1 << nmb << std::endl;
      *ofs << "400 1 0 4 255 0 0 255" << std::endl;
      *ofs << nmb << std::endl;
      for (int ka=0; ka<nreg; ++ka)
	{
	  RevEngRegion *reg =  surfaces_[ki]->getRegion(ka);
	  int nmb2 = reg->numPoints();;
	  for (int kb=0; kb<nmb2; ++kb)
	    {
	      of1 << reg->getPoint(kb)->getPoint() << std::endl;
	      *ofs << reg->getPoint(kb)->getPoint() << std::endl;
	    }
	}
    }

  int stop_break2 = 1;
}

//===========================================================================
shared_ptr<HedgeSurface> RevEng::doMerge(vector<size_t>& cand_ix)
//===========================================================================
{
  size_t candsize = cand_ix.size();
  shared_ptr<HedgeSurface> merged_surf;

  std::ofstream ofp("all_merge_points.g2");
  for (size_t kh=0; kh<cand_ix.size(); ++kh)
    {
      // Try merge current candidate set
      // Collect regions and point clouds
      vector<RevEngRegion*> regions;
      vector<pair<vector<RevEngPoint*>::iterator,
		  vector<RevEngPoint*>::iterator> > points;
      BoundingBox bbox(3);
      vector<int> nmbpts;
      for (size_t ki=0; ki<cand_ix.size(); ++ki)
	{
	  HedgeSurface* surf = surfaces_[cand_ix[ki]].get();
	  vector<RevEngRegion*> reg = surf->getRegions();
	  regions.insert(regions.end(), reg.begin(), reg.end());
	  for (size_t kj=0; kj<reg.size(); ++kj)
	    {
	      nmbpts.push_back(reg[kj]->numPoints());
	      points.push_back(std::make_pair(reg[kj]->pointsBegin(),
					      reg[kj]->pointsEnd()));
	      bbox.addUnionWith(reg[kj]->boundingBox());
	      ofp << "400 1 0 0" << std::endl;
	      int npt = reg[kj]->numPoints();
	      ofp << npt << std::endl;
	      for (int ka=0; ka<npt; ++ka)
		ofp << reg[kj]->getPoint(ka)->getPoint() << std::endl;
	    }
	}

      // Distinguish between the different surface types
      shared_ptr<ParamSurface> surf;
      if (surfaces_[cand_ix[kh]]->isPlane())
	{
	  surf = doMergePlanes(points, bbox, nmbpts);
	}
      else if (surfaces_[cand_ix[kh]]->isCylinder())
	{
	  surf = doMergeCylinders(points, bbox, nmbpts);
	}
      else if (surfaces_[cand_ix[kh]]->isTorus())
	{
	  surf = doMergeTorus(points, bbox, nmbpts);
	}
      if (!surf.get())
	return merged_surf;

      std::ofstream of0("merge_surface.g2");
      surf->writeStandardHeader(of0);
      surf->write(of0);

      // Check accuracy
      double dfac = 5.0;
      std::ofstream of1("regions_merge.g2");
      std::ofstream of2("in_out_merge.g2");
      for (size_t ki=0; ki<regions.size(); )
	{
	  regions[ki]->writeRegionInfo(of1);
	  double maxd, avd;
	  int num2;
	  vector<RevEngPoint*> in, out;
	  vector<pair<double,double> > distang;
	  RevEngUtils::distToSurf(points[ki].first, points[ki].second,
				  surf, approx_tol_, maxd, avd, num2, in, out,
				  distang);

	  of2 << "400 1 0 4 155 50 50 255" << std::endl;
	  of2 << in.size() << std::endl;
	  for (size_t kr=0; kr<in.size(); ++kr)
	    of2 << in[kr]->getPoint() << std::endl;
	  of2 << "400 1 0 4 50 155 50 255" << std::endl;
	  of2 << out.size() << std::endl;
	  for (size_t kr=0; kr<out.size(); ++kr)
	    of2 << out[kr]->getPoint() << std::endl;
	  
	  double maxd_init, avd_init;
	  int num2_init;
	  regions[ki]->getAccuracy(maxd_init, avd_init, num2_init);
	  int num = regions[ki]->numPoints();

	  if (num2 < num/2 /*|| maxd > dfac*maxd_init || avd > dfac*avd_init */
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
	      regions.erase(regions.begin()+ki);
	      points.erase(points.begin()+ki);
	      nmbpts.erase(nmbpts.begin()+ki);
	    }
	  else
	    ++ki;
	  int stop_break = 1;
	}
      if (cand_ix.size() <= 1)
	break;
      if (cand_ix.size() == candsize)
	{
	  merged_surf =
	    shared_ptr<HedgeSurface>(new HedgeSurface(surf, regions));
	  for (size_t kj=0; kj<regions.size(); ++kj)
	    regions[kj]->setHedge(merged_surf.get());
	  break;
	}
      candsize = cand_ix.size();
    }

  return merged_surf;
}

//===========================================================================
void RevEng::recognizePlanes()
//===========================================================================
{
  std::ofstream planeout("plane.g2");
  double anglim = 0.1;
  size_t nmbsfs = surfaces_.size();
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      // Check type
      bool planar = regions_[ki]->planartype();
      DirectionCone normalcone = regions_[ki]->getNormalCone();
      if (!planar /*&& (normalcone.angle() > anglim ||
		    normalcone.greaterThanPi())*/)
	continue;  // Not a probable plane

      // Try to fit the point cloud with a plane
      vector<shared_ptr<HedgeSurface> > plane_sfs;
      vector<HedgeSurface*> prev_surfs;
      bool found = regions_[ki]->extractPlane(approx_tol_, min_point_region_,
					      plane_sfs, prev_surfs, planeout);
      for (size_t kr=0; kr<prev_surfs.size(); ++kr)
	{
	  size_t kj;
	  for (kj=0; kj<surfaces_.size(); ++kj)
	    if (surfaces_[kj].get() == prev_surfs[kr])
	      break;
	  if (kj < surfaces_.size())
	    surfaces_.erase(surfaces_.begin()+kj);
	}
      if (plane_sfs.size() > 0)
	{
	  surfaces_.insert(surfaces_.end(), plane_sfs.begin(), plane_sfs.end());

	  vector<RevEngRegion*> grown_regions;
	  int min_nmb = 5*min_point_region_;  // Should be set from distribution of how many
	  // points the regions have
	  vector<HedgeSurface*> adj_surfs;
	  regions_[ki]->growWithSurf(min_nmb, approx_tol_, grown_regions, adj_surfs);
	  if (grown_regions.size() > 0)
	    {
	      for (size_t kr=0; kr<grown_regions.size(); ++kr)
		{
		  size_t kj;
		  for (kj=0; kj<regions_.size(); )
		    {
		      if (kj == ki)
			{
			  ++kj;
			  continue;
			}

		      if (grown_regions[kr] == regions_[kj].get())
			{
			  regions_.erase(regions_.begin()+kj);
			  if (kj < ki)
			    --ki;
			}
		      else
			++kj;
		    }
		}
	      for (size_t kr=0; kr<adj_surfs.size(); ++kr)
		{
		  size_t kj;
		  for (kj=0; kj<surfaces_.size(); ++kj)
		    if (surfaces_[kj].get() == adj_surfs[kr])
		      break;
		  if (kj < surfaces_.size())
		    surfaces_.erase(surfaces_.begin()+kj);
		}
	    }
	}

      for (size_t kj=0; kj<regions_.size(); ++kj)
	regions_[kj]->setVisited(false);


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
  
  std::ofstream ofp("planes2.g2");
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
      surf->writeStandardHeader(ofp);
      surf->write(ofp);
    }

  if (regions_.size() > 0)
    {
      std::ofstream ofr("regions3.g2");
      for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  // BoundingBox bbox = regions_[kr]->boundingBox();
	  // if (bbox.low().dist(bbox.high()) < 0.1)
	  //   std::cout << "Small bounding box" << std::endl;
	  if (regions_[kr]->numPoints() < 5)
	    continue;
	  ofr << "400 1 0 0" << std::endl;
	  int nmb = regions_[kr]->numPoints();
	  ofr << nmb << std::endl;
	  for (int ki=0; ki<nmb; ++ki)
	    {
	      ofr << regions_[kr]->getPoint(ki)->getPoint() << std::endl;
	    }
	}
    }

  if (surfaces_.size() - nmbsfs > 1)
    mergePlanes(nmbsfs, surfaces_.size());

  std::ofstream of("surfaces0.g2");
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
      surf->writeStandardHeader(of);
      surf->write(of);
    }
}

//===========================================================================
void RevEng::recognizeCylinders()
//===========================================================================
{
  std::ofstream cylout("cylinder.g2");
  int nmbsfs = surfaces_.size();
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      // Check type
      bool cyl = regions_[ki]->cylindertype();
      // if (!cyl)
      // 	continue;
      
      // So far, try to regognize cylinders
      vector<shared_ptr<HedgeSurface> > cyl_sfs;
      vector<HedgeSurface*> prev_surfs;
      bool found = regions_[ki]->extractCylinder(approx_tol_, min_point_region_,
						 mean_edge_len_, cyl_sfs, prev_surfs,
						 cylout);
      for (size_t kr=0; kr<prev_surfs.size(); ++kr)
	{
	  size_t kj;
	  for (kj=0; kj<surfaces_.size(); ++kj)
	    if (surfaces_[kj].get() == prev_surfs[kr])
	      break;
	  if (kj < surfaces_.size())
	    {
	      surfaces_.erase(surfaces_.begin()+kj);
	      nmbsfs--;
	    }
	}
      if (cyl_sfs.size() > 0)
	{
	  surfaces_.insert(surfaces_.end(), cyl_sfs.begin(), cyl_sfs.end());
      
	  vector<RevEngRegion*> grown_regions;
	  int min_nmb = 10*min_point_region_;  // Should be set from distribution of how many
	  // points the regions have
	  vector<HedgeSurface*> adj_surfs;
	  regions_[ki]->growWithSurf(min_nmb, approx_tol_, grown_regions, adj_surfs);
	  if (grown_regions.size() > 0)
	    {
	      for (size_t kr=0; kr<grown_regions.size(); ++kr)
		{
		  size_t kj;
		  for (kj=0; kj<regions_.size(); )
		    {
		      if (kj == ki)
			{
			  ++kj;
			  continue;
			}

		      if (grown_regions[kr] == regions_[kj].get())
			{
			  regions_.erase(regions_.begin()+kj);
			  if (kj < ki)
			    ki--;
			}
		      else
			++kj;
		    }
		}
	      for (size_t kr=0; kr<adj_surfs.size(); ++kr)
		{
		  size_t kj;
		  for (kj=0; kj<surfaces_.size(); ++kj)
		    if (surfaces_[kj].get() == adj_surfs[kr])
		      break;
		  if (kj < surfaces_.size())
		    {
		      surfaces_.erase(surfaces_.begin()+kj);
		      if (kj < nmbsfs)
			nmbsfs--;
		    }
		}
	    }
	}

      for (size_t kj=0; kj<regions_.size(); ++kj)
	regions_[kj]->setVisited(false);
      
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

  
  std::ofstream ofp("cylinders2.g2");
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
      surf->writeStandardHeader(ofp);
      surf->write(ofp);
    }

  if (regions_.size() > 0)
    {
      std::ofstream ofr("regions4.g2");
      for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  // BoundingBox bbox = regions_[kr]->boundingBox();
	  // if (bbox.low().dist(bbox.high()) < 0.1)
	  //   std::cout << "Small bounding box" << std::endl;
	  if (regions_[kr]->numPoints() < 5)
	    continue;
	  ofr << "400 1 0 0" << std::endl;
	  int nmb = regions_[kr]->numPoints();
	  ofr << nmb << std::endl;
	  for (int ki=0; ki<nmb; ++ki)
	    {
	      ofr << regions_[kr]->getPoint(ki)->getPoint() << std::endl;
	    }
	}
    }

  // Now the surfaces between nmbsfs and surfaces_.size() are cylinder surfaces
  // and can be merged if they represent the same cylinder
  if (surfaces_.size() + nmbsfs > 1)
    mergeCylinders(nmbsfs, surfaces_.size());
  //    mergeCylinders(0, surfaces_.size());

  std::ofstream of("surfaces.g2");
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
      surf->writeStandardHeader(of);
      surf->write(of);
    }
}

//===========================================================================
void RevEng::mergePlanes(size_t first, size_t last)
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

  std::ofstream of("planes0.g2");
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
      surf->writeStandardHeader(of);
      surf->write(of);
    }
  
  // Find similar planes
  for (size_t ki=first; ki<last; ++ki)
    {
      vector<size_t> cand_ix;
      int code;
      ClassType type1 = surfaces_[ki]->instanceType(code);
      if (type1 != Class_Plane && type1 != Class_BoundedSurface)
	continue;

      shared_ptr<ParamSurface> surf1 = surfaces_[ki]->surface();
      shared_ptr<Plane> pla1 =
	dynamic_pointer_cast<Plane,ParamSurface>(surf1);
      if (!pla1.get())
	{
	  shared_ptr<BoundedSurface> bdsf1 =
	    dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf1);
	  if (bdsf1.get())
	    {
	      surf1 = bdsf1->underlyingSurface();
	      pla1 = dynamic_pointer_cast<Plane,ParamSurface>(surf1);
	    }
	}
      if (!pla1.get())
	continue; 

      cand_ix.push_back(ki); 
      Point norm1 = pla1->getNormal();
      Point pos1 = pla1->getPoint();
      for (size_t kj=ki+1; kj<last; ++kj)
	{
	  ClassType type2 = surfaces_[kj]->instanceType(code);
	  if (type2 != Class_Plane && type2 != Class_BoundedSurface)
	    continue;
	  
	  shared_ptr<ParamSurface> surf2 = surfaces_[kj]->surface();
	  shared_ptr<Plane> pla2 =
	    dynamic_pointer_cast<Plane,ParamSurface>(surf2);
	  if (!pla2.get())
	    {
	      shared_ptr<BoundedSurface> bdsf2 =
		dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf2);
	      if (bdsf2.get())
		{
		  surf2 = bdsf2->underlyingSurface();
		  pla2 = dynamic_pointer_cast<Plane,ParamSurface>(surf2);
		}
	    }
 	  if (!pla2.get())
	    continue;  // Should not happen
	  Point norm2 = pla2->getNormal();
	  Point pos2 = pla2->getPoint();

	  double ang = norm1.angle(norm2);
	  ang = std::min(ang, M_PI-ang);
	  Point pos2_0 = pos2 - ((pos2-pos1)*norm1)*norm1;
	  Point pos1_0 = pos1 - ((pos2-pos1)*norm2)*norm2;
	  double pdist1 = pos2.dist(pos2_0);
	  double pdist2 = pos1.dist(pos1_0);

	  if (ang > 2.0*anglim_ || pdist1 > 5.0*approx_tol_ || pdist2 > 5.0*approx_tol_)
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
	  shared_ptr<HedgeSurface> merged_surf;// = doMergePlanes(cand_ix);
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
shared_ptr<ParamSurface> RevEng::doMergePlanes(vector<pair<vector<RevEngPoint*>::iterator,
					       vector<RevEngPoint*>::iterator> > points,
					       const BoundingBox& bbox,
					       vector<int>& nmbpts)
//===========================================================================
{
  int totnmb = 0;
  for (size_t kh=0; kh<nmbpts.size(); ++kh)
    totnmb += nmbpts[kh];
  
  Point pos(0.0, 0.0, 0.0);
  Point norm(0.0, 0.0, 0.0);
  vector<RevEngPoint*> all_pts;
  all_pts.reserve(totnmb);
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      double wgt = 1.0/(double)totnmb;

      for (auto it=points[ki].first; it!=points[ki].second; ++it)
	{
	  Point curr = (*it)->getPCANormal();
	  Vector3D xyz = (*it)->getPoint();
	  pos +=  wgt*Point(xyz[0], xyz[1], xyz[2]);
	  norm += wgt*curr;
	  all_pts.push_back(*it);
	}
    }
  
  // Perform approximation with combined point set
  shared_ptr<ImplicitApprox> impl(new ImplicitApprox());
  impl->approx(points, 1);
  Point pos2, normal2;
  impl->projectPoint(pos, norm, pos2, normal2);
  std::ofstream outviz("implsf_merge.g2");
  impl->visualize(all_pts, outviz);
 
  shared_ptr<Plane> surf(new Plane(pos2, normal2));
  Point low = bbox.low();
  Point high = bbox.high();
  double len = low.dist(high);
  surf->setParameterBounds(-0.5*len, -0.5*len, 0.5*len, 0.5*len);

  return surf;
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

  std::ofstream cyl("sorted_cylinders.g2");
  for (size_t ki=first; ki<last; ++ki)
    {
      shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
      surf->writeStandardHeader(cyl);
      surf->write(cyl);
    }
  
  // Find similar cylinders
  for (size_t ki=first; ki<last; ++ki)
    {
      vector<size_t> cand_ix;
      int code;
      ClassType type1 = surfaces_[ki]->instanceType(code);
      if (type1 != Class_Cylinder && type1 != Class_BoundedSurface)
	continue;

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

      cand_ix.push_back(ki); 
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
	  shared_ptr<HedgeSurface> merged_surf;// = doMergeCylinders(cand_ix);
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
shared_ptr<ParamSurface> RevEng::doMergeCylinders(vector<pair<vector<RevEngPoint*>::iterator,
						  vector<RevEngPoint*>::iterator> > points,
						  const BoundingBox& bbox,
						  vector<int>& nmbpts)
//===========================================================================
{
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
  shared_ptr<Cylinder> surf(new Cylinder(rad, pnt, axis, Cy));
  double len = low.dist(high);
  surf->setParamBoundsV(-len, len);

  return surf;
}

//===========================================================================
shared_ptr<ParamSurface> RevEng::doMergeTorus(vector<pair<vector<RevEngPoint*>::iterator,
					      vector<RevEngPoint*>::iterator> > points,
					      const BoundingBox& bbox,
					      vector<int>& nmbpts)
//===========================================================================
{
  shared_ptr<ParamSurface> dummy;
  
  int totnmb = 0;
  for (size_t kh=0; kh<nmbpts.size(); ++kh)
    totnmb += nmbpts[kh];
  
  // Compute mean curvature and initial point in plane
  double k2mean = 0.0;
  double wgt = 1.0/(double)totnmb;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      for (auto it=points[ki].first; it!=points[ki].second; ++it)
	{
	  double kmax = (*it)->maxPrincipalCurvature();
	  k2mean += wgt*kmax;
	}
    }
  double rd = 1.0/k2mean;
  
  vector<Point> centr(totnmb);
  Point mid(0.0, 0.0, 0.0);
  size_t kr = 0;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      for (auto it=points[ki].first; it!=points[ki].second; ++it, ++kr)
	{
	  Point norm = (*it)->getPCANormal();
	  Vector3D xyz = (*it)->getPoint();
	  Point xyz2(xyz[0], xyz[1], xyz[2]);
	  centr[kr] = xyz2 + rd*norm;
	  mid += wgt*centr[kr];
	}
    }
  
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
    return dummy;
  
  Point Cx = centr[0] - mid;
  Cx -= (Cx*normal)*normal;
  Cx.normalize();
  Point Cy = Cx.cross(normal);
  
  double rad;
  Point pnt;
  RevEngUtils::computeCircPosRadius(centr, normal, Cx, Cy, pnt, rad);
  pnt -= ((pnt - pos)*normal)*normal;

  vector<Point> rotated;
  RevEngUtils::rotateToPlane(points, Cx, normal, pnt, rotated);
  Point cpos;
  double crad;
  RevEngUtils::computeCircPosRadius(rotated, Cy, Cx, normal, cpos, crad);
  Point cvec = cpos - pnt;
  double R1 = (cvec - (cvec*normal)*normal).length();
  double R2 = (cvec*normal)*normal.length();
 
  shared_ptr<Torus> surf(new Torus(R1, crad, pnt+R2*normal, normal, Cy));

  return surf;
}

//===========================================================================
void RevEng::adaptToMainAxis()
//===========================================================================
{
  vector<SurfaceProperties> sfprop;
  collectAxis(sfprop);

  // Sort surfaces into groups with roughly the same axis
  double epsang = 0.05;
  vector<DirectionCone> axis_cone;
  vector<vector<size_t> > group_ixs;
  for (size_t ki=0; ki<sfprop.size(); ++ki)
    {
      size_t kr;
      for (kr=0; kr<axis_cone.size(); ++kr)
	{
	  Point axis = sfprop[ki].dir_;
	  Point centre = axis_cone[kr].centre();
	  if (axis*centre < 0.0)
	    axis *= -1;
	  double angle = axis_cone[kr].unionAngle(axis);
	  if (angle <= epsang)
	    {
	      axis_cone[kr].addUnionWith(axis);
	      group_ixs[kr].push_back(ki);
	      break;
	    }
	}
      if (kr == axis_cone.size())
	{
	  DirectionCone cone(sfprop[ki].dir_);
	  axis_cone.push_back(cone);
	  vector<size_t> ixs;
	  ixs.push_back(ki);
	  group_ixs.push_back(ixs);
	}
    }
  int stop_break = 1;
}

//===========================================================================
void RevEng::collectAxis(vector<SurfaceProperties>& sfprop)
//===========================================================================
{
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      int code;
      ClassType type = surfaces_[ki]->instanceType(code);
      shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
      shared_ptr<ElementarySurface> elemsf =
	dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf);
      if (!elemsf.get())
	continue;
      
      Point loc, dir;
      double rad1, rad2;
      loc = elemsf->location();
      dir = elemsf->direction();
      rad1 = elemsf->radius(0.0, 0.0);   // Not good enough for cones
      rad2 = elemsf->radius2(0.0, 0.0);   // Not good enough for cones
      SurfaceProperties currprop(type, dir, loc, rad1, rad2);
      sfprop.push_back(currprop);
    }
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
  max_next_ = 80; //500;
  rfac_ = 3.0;  // Factor for radius in which to search for neighbouring points
  cfac_ = 8.0;  // Edge points from curvature is given by
  // cfac_ times the average length of triangulation edges in a vertex
  norm_plane_lim_= 0.005; // Limit for when the cone angle corresponding
  // to triangle normals indicate an edge
  zero_H_ = 0.001; //0.0001;  // When mean curvature is considered zero
  zero_K_ = 0.001; //0.0001;  // When Gauss curvature is considered zero
  zero_si_ = 0.0075; //0.001; // When shape index is considered zero
  norm_ang_lim_ = 0.1*M_PI; // Limit for when the cone angle corresponding
    // to triangle normals indicate an edge
  pca_lim_ = cness_lim_ = std::numeric_limits<double>::max();
  min_point_region_ = 200; //50; //10;  // Should be updated with regard to the total
  // number of points
  approx_tol_ = 0.001;  // Very preliminary
  anglim_ = 0.01;
  max_nmb_outlier_ = 3;
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
void RevEng::storeGrownRegions(ostream& os) const
//===========================================================================
{
  storeClassified(os);
  os << regions_.size() << std::endl;
  for (size_t ki=0; ki<regions_.size(); ++ki)
    regions_[ki]->store(os);
}

 //===========================================================================
void RevEng::readGrownRegions(istream& is)
//===========================================================================
{
  readClassified(is);
  curvatureFilter();
  int num_regions;
  is >> num_regions;
  regions_.resize(num_regions);
  for (int ki=0; ki<num_regions; ++ki)
    {
      regions_[ki] = shared_ptr<RevEngRegion>(new RevEngRegion());
      regions_[ki]->read(is, tri_sf_);
    }

  for (int ki=0; ki<num_regions; ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }
}
  

 //===========================================================================
void RevEng::storeParams(ostream& os) const
//===========================================================================
{
  os << mean_edge_len_ << " " << min_next_ << " " << rfac_ << " " << cfac_;
  os << " " << pca_lim_ << " " << cness_lim_ << " " << norm_ang_lim_;
  os << " " << norm_plane_lim_ << " " << zero_H_ << " " << zero_K_;
  os << " " << zero_si_ << " " << min_point_region_ << " " << approx_tol_ ;
  os << " " << anglim_ << " " << max_nmb_outlier_ << std::endl;
}

 //===========================================================================
void RevEng::readParams(istream& is)
//===========================================================================
{
  is >> mean_edge_len_ >> min_next_ >> rfac_ >> cfac_ >> pca_lim_ >> cness_lim_;
  is >> norm_ang_lim_ >> norm_plane_lim_ >> zero_H_ >> zero_K_ >> zero_si_;
  is >> min_point_region_ >> approx_tol_ >> anglim_ >> max_nmb_outlier_;
}
