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
#include "GoTools/compositemodel/SurfaceModelUtils.h"
#include "GoTools/utils/DirectionCone.h"
#include "GoTools/utils/MatrixXD.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/Plane.h"
#include "GoTools/geometry/Torus.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/geometry/SISLconversion.h"
#include "sislP.h"
#include <vector>
#include <fstream>
#include <iostream> // @@ debug

using namespace Go;
using std::vector;
using std::pair;
using std::istream;
using std::ostream;

typedef MatrixXD<double, 3> Matrix3D;

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
  max_next_ = std::min(80, tri_sf_->size()/200);
  max_next_ = std::max(2*min_next_, max_next_);
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


int compare_x_par(const RevEngPoint* p1, const RevEngPoint* p2)
{
  return (p1->getPoint()[0] < p2->getPoint()[0]);
  // if (p1->getPoint()[0] < p2->getPoint()[0])
  //   return -1;
  // else if (p1->getPoint()[0] > p2->getPoint()[0])
  //   return 1;
  // else
  //   return 0;
}

int compare_y_par(const RevEngPoint* p1, const RevEngPoint* p2)
{
  return (p1->getPoint()[1] < p2->getPoint()[1]);
  // if (p1->getPoint()[1] < p2->getPoint()[1])
  //   return -1;
  // else if (p1->getPoint()[1] > p2->getPoint()[1])
  //   return 1;
  // else
  //   return 0;
}

int compare_z_par(const RevEngPoint* p1, const RevEngPoint* p2)
{
  return (p1->getPoint()[2] < p2->getPoint()[2]);
  // if (p1->getPoint()[2] < p2->getPoint()[2])
  //   return -1;
  // else if (p1->getPoint()[2] > p2->getPoint()[2])
  //   return 1;
  // else
  //   return 0;
}

//===========================================================================
void RevEng::enhancePoints2()
//===========================================================================
{
  int nmbpt = tri_sf_->size();
  vector<RevEngPoint*> all_pts(nmbpt);
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint* pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      Vector3D xyz = pt->getPoint();
      Point xyz2 = Point(xyz[0], xyz[1], xyz[2]);
      all_pts[ki] = pt;
      if (ki == 0)
	bbox_ = BoundingBox(xyz2, xyz2);
      else
	bbox_.addUnionWith(xyz2);
    }

}

//===========================================================================
void RevEng::enhancePoints()
//===========================================================================
{
  enhancePoints2();
  std::cout << "Bounding box, min: " << bbox_.low() << ", max: " << bbox_.high() << std::endl;
  
  int writepoints = 0;
  int nmbpt = tri_sf_->size();
  std::ofstream of01("minc1.g2");
  //std::ofstream of02("minc2.g2");
  std::ofstream of03("maxc1.g2");
   if (writepoints) {
      //std::ofstream of04("maxc2.g2");
      of01 << "410 1 0 4 200 50 0 255" << std::endl;
      of01 << nmbpt << std::endl;
      // of02 << "410 1 0 4 200 50 0 255" << std::endl;
      // of02 << nmbpt << std::endl;
      of03 << "410 1 0 4 0 50 200 255" << std::endl;
      of03 << nmbpt << std::endl;
      // of04 << "410 1 0 4 0 50 200 255" << std::endl;
      // of04 << nmbpt << std::endl;
  }

  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);

      // Compute surface normal from triangulation
      pt->computeTriangNormal(100.0*mean_edge_len_);
      double avlen = pt->getMeanEdgLen();
    }
  
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      if (pt->nmbMonge() > 0)
	continue;  // Already enhanced

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
      vector<RevEngPoint*> nearpts;
      double local_len = pt->getMeanEdgLen(10.0*mean_edge_len_);
      double radius = rfac_*(local_len + mean_edge_len_);
      double radius2 = 0.5*radius;
      radius = std::min(radius, 20.0*mean_edge_len_);
      //radius *= 1.5; // TEST 
      //double radius = 0.5*rfac_*(local_len + mean_edge_len_);
      pt->fetchClosePoints2(radius, min_next_, max_next_, nearpts);

      Point mincvec(0.0, 0.0, 0.0), maxcvec(0.0, 0.0, 0.0);
      //      Point mincvec2(0.0, 0.0, 0.0), maxcvec2(0.0, 0.0, 0.0);

      if (nearpts.size() >= 3)
	{

	  std::ofstream of("nearpts.g2");
      if (writepoints)
	{
	  of << "400 1 0 4 255 0 0 255" << std::endl;
	  of << "1" << std::endl;
	  of << pt->getPoint() << std::endl << std::endl;
	  of << "400 1 0 4 0 255 0 255" << std::endl;
	  of << nearpts.size() << std::endl;
	  for (size_t kr=0; kr<nearpts.size(); ++kr)
	    of << nearpts[kr]->getPoint() << std::endl;
	}
      
      // Compute eigenvectors and values of covariance matrix
      nearpts.push_back(pt);
      double lambda[3];
      double eigenvec[3][3];
      RevEngUtils::principalAnalysis(nearpts, lambda, eigenvec);
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

      for (size_t kr=0; kr<nearpts.size(); ++kr)
	{
	  if (pt->pntDist(nearpts[kr]) <= radius2)
	    nearpts[kr]->addCovarianceEigen(eigen1, lambda[0], eigen2, lambda[1],
					    eigen3, lambda[2]);
	}
      
      if (writepoints)
	{
	  for (int ki=0; ki<3; ++ki)
	    {
	      Vector3D vec(eigenvec[ki][0], eigenvec[ki][1], eigenvec[ki][2]);
	      of << "410 1 0 4 0 0 0 255" << std::endl;
	      of << "1" << std::endl;
	      Vector3D curr = pt->getPoint();
	      of << curr << " " << curr+0.1*vec << std::endl;
	    }
	}
      // Compute normal and curvature using Monge patch
      // Point normal;//, mincvec, maxcvec;
      // double minc, maxc;
      // double currdist, avdist;
      // RevEngUtils::computeMonge(curr, nearpts, eigen1, eigen3, normal, mincvec, minc,
      // 				maxcvec, maxc, currdist, avdist);
      computeMonge(pt, nearpts, eigen1, eigen3, radius2);
      // Orient vectors with respect to triangulation normal
      // The normal vectors should be OK. Curvature vectors are not necessarily
      // consistent with regard to orientation

      // double minc2, maxc2;
      // RevEngUtils::TaubinCurvature(curr, nearpts, eigen1, eigen3, mincvec2, minc2,
      // 				   maxcvec2, maxc2);
      
      // if (normal*tnorm < 0.0)
      // 	normal *= -1;
      // pt->addMongeInfo(normal, mincvec, minc, maxcvec, maxc, currdist, avdist,
      // 		       zero_si_);
      // 	}
      // of01 << curr << " " << curr+mincvec << std::endl;
      // //of02 << curr << " " << curr+mincvec2 << std::endl;
      // of03 << curr << " " << curr+maxcvec << std::endl;
      //of04 << curr << " " << curr+maxcvec2 << std::endl;

      int stop_break = 1;
    }
    }

  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]); 
      if (pt->isOutlier())
	continue;
      double rp[2];
      setRp(pt, rp);
      pt->setRp(rp);
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
  
  setClassificationParams();
  vector<Vector3D> triangcorners;
  vector<Vector3D> triangplane;
  vector<Vector3D> curvaturecorners;
  vector<Vector3D> smoothnesscorners;
  vector<Vector3D> PCAcorners;
  vector<Vector3D> Rpcorners;
  std::ofstream of("triangnorm.g2");
  std::ofstream ofM("Mongenorm.g2");
  std::ofstream ofP("PCAnorm.g2");
  int write_file = 0;
  if (write_file) {
      of << "410 1 0 4 0 0 255 255" << std::endl;
      of << nmbpt << std::endl;
      ofM << "410 1 0 4 255 0 0 255" << std::endl;
      ofM << nmbpt << std::endl;
      ofP << "410 1 0 4 0 200 55 255" << std::endl;
      ofP << nmbpt << std::endl;
  }
  // of01 << "410 1 0 4 0 0 255 255" << std::endl;
  // of01 << nmbpt << std::endl;
  // of03 << "410 1 0 4 0 255 0 255" << std::endl;
  // of03 << nmbpt << std::endl;
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      Point minc = pt->minCurvatureVec();
      Point maxc = pt->maxCurvatureVec();
      Vector3D xyz = pt->getPoint();
      Point xyz2(xyz[0], xyz[1], xyz[2]);
      double avlen = pt->getMeanEdgLen();

      double fac = (pt->nmbMonge() == 0) ? 0.0 : 5.0;
      if (writepoints) {
          of01 << xyz2 << " " << xyz2+fac*avlen*minc << std::endl;
          of03 << xyz2 << " " << xyz2+fac*avlen*maxc << std::endl;
      }

      Point norm = pt->getTriangNormal();
      double ang = pt->getTriangAngle();
      if (write_file) {
          of << xyz2 << " " << xyz2+5.0*avlen*norm << std::endl;
      }
      if (pt->getTriangAngle() > norm_ang_lim_)
	triangcorners.push_back(xyz);
      Point Mnorm = pt->getMongeNormal();
      if (write_file) {
          ofM << xyz2 << " " << xyz2+5.0*avlen*Mnorm << std::endl;
      }

      Point Pnorm = pt->getPCANormal();
      if (write_file) {
          ofP << xyz2 << " " << xyz2+5.0*avlen*Pnorm << std::endl;
      }

      if (pt->isOutlier())
	continue;
      
      if (ang <= norm_plane_lim_)
	triangplane.push_back(xyz);
      double maxpc = std::max(fabs(pt->maxPrincipalCurvature()),
			      fabs(pt->minPrincipalCurvature()));
      double crvrad = 1.0/maxpc; //fabs(maxpc);
      if (crvrad < cfac_*avlen)
	curvaturecorners.push_back(xyz);

      double curved = pt->getCurvedness();
      if (curved > cness_lim_)
	smoothnesscorners.push_back(xyz);

      double sfvar = pt->getSurfaceVariation();
      if (sfvar > pca_lim_)
	PCAcorners.push_back(xyz);

      double rp = pt->getRp(1);
      if (rp > 0.03)
	Rpcorners.push_back(xyz);
    }

  if (write_file) {
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

      std::ofstream of5("smoothnesscorners.g2");
      of5 << "400 1 0 4 10 10 10 255" << std::endl;
      of5 << smoothnesscorners.size() << std::endl;
      for (size_t kj=0; kj<smoothnesscorners.size(); ++kj)
          of5 << smoothnesscorners[kj] << std::endl;

      std::ofstream of6("PCAcorners.g2");
      of6 << "400 1 0 4 10 10 10 255" << std::endl;
      of6 << PCAcorners.size() << std::endl;
      for (size_t kj=0; kj<PCAcorners.size(); ++kj)
          of6 << PCAcorners[kj] << std::endl;

      std::ofstream of7("Rpcorners.g2");
      of7 << "400 1 0 4 10 10 10 255" << std::endl;
      of7 << Rpcorners.size() << std::endl;
      for (size_t kj=0; kj<Rpcorners.size(); ++kj)
          of7 << Rpcorners[kj] << std::endl;

      std::ofstream of4("triangplane.g2");
      of4 << "400 1 0 4 200 0 200 255" << std::endl;
      of4 << triangplane.size() << std::endl;
      for (size_t kj=0; kj<triangplane.size(); ++kj)
          of4 << triangplane[kj] << std::endl;
  }

  double minptdist = std::numeric_limits<double>::max();
  double maxptdist = 0.0;
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      double ptdist = pt->getPointDistance();
      minptdist = std::min(minptdist, ptdist);
      maxptdist = std::max(maxptdist, ptdist);
    }
  vector<vector<Vector3D> > ptd(12);
  double ptdel = (maxptdist - minptdist)/(double)12;
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      double ptdist = pt->getPointDistance();
      int ix = ptdist/ptdel;
      ix = std::min(ix, 11);
      ptd[ix].push_back(pt->getPoint());
    }
  std::cout << "Min ptdist: " << minptdist << ", maxptdist: " << maxptdist << std::endl;
  for (int ka=0; ka<12; ++ka)
    std::cout << ptd[ka].size() << ", ";
  std::cout << std::endl;
  
  if (write_file) {
      std::ofstream ofd("ptdist.g2");
      for (int ka=0; ka<12; ++ka)
      {
          if (ptd[ka].size() > 0)
          {
              ofd << "400 1 0 4 ";
              for (int kb=0; kb<3; ++kb)
                  ofd << colors[ka][kb] << " ";
              ofd << " 255" << std::endl;
              ofd << ptd[ka].size() << std::endl;
              for (size_t kh=0; kh<ptd[ka].size(); ++kh)
                  ofd << ptd[ka][kh] << std::endl;
          }
      }
  }
  
  std::cout << "Start curvature filter" << std::endl;
  curvatureFilter();
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      pt->resetPointAssociation();
    }
  
  std::cout << "Finish curvature filter" << std::endl;

}

//===========================================================================
void RevEng::computeMonge(RevEngPoint* pt, std::vector<RevEngPoint*>& points,
			  Point& vec1, Point& vec2, double radius)
//===========================================================================
{
  // Transform points to coordinate system given by vec1 (x-axis) and vec2 (y-axis)
  Matrix3D mat1, mat2, rotmat;
  Vector3D vec1_2(vec1[0], vec1[1], vec1[2]);
  Vector3D vec2_2(vec2[0], vec2[1], vec2[2]);
  Vector3D xaxis(1, 0, 0);
  Vector3D zaxis(0, 0, 1);
  mat1.setToRotation(vec1_2, xaxis);
  Vector3D v1 = mat1*vec1_2;
  Vector3D vec2_3 = mat1*vec2_2;
  mat2.setToRotation(vec2_3, zaxis);
  Vector3D v2 = mat2*vec2_3;
  rotmat = mat2*mat1;
  //rotmat.identity();

  // Perform rotation and sort parameter values and z-value
  int nmbpts = (int)points.size();
  vector<double> par(2*nmbpts);
  vector<double> zval(nmbpts);
  Vector3D curr = pt->getPoint();
  for (int ki=0; ki<nmbpts; ++ki)
    {
      Vector3D dv = points[ki]->getPoint() - curr;
      Vector3D dvrot = rotmat*dv;
      //Vector3D dvrot = mat2*dvrot0;
      par[2*ki] = curr[0] + dvrot[0];
      par[2*ki+1] = curr[1] + dvrot[1];
      zval[ki] = curr[2] + dvrot[2];
    }

  // Approximate z-component by biquadratic Bezier function in x and y
  int order = 3;
  shared_ptr<SplineSurface> mongesf = RevEngUtils::surfApprox(zval, 1, par, order,
							      order, order, order);

  vector<double> coefs2(3*order*order);
  std::vector<double>::iterator cf = mongesf->coefs_begin();
  for (int ka=0; ka<order; ++ka)
    {
      double vpar = mongesf->basis_v().grevilleParameter(ka);
      for (int kb=0; kb<order; ++kb, ++cf)
	{
	  double upar = mongesf->basis_u().grevilleParameter(kb);
	  coefs2[(ka*order+kb)*3] = upar;
	  coefs2[(ka*order+kb)*3+1] = vpar;
	  coefs2[(ka*order+kb)*3+2] = *cf;
	}
    }
  shared_ptr<SplineSurface> tmp(new SplineSurface(order, order, order, order, 
						  mongesf->basis_u().begin(),
						  mongesf->basis_v().begin(), &coefs2[0], 3));
  int writesurface = 0;
  if (writesurface)
    {
      std::ofstream of("approx_sf.g2");
      tmp->writeStandardHeader(of);
      tmp->write(of);
      of << "400 1 0 4 0 255 0 255" << std::endl;
      of << 1 << std::endl;
      of << curr << std::endl;
      of << "400 1 0 4 255 0 0 255" << std::endl;
      of << nmbpts << std::endl;
      for (int ka=0; ka<nmbpts; ++ka)
	{
	  Point tmppt(par[2*ka], par[2*ka+1], zval[ka]);
	  of << tmppt << std::endl;
	}
    }
  
  // Compute surface normal 
  double avdist = 0.0;
  for (int ki=0; ki<nmbpts; ++ki)
    {
      Point pos;
      mongesf->point(pos, par[2*ki], par[2*ki+1]);
      avdist += fabs(zval[ki] - pos[0]);
    }
  avdist /= (double)nmbpts;

  std::ofstream of2("Monge_curvature.g2");
  std::ofstream of3("Monge_curvature2.g2");
  vector<Point> monge1, monge2, monge3, monge4;
   for (size_t kr=0; kr<points.size(); ++kr)
    {
      if (pt->pntDist(points[kr]) > radius)
	continue;

      Point triang_norm = points[kr]->getTriangNormal();
      vector<Point> der(3);
      mongesf->point(der, par[2*kr], par[2*kr+1], 1);
      Vector3D norm(-der[1][0], -der[2][0], 1.0);
      norm.normalize();

      // Accuracy of approximation
      double currdist = fabs(zval[kr] - der[0][0]);
  
      // Compute principal curvatures in curr
      shared_ptr<SISLSurf> sislsf(GoSurf2SISL(*mongesf, false));
      int left1 = 0, left2 = 0;
      int stat = 0;
      double minc, maxc;
      double d1[2], d2[2];
      s2542(sislsf.get(), 0, 0, 0, &par[0], &left1, &left2, &minc, &maxc, d1, d2, &stat);
      Vector3D du(1.0, 0.0, der[1][0]);
      Vector3D dv(0.0, 1.0, der[2][0]);
      Vector3D cvec1 = d1[0]*du + d1[1]*dv;
      Vector3D cvec2 = d2[0]*du + d2[1]*dv;

      // Vector3D origin(par[0], par[1], zval[0]);
      // of << "410 1 0 4 0 0 0 255" << std::endl;
      // of << "1" << std::endl;
      // of << origin << " " << origin+norm << std::endl;

      // of << "410 1 0 4 0 55 155 255" << std::endl;
      // of << "1" << std::endl;
      // of << origin << " " << origin+cvec1 << std::endl;

  
      // of << "410 1 0 4 155 55 0 255" << std::endl;
      // of << "1" << std::endl;
      // of << origin << " " << origin+cvec2 << std::endl;

  
  
      // Transform results to original coordinate system
      Matrix3D mat3, mat4, rotmat2;
      mat4.setToRotation(zaxis, vec2_3);
      mat3.setToRotation(xaxis, vec1_2);
      rotmat2 = mat3*mat4;
      //rotmat2.identity();
      //Vector3D norm0 = mat4*norm;
      Vector3D norm2 = rotmat2*norm;
      Point normal = Point(norm2[0], norm2[1], norm2[2]);
      if (triang_norm.length() > 1.0e-10 && normal*triang_norm < 0.0)
	normal *= -1;
  
      Vector3D cvec3 = rotmat2*cvec1;
      Point mincvec = Point(cvec3[0], cvec3[1],cvec3[2]); 
      Vector3D cvec4 = rotmat2*cvec2;
      Point maxcvec = Point(cvec4[0], cvec4[1],cvec4[2]);
      points[kr]->addMongeInfo(normal, mincvec, minc, maxcvec, maxc, currdist, avdist,
      		       zero_si_);

      Vector3D xyz = points[kr]->getPoint();
      Point xyz2 = Point(xyz[0], xyz[1], xyz[2]);
      Vector3D der2(der[0][0], der[0][1], der[0][2]);
      Vector3D der3 = rotmat2*der2;
      Point der4(der3[0], der3[1], der3[2]);
      monge1.push_back(xyz2);
      monge1.push_back(xyz2+mincvec);
      monge2.push_back(xyz2);
      monge2.push_back(xyz2+maxcvec);
      monge3.push_back(der4);
      monge3.push_back(der4+mincvec);
      monge4.push_back(der4);
      monge4.push_back(der4+maxcvec);
    }

   int writeMonge = 0;

   if (writeMonge)
     {
       of2 << "410 1 0 4 255 0 0 255" << std::endl;
       of2 << monge1.size()/2 << std::endl;
       for (size_t kr=0; kr<monge1.size(); kr+=2)
	 of2 << monge1[kr] << " " << monge1[kr+1] << std::endl;
       of2 << "410 1 0 4 0 0 255 255" << std::endl;
       of2 << monge2.size()/2 << std::endl;
       for (size_t kr=0; kr<monge2.size(); kr+=2)
	 of2 << monge2[kr] << " " << monge2[kr+1] << std::endl;
   
       of3 << "400 1 0 4 0 255 0 255" << std::endl;
       of3 << monge3.size()/2 << std::endl;
       for (size_t kr=0; kr<monge3.size(); kr+=2)
	 of3 << monge3[kr] << std::endl;
       of3 << "410 1 0 4 255 0 0 255" << std::endl;
       of3 << monge3.size()/2 << std::endl;
       for (size_t kr=0; kr<monge3.size(); kr+=2)
	 of3 << monge3[kr] << " " << monge3[kr+1] << std::endl;
       of3 << "410 1 0 4 0 0 255 255" << std::endl;
       of3 << monge4.size()/2 << std::endl;
       for (size_t kr=0; kr<monge4.size(); kr+=2)
	 of3 << monge4[kr] << " " << monge4[kr+1] << std::endl;
     }
  int stop_break = 1;
}
 
//===========================================================================
void RevEng::setRp(RevEngPoint* first, double rp[2])
//===========================================================================
{
  if (first->nmbMonge() == 0)
    return;

  // Fetch associated triangles
  vector<RevEngPoint*> next;
  first->setVisited();
  next.push_back(first);

  vector<vector<int> > tri;
  first->getAttachedTriangles(tri);
  rp[0] = computeRp(first, tri);
  size_t ki=0;
  for (int kb=1; kb<2; ++kb)
    {
      size_t nn = next.size();
      for (; ki<nn; ++ki)
	{
	  vector<ftSamplePoint*> curr_next = next[ki]->getNeighbours();
	  for (size_t kj=0; kj<curr_next.size(); ++kj)
	    {
	      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(curr_next[kj]);
	      if (pt->visited())
		continue;
	      next.push_back(pt);
	      vector<vector<int> > tri2;
	      pt->getAttachedTriangles(tri2);
	      for (size_t kr=0; kr<tri2.size(); ++kr)
		{
		  auto it = std::find(tri.begin(), tri.end(), tri2[kr]);
		  if (it == tri.end())
		    tri.push_back(tri2[kr]);
		}
	    }
	}
      rp[kb] = computeRp(first, tri);
    }

  for (size_t kr=0; kr<next.size(); ++kr)
    next[kr]->unsetVisited();
}

//===========================================================================
double RevEng::computeRp(RevEngPoint* first, vector<vector<int> >& tri)
//===========================================================================
{
  double rp = 0.0;
  Point norm = first->getMongeNormal();
  Vector3D pnt0 = first->getPoint();
  Point pnt(pnt0[0], pnt0[1], pnt0[2]);
  size_t ntri = tri.size();
  Point cvecmax = first->maxCurvatureVec();
  double kmax = first->maxPrincipalCurvature();
  Point cvecmin = first->minCurvatureVec();
  double kmin = first->minPrincipalCurvature();
  for (size_t ki=0; ki<ntri; ++ki)
    {
      RevEngPoint* pt1 = dynamic_cast<RevEngPoint*>((*tri_sf_)[tri[ki][0]]);
      RevEngPoint* pt2 = dynamic_cast<RevEngPoint*>((*tri_sf_)[tri[ki][1]]);
      RevEngPoint* pt3 = dynamic_cast<RevEngPoint*>((*tri_sf_)[tri[ki][2]]);
      Point tnorm0 = pt1->getTriangNormal() + pt2->getTriangNormal() + pt3->getTriangNormal();
      Vector3D xyz1 = pt1->getPoint();
      Point pos1(xyz1[0], xyz1[1], xyz1[2]);
      Vector3D xyz2 = pt2->getPoint();
      Point pos2(xyz2[0], xyz2[1], xyz2[2]);
      Vector3D xyz3 = pt3->getPoint();
      Point pos3(xyz3[0], xyz3[1], xyz3[2]);
      Point btri = (pos1 + pos2 + pos3)/3.0;
      Point vec1 = pos2 - pos1;
      Point vec2 = pos3 - pos1;
      Point tnorm = vec1.cross(vec2);
      if (tnorm0*tnorm < 0.0)
	tnorm *= -1;
      tnorm.normalize_checked();
      Point vec3 = btri - pnt;
      vec3.normalize_checked();
      double fac1 = vec3*cvecmax;
      double fac2 = vec3*cvecmin;
      double curv = fac1*kmax + fac2*kmin;
      double div = 1.0 - pnt.dist2(btri)*curv*curv;
      double tmp = (div > 1.0e-12) ? 1.0 - ((tnorm*norm)/sqrt(div)) : 0.0;
      rp += fabs(tmp);
    }
  rp /= (double)ntri;
  return rp;
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
  double hfac = 14.0/(15.0*nmbpts); //6.0/(nmbpts*7.0);
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
  std::cout << "cness_lim_: " << cness_lim_ << ", pca_lim_: " << pca_lim_ << std::endl;
  // std::cout << "Give rpix: " << std::endl;
  // std::cin >> rpix_;
  rpix_ = 1;

#if 0
  std::cout << "Give rpfac: " << std::endl;
  std::cin >> rpfac_;
  std::cout << "Give ffac: " << std::endl;
  std::cin >> ffac_;
  std::cout << "Give sfac: " << std::endl;
  std::cin >> sfac_;
#else
  std::cout << "Turned off setting rpfac_, ffac_ & sfac_ by the user (command line input)." << std::endl;
#endif

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
  //setClassificationParams();

  vector<vector<Vector3D> > class_pts(9);
  vector<vector<Vector3D> > class_shape(10);
  vector<vector<Vector3D> > class_pfs(10);
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

      int ps_type = PS_UNDEF;
      double fpa = pt->getfpa();
      double spa = pt->getspa();
      if (fpa <= ffac_)
	{
	  ps_type = PS_PLANE;
	  class_pfs[0].push_back(xyz);
	}
      else if (spa <= -1+sfac_)
	{
	  ps_type = PS_UC;
	  class_pfs[1].push_back(xyz);
	}
      else if (spa < -0.5-sfac_)
	{
	  ps_type = PS_EC;
	  class_pfs[2].push_back(xyz);
	}
      else if (spa <= -0.5+sfac_)
	{
	  ps_type = PS_PC;
	  class_pfs[3].push_back(xyz);
	}
      else if (spa < -sfac_)
	{
	  ps_type = PS_HC;
	  class_pfs[4].push_back(xyz);
	}
      else if (spa <= sfac_)
	{
	  ps_type = PS_HS;
	  class_pfs[5].push_back(xyz);
	}
      else if (spa < 0.5-sfac_)
	{
	  ps_type = PS_HX; 
	  class_pfs[6].push_back(xyz);
	}
      else if (spa <= 0.5+sfac_)
	{
	  ps_type = PS_PX;
	  class_pfs[7].push_back(xyz);
	}
      else if (spa < 1-sfac_)
	{
	  ps_type = PS_EX;
	  class_pfs[8].push_back(xyz);
	}
      else
	{
	  ps_type = PS_UX;
	  class_pfs[9].push_back(xyz);
	}

      // Edge classification with curvedness
      double curved = pt->getCurvedness();
      int c2_edge = (curved > cness_lim_) ? C2_EDGE : C2_NOT_EDGE;

      // Edge classification with surface variation
      double var = pt->getSurfaceVariation();
      int pca_edge = (var > pca_lim_) ? PCA_EDGE : PCA_NOT_EDGE;

      double rp = pt->getRp(rpix_);
      int rp_edge =  (rp >= rpfac_) ? RP_EDGE : RP_NOT_EDGE;

      // Store classification in point
      pt->setClassification(ctype, c1_edge, si_type, c2_edge, pca_edge, ps_type, rp_edge);
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
  
   // Point classification association
  std::ofstream of4("ptsclass_segments.g2");
  for (int ka=0; ka<10; ++ka)
    {
      of4 << "400 1 0 4 ";
      for (int kc=0; kc<3; ++kc)
	of4 << colors[ka][kc] << " ";
      of4 << "255" << std::endl;
      of4 << class_pfs[ka].size() << std::endl;
      for (size_t kr=0; kr<class_pfs[ka].size(); ++kr)
	of4 << class_pfs[ka][kr] << std::endl;
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

#if 0
  std::cout << "New tolerance: " << std::endl;
  std::cin >> approx_tol_;
#else
  std::cout << "Turned off setting approx_tol_ by the user (command line input)." << std::endl;
#endif

  std::ofstream ofd("ptdist.g2");
  double ptd1 = pointdist[0];
  double ptd2 = pointdist[dix];
  int nmbd = 12;
  double pdel = (ptd2 - ptd1)/(double)nmbd;
  vector<vector<Vector3D> > ptrs(nmbd);
  for (int ki=0; ki<nmbpts; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      Vector3D xyz = pt->getPoint();
      double dd = pt->getPointDistance();
      int ix = (int)((dd-ptd1)/pdel);
      ix = std::min(ix, nmbd-1);
      ptrs[ix].push_back(xyz);
    }

  for (int ki=0; ki<nmbd; ++ki)
    {
      if (ptrs[ki].size() > 0)
	{
	  ofd << "400 1 0 4 " << colors[ki][0] << " " << colors[ki][1] << " ";
	  ofd << colors[ki][2] << " 255"  << std::endl;
	  ofd << ptrs[ki].size();
	  for (size_t kr=0; kr<ptrs[ki].size(); ++kr)
	    ofd << ptrs[ki][kr] << std::endl;
	}
    }

  std::sort(regions_.begin(), regions_.end(), sort_region);
  if (regions_.size() > 0)
    {
      std::cout << "Regions 1, size: " << regions_.size() << std::endl;
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

  bool added_phase= false;
  // Sort regions according to number of points
  std::sort(regions_.begin(), regions_.end(), sort_region);

  min_point_region_ = setSmallRegionNumber();
  std::cout << "Min point region: " << min_point_region_ << std::endl;
  
  if (added_phase)
    {
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
      std::ofstream of1_1("regions1_1.g2");
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
	      of1_1 << "400 1 0 0" << std::endl;
	      int nmb = regions_[kr]->numPoints();
	      of1_1 << nmb << std::endl;
	      for (int ki=0; ki<nmb; ++ki)
		{
		  of1_1 << regions_[kr]->getPoint(ki)->getPoint() << std::endl;
		}
	    }
	}
      std::ofstream ofs1_1("small_regions1_1.g2");
      ofs1_1 << "400 1 0 4 0 0 0 255" << std::endl;
      ofs1_1 << small.size() << std::endl;
      for (size_t kr=0; kr<small.size(); ++kr)
	ofs1_1 << small[kr] << std::endl;
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
      std::ofstream of1_2("regions1_2.g2");
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
	      of1_2 << "400 1 0 0" << std::endl;
	      int nmb = regions_[kr]->numPoints();
	      of1_2 << nmb << std::endl;
	      for (int ki=0; ki<nmb; ++ki)
		{
		  of1_2 << regions_[kr]->getPoint(ki)->getPoint() << std::endl;
		}
	    }
	}
      std::ofstream ofs1_2("small_regions1_2.g2");
      ofs1_2 << "400 1 0 4 0 0 0 255" << std::endl;
      ofs1_2 << small.size() << std::endl;
      for (size_t kr=0; kr<small.size(); ++kr)
	ofs1_2 << small[kr] << std::endl;
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
      vector<shared_ptr<RevEngRegion> > added_groups;
      regions_[kr]->splitComposedRegions(classtype, added_groups);
      regions_.insert(regions_.end(), added_groups.begin(), added_groups.end());
    }
  
  std::sort(regions_.begin(), regions_.end(), sort_region);
  
  if (regions_.size() > 0)
    {
      std::cout << "Regions 2, size: " << regions_.size() << std::endl;
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
  
  // bool thirdpass = false;
  // if (thirdpass)
  //   {
  //     // Third pass: Verify/update regions with surface approximation
  //     std::ofstream of2("seed_points.g2");
  //     for (size_t kr=0; kr<regions_.size(); ++kr)
  // 	{
  // 	  if (regions_[kr]->numPoints() < min_point_region_)
  // 	    break;

  // 	  // Fetch seed point
  // 	  RevEngPoint *seed = regions_[kr]->seedPoint();
  // 	  if (seed)
  // 	    {
  // 	      of2 << "400 1 0 4 255 0 0 255" << std::endl;
  // 	      of2 << "1" << std::endl;
  // 	      of2 << seed->getPoint() << std::endl;
  // 	    }

  // 	  vector<RevEngPoint*> outpts;
  // 	  double local_len = seed->getMeanEdgLen();
  // 	  double radius = 2.0*rfac_*local_len;
  // 	  int mnext = std::max(min_next_, regions_[kr]->numPoints()/100);
  // 	  regions_[kr]->growLocal(seed, approx_tol_, radius, mnext, outpts);
  // 	  int stop_break = 1;
  // 	}
  
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
    // }

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

  bool joinreg = true;
  if (joinreg)
    {
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
  std::ofstream sphout("sphere.g2");
  std::ofstream coneout("cone.g2");
  std::ofstream torout("torus.g2");
  std::ofstream splout("spline.g2");
  std::ofstream linsweep_out("linearswept.g2");
  double frac = 0.75;   // Fraction of points with a certain property
  double angfac = 10.0;
  double angtol = -1; //angfac*anglim_;
  
  //min_point_region_ = 200;  // Testing
  int min_point_in = 10; //20;
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      if (regions_[ki]->numPoints() < min_point_region_)
	continue;
      int classtype = regions_[ki]->getClassification();
      
      std::ofstream of1("region.g2");
      regions_[ki]->writeRegionInfo(of1);
      std::ofstream of2("unitsphere.g2");
      regions_[ki]->writeUnitSphereInfo(of2);

      bool found1 = false, found2 = false, found3 = false, found4 = false;
      bool found5 = false, found6 = false, foundls = false;
      if (true) //regions_[ki]->possiblePlane(2.0*angfac*anglim_, frac))
      {
	vector<shared_ptr<HedgeSurface> > plane_sfs;
	vector<HedgeSurface*> prev_surfs;
	vector<vector<RevEngPoint*> > out_groups;
	found1 = regions_[ki]->extractPlane(approx_tol_, angtol, min_point_in,
					    plane_sfs, prev_surfs, out_groups,
					    planeout);

	for (size_t kr=0; kr<out_groups.size(); ++kr)
	  {
	    shared_ptr<RevEngRegion> reg(new RevEngRegion(classtype,
							  out_groups[kr]));
	    regions_.push_back(reg);
	  }	  
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
	if (found1 && regions_[ki]->getMaxSfDist() < 0.5*approx_tol_)
	  continue;
      }

      
      if (true) //regions_[ki]->possibleCylinder(angfac*anglim_, frac))
	{
	  vector<shared_ptr<HedgeSurface> > cyl_sfs;
	  vector<HedgeSurface*> prev_surfs;
	  vector<vector<RevEngPoint*> > out_groups;
	   found2 = regions_[ki]->extractCylinder(approx_tol_, min_point_in,
						  angtol, mean_edge_len_,
						  cyl_sfs, prev_surfs,
						  out_groups, cylout);

      
	for (size_t kr=0; kr<out_groups.size(); ++kr)
	  {
	    shared_ptr<RevEngRegion> reg(new RevEngRegion(classtype,
							  out_groups[kr]));
	    regions_.push_back(reg);
	  }	  
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
	if (found2 && regions_[ki]->getMaxSfDist() < 0.5*approx_tol_)
	  continue;
	 }


      if (regions_[ki]->hasSweepInfo() && regions_[ki]->sweepType() == 1)
	{
	  vector<shared_ptr<HedgeSurface> > linsweep_sfs;
	  vector<HedgeSurface*> prev_surfs;
	  foundls = regions_[ki]->extractLinearSweep(approx_tol_, min_point_in,
						     angtol, linsweep_sfs, 
						     prev_surfs, linsweep_out);
	   for (size_t kr=0; kr<prev_surfs.size(); ++kr)
	     {
	       size_t kj;
	       for (kj=0; kj<surfaces_.size(); ++kj)
		 if (surfaces_[kj].get() == prev_surfs[kr])
		   break;
	       if (kj < surfaces_.size())
		 surfaces_.erase(surfaces_.begin()+kj);
	     }

	   if (linsweep_sfs.size() > 0)
	     surfaces_.insert(surfaces_.end(), linsweep_sfs.begin(), linsweep_sfs.end());
	   
	if (foundls && regions_[ki]->getMaxSfDist() < 0.5*approx_tol_)
	  continue;
	}
      
      if (true) //regions_[ki]->possibleSphere(angfac*anglim_, frac))
	 {
	   vector<shared_ptr<HedgeSurface> > sph_sfs;
	   vector<HedgeSurface*> prev_surfs;
	  vector<vector<RevEngPoint*> > out_groups;
	   found3 = regions_[ki]->extractSphere(approx_tol_, min_point_in,
						angtol, mean_edge_len_,
						sph_sfs, prev_surfs,
						out_groups, sphout);

      
	for (size_t kr=0; kr<out_groups.size(); ++kr)
	  {
	    shared_ptr<RevEngRegion> reg(new RevEngRegion(classtype,
							  out_groups[kr]));
	    regions_.push_back(reg);
	  }	  
	   for (size_t kr=0; kr<prev_surfs.size(); ++kr)
	     {
	       size_t kj;
	       for (kj=0; kj<surfaces_.size(); ++kj)
		 if (surfaces_[kj].get() == prev_surfs[kr])
		   break;
	       if (kj < surfaces_.size())
		 surfaces_.erase(surfaces_.begin()+kj);
	     }
	   if (sph_sfs.size() > 0)
	     surfaces_.insert(surfaces_.end(), sph_sfs.begin(), sph_sfs.end());
	if (found3 && regions_[ki]->getMaxSfDist() < 0.5*approx_tol_)
	  continue;
	 }

       
      if (regions_[ki]->possibleCone(angfac*anglim_, frac))
	 {
	   vector<shared_ptr<HedgeSurface> > cone_sfs;
	   vector<HedgeSurface*> prev_surfs;
	  vector<vector<RevEngPoint*> > out_groups;
	   found4 = regions_[ki]->extractCone(approx_tol_, min_point_in,
					      angtol, mean_edge_len_,
					      cone_sfs, prev_surfs,
					      out_groups, coneout);

      
	for (size_t kr=0; kr<out_groups.size(); ++kr)
	  {
	    shared_ptr<RevEngRegion> reg(new RevEngRegion(classtype,
							  out_groups[kr]));
	    regions_.push_back(reg);
	  }	  
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
	if (found4 && regions_[ki]->getMaxSfDist() < 0.5*approx_tol_)
	  continue;
	 }

       
       if (regions_[ki]->possibleTorus(angfac*anglim_, frac))  // Not the final parameters
	{
	  vector<shared_ptr<HedgeSurface> > tor_sfs;
	  vector<HedgeSurface*> prev_surfs;
	  vector<vector<RevEngPoint*> > out_groups;
	  found5 = regions_[ki]->extractTorus(approx_tol_, min_point_in,
					      angtol, mean_edge_len_,
					      tor_sfs, prev_surfs,
					      out_groups, torout);

	for (size_t kr=0; kr<out_groups.size(); ++kr)
	  {
	    shared_ptr<RevEngRegion> reg(new RevEngRegion(classtype,
							  out_groups[kr]));
	    regions_.push_back(reg);
	  }	  
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
	if (found5 && regions_[ki]->getMaxSfDist() < 0.5*approx_tol_)
	  continue;
	 }

       if (true)
	 {
	  vector<shared_ptr<HedgeSurface> > spl_sfs;
	  vector<HedgeSurface*> prev_surfs;
	  vector<vector<RevEngPoint*> > out_groups;
	  found6 = regions_[ki]->extractFreeform(approx_tol_, min_point_in,
						 angtol, mean_edge_len_,
						 spl_sfs, prev_surfs,
						 out_groups, splout);

	for (size_t kr=0; kr<out_groups.size(); ++kr)
	  {
	    shared_ptr<RevEngRegion> reg(new RevEngRegion(classtype,
							  out_groups[kr]));
	    regions_.push_back(reg);
	  }	  
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
       
       if ((!found1) && (!found2) && (!found3) && (!found4) &&
	   (!found5) && (!found6) &&  (!foundls))
	 {
	   int write_sub_triang = 0;
	   if (write_sub_triang)
	     {
	       std::ofstream ofreg("sub_triang.txt");
	       regions_[ki]->writeSubTriangulation(ofreg);
	     }

	   int class_type = regions_[ki]->getClassification();
	   if (class_type != CLASSIFICATION_SHAPEINDEX)
	     {
	       // Split regions with shape index
	       vector<shared_ptr<RevEngRegion> > updated_regs;
	       regions_[ki]->splitWithShapeIndex(updated_regs);
	       if (updated_regs.size() > 1)
		 {
		   std::sort(updated_regs.begin(), updated_regs.end(), sort_region);
		   
		   std::ofstream ofupd("updated_subreg.g2");
		   vector<Vector3D> small_upd;
		   for (size_t kr=0; kr<updated_regs.size(); ++kr)
		     {
		       int nmb = updated_regs[kr]->numPoints();
		       if (nmb < 50)
			 {
			   for (int ki=0; ki<nmb; ++ki)
			     small_upd.push_back(updated_regs[kr]->getPoint(ki)->getPoint());
			 }
		       else
			 {
			   ofupd << "400 1 0 0" << std::endl;
			   ofupd << nmb << std::endl;
			   for (int ka=0; ka<nmb; ++ka)
			     {
			       ofupd << updated_regs[kr]->getPoint(ka)->getPoint() << std::endl;
			     }
			 }
		     }
		   std::ofstream ofupds("small_upd.g2");
		   ofupds << "400 1 0 4 0 0 0 255" << std::endl;
		   ofupds << small_upd.size() << std::endl;
		   for (size_t kr=0; kr<small_upd.size(); ++kr)
		     ofupds << small_upd[kr] << std::endl;
      
		   std::swap(regions_[ki], updated_regs[0]);
		   regions_.insert(regions_.end(), updated_regs.begin()+1, updated_regs.end());
		   --ki; // Repeat surface generation for the reduced point group
		 }
	     }
	   //regions_[ki]->implicitizeSplit();
	   int stop_no_surf = 1;
	 }
       
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
		  regions_[kr]->writeSurface(of);
		  // HedgeSurface *sf = regions_[kr]->getSurface(0);
		  // shared_ptr<ParamSurface> sf2 = sf->surface();
		  // sf2->writeStandardHeader(of);
		  // sf2->write(of);
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
      std::cout << "Number of surfaces: " << surfaces_.size() << std::endl;

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
	growSurface(ki);
    }
      std::cout << "Number of regions, post grow with surf: " << regions_.size() << std::endl;
      std::cout << "Number of surfaces: " << surfaces_.size() << std::endl;

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
		  regions_[kr]->writeSurface(of2);
		  // HedgeSurface *sf = regions_[kr]->getSurface(0);
		  // shared_ptr<ParamSurface> sf2 = sf->surface();
		  // sf2->writeStandardHeader(of2);
		  // sf2->write(of2);
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

  bool adjust_sf = true; //false;
  if (adjust_sf)
    {
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      if (regions_[ki]->hasSurface())
	regions_[ki]->adjustWithSurf(approx_tol_, 10.0*anglim_);
    }
      std::cout << "Number of regions, post grow with surf: " << regions_.size() << std::endl;
      std::cout << "Number of surfaces: " << surfaces_.size() << std::endl;

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
		  regions_[kr]->writeSurface(of5);
		  // HedgeSurface *sf = regions_[kr]->getSurface(0);
		  // shared_ptr<ParamSurface> sf2 = sf->surface();
		  // sf2->writeStandardHeader(of5);
		  // sf2->write(of5);
		}
	    }
	}
      std::ofstream ofs5("small_regions5.g2");
      ofs5 << "400 1 0 4 0 0 0 255" << std::endl;
      ofs5 << small.size() << std::endl;
      for (size_t kr=0; kr<small.size(); ++kr)
	ofs5 << small[kr] << std::endl;
    }
    }
      std::cout << "Number of regions, pre grow with surf: " << regions_.size() << std::endl;
      std::cout << "Number of surfaces: " << surfaces_.size() << std::endl;
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      if (regions_[ki]->hasSurface())
	growSurface(ki);
    }
      std::cout << "Number of regions, post grow with surf: " << regions_.size() << std::endl;
      std::cout << "Number of surfaces: " << surfaces_.size() << std::endl;

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
		  regions_[kr]->writeSurface(of2);
		  // HedgeSurface *sf = regions_[kr]->getSurface(0);
		  // shared_ptr<ParamSurface> sf2 = sf->surface();
		  // sf2->writeStandardHeader(of2);
		  // sf2->write(of2);
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

  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

  std::cout << "Pre merge2: " << surfaces_.size() << " surfaces" << std::endl;
  mergeSplineSurfaces();
  std::cout << "Post merge2: " << surfaces_.size() << " surfaces" << std::endl;

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
      vector<double> cand_score;
      cand_ix.push_back(ki);
      cand_score.push_back(0.0);
      ClassType type;
      for (size_t kj=ki+1; kj<surfaces_.size(); ++kj)
	{
	  double score;
	  if (surfaces_[ki]->isCompatible(surfaces_[kj].get(), anglim_,
					  approx_tol_, type, score))
	    {
	      cand_ix.push_back(kj);
	      cand_score.push_back(score);
	    }
	}

     if (cand_ix.size() > 1)
	{
	  // Sort accorading to compability
	  for (size_t kr=1; kr<cand_ix.size(); ++kr)
	    for (size_t kh=kr+1; kh<cand_ix.size(); ++kh)
	      {
		if (cand_score[kh] < cand_score[kr])
		  {
		    std::swap(cand_ix[kr], cand_ix[kh]);
		    std::swap(cand_score[kr], cand_score[kh]);
		  }
	      }
	  
	  std::ofstream pre("pre_merge.g2");
	  for (size_t kr=0; kr<cand_ix.size(); ++kr)
	    {
	      shared_ptr<ParamSurface> surf = surfaces_[cand_ix[kr]]->surface();
	      surf->writeStandardHeader(pre);
	      surf->write(pre);
	    }
	  shared_ptr<HedgeSurface> merged_surf = doMerge(cand_ix, cand_score,
							 type);
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
	     

	      std::sort(cand_ix.begin(), cand_ix.end());
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

  // Limit primary surfaces
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      surfaces_[ki]->limitSurf();
    }
  
  std::ofstream of1("surfs1.g2");
  shared_ptr<std::ofstream> ofp(new std::ofstream("planes_2.g2"));
  shared_ptr<std::ofstream> ofc1(new std::ofstream("cylinders_2.g2"));
  shared_ptr<std::ofstream> of0(new std::ofstream("spheres_2.g2"));
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
       else if (surfaces_[ki]->isSphere())
	ofs = of0;
     else if (surfaces_[ki]->isCone())
	ofs = ofc2;
      else if (surfaces_[ki]->isTorus())
	ofs = oft;
      else if (surfaces_[ki]->isSpline())
	ofs = ofsp;
      shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
      surf->writeStandardHeader(of1);
      surf->write(of1);
      if (ofs.get())
	{
	  surf->writeStandardHeader(*ofs);
	  surf->write(*ofs);
	}
      
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
void RevEng::mergeSplineSurfaces()
//===========================================================================
{
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      int code1;
      int type1 = surfaces_[ki]->instanceType(code1);
      if (type1 != Class_SplineSurface)
	continue;
      
      std::ofstream ofn("surfsn.g2");
      for (size_t kh=0; kh<surfaces_.size(); ++kh)
	{
	  shared_ptr<ParamSurface> surf = surfaces_[kh]->surface();
	  surf->writeStandardHeader(ofn);
	  surf->write(ofn);
	}

      vector<RevEngRegion*> regions1 = surfaces_[ki]->getRegions();
      DirectionCone cone1 = surfaces_[ki]->surface()->normalCone();
      
      // Identify possible merge candidates
      for (size_t kj=ki+1; kj<surfaces_.size(); ++kj)
	{
	  int code2;
	  int type2 = surfaces_[kj]->instanceType(code2);
	  if (type2 != Class_SplineSurface)
	    continue;

	  // Check adjacency
	  vector<RevEngRegion*> regions2 = surfaces_[kj]->getRegions();
	  DirectionCone cone2 = surfaces_[kj]->surface()->normalCone();
	  DirectionCone cone3 = cone2;
	  cone3.addUnionWith(cone1);
	  // if (cone3.greaterThanPi())
	  //   continue;
	  bool can_merge = false;
	  int num_edge_between = 0;
	  for (size_t kr=0; kr<regions1.size(); ++kr)
	    for (size_t kh=0; kh<regions2.size(); ++kh)
	      {
		bool adjacent = regions1[kr]->hasAdjacentRegion(regions2[kh]);
		if (adjacent)
		  {
		    can_merge = true;
		    bool edge_between =
		      regions1[kr]->hasEdgeBetween(regions2[kh]);
		    if (edge_between)
		      ++num_edge_between;
		    int stop_break = 1;
		  }
	      }
	  if (can_merge && num_edge_between == 0)
	    {
	      std::ofstream pre("pre_merge_spline.g2");
	      shared_ptr<ParamSurface> surf1 = surfaces_[ki]->surface();
	      surf1->writeStandardHeader(pre);
	      surf1->write(pre);
	      shared_ptr<ParamSurface> surf2 = surfaces_[kj]->surface();
	      surf2->writeStandardHeader(pre);
	      surf2->write(pre);

	      shared_ptr<HedgeSurface> merged =
		doMergeSpline(surfaces_[ki], cone1, surfaces_[kj], cone2);
	      if (merged)
		{
		  std::ofstream post("post_merge_spline.g2");
		  shared_ptr<ParamSurface> surfm = merged->surface();
		  surfm->writeStandardHeader(post);
		  surfm->write(post);
		  int nreg = merged->numRegions();
		  for (int ka=0; ka<nreg; ++ka)
		    {
		      RevEngRegion *reg =  merged->getRegion(ka);
		      int nmb = reg->numPoints();
		      post << "400 1 0 4 255 0 0 255" << std::endl;
		      post << nmb << std::endl;
		      for (int kb=0; kb<nmb; ++kb)
			{
			  post << reg->getPoint(kb)->getPoint() << std::endl;
			}
		    }
		  std::swap(surfaces_[ki], merged);
		  surfaces_.erase(surfaces_.begin()+kj);
		  kj--;
		}
	      int stop_break2 = 1;
	    }
	}
    }
  std::ofstream of1("surfs2.g2");
  std::ofstream of1n("surfs2n.g2");
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
      surf->writeStandardHeader(of1);
      surf->write(of1);
      surf->writeStandardHeader(of1n);
      surf->write(of1n);
      int nreg = surfaces_[ki]->numRegions();
      int nmb = 0;
      for (int ka=0; ka<nreg; ++ka)
	{
	  RevEngRegion *reg =  surfaces_[ki]->getRegion(ka);
	  nmb += reg->numPoints();
	}
      of1 << "400 1 0 4 255 0 0 255" << std::endl;
      of1 << nmb << std::endl;
      for (int ka=0; ka<nreg; ++ka)
	{
	  RevEngRegion *reg =  surfaces_[ki]->getRegion(ka);
	  int nmb2 = reg->numPoints();;
	  for (int kb=0; kb<nmb2; ++kb)
	    {
	      of1 << reg->getPoint(kb)->getPoint() << std::endl;
	    }
	}
    }
}

// //===========================================================================
// shared_ptr<HedgeSurface> RevEng::doMerge(vector<size_t>& cand_ix,
// 					 vector<double>& cand_score,
// 					 ClassType type)
// //===========================================================================
// {
//   size_t candsize = cand_ix.size();
//   shared_ptr<HedgeSurface> merged_surf;
//   double delscore = (cand_score[candsize-1]-cand_score[0])/(double)(candsize-1);

//   std::ofstream ofp("all_merge_points.g2");

//   // Collect regions and point clouds
//   vector<RevEngRegion*> regions;
//   vector<pair<vector<RevEngPoint*>::iterator,
// 	      vector<RevEngPoint*>::iterator> > points;
//   BoundingBox bbox(3);
//   vector<int> nmbpts;
//   vector<int> nmbreg;
//   for (size_t ki=0; ki<cand_ix.size(); ++ki)
//     {
//       HedgeSurface* surf = surfaces_[cand_ix[ki]].get();
//       vector<RevEngRegion*> reg = surf->getRegions();
//       regions.insert(regions.end(), reg.begin(), reg.end());
//       nmbreg.push_back((int)reg.size());
//       for (size_t kj=0; kj<reg.size(); ++kj)
// 	{
// 	  nmbpts.push_back(reg[kj]->numPoints());
// 	  points.push_back(std::make_pair(reg[kj]->pointsBegin(),
// 					  reg[kj]->pointsEnd()));
// 	  bbox.addUnionWith(reg[kj]->boundingBox());
// 	  ofp << "400 1 0 0" << std::endl;
// 	  int npt = reg[kj]->numPoints();
// 	  ofp << npt << std::endl;
// 	  for (int ka=0; ka<npt; ++ka)
// 	    ofp << reg[kj]->getPoint(ka)->getPoint() << std::endl;
// 	}
//     }

//   for (size_t kh=0; kh<cand_ix.size(); ++kh)
//     {
//       // Extract appropriate point clouds
//       vector<pair<vector<RevEngPoint*>::iterator,
// 		  vector<RevEngPoint*>::iterator> > curr_points;
//       vector<int> curr_nmbpts;
//       size_t ki, kj, kr;
//       for (ki=0, kr=0; ki<cand_ix.size(); ++ki)
// 	{
// 	  if (ki > 0 && cand_score[ki]-cand_score[ki-1] > delscore)
// 	    break;
// 	  for (kj=0; kj<nmbreg[ki]; ++kj, ++kr)
// 	    {
// 	      curr_points.push_back(points[kr]);
// 	      curr_nmbpts.push_back(nmbpts[kr]);
// 	    }
// 	}
//       size_t appn = ki;
      
//       // Try merge current candidate set
//       // Distinguish between the different surface types
//       shared_ptr<ParamSurface> surf;
//       if (type == Class_Plane)
// 	{
// 	  surf = doMergePlanes(curr_points, bbox, curr_nmbpts);
// 	}
//       else if (type == Class_Cylinder)
// 	{
// 	  surf = doMergeCylinders(curr_points, bbox, curr_nmbpts);
// 	}
//       else if (type == Class_Torus)
// 	{
// 	  surf = doMergeTorus(curr_points, bbox, curr_nmbpts);
// 	}
//       if (!surf.get())
// 	return merged_surf;

//       std::ofstream of0("merge_surface.g2");
//       surf->writeStandardHeader(of0);
//       surf->write(of0);

//       // Check accuracy
//       double dfac = 5.0;
//       std::ofstream of1("regions_merge.g2");
//       std::ofstream of2("in_out_merge.g2");
//       vector<vector<RevEngPoint*> > all_in;
//       double all_maxd = 0.0, all_avd = 0.0;
//       int all_inside = 0;
//       for (ki=0, kr=0; ki<cand_ix.size(); ++ki)
// 	{
// 	  for (kj=0; kj<nmbreg[ki]; ++kj, ++kr)
// 	    {
// 	      regions[kr]->writeRegionInfo(of1);
	      
// 	      double maxd, avd;
// 	      int num2;
// 	      vector<RevEngPoint*> in, out;
// 	      vector<pair<double,double> > distang;
// 	      RevEngUtils::distToSurf(points[kr].first, points[kr].second,
// 				      surf, approx_tol_, maxd, avd, num2, 
// 				      in, out, distang);

// 	      of2 << "400 1 0 4 155 50 50 255" << std::endl;
// 	      of2 << in.size() << std::endl;
// 	      for (size_t kn=0; kn<in.size(); ++kn)
// 		of2 << in[kn]->getPoint() << std::endl;
// 	      of2 << "400 1 0 4 50 155 50 255" << std::endl;
// 	      of2 << out.size() << std::endl;
// 	      for (size_t kn=0; kn<out.size(); ++kn)
// 		of2 << out[kn]->getPoint() << std::endl;
	  
// 	      double maxd_init, avd_init;
// 	      int num2_init;
// 	      regions[kr]->getAccuracy(maxd_init, avd_init, num2_init);
// 	      int num = regions[kr]->numPoints();

// 	      if (num2 < num/2 || avd > approx_tol_)
// 		{
// 		  if (ki < appn)
// 		    {
// 		      cand_ix.erase(cand_ix.begin() + ki);
// 		      cand_score.erase(cand_score.begin() + ki);

// 		      for ()
// 			{
// 			  regions.erase(regions.begin()+ki);
// 			  points.erase(points.begin()+ki);
// 			  nmbpts.erase(nmbpts.begin()+ki);
// 			}
// 		    }
// 		  else
// 		    {
// 		    }
// 		}
// 	      else
// 		{
// 		  all_in.push_back(in);
// 		  all_maxd = std::max(all_maxd, maxd);
// 		  all_avd += num*avd;
// 		  all_inside += num2;
// 		  ++ki;
// 		}
// 	  int stop_break = 1;
// 	}
//       if (cand_ix.size() <= 1)
// 	break;
//       candsize = cand_ix.size();
//     }
//   if (cand_ix.size() == candsize)
//     {
//       merged_surf =
// 	shared_ptr<HedgeSurface>(new HedgeSurface(surf, regions));
//       for (size_t kj=0; kj<regions.size(); ++kj)
// 	regions[kj]->setHedge(merged_surf.get());
//     }

//   return merged_surf;
// }

//===========================================================================
shared_ptr<HedgeSurface> RevEng::doMerge(vector<size_t>& cand_ix,
					 vector<double>& cand_score,
					 ClassType type)
//===========================================================================
{
  size_t candsize = cand_ix.size();
  shared_ptr<HedgeSurface> merged_surf;
  if (cand_ix.size() <= 1)
    return merged_surf;
  double delscore = 2.0*(cand_score[1]-cand_score[0]);

  vector<size_t> select_ix;
  select_ix.push_back(0);
  for (size_t ki=1; ki<cand_ix.size(); ++ki)
    {
      if (cand_score[ki]-cand_score[ki-1] > delscore)
	break;
      select_ix.push_back(ki);
    }
  size_t init_select = select_ix.size();

  shared_ptr<ParamSurface> prev_surf, surf;
  vector<vector<RevEngPoint*> > all_in, all_in2;
  BoundingBox bbox(3);
  double avd_all = 0.0, maxd_all = 0.0;
  int num_in_all = 0, num_all = 0, num_all2 = 0;
  vector<RevEngRegion*> regions, regions2;
  double tolfac = 2.0;
  vector<vector<pair<double,double> > > distang;
  vector<vector<double> > parvals;
  for (size_t kh=0; kh<init_select; ++kh)
    {
      surf = approxMergeSet(cand_ix, select_ix, type);
      if (!surf.get())
	break;

      std::ofstream of1("merge_sf.g2");
      std::ofstream of2("in_out_merge.g2");
      surf->writeStandardHeader(of1);
      surf->write(of1);
      
      // Test accuracy
      all_in.clear();
      regions.clear();
      all_in2.clear();
      regions2.clear();
      distang.clear();
      parvals.clear();
      num_in_all = num_all = num_all2 = 0;
       vector<size_t> select_ix2;
      BoundingBox bb(3);
      size_t ki, kj;
      for (ki=0; ki<cand_ix.size(); ++ki)
	{
	  HedgeSurface* hsurf = surfaces_[cand_ix[ki]].get();
	  vector<RevEngRegion*> reg = hsurf->getRegions();
	  vector<RevEngPoint*> in, out;
	  double avd_sf = 0.0, maxd_sf = 0.0;
	  int num_in_sf = 0, num_sf = hsurf->numPoints();
	  vector<vector<pair<double,double> > > curr_distang(reg.size());
	  vector<vector<double> > curr_parvals(reg.size());
	  for (kj=0; kj<reg.size(); ++kj)
	    {
 	      double maxd, avd;
	      int num_in;
	      int num = reg[kj]->numPoints();
	      RevEngUtils::distToSurf(reg[kj]->pointsBegin(),
				      reg[kj]->pointsEnd(),
				      surf, approx_tol_, maxd, avd, num_in, 
				      in, out, curr_parvals[kj], curr_distang[kj]);
	      maxd_sf = std::max(maxd_sf, maxd);
	      avd_sf += num*avd/(double)num_sf;
	      num_in_sf += num_in;
	    }
	    of2 << "400 1 0 4 155 50 50 255" << std::endl;
	    of2 << in.size() << std::endl;
	    for (size_t kn=0; kn<in.size(); ++kn)
	      of2 << in[kn]->getPoint() << std::endl;
	    of2 << "400 1 0 4 50 155 50 255" << std::endl;
	    of2 << out.size() << std::endl;
	    for (size_t kn=0; kn<out.size(); ++kn)
	      of2 << out[kn]->getPoint() << std::endl;

	    if (num_in_sf > num_sf/2 && avd_sf < approx_tol_)
	      {
		all_in.push_back(in);
		select_ix2.push_back(ki);
		bb.addUnionWith(hsurf->regionsBox());
		avd_all += num_sf*avd_sf;
		maxd_all = std::max(maxd_all, maxd_sf);
		num_in_all += num_in_sf;
		num_all += num_sf;
		regions.insert(regions.end(), reg.begin(), reg.end());
		distang.insert(distang.end(), curr_distang.begin(), curr_distang.end());
		parvals.insert(parvals.end(), curr_parvals.begin(), curr_parvals.end());
	      }
	    else
	      {
		if (num_in_sf > num_sf/2 && avd_sf < tolfac*approx_tol_)
		  {
		    all_in2.push_back(in);
		    regions2.insert(regions2.end(), reg.begin(), reg.end());
		    num_all2 += num_sf;
		  }
		if (ki < init_select)
		  init_select--;
		else
		  break;
	      }
	}
      
      avd_all /= (double)num_all;
      bbox = bb;
      if (select_ix == select_ix2)
	break;
      select_ix = select_ix2;
      
      if (ki < cand_ix.size())
	break;
      prev_surf = surf;
    }
  
  for (size_t ki=0; ki<cand_ix.size(); )
    {
      vector<size_t>::iterator found =
	std::find(select_ix.begin(), select_ix.end(), ki);
      if (found == select_ix.end())
	cand_ix.erase(cand_ix.begin()+ki);
      else
	++ki;
    }

  if (cand_ix.size() <= 1 || all_in.size() + all_in2.size() <= 1)
    return merged_surf;

  
  vector<RevEngRegion*> regions3 = regions;
  if (all_in2.size() > 0)
    {
      all_in.insert(all_in.end(), all_in2.begin(), all_in2.end());
      regions3.insert(regions3.end(), regions2.begin(), regions2.end());
      num_all2 += num_all;
    }
  
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > points(all_in.size());
  vector<int> nmbpts(all_in.size());
  for (size_t kj=0; kj<all_in.size(); ++kj)
    {
      points[kj] = std::make_pair(all_in[kj].begin(),
				  all_in[kj].end());
      nmbpts[kj] = (int)all_in[kj].size();
    }

  BoundingBox bbox2(3);
  for (size_t kj=0; kj<regions3.size(); ++kj)
    {
      bbox2.addUnionWith(regions3[kj]->boundingBox());
    }

  shared_ptr<ParamSurface> surf2;
  if (type == Class_Plane)
    {
      surf2 = doMergePlanes(points, bbox2, nmbpts);
    }
  else if (type == Class_Cylinder)
    {
      surf2 = doMergeCylinders(points, bbox2, nmbpts);
    }
  else if (type == Class_Torus)
    {
      surf2 = doMergeTorus(points, bbox2, nmbpts);
    }
  
  std::ofstream of3("merge_sf2.g2");
  std::ofstream of4("in_out_merge2.g2");
  surf2->writeStandardHeader(of3);
  surf2->write(of3);
  vector<RevEngPoint*> in2, out2;
  vector<vector<pair<double,double> > > distang2(regions3.size());
  vector<vector<double> > parvals2(regions3.size());
  double avd_all2 = 0.0, maxd_all2 = 0.0;
  int num_in_all2 = 0;
  for (size_t kj=0; kj<regions3.size(); ++kj)
    {
      double maxd, avd;
      int num_in;
      RevEngUtils::distToSurf(regions3[kj]->pointsBegin(),
			      regions3[kj]->pointsEnd(),
			      surf2, approx_tol_, maxd, avd, num_in, 
			      in2, out2, parvals2[kj], distang2[kj]);
      maxd_all2 = std::max(maxd_all2, maxd);
      avd_all2 += regions3[kj]->numPoints()*avd/(double)num_all2;
      num_in_all2 += num_in;
    }
  of4 << "400 1 0 4 155 50 50 255" << std::endl;
  of4 << in2.size() << std::endl;
  for (size_t kn=0; kn<in2.size(); ++kn)
    of4 << in2[kn]->getPoint() << std::endl;
  of4 << "400 1 0 4 50 155 50 255" << std::endl;
  of4 << out2.size() << std::endl;
  for (size_t kn=0; kn<out2.size(); ++kn)
    of4 << out2[kn]->getPoint() << std::endl;

  if (surf2.get() && num_in_all2 >= num_in_all && avd_all2 < avd_all &&
      num_in_all2 > num_all2/2 && avd_all2 < approx_tol_)
    {
      merged_surf =
	shared_ptr<HedgeSurface>(new HedgeSurface(surf2, regions3));
      for (size_t kj=0; kj<regions3.size(); ++kj)
	{
	  regions3[kj]->setHedge(merged_surf.get());
	  int numpt = regions3[kj]->numPoints();
	  for (int ka=0; ka<numpt; ++ka)
	    {
	      RevEngPoint *pt = regions3[kj]->getPoint(ka);
	      pt->setPar(Vector2D(parvals2[kj][2*ka],parvals2[kj][2*ka+1]));
	      pt->setSurfaceDist(distang2[kj][ka].first, distang2[kj][ka].second);
	    }
	}
    }
  else if (surf.get() && num_in_all > num_all/2 && avd_all < approx_tol_)
    {
      merged_surf =
	shared_ptr<HedgeSurface>(new HedgeSurface(surf, regions));
      for (size_t kj=0; kj<regions.size(); ++kj)
	{
	  regions[kj]->setHedge(merged_surf.get());
	  int numpt = regions[kj]->numPoints();
	  for (int ka=0; ka<numpt; ++ka)
	    {
	      RevEngPoint *pt = regions[kj]->getPoint(ka);
	      pt->setPar(Vector2D(parvals[kj][2*ka],parvals[kj][2*ka+1]));
	      pt->setSurfaceDist(distang[kj][ka].first, distang[kj][ka].second);
	    }
	}
    }
  return merged_surf;
}
 


//===========================================================================
  shared_ptr<ParamSurface> RevEng::approxMergeSet(vector<size_t>& cand_ix,
						  vector<size_t>& select_ix,
						  ClassType type)
						  
//===========================================================================
  {
  // Collect point clouds
  std::ofstream ofp("all_merge_points.g2");
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > points;
  BoundingBox bbox(3);
  vector<int> nmbpts;
  for (size_t ki=0; ki<select_ix.size(); ++ki)
    {
      HedgeSurface* surf = surfaces_[cand_ix[select_ix[ki]]].get();
      vector<RevEngRegion*> reg = surf->getRegions();
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
  
  shared_ptr<ParamSurface> surf;
  if (type == Class_Plane)
    {
      surf = doMergePlanes(points, bbox, nmbpts);
    }
  else if (type == Class_Cylinder)
    {
      surf = doMergeCylinders(points, bbox, nmbpts);
    }
  else if (type == Class_Torus)
    {
      surf = doMergeTorus(points, bbox, nmbpts);
    }
  return surf;
  }
  
//===========================================================================
shared_ptr<HedgeSurface> RevEng::doMergeSpline(shared_ptr<HedgeSurface>& surf1,
					       DirectionCone& cone1,
					       shared_ptr<HedgeSurface>& surf2,
					       DirectionCone& cone2)
//===========================================================================
{
  ClassType type1 = Class_Unknown, type2 = Class_Unknown;

  // Parameterize
  // Collect all points
  vector<RevEngPoint*> all_pts;
  vector<RevEngRegion*> regions;
  int dim = surf1->surface()->dimension();
  BoundingBox bbox(dim);
  int num1 = surf1->numRegions();
  for (int ka=0; ka<num1; ++ka)
    {
      RevEngRegion *reg = surf1->getRegion(ka);
      regions.push_back(reg);
      bbox.addUnionWith(reg->getBbox());
      vector<RevEngPoint*> pts = reg->getPoints();
      all_pts.insert(all_pts.end(), pts.begin(), pts.end());
      if (reg->hasPrimary())
	{
	  ClassType ctype = reg->getPrimary()->instanceType();
	  if (type1 == Class_Unknown)
	    type1 = ctype;
	  else if (type1 != ctype)
	    type1 = Class_Unknown;  // Can be a problem with several regions
	}
    }
  int num2 = surf2->numRegions();
  for (int ka=0; ka<num2; ++ka)
    {
      RevEngRegion *reg = surf2->getRegion(ka);
      regions.push_back(reg);
      bbox.addUnionWith(reg->getBbox());
      vector<RevEngPoint*> pts = reg->getPoints();
      all_pts.insert(all_pts.end(), pts.begin(), pts.end());
      if (reg->hasPrimary())
	{
	  ClassType ctype = reg->getPrimary()->instanceType();
	  if (type2 == Class_Unknown)
	    type2 = ctype;
	  else if (type1 != ctype)
	    type2 = Class_Unknown;
	}
    }
  
  double avd_all = 0.0, maxd_all = 0.0;
  int num_in_all = 0;
  for (size_t ki=0; ki<regions.size(); ++ki)
    {
      double avdr, maxdr;
      int num_inr;
      regions[ki]->getAccuracy(maxdr, avdr, num_inr);
      maxd_all = std::max(maxd_all, maxdr);
      double wgt = (double)regions[ki]->numPoints()/(double)all_pts.size();
      avd_all += wgt*avdr;
      num_in_all += num_inr;
    }
  
  // Define base surface
  ClassType type = Class_Unknown;
  if (type1 == type2 && type1 != Class_Unknown)
    {
      type = type1;
    }
  else if (type1 != Class_Unknown && type2 != Class_Unknown)
    {
    }
  else if (type1 != Class_Unknown)
    {
      type = type1;
    }
  else if (type2 != Class_Unknown)
    {
      type = type2;
    }

  shared_ptr<ParamSurface> projectsf;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > all_pts2(1);
  all_pts2[0] = std::make_pair(all_pts.begin(),all_pts.end());
  vector<int> nmbpts(1, (int)all_pts.size());
  if (type == Class_Plane)
    {
      projectsf = doMergePlanes(all_pts2, bbox, nmbpts, false);
    }
  else if (type == Class_Cylinder)
    {
      projectsf = doMergeCylinders(all_pts2, bbox, nmbpts, false);
    }
  else if (type == Class_Torus)
    {
      projectsf = doMergeTorus(all_pts2, bbox, nmbpts);
    }

  if (projectsf.get())
    {
      std::ofstream of1("projectsf.g2");
      projectsf->writeStandardHeader(of1);
      projectsf->write(of1);
    }
  
  // Parameterize on surface
  vector<double> data;
  vector<double> param;
  int inner_u = 0, inner_v = 0;
  bool close1 = false, close2 = false;
  if (projectsf.get())
    RevEngUtils::parameterizeOnPrimary(all_pts, projectsf, data, param,
				       inner_u, inner_v, close1, close2);
  else
    {
     double lambda[3];
     double eigenvec[3][3];
     RevEngUtils::principalAnalysis(all_pts, lambda, eigenvec);
     Point xAxis(eigenvec[0][0], eigenvec[0][1], eigenvec[0][2]);
     Point yAxis(eigenvec[1][0], eigenvec[1][1], eigenvec[1][2]);
     RevEngUtils::parameterizeWithPlane(all_pts, bbox, xAxis, yAxis,
					data, param);
    }
  
  // Approximate
  double maxd, avd;
  int num_out;
  int order = 4;
  int ncoef1 = order + inner_u;
  int ncoef2 = order + inner_v;
  int max_iter = 6;
  double belt = 0.01;
  shared_ptr<HedgeSurface> merged;
  shared_ptr<SplineSurface> surf;
  try {
    surf =
    RevEngUtils::surfApprox(data, dim, param, order, order, ncoef1, ncoef2,
			    close1, close2, max_iter, approx_tol_, maxd, avd,
			    num_out, belt);
  }
  catch (...)
    {
      return merged;
    }
  int num_in = (int)all_pts.size() - num_out;

  if (surf.get())
    {
      std::ofstream of2("mergedsf.g2");
      surf->writeStandardHeader(of2);
      surf->write(of2);
      for (size_t ki=0; ki<regions.size(); ++ki)
	regions[ki]->writeRegionInfo(of2);
    }

  BoundingBox bb1 = surf1->boundingBox();
  BoundingBox bb2 = surf2->boundingBox();
  bb1.addUnionWith(bb2);
  BoundingBox bb3 = surf->boundingBox();
  double diag1 = bb1.low().dist(bb1.high());
  double diag3 = bb3.low().dist(bb3.high());
  double frac1 = (double)num_in_all/(double)all_pts.size();
  double frac2 = (double)num_in/(double)all_pts.size();
  if (frac2 > 2.0*frac1/3.0 && frac2 > 0.5 && avd <= approx_tol_ &&
      diag3 < 1.5*diag1)
    {
      merged = shared_ptr<HedgeSurface>(new HedgeSurface(surf, regions));
      for (size_t kj=0; kj<regions.size(); ++kj)
	regions[kj]->setHedge(merged.get());
    }
  return merged;
}

// //===========================================================================
// void RevEng::recognizePlanes()
// //===========================================================================
// {
//   std::ofstream planeout("plane.g2");
//   double anglim = 0.1;
//   size_t nmbsfs = surfaces_.size();
//   for (size_t ki=0; ki<regions_.size(); ++ki)
//     {
//       // Check type
//       bool planar = regions_[ki]->planartype();
//       DirectionCone normalcone = regions_[ki]->getNormalCone();
//       if (!planar /*&& (normalcone.angle() > anglim ||
// 		    normalcone.greaterThanPi())*/)
// 	continue;  // Not a probable plane

//       // Try to fit the point cloud with a plane
//       vector<shared_ptr<HedgeSurface> > plane_sfs;
//       vector<HedgeSurface*> prev_surfs;
//       bool found = regions_[ki]->extractPlane(approx_tol_, min_point_region_,
// 					      10.0*anglim_, plane_sfs, prev_surfs, planeout);
//       for (size_t kr=0; kr<prev_surfs.size(); ++kr)
// 	{
// 	  size_t kj;
// 	  for (kj=0; kj<surfaces_.size(); ++kj)
// 	    if (surfaces_[kj].get() == prev_surfs[kr])
// 	      break;
// 	  if (kj < surfaces_.size())
// 	    surfaces_.erase(surfaces_.begin()+kj);
// 	}
//       if (plane_sfs.size() > 0)
// 	{
// 	  surfaces_.insert(surfaces_.end(), plane_sfs.begin(), plane_sfs.end());

// 	  vector<RevEngRegion*> grown_regions;
// 	  int min_nmb = 5*min_point_region_;  // Should be set from distribution of how many
// 	  // points the regions have
// 	  vector<HedgeSurface*> adj_surfs;
// 	  regions_[ki]->growWithSurf(min_nmb, approx_tol_, grown_regions, adj_surfs);
// 	  if (grown_regions.size() > 0)
// 	    {
// 	      for (size_t kr=0; kr<grown_regions.size(); ++kr)
// 		{
// 		  size_t kj;
// 		  for (kj=0; kj<regions_.size(); )
// 		    {
// 		      if (kj == ki)
// 			{
// 			  ++kj;
// 			  continue;
// 			}

// 		      if (grown_regions[kr] == regions_[kj].get())
// 			{
// 			  regions_.erase(regions_.begin()+kj);
// 			  if (kj < ki)
// 			    --ki;
// 			}
// 		      else
// 			++kj;
// 		    }
// 		}
// 	      for (size_t kr=0; kr<adj_surfs.size(); ++kr)
// 		{
// 		  size_t kj;
// 		  for (kj=0; kj<surfaces_.size(); ++kj)
// 		    if (surfaces_[kj].get() == adj_surfs[kr])
// 		      break;
// 		  if (kj < surfaces_.size())
// 		    surfaces_.erase(surfaces_.begin()+kj);
// 		}
// 	    }
// 	}

//       for (size_t kj=0; kj<regions_.size(); ++kj)
// 	regions_[kj]->setVisited(false);


//       // if not planar
//       // continue;

//       // apply method for recognition of plane
//       // Should RevEngPoint be given as input or should the method take a
//       // vector of Points as input to enforce independence on a triangulation?
//       // The class HedgeSurface must lie in compositemodel due to the
//       // connection to RevEngRegion
//       // Should the computation take place in HedgeSurface? Or in RevEngRegion?
//       // Or in this class?
      
//       // some points may be disassembled

//       // check if the number of points in the plane is large enough
//       // if not, should the region be removed?
//       int stop_break = 1;
//     }
  
//   std::ofstream ofp("planes2.g2");
//   for (size_t ki=0; ki<surfaces_.size(); ++ki)
//     {
//       shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
//       surf->writeStandardHeader(ofp);
//       surf->write(ofp);
//     }

//   if (regions_.size() > 0)
//     {
//       std::ofstream ofr("regions3.g2");
//       for (size_t kr=0; kr<regions_.size(); ++kr)
// 	{
// 	  // BoundingBox bbox = regions_[kr]->boundingBox();
// 	  // if (bbox.low().dist(bbox.high()) < 0.1)
// 	  //   std::cout << "Small bounding box" << std::endl;
// 	  if (regions_[kr]->numPoints() < 5)
// 	    continue;
// 	  ofr << "400 1 0 0" << std::endl;
// 	  int nmb = regions_[kr]->numPoints();
// 	  ofr << nmb << std::endl;
// 	  for (int ki=0; ki<nmb; ++ki)
// 	    {
// 	      ofr << regions_[kr]->getPoint(ki)->getPoint() << std::endl;
// 	    }
// 	}
//     }

//   if (surfaces_.size() - nmbsfs > 1)
//     mergePlanes(nmbsfs, surfaces_.size());

//   std::ofstream of("surfaces0.g2");
//   for (size_t ki=0; ki<surfaces_.size(); ++ki)
//     {
//       shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
//       surf->writeStandardHeader(of);
//       surf->write(of);
//     }
// }

// //===========================================================================
// void RevEng::recognizeCylinders()
// //===========================================================================
// {
//   std::ofstream cylout("cylinder.g2");
//   int nmbsfs = surfaces_.size();
//   for (size_t ki=0; ki<regions_.size(); ++ki)
//     {
//       // Check type
//       bool cyl = regions_[ki]->cylindertype();
//       // if (!cyl)
//       // 	continue;
      
//       // So far, try to regognize cylinders
//       vector<shared_ptr<HedgeSurface> > cyl_sfs;
//       vector<HedgeSurface*> prev_surfs;
//       bool found = regions_[ki]->extractCylinder(approx_tol_, min_point_region_,
// 						 10.0*anglim_,
// 						 mean_edge_len_, cyl_sfs, prev_surfs,
// 						 cylout);
//       for (size_t kr=0; kr<prev_surfs.size(); ++kr)
// 	{
// 	  size_t kj;
// 	  for (kj=0; kj<surfaces_.size(); ++kj)
// 	    if (surfaces_[kj].get() == prev_surfs[kr])
// 	      break;
// 	  if (kj < surfaces_.size())
// 	    {
// 	      surfaces_.erase(surfaces_.begin()+kj);
// 	      nmbsfs--;
// 	    }
// 	}
//       if (cyl_sfs.size() > 0)
// 	{
// 	  surfaces_.insert(surfaces_.end(), cyl_sfs.begin(), cyl_sfs.end());
      
// 	  vector<RevEngRegion*> grown_regions;
// 	  int min_nmb = 10*min_point_region_;  // Should be set from distribution of how many
// 	  // points the regions have
// 	  vector<HedgeSurface*> adj_surfs;
// 	  regions_[ki]->growWithSurf(min_nmb, approx_tol_, grown_regions, adj_surfs);
// 	  if (grown_regions.size() > 0)
// 	    {
// 	      for (size_t kr=0; kr<grown_regions.size(); ++kr)
// 		{
// 		  size_t kj;
// 		  for (kj=0; kj<regions_.size(); )
// 		    {
// 		      if (kj == ki)
// 			{
// 			  ++kj;
// 			  continue;
// 			}

// 		      if (grown_regions[kr] == regions_[kj].get())
// 			{
// 			  regions_.erase(regions_.begin()+kj);
// 			  if (kj < ki)
// 			    ki--;
// 			}
// 		      else
// 			++kj;
// 		    }
// 		}
// 	      for (size_t kr=0; kr<adj_surfs.size(); ++kr)
// 		{
// 		  size_t kj;
// 		  for (kj=0; kj<surfaces_.size(); ++kj)
// 		    if (surfaces_[kj].get() == adj_surfs[kr])
// 		      break;
// 		  if (kj < surfaces_.size())
// 		    {
// 		      surfaces_.erase(surfaces_.begin()+kj);
// 		      if (kj < nmbsfs)
// 			nmbsfs--;
// 		    }
// 		}
// 	    }
// 	}

//       for (size_t kj=0; kj<regions_.size(); ++kj)
// 	regions_[kj]->setVisited(false);
      
//       // if not planar
//       // continue;

//       // apply method for recognition of plane
//       // Should RevEngPoint be given as input or should the method take a
//       // vector of Points as input to enforce independence on a triangulation?
//       // The class HedgeSurface must lie in compositemodel due to the
//       // connection to RevEngRegion
//       // Should the computation take place in HedgeSurface? Or in RevEngRegion?
//       // Or in this class?
      
//       // some points may be disassembled

//       // check if the number of points in the plane is large enough
//       // if not, should the region be removed?
//     }

  
//   std::ofstream ofp("cylinders2.g2");
//   for (size_t ki=0; ki<surfaces_.size(); ++ki)
//     {
//       shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
//       surf->writeStandardHeader(ofp);
//       surf->write(ofp);
//     }

//   if (regions_.size() > 0)
//     {
//       std::ofstream ofr("regions4.g2");
//       for (size_t kr=0; kr<regions_.size(); ++kr)
// 	{
// 	  // BoundingBox bbox = regions_[kr]->boundingBox();
// 	  // if (bbox.low().dist(bbox.high()) < 0.1)
// 	  //   std::cout << "Small bounding box" << std::endl;
// 	  if (regions_[kr]->numPoints() < 5)
// 	    continue;
// 	  ofr << "400 1 0 0" << std::endl;
// 	  int nmb = regions_[kr]->numPoints();
// 	  ofr << nmb << std::endl;
// 	  for (int ki=0; ki<nmb; ++ki)
// 	    {
// 	      ofr << regions_[kr]->getPoint(ki)->getPoint() << std::endl;
// 	    }
// 	}
//     }

//   // Now the surfaces between nmbsfs and surfaces_.size() are cylinder surfaces
//   // and can be merged if they represent the same cylinder
//   if (surfaces_.size() + nmbsfs > 1)
//     mergeCylinders(nmbsfs, surfaces_.size());
//   //    mergeCylinders(0, surfaces_.size());

//   std::ofstream of("surfaces.g2");
//   for (size_t ki=0; ki<surfaces_.size(); ++ki)
//     {
//       shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
//       surf->writeStandardHeader(of);
//       surf->write(of);
//     }
// }

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
					       vector<int>& nmbpts,
					       bool set_bound)
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
	  Point curr = (*it)->getMongeNormal();
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
  if (set_bound)
    {
      double len = low.dist(high);
      surf->setParameterBounds(-0.5*len, -0.5*len, 0.5*len, 0.5*len);
    }

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
						  vector<int>& nmbpts,
						  bool set_bound)
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
  if (set_bound)
    {
      double len = low.dist(high);
      surf->setParamBoundsV(-len, len);
    }

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
	  Point norm = (*it)->getMongeNormal();
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
  vector<int> num_pts;
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
	      num_pts[kr] += sfprop[ki].num_points_;
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
	  num_pts.push_back(sfprop[ki].num_points_);
	}
    }

  // Sort according to number of points assiciated to the surface
  for (size_t ki=0; ki<num_pts.size(); ++ki)
    for (size_t kj=ki+1; kj<num_pts.size(); ++kj)
      {
	if (num_pts[kj] > num_pts[ki])
	  {
	    std::swap(num_pts[ki], num_pts[kj]);
	    std::swap(axis_cone[ki], axis_cone[kj]);
	    std::swap(group_ixs[ki], group_ixs[kj]);
	  }
      }

  for (size_t ki=0; ki<num_pts.size(); ++ki)
    {
      // Identify cylinders with the "same" axis including location
      if (group_ixs[ki].size() < 2)
	continue;
      for (size_t kj=0; kj<group_ixs[ki].size(); ++kj)
	{
	  if (sfprop[group_ixs[ki][kj]].type_ != Class_Cylinder &&
	      sfprop[group_ixs[ki][kj]].prev_type_  != Class_Cylinder)
	    continue;
	  Point loc1 = sfprop[group_ixs[ki][kj]].loc_;
	  Point vec1 = sfprop[group_ixs[ki][kj]].dir_;
	  double rad1 = sfprop[group_ixs[ki][kj]].rad1_;
	  vector<int> adapt_ix;
	  adapt_ix.push_back(sfprop[group_ixs[ki][kj]].sfix_);
	  for (size_t kr=kj+1; kr<group_ixs[ki].size(); ++kr)
	    {
	      if (sfprop[group_ixs[ki][kr]].type_ != Class_Cylinder &&
		  sfprop[group_ixs[ki][kr]].prev_type_  != Class_Cylinder)
		continue;
	      Point loc2 = sfprop[group_ixs[ki][kr]].loc_;
	      Point loc2_0 = loc1 + ((loc2-loc1)*vec1)*vec1;
	      double rad2 = sfprop[group_ixs[ki][kr]].rad1_;
	      double pdist = loc2.dist(loc2_0);
	      double dlim = std::max(0.05*std::max(rad1,rad2), approx_tol_);
	      if (pdist < dlim)
		{
		  adapt_ix.push_back(sfprop[group_ixs[ki][kr]].sfix_);
		}
	    }

	  if (adapt_ix.size() > 1)
	    {
	      // Try simultanous fitting of cylinders
	      cylinderFit(adapt_ix, axis_cone[ki].centre());
	      int stop_break = 1;
	    }
	}
    }
  int stop_break = 1;
}

//===========================================================================
void RevEng::cylinderFit(vector<int>& sf_ix, Point normal)
//===========================================================================
{
  // Collect data points
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > points;
  BoundingBox bbox(3);
  for (size_t ki=0; ki<sf_ix.size(); ++ki)
    {
      HedgeSurface* surf = surfaces_[sf_ix[ki]].get();
      vector<RevEngRegion*> reg = surf->getRegions();
      for (size_t kj=0; kj<reg.size(); ++kj)
	{
	  points.push_back(std::make_pair(reg[kj]->pointsBegin(),
					  reg[kj]->pointsEnd()));
	  bbox.addUnionWith(reg[kj]->boundingBox());
	}
    }

  Point axis, Cx, Cy;
  RevEngUtils::computeAxis(points, axis, Cx, Cy);

  Point low = bbox.low();
  Point high = bbox.high();
  Point pos;
  double radius;
  RevEngUtils::computeCylPosRadius(points, low, high, axis, Cx, Cy, pos,
				   radius);

  std::ofstream of("cylinderfit.g2");
  
  for (size_t ki=0; ki<sf_ix.size(); ++ki)
    {
      HedgeSurface* surf = surfaces_[sf_ix[ki]].get();
      vector<RevEngRegion*> reg = surf->getRegions();
      vector<pair<vector<RevEngPoint*>::iterator,
		  vector<RevEngPoint*>::iterator> > points0;
      BoundingBox bbox0(3);
      for (size_t kj=0; kj<reg.size(); ++kj)
	{
	  points0.push_back(std::make_pair(reg[kj]->pointsBegin(),
					  reg[kj]->pointsEnd()));
	  bbox0.addUnionWith(reg[kj]->boundingBox());
	}
      Point low0 = bbox0.low();
      Point high0 = bbox0.high();
      Point pos0;
      double radius0;
      RevEngUtils::computeCylPosRadius(points0, low0, high0, axis, Cx, Cy, 
				      pos0, radius0);
      double radius2 = computeCylRadius(points0, pos, Cx, Cy);
      
      shared_ptr<Cylinder> cyl(new Cylinder(radius0, pos0, axis, Cy));
      cyl->writeStandardHeader(of);
      cyl->write(of);
      shared_ptr<Cylinder> cyl2(new Cylinder(radius2, pos, axis, Cy));
      cyl2->writeStandardHeader(of);
      cyl2->write(of);

      shared_ptr<Cylinder> cyl3(new Cylinder(radius0, pos, axis, Cy));
      cyl2->writeStandardHeader(of);
      cyl2->write(of);

      for (size_t kj=0; kj<reg.size(); ++kj)
	{
	  double maxd, avd;
	  int num_in;
	  vector<RevEngPoint*> in, out;
	  vector<pair<double,double> > distang;
	  vector<double> parvals;
	  RevEngUtils::distToSurf(reg[kj]->pointsBegin(), reg[kj]->pointsEnd(),
				  cyl, approx_tol_, maxd, avd, num_in, in, out,
				  parvals, distang);
	  of << "400 1 0 4 155 50 50 255" << std::endl;
	  of << in.size() << std::endl;
	  for (size_t kr=0; kr<in.size(); ++kr)
	    of << in[kr]->getPoint() << std::endl;
	  of << "400 1 0 4 50 155 50 255" << std::endl;
	  of << out.size() << std::endl;
	  for (size_t kr=0; kr<out.size(); ++kr)
	    of << out[kr]->getPoint() << std::endl;
	  
	  double maxd2, avd2;
	  int num_in2;
	  vector<RevEngPoint*> in2, out2;
	  vector<pair<double,double> > distang2;
	  vector<double> parvals2;
	  RevEngUtils::distToSurf(reg[kj]->pointsBegin(), reg[kj]->pointsEnd(),
				  cyl2, approx_tol_, maxd2, avd2, num_in2, in2,
				  out2, parvals2, distang2);
	  of << "400 1 0 4 155 50 50 255" << std::endl;
	  of << in2.size() << std::endl;
	  for (size_t kr=0; kr<in2.size(); ++kr)
	    of << in2[kr]->getPoint() << std::endl;
	  of << "400 1 0 4 50 155 50 255" << std::endl;
	  of << out2.size() << std::endl;
	  for (size_t kr=0; kr<out2.size(); ++kr)
	    of << out2[kr]->getPoint() << std::endl;
	  
	  double maxd3, avd3;
	  int num_in3;
	  vector<RevEngPoint*> in3, out3;
	  vector<pair<double,double> > distang3;
	  vector<double> parvals3;
	  RevEngUtils::distToSurf(reg[kj]->pointsBegin(), reg[kj]->pointsEnd(),
				  cyl3, approx_tol_, maxd3, avd3, num_in3, in3,
				  out3, parvals3, distang3);
	  of << "400 1 0 4 155 50 50 255" << std::endl;
	  of << in3.size() << std::endl;
	  for (size_t kr=0; kr<in3.size(); ++kr)
	    of << in3[kr]->getPoint() << std::endl;
	  of << "400 1 0 4 50 155 50 255" << std::endl;
	  of << out3.size() << std::endl;
	  for (size_t kr=0; kr<out3.size(); ++kr)
	    of << out3[kr]->getPoint() << std::endl;
	  
	  int stop_break = 1;
	}
    }
}

//===========================================================================
double RevEng::computeCylRadius(vector<pair<vector<RevEngPoint*>::iterator,
				vector<RevEngPoint*>::iterator> >& points,
				Point mid, Point vec1, Point vec2)
//===========================================================================
{
  int nmb = 0;
  for (size_t ki=0; ki<points.size(); ++ki)
    nmb += (int)(points[ki].second - points[ki].first);

  double wgt = 1.0/(double)nmb;

  // Make tranformation matrix
  Matrix3D mat1, mat2, rotmat;
  Vector3D vec1_2(vec1[0], vec1[1], vec1[2]);
  Vector3D vec2_2(vec2[0], vec2[1], vec2[2]);
  Vector3D xaxis(1, 0, 0);
  Vector3D yaxis(0, 1, 0);
  mat1.setToRotation(vec1_2, xaxis);
  Vector3D v1 = mat1*vec1_2;
  Vector3D vec2_3 = mat1*vec2_2;
  mat2.setToRotation(vec2_3, yaxis);
  Vector3D v2 = mat2*vec2_3;
  rotmat = mat2*mat1;

  // Rotate points and add to circle decription
  double r2 = 0.0;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      vector<RevEngPoint*>::iterator start = points[ki].first;
      vector<RevEngPoint*>::iterator end = points[ki].second;
      for (auto it=start; it!=end; ++it)
	{
	  RevEngPoint *pt = *it;
	  Vector3D pnt = pt->getPoint();
	  Vector3D rpnt = rotmat*pnt;
	  double curr = rpnt[0]*rpnt[0] + rpnt[1]*rpnt[1];
	  r2 += wgt*curr;
	}
    }

  double radius = sqrt(r2);
  return radius;
}

//===========================================================================
void RevEng::collectAxis(vector<SurfaceProperties>& sfprop)
//===========================================================================
{
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      int code;
      ClassType type = surfaces_[ki]->instanceType(code);
      ClassType type2 = Class_Unknown;
      shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
      shared_ptr<ElementarySurface> elemsf =
	dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf);

      int nreg = surfaces_[ki]->numRegions();
      int num_pts = 0;
      shared_ptr<ParamSurface> primary;
      int num_pt_primary = 0;
      for (int ka=0; ka<nreg; ++ka)
	{
	  RevEngRegion *reg = surfaces_[ki]->getRegion(ka);
	  int num = reg->numPoints();
	  num_pts += num;
	  if (reg->hasPrimary() && num > num_pt_primary)
	    {
	      double maxdp, avdp;
	      int num_inp;
	      reg->getPrimaryInfo(maxdp, avdp, num_inp);
	      if (num_inp > num/2 && avdp < approx_tol_)
		{
		  primary = reg->getPrimary();
		  num_pt_primary = num;
		}
	    }
	}
      if (primary.get())
	type2 = primary->instanceType();
      
      if (!elemsf.get())
	{
	  if (primary.get())
	    elemsf =
	      dynamic_pointer_cast<ElementarySurface,ParamSurface>(primary);
	}
      if (!elemsf.get())
	continue;
      
      Point loc, dir;
      double rad1, rad2;
      loc = elemsf->location();
      dir = elemsf->direction();
      rad1 = elemsf->radius(0.0, 0.0);   // Not good enough for cones
      rad2 = elemsf->radius2(0.0, 0.0);   // Not good enough for cones
      SurfaceProperties currprop(ki, type, num_pts, dir, loc, type2,
				 rad1, rad2);
      sfprop.push_back(currprop);
    }
}


//===========================================================================
void RevEng::trimSurfaces()
//===========================================================================
{
  std::ofstream of1("surfbd.g2");
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      // Restrict unbounded surfaces
      surfaces_[ki]->ensureSurfaceBounded();
      surfaces_[ki]->surface()->writeStandardHeader(of1);
      surfaces_[ki]->surface()->write(of1);
    }

  vector<shared_ptr<ParamSurface> > sub_sfs;
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
      vector<shared_ptr<ParamSurface> > splitsfs =
	SurfaceModelUtils::checkClosedFaces(surf, 10.0*int_tol_);
      sub_sfs.insert(sub_sfs.end(), splitsfs.begin(), splitsfs.end());
    }
  
  vector<vector<shared_ptr<CurveOnSurface> > > all_int_cvs(surfaces_.size());
  vector<shared_ptr<BoundedSurface> > bd_sfs(surfaces_.size());
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      shared_ptr<ParamSurface> sf1 = surfaces_[ki]->surface();
      for (size_t kj=ki+1; kj<surfaces_.size(); ++kj)
  	{
	  if (surfaces_[ki]->isTangential(surfaces_[kj].get()))
	    continue;
	  shared_ptr<ParamSurface> sf2 = surfaces_[kj]->surface();

  	  // Intersect sf1 and sf2
  	  // Remember intersection curves
	  shared_ptr<BoundedSurface> bd1, bd2;
	  vector<shared_ptr<CurveOnSurface> > int_cvs1, int_cvs2;
	  BoundedUtils::getSurfaceIntersections(sf1, sf2, int_tol_,
						int_cvs1, bd1, int_cvs2, bd2);
	  bd_sfs[ki] = bd1;
	  bd_sfs[kj] = bd2;
	  if (int_cvs1.size() > 0)
	    all_int_cvs[ki].insert(all_int_cvs[ki].end(), int_cvs1.begin(), int_cvs1.end());
	  if (int_cvs2.size() > 0)
	    all_int_cvs[kj].insert(all_int_cvs[kj].end(), int_cvs2.begin(), int_cvs2.end());
  	}
    }
  
  std::ofstream of2("intcvs.g2");
  for (size_t ki=0; ki<all_int_cvs.size(); ++ki)
    for (size_t kj=0; kj<all_int_cvs[ki].size(); ++kj)
      {
	shared_ptr<ParamCurve> cv = all_int_cvs[ki][kj]->spaceCurve();
	cv->writeStandardHeader(of2);
	cv->write(of2);
      }

  size_t nmb_sfs = surfaces_.size();
  for (size_t ki=0; ki<nmb_sfs; ++ki)
    {
      vector<shared_ptr<HedgeSurface> > added_sfs;
      surfaces_[ki]->doTrim(all_int_cvs[ki], bd_sfs[ki], int_tol_, added_sfs);
      if (added_sfs.size() > 0)
	surfaces_.insert(surfaces_.end(), added_sfs.begin(), added_sfs.end());
    }
  
  std::ofstream of3("trimsurfs.g2");
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      surfaces_[ki]->surface()->writeStandardHeader(of3);
      surfaces_[ki]->surface()->write(of3);
    }
  int stop_break = 1;
}

 //===========================================================================
shared_ptr<SurfaceModel> RevEng::createModel()
//===========================================================================
{
  vector<shared_ptr<ftSurface> > tmpsfs(surfaces_.begin(), surfaces_.end());
  sfmodel_ = shared_ptr<SurfaceModel>(new SurfaceModel(approx_tol_, 10.0*int_tol_,
						       100*int_tol_, anglim_, 10*anglim_,
						       tmpsfs));
  return sfmodel_;
}

 //===========================================================================
void RevEng::initParameters()
//===========================================================================
{
  // Set default parameters
  min_next_ = 10;  // Minimum number of neighbouring points
  max_next_ = 500; //std::min(80, tri_sf_->size()/200); //500;
  rfac_ = 3.0;  // Factor for radius in which to search for neighbouring points
  cfac_ = 8.0;  // Edge points from curvature is given by
  // cfac_ times the average length of triangulation edges in a vertex
  norm_plane_lim_= 0.005; // Limit for when the cone angle corresponding
  // to triangle normals indicate an edge
  zero_H_ = 0.005; //0.001; //0.0001;  // When mean curvature is considered zero
  zero_K_ = 0.005; //0.001; //0.0001;  // When Gauss curvature is considered zero
  zero_si_ = 0.0075; //0.001; // When shape index is considered zero
  norm_ang_lim_ = 0.1*M_PI; // Limit for when the cone angle corresponding
    // to triangle normals indicate an edge
  pca_lim_ = cness_lim_ = std::numeric_limits<double>::max();
  min_point_region_ = 200; //50; //10;  // Should be updated with regard to the total
  // number of points
  approx_tol_ = 0.001;  // Very preliminary
  int_tol_ = 1.0e-6;
  anglim_ = 0.01;
  max_nmb_outlier_ = 3;
  rpix_ = 1;
  rpfac_ = 0.05;
  ffac_ = 0.01;
  sfac_ = 0.05;
}

 //===========================================================================
int RevEng::setSmallRegionNumber()
//===========================================================================
{
  vector<int> nmb_pt_reg(regions_.size());
  for (size_t ki=0; ki<regions_.size(); ++ki)
    nmb_pt_reg[ki] = regions_[ki]->numPoints();

  std::sort(nmb_pt_reg.begin(), nmb_pt_reg.end());
  int tot_num = tri_sf_->size();
  int num_reg = (int)regions_.size();
  int idel = tot_num/num_reg;
  int min_num = std::min(10, tot_num);
  int ixmax = (int)(0.99*num_reg);
  int max_num = std::max(min_num, nmb_pt_reg[ixmax]);
  max_num = std::min(max_num, 10*idel);
  int ixdel = num_reg/100;
  int prev = nmb_pt_reg[0], prev0 = 0;
  int ix;
  int fac = 2;
  for (ix=ixdel; ix<num_reg; ix+=ixdel)
    {
      int diff = nmb_pt_reg[ix] - prev;
      if (diff > fac*(prev-prev0))
	break;
      if (diff > 0)
	prev0 = prev;
      prev = nmb_pt_reg[ix];
    }
  ix = std::min(ix, ixmax);
      
  int num = std::max(min_num, std::min(nmb_pt_reg[ix], max_num));
  return num;
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
  
  for (int ki=0; ki<nmbpts; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]); 
      double rp[3];
      setRp(pt, rp);
      pt->setRp(rp);
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
