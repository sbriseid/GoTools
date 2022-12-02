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

#include "GoTools/compositemodel/RevEngPoint.h"
#include "GoTools/compositemodel/RevEngRegion.h"
#include "GoTools/utils/Array.h"
#include "GoTools/utils/Point.h"
#include <vector>

using namespace Go;
using std::vector;
  
//===========================================================================
RevEngPoint::RevEngPoint()
  : ftSamplePoint()
//===========================================================================
{
  int dim = 3;
  eigen1_ = Point(dim);
  eigen2_ = Point(dim);
  eigen3_ = Point(dim);
  Mongenormal_ = Point(dim);
  kvecmin_ = Point(dim);
  kvecmax_ = Point(dim);
  lambda1_ = lambda2_ = lambda3_ = -1.0;
  kmin_ = kmax_ = 0.0;
  ptdist_ = avdist_ = 0.0;
  Point dummy(0.0, 0.0, 0.0);
  // normalcone_.setFromArray(dummy.begin(), dummy.end(), 3);
  normalcone_ = DirectionCone(dummy);
  avedglen_ = -1.0;
  region_ = 0;
  visited_ = 0;
  moved_ = 0;
  outlier_ = false;
  sfdist_ = -1.0;
  sfang_ = -1.0;
}

//===========================================================================
RevEngPoint::RevEngPoint(Vector3D xyz, int bnd)
  : ftSamplePoint(xyz, bnd)
//===========================================================================
{
  int dim = 3;
  eigen1_ = Point(dim);
  eigen2_ = Point(dim);
  eigen3_ = Point(dim);
  Mongenormal_ = Point(dim);
  kvecmin_ = Point(dim);
  kvecmax_ = Point(dim);
  lambda1_ = lambda2_ = lambda3_ = -1.0;
  kmin_ = kmax_ = 0.0;
  Point dummy(0.0, 0.0, 0.0);
  // normalcone_.setFromArray(dummy.begin(), dummy.end(), 3);
  normalcone_ = DirectionCone(dummy);
  avedglen_ = -1.0;
  region_ = 0;
  visited_ = 0;
  outlier_ = false;
  sfdist_ = -1.0;
  sfang_ = -1.0;

}

//===========================================================================
RevEngPoint::~RevEngPoint()
//===========================================================================
{
}

//===========================================================================
double RevEngPoint::getMeanEdgLen()
//===========================================================================
{
  if (avedglen_ < 0.0)
    {
      double len = 0.0;
      for (size_t ki=0; ki<next_.size(); ++ki)
	{
	  double currlen = xyz_.dist(next_[ki]->getPoint());
	  len += currlen;
	}
      if (next_.size() > 0)
	len /= (double)next_.size();
      avedglen_ = len;
    }
  return avedglen_;

}

//===========================================================================
double RevEngPoint::getMeanEdgLen(double maxlen)
//===========================================================================
{
  int nmb = 0;
  double len = 0.0;
  for (size_t ki=0; ki<next_.size(); ++ki)
    {
      double currlen = xyz_.dist(next_[ki]->getPoint());
      if (currlen > maxlen)
	continue;
	  len += currlen;
      ++nmb;
    }
  if (nmb > 0)
    len /= (double)nmb;
  return len;

}

//===========================================================================
void RevEngPoint::computeTriangNormal(double lim)
//===========================================================================
{
  if (next_.size() == 0)
    return;
  double eps = 1.0e-10;
  size_t prev = next_.size()-1;
  Vector3D vec1 = next_[prev]->getPoint() - xyz_;

  size_t ki, kj;
  for (ki=0, kj=0; ki<next_.size(); prev=ki, ++ki)
    {
      Vector3D vec2 = next_[ki]->getPoint() - xyz_;
      Vector3D norm = vec1 % vec2;
      bool neighbour = next_[ki]->isNeighbour(next_[prev]);
  
      if (neighbour && vec1.length() <= lim && vec2.length() <= lim &&
	  norm.length() > eps)
	{
	  if (kj == 0)
	    normalcone_.setFromArray(norm.begin(), norm.end(), 3);
	  else
	    {
	      Point norm2(norm[0], norm[1], norm[2]);
	      normalcone_.addUnionWith(norm2);
	    }
	  ++kj;
	}
      else
	{
	  int stop_break = 1;
	}
      vec1 = vec2;
    }
  if (kj == 0)
    setOutlier();
}

//===========================================================================
int RevEngPoint::surfaceClassification(int classification_type)
//===========================================================================
{
  if (classification_type == CLASSIFICATION_CURVATURE)
    return surf_[1];
  else if (classification_type == CLASSIFICATION_SHAPEINDEX)
    return surf_[2];
  else
    return CLASSIFICATION_UNDEF;
 }

//===========================================================================
Point RevEngPoint::fetchClosePoints(double radius, int min_nmb, int max_nmb,
				    vector<Point>& nearpts)
//===========================================================================
{
  int nmb_iter = 0;
  int max_iter = 5;
  while ((int)nearpts.size() < min_nmb)
    {
      setVisited();
      vector<RevEngPoint*> near;
      for (size_t ki=0; ki<next_.size(); ++ki)
	{
	  RevEngPoint* curr = dynamic_cast<RevEngPoint*>(next_[ki]);
	  if (curr->visited())
	    continue;
	  if (xyz_.dist(curr->getPoint()) <= radius)
	    {
	      curr->setVisited();
	      near.push_back(curr);
	      curr->getNearby(xyz_, radius, max_nmb, near);
	    }
	}
      
      unsetVisited();
      for (size_t ki=0; ki<near.size(); ++ki)
	{
	  Vector3D vx = near[ki]->getPoint();
	  nearpts.push_back(Point(vx[0], vx[1], vx[2]));
	  near[ki]->unsetVisited();
	}

      if (nmb_iter > max_iter)
	break;
      
      if (nearpts.size() < min_nmb)
	{
	  radius *= std::max(1.1, (double)min_nmb/(double)nearpts.size());
	  nearpts.clear();
	}
      ++nmb_iter;
    }
  return Point(xyz_[0], xyz_[1], xyz_[2]);
}

//===========================================================================
void RevEngPoint::fetchClosePoints2(double radius, int min_nmb, int max_nmb,
				    vector<RevEngPoint*>& nearpts,
				    RevEngRegion *region)
//===========================================================================
{
  int nmb_iter = 0;
  int max_iter = 5;
  while ((int)nearpts.size() < min_nmb)
    {
      setVisited();
      vector<RevEngPoint*> near;
      for (size_t ki=0; ki<next_.size(); ++ki)
	{
	  RevEngPoint* curr = dynamic_cast<RevEngPoint*>(next_[ki]);
	  if (curr->visited())
	    continue;
	  if (region && curr->region() != region)
	    continue;
	  if (xyz_.dist(curr->getPoint()) <= radius)
	    {
	      curr->setVisited();
	      near.push_back(curr);
	      curr->getNearby(xyz_, radius, max_nmb, near, region);
	    }
	}
      
      unsetVisited();
      for (size_t ki=0; ki<near.size(); ++ki)
	{
	  nearpts.push_back(near[ki]);
	  near[ki]->unsetVisited();
	}

      if (nmb_iter > max_iter)
	break;
      
      if (nearpts.size() < min_nmb)
	{
	  radius *= std::max(1.1, (double)min_nmb/(double)nearpts.size());
	  nearpts.clear();
	}
      ++nmb_iter;
    }
}

//===========================================================================
void RevEngPoint::fetchConnected(RevEngRegion *region, int max_nmb,
				 vector<RevEngPoint*>& group)
//===========================================================================
{
  double radius = std::numeric_limits<double>::max();

  setVisited();
  vector<RevEngPoint*> connected;
  for (size_t ki=0; ki<next_.size(); ++ki)
    {
      RevEngPoint* curr = dynamic_cast<RevEngPoint*>(next_[ki]);
      if (curr->visited())
	continue;
      if (curr->region() != region)
	continue;
      curr->setVisited();
      connected.push_back(curr);
      curr->getNearby(xyz_, radius, max_nmb, connected, region);
    }
  group.push_back(this);
  group.insert(group.end(), connected.begin(), connected.end());
}

//===========================================================================
void RevEngPoint::getNearby(Vector3D xyz, double radius, int max_nmb,
			    vector<RevEngPoint*>& near, RevEngRegion* region)
//===========================================================================
{
  for (size_t ki=0; ki<next_.size(); ++ki)
    {
      RevEngPoint* curr = dynamic_cast<RevEngPoint*>(next_[ki]);
      if (curr->visited())
	continue;
      if (region && curr->region() != region)
	continue;
      if (xyz.dist(curr->getPoint()) <= radius)
	{
	  curr->setVisited();
	  near.push_back(curr);
	  if (near.size() < max_nmb)
	    curr->getNearby(xyz, radius, max_nmb, near, region);
	}
    }
}
 
//===========================================================================
void
RevEngPoint::addCovarianceEigen(Point& eigen1, double lambda1, Point& eigen2,
				double lambda2, Point& eigen3, double lambda3)
//===========================================================================
{
  eigen1_ = eigen1;
  eigen2_ = eigen2;
  eigen3_ = eigen3;
  lambda1_ = lambda1;
  lambda2_ = lambda2;
  lambda3_ = lambda3;
  sfvariation_ = lambda3_/(lambda1_ + lambda2_ + lambda3_);
}

//===========================================================================
void RevEngPoint::addMongeInfo(Point& norm, Point& mincvec, double minc, Point& maxcvec,
			       double maxc, double currdist, double avdist,
			       double eps)
//===========================================================================
{
  Mongenormal_ = norm;
  kvecmin_ = mincvec;
  kvecmax_ = maxcvec;
  kmin_ = minc;
  kmax_ = maxc;
  ptdist_ = currdist;
  avdist_ = avdist;

  meancurv0_ = meancurv_ = 0.5*(kmin_ + kmax_);
  gausscurv0_ = gausscurv_ = kmin_*kmax_;
  curvedness_ = sqrt(0.5*(kmin_*kmin_ + kmax_*kmax_));
  if (fabs(kmin_) < eps && fabs(kmax_) < eps)
    shapeindex_ = 2;
  else if (fabs(kmax_ - kmin_) < eps)
    shapeindex_ = (kmax_ + kmin_ > 0.0) ? -1.0 : 1.0;
  else
    shapeindex_ = -2.0*atan((kmax_ + kmin_)/(kmax_ - kmin_))/M_PI;
}

//===========================================================================
bool RevEngPoint::isolatedEdge(int nmb, bool close)
//===========================================================================
{
  if (notEdge())
    return false;

  int nn = 0;
  for (size_t ki=0; ki<next_.size(); ++ki)
    {
      RevEngPoint* curr = dynamic_cast<RevEngPoint*>(next_[ki]);
      bool found = (close) ? curr->closeEdge() : curr->isEdge();
      if (found)
	++nn;
    }
  return (nn <= nmb);
}


//===========================================================================
void RevEngPoint::adjustWithTriangNorm(double anglim)
//===========================================================================
{
  double ang = getTriangAngle();
  if (ang < anglim)
    {
      if (edge_[0] == PCA_EDGE)
	edge_[0] = PCA_CLOSE_EDGE;
     if (edge_[1] == C1_EDGE)
	edge_[1] = C1_CLOSE_EDGE;
     if (edge_[2] == C2_EDGE)
	edge_[2] = C2_CLOSE_EDGE;
    }
}

//===========================================================================
bool RevEngPoint::isNeighbour(RevEngRegion* reg) const
//===========================================================================
{
  for (size_t ki=0; ki<next_.size(); ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(next_[ki]);
      if (pt->region() == reg)
	return true;
    }
  return false;
}

//===========================================================================
void RevEngPoint::store(std::ostream& os) const
//===========================================================================
{
  os << index_ << " " << xyz_ << std::endl;
  os << next_.size();
  for (size_t ki=0; ki<next_.size(); ++ki)
    os << " " << next_[ki]->getIndex();
  os << std::endl;
  os << avedglen_ << " " << eigen1_ << " " << lambda1_ << " " << eigen2_;
  os << " " << lambda2_ << " " << eigen3_ << " " << lambda3_ << std::endl;
  os << Mongenormal_ << " " << kvecmin_ << " " << kmin_ << " " << kvecmax_;
  os << " " << kmax_ << std::endl;
  os << ptdist_ << " " << avdist_ << std::endl;
  normalcone_.write(os);
  for (int ka=0; ka<3; ++ka)
    os << " " << edge_[ka];
  for (int ka=0; ka<3; ++ka)
    os << " " << surf_[ka];
  os << " " << outlier_ << std::endl;
  os << std::endl;
}

//===========================================================================
void RevEngPoint::read(std::istream& is, double eps, vector<int>& next_ix) 
//===========================================================================
{
  is >> index_ >> xyz_;
  int nmb_next;
  is >> nmb_next;
  if (nmb_next > 0)
    next_ix.resize(nmb_next);
  for (int ki=0; ki<nmb_next; ++ki)
    is >> next_ix[ki];
  is >> avedglen_ >> eigen1_ >> lambda1_ >> eigen2_ >> lambda2_;
  is >> eigen3_ >> lambda3_ >> Mongenormal_ >> kvecmin_ >> kmin_;
  is >> kvecmax_ >> kmax_ >> ptdist_ >> avdist_;
  Point dummy(3);
  //normalcone_ = DirectionCone(dummy);
  normalcone_.read(is);
    meancurv0_ = meancurv_ = 0.5*(kmin_ + kmax_);
  gausscurv0_ = gausscurv_ = kmin_*kmax_;
  curvedness_ = sqrt(0.5*(kmin_*kmin_ + kmax_*kmax_));
  if (fabs(kmin_) < eps && fabs(kmax_) < eps)
    shapeindex_ = 2;
  else if (fabs(kmax_ - kmin_) < eps)
    shapeindex_ = (kmax_ + kmin_ > 0.0) ? -1.0 : 1.0;
  else
    shapeindex_ = -2.0*atan((kmax_ + kmin_)/(kmax_ - kmin_))/M_PI;
  for (int ka=0; ka<3; ++ka)
    is >> edge_[ka];
  for (int ka=0; ka<3; ++ka)
    is >> surf_[ka];
  is >> outlier_;
  sfdist_ = -1.0;
  sfang_ = -1.0;

}


