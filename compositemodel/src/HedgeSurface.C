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

#include "GoTools/compositemodel/HedgeSurface.h"
#include "GoTools/compositemodel/RevEngRegion.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/ElementarySurface.h"

using namespace Go;
using std::vector;

// //===========================================================================
// HedgeSurface::HedgeSurface()
//   : ftSurface()
// //===========================================================================
// {
// }

//===========================================================================
HedgeSurface::HedgeSurface(shared_ptr<ParamSurface> sf, RevEngRegion *region)
  : ftSurface(sf, -1)
//===========================================================================
{
  regions_.push_back(region);
  bbox_ = region->boundingBox();
}

//===========================================================================
HedgeSurface::HedgeSurface(shared_ptr<ParamSurface> sf,
			   vector<RevEngRegion*>& region)
  : ftSurface(sf, -1), regions_(region)
//===========================================================================
{
  if (region.size() > 0)
    {
      bbox_ = region[0]->boundingBox();
      for (size_t ki=1; ki<region.size(); ++ki)
	bbox_.addUnionWith(region[ki]->boundingBox());
    }
}

//===========================================================================
HedgeSurface::~HedgeSurface()
//===========================================================================
{
}

//===========================================================================
int HedgeSurface::numPoints()
//===========================================================================
{
  int num = 0;
  for (size_t ki=0; ki<regions_.size(); ++ki)
    num += regions_[ki]->numPoints();
  return num;
 }

//===========================================================================
ClassType HedgeSurface::instanceType(int& code)
//===========================================================================
{
  code = 0;  // For later use
  ClassType type = surface()->instanceType();
  if (type == Class_BoundedSurface)
    {
      shared_ptr<ParamSurface> surf = surface();
      shared_ptr<BoundedSurface> bdsurf =
	dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf);
      type = bdsurf->underlyingSurface()->instanceType();
    }
  return type;
}


//===========================================================================
bool HedgeSurface::isCompatible(HedgeSurface* other, double angtol, double approx_tol, ClassType& type, double& score)
//===========================================================================
{
  score = std::numeric_limits<double>::max();
  int code1 = -1, code2 = -1;
  ClassType type1 = instanceType(code1);
  ClassType type2 = other->instanceType(code2);

  shared_ptr<ParamSurface> surf1 = surface();
  shared_ptr<ElementarySurface> psurf1 =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf1);
  if (!psurf1.get())
    {
      shared_ptr<BoundedSurface> bdsf1 =
	dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf1);
      if (bdsf1.get())
	{
	  surf1 = bdsf1->underlyingSurface();
	  psurf1 = dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf1);
	}
    }
  
  shared_ptr<ParamSurface> surf2 = other->surface();
  shared_ptr<ElementarySurface> psurf2 =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf2);
  if (!psurf2.get())
    {
      shared_ptr<BoundedSurface> bdsf2 =
	dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf2);
      if (bdsf2.get())
	{
	  surf2 = bdsf2->underlyingSurface();
	  psurf2 = dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf2);
	}
    }

  if (!psurf1.get())
    {
      RevEngRegion *preg = 0;
      int numpt = 0;
      int nreg = numRegions();
      for (int ka=0; ka<nreg; ++ka)
	{
	  RevEngRegion *reg =  getRegion(ka);
	  if (reg->hasPrimary())
	    {
	      double maxdp, avdp;
	      int num_inp;
	      int curr_numpt = reg->numPoints();
	      reg->getPrimaryInfo(maxdp, avdp, num_inp);
	      if (avdp < approx_tol && num_inp > curr_numpt/2 && curr_numpt > numpt)
		{
		  preg = reg;
		  numpt = curr_numpt;
		}
	    }
	}
      if (preg)
	{
	  psurf1 =
	    dynamic_pointer_cast<ElementarySurface,ParamSurface>(preg->getPrimary());
	  type1 = preg->getPrimary()->instanceType();
	}
    }
  
  if (!psurf2.get())
    {
      RevEngRegion *preg = 0;
      int numpt = 0;
      int nreg = other->numRegions();
      for (int ka=0; ka<nreg; ++ka)
	{
	  RevEngRegion *reg =  other->getRegion(ka);
	  if (reg->hasPrimary())
	    {
	      double maxdp, avdp;
	      int num_inp;
	      int curr_numpt = reg->numPoints();
	      reg->getPrimaryInfo(maxdp, avdp, num_inp);
	      if (avdp < approx_tol && num_inp > curr_numpt/2 && curr_numpt > numpt)
		{
		  preg = reg;
		  numpt = curr_numpt;
		}
	    }
	}
      if (preg)
	{
	  psurf2 =
	    dynamic_pointer_cast<ElementarySurface,ParamSurface>(preg->getPrimary());
	  type2 = preg->getPrimary()->instanceType();
	}
    }
  
  if (!psurf1.get() || !psurf2.get())
    return false;
  if (type1 != type2 || code1 != code2)
    return false;  // Not the same surface type
  type = type1;

  Point loc1 = psurf1->location();
  Point loc2 = psurf2->location();
  Point vec1 = psurf1->direction();
  Point vec2 = psurf2->direction();
  double rad1 = psurf1->radius(0.0, 0.0);  // Parameters must be set properly for cone
  double rad2 = psurf2->radius(0.0, 0.0);
  double smallrad1 = psurf1->radius2(0.0, 0.0);
  double smallrad2 = psurf2->radius2(0.0, 0.0);
  vec1.normalize_checked();
  vec2.normalize_checked();

  int sgn = (type1 == Class_Plane) ? -1 : 1;
  double ang = vec1.angle(vec2);
  ang = std::min(ang, M_PI-ang);

  double dlim = (rad1 < 0.0) ? approx_tol : std::max(0.05*rad1, approx_tol);
  double anglim = 10.0*angtol;
  double eps = 1.0e-8;
  if (ang > anglim)
    return false;
  // if (fabs(rad2-rad1) > dlim && fabs(smallrad2-smallrad1) < eps)
  //   return false;
  if (fabs(rad2-rad1) > std::min(rad1, rad2) && fabs(smallrad2-smallrad1) < eps)
    return false;
  else if (smallrad1 > 0.0 &&
	   (rad1 < rad2-smallrad2 || rad1 > rad2+smallrad2 ||
	    rad2 < rad1-smallrad1 || rad2 > rad1+smallrad1))
    return false;
    
  double pdist1 = 0.0, pdist2 = 0.0;
  if (type1 == Class_Plane)
    {
      Point loc2_0 = loc2 - ((loc2-loc1)*vec1)*vec1;
      Point loc1_0 = loc1 - ((loc2-loc1)*vec2)*vec2;
      pdist1 = loc2.dist(loc2_0);
      pdist2 = loc1.dist(loc1_0);
      pdist1 = pdist2 = std::min(pdist1, pdist2);
      if (pdist1 > 2.0*dlim || pdist2 > 2.0*dlim)    
	return false;
    }
  else if (type1 == Class_Cylinder)
    {
      Point loc2_0 = loc1 + ((loc2-loc1)*vec1)*vec1;
      pdist1 = loc2.dist(loc2_0);
      if (pdist1 + fabs(rad2-rad1) > std::min(rad1, rad2))
      // if (pdist1 > dlim)
	return false;
    }
  else if (type1 == Class_Torus)
    {
      pdist1 = loc1.dist(loc2);
      if (pdist1 > 2.0*dlim || pdist2 > 2.0*dlim)    
	return false;
    }

  score = 10.0*ang + fabs(rad2-rad1) + fabs(smallrad2-smallrad1) +
    pdist1 + pdist2;
  return true;
}

//===========================================================================
bool HedgeSurface::hasPrimary()
//===========================================================================
{
  for (size_t ki=0; ki<regions_.size(); ++ki)
    if (regions_[ki]->hasPrimary())
      return true;
  return false;
}
