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

#include "GoTools/compositemodel/ImplicitApprox.h"
#include "GoTools/implicitization/ImplicitizePointCloudAlgo.h"
#include "GoTools/geometry/PointCloud.h"

using namespace Go;
using std::vector;

//===========================================================================
ImplicitApprox::ImplicitApprox()
  : eps_(1.0e-12)
//===========================================================================
{
}

//===========================================================================
ImplicitApprox::~ImplicitApprox()
//===========================================================================
{
}

//===========================================================================
void ImplicitApprox::approx(vector<RevEngPoint*> points, int degree)
//===========================================================================
{
  // Extract xyz values
  vector<Vector3D> xyz(points.size());
  for (size_t ki=0; ki<points.size(); ++ki)
    xyz[ki] = points[ki]->getPoint();

  PointCloud3D pointset(xyz);

  // Implicitize
  degree_ = degree;
  ImplicitizePointCloudAlgo implicitize(pointset, degree);
  implicitize.perform();
  
  // Get result
  implicitize.getResultData(implicit_, bc_, sigma_min_);
  
  // Differentiate
  Vector4D bdir1(1.0, 0.0, 0.0, 0.0);
  Vector4D bdir2(0.0, 1.0, 0.0, 0.0);
  Vector4D bdir3(0.0, 0.0, 1.0, 0.0);
  Vector4D bdir4(0.0, 0.0, 0.0, 1.0);
  implicit.deriv(1, bdir1, deriv1_);
  implicit.deriv(1, bdir2, deriv2_);
  implicit.deriv(1, bdir3, deriv3_);
  implicit.deriv(1, bdir4, deriv4_);

}

//===========================================================================
double ImplicitApprox::estimateDist(RevEngPoint* pt)
//===========================================================================
{
  Vector3D xyz = pt->getPoint();
  Vector4D bary = bc.cartToBary(xyz);
  double dist0 = implicit(bary);
  double d1 = deriv1_(bary);
  double d2 = deriv2_(bary);
  double d3 = deriv3_(bary);
  double d4 = deriv4_(bary);
  Vector4D dv(d1,d2,d3,d4);
  Vector4D bary2 = bary+dv;
  Vector3D pt2 = bc.baryToCart(bary2);
  Vector3D grad = pt2 - curr;
  double len = grad.length();
  return (len > eps_) ? dist/len : dist;
}

