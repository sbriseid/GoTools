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

#ifndef _REVENGUTILS_H
#define _REVENGUTILS_H

#include "GoTools/utils/Point.h"
#include "GoTools/compositemodel/RevEngPoint.h"
#include "GoTools/geometry/ParamSurface.h"
#include <vector>

namespace Go {

  class RevEngPoint;
  
  namespace RevEngUtils
  {
    void principalAnalysis(Point& curr, std::vector<Point>& points, 
			   double lambda[3], double eigenvec[3][3]);

    void computeMonge(Point& curr, std::vector<Point>& points,
		      Point& vec1, Point& vec2, Point& normal, Point& mincvec,
		      double& minc, Point& maxcvec, double& maxc,
		      double& currdist, double& avdist);

    void computeAxis(std::vector<Point>& points,
		     Point& axis, Point& Cx, Point& Cy);

    void computeAxis(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
		     std::vector<RevEngPoint*>::iterator> >& points,
		     Point& axis, Point& Cx, Point& Cy);

    void coneAxis(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
		     std::vector<RevEngPoint*>::iterator> >& points,
		     Point& axis, Point& Cx, Point& Cy);

    void coneApex(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
		  std::vector<RevEngPoint*>::iterator> >& points,
		  Point axis, Point& apex, double& phi);

    void computeCylPosRadius(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
			     std::vector<RevEngPoint*>::iterator> >& points,
			     Point& low, Point& high, Point& axis, Point& Cx, 
			     Point& Cy, Point& pos, double& radius);

    void computeCircPosRadius(std::vector<Point>& points,
			      Point& axis, Point& Cx, Point& Cy,
			      Point& pos, double& radius);
    
    void computeRadius(std::vector<Point>& points,
		       Point& axis, Point& Cx, Point& Cy, double& radius);
    
    void computePlane(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
		      std::vector<RevEngPoint*>::iterator> >& points,
		      Point& pos, Point& norm);

    void rotateToPlane(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
		       std::vector<RevEngPoint*>::iterator> >& points,
		       Point& xvec, Point& axis, Point& mid, std::vector<Point>& rotated);

    void projectToPlane(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
			std::vector<RevEngPoint*>::iterator> >& points,
			Point& axis, Point& mid, std::vector<Point>& projected,
			double& maxdist, double& avdist);

     void rotateToPlane(std::vector<Point>& points,
			Point& xvec, Point& axis, Point& mid, std::vector<Point>& rotated);

   void distToSurf(std::vector<RevEngPoint*>::iterator start,
		    std::vector<RevEngPoint*>::iterator end,
		    shared_ptr<ParamSurface> surf, double tol,
		   double& maxdist, double& avdist, int& inside,
		   std::vector<RevEngPoint*>& in,
		   std::vector<RevEngPoint*>& out);
    
    void distToSurf(std::vector<Point>& points,
		    shared_ptr<ParamSurface> surf, double tol,
		    double& maxdist, double& avdist, int& inside);
  }
  
} // namespace Go


#endif // _REVENGUTILS_H

