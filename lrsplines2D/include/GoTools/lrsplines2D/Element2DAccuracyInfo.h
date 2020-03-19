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

#ifndef _ELEMENT2DACCURACYINFO_H
#define _ELEMENT2DACCURACYINFO_H

#include <vector>

namespace Go
{
  class Element2D;

  class Element2DAccuracyInfo
  {
  public:
    Element2DAccuracyInfo(double umin, double umax, double vmin, double vmax)
      : umin_(umin), umax_(umax), vmin_(vmin), vmax_(vmax), info_set_(false),
      maxerr_(0.0), averr_(0.0), nmb_pts_(0), nmb_out_(0), prev_element_(0)
    {
    }

    void addNmbPoints(int nmb)
    {
      nmb_pts_ += nmb;
    }

    void updateAccuracyPtr(Element2DAccuracyInfo* prev)
    {
      if (prev)
	prev->next_element_.push_back(this);
      prev_element_ = prev;
    }

    void setAccuracyInfo(double max_err, double av_err, int nmb_out)
    {
      maxerr_ = max_err;
      averr_ = av_err;
      nmb_out_ = nmb_out;
      info_set_ = true;
    }

    void resetElementInfo(Element2D* element);

    Element2DAccuracyInfo* getPrevious()
    {
      return prev_element_;
    }

    void getAccuracyInfo(int& nmb_pts, int& nmb_out, 
			 double& max_err, double& av_err)
    {
      nmb_pts = nmb_pts_;
      nmb_out = nmb_out_;
      max_err = maxerr_;
      av_err = averr_;
    }

    void getElementDomain(double& umin, double& umax, 
			  double& vmin, double& vmax)
    {
      umin = umin_;
      umax = umax_;
      vmin = vmin_;
      vmax = vmax_;
    }

    std::vector<Element2DAccuracyInfo*> getNextElements()
      {
	return next_element_;
      }

  private:
    double umin_, umax_, vmin_, vmax_;  // Element boundaries
    bool info_set_;  // Whether or not the accuracy information is updated
    double maxerr_;  // Maximum approximation error
    double averr_;   // Average approximation error;
    int nmb_pts_;    // Number of data points in element
    int nmb_out_;    // Number of data points outside tolerance threshold
    Element2DAccuracyInfo* prev_element_;  // Element at previous iteration level;
    std::vector<Element2DAccuracyInfo*> next_element_;  // Elements at next level
  }; // end class Element2DAccuracyInfo
};  // end namespace Go

#endif // _ELEMENT2DACCURACYINFO_H
