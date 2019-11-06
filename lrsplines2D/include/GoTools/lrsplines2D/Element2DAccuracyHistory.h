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

#ifndef _ELEMENT2DACCURACYHISTORY_H
#define _ELEMENT2DACCURACYHISTORY_H

#include <vector>
#include "GoTools/lrsplines2D/Element2DAccuracyInfo.h"
#include "GoTools/lrsplines2D/Element2D.h"
#include <fstream>
#include <iostream> // @@ debug

namespace Go
{
  class Element2DAccuracyHistory
  {
  public:
    // Constructors
    Element2DAccuracyHistory()
      : curr_level_(0)
      {
      }

    Element2DAccuracyHistory(int max_iter)
      : curr_level_(0)
      {
	elementhist_.resize(max_iter+1);
      }

    void setMaxIter(int max_iter)
      {
	elementhist_.resize(max_iter+1);
      }

    int getMaxIter()
    {
      return (int)elementhist_.size() - 1;
    }

    void setCurrentLevel(int curr_level)
    {
      if (curr_level < (int)elementhist_.size())
	curr_level_ = curr_level;
    }

    Element2DAccuracyInfo* addElementInfo(Element2D* element)
    {
      size_t nmb = elementhist_[curr_level_].size();
      elementhist_[curr_level_].push_back(std::unique_ptr<Element2DAccuracyInfo>
				   (new Element2DAccuracyInfo(element->umin(),
							      element->umax(),
							      element->vmin(),
							      element->vmax())));
      if (element->hasElementAccuracyInfo())
	  elementhist_[curr_level_][nmb]->updateAccuracyPtr(element->getElementAccuracyInfo());
      int nmb_pts = element->nmbDataPoints() + element->nmbSignificantPoints();
      elementhist_[curr_level_][nmb]->addNmbPoints(nmb_pts);
      element->setElementAccuracyInfo(elementhist_[curr_level_][nmb].get());
      return elementhist_[curr_level_][nmb].get();
    }

    //DEBUG
    void checkAccuracyChange(int level, char* outfile, char* outfile2,
			     double frac);

  private:
    // Contains accuracy information distributed with regard to
    // iteration level and element
    int curr_level_;
    std::vector<std::vector<std::unique_ptr<Element2DAccuracyInfo> > > elementhist_;

  }; // end class Element2DAccuracyHistory
};  // end namespace Go

#endif // _ELEMENT2DACCURACYHISTORY_H
