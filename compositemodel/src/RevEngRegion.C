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

#include "GoTools/compositemodel/RevEngRegion.h"
#include "GoTools/compositemodel/RevEngPoint.h"
//#include "GoTools/utils/Array.h"
//#include "GoTools/utils/Point.h"
#include <vector>

using namespace Go;
using std::vector;
  
//===========================================================================
RevEngRegion::RevEngRegion()
//===========================================================================
  : classification_type_(CLASSIFICATION_UNDEF), associated_sf_(0)
{
}

//===========================================================================
RevEngRegion::RevEngRegion(int classification_type)
//===========================================================================
  : classification_type_(classification_type), associated_sf_(0)
{
}

//===========================================================================
void RevEngRegion::collect(RevEngPoint *pt)
//===========================================================================
{
  if (pt->hasRegion() && pt->region() != this)
    return;  // Cannot grow
  group_points_.push_back(pt);
  if (classification_type_ == CLASSIFICATION_UNDEF)
    return; // Cannot grow
  int type = pt->surfaceClassification(classification_type_);
  if (type == C1_UNDEF)  // SI_UNDEF == C1_UNDEF
    return; // Cannot grow

  vector<RevEngPoint*> grouped;
  grouped.push_back(pt);
  for (size_t kj=0; kj<grouped.size(); ++kj)
    {
      vector<ftSamplePoint*> next = grouped[kj]->getNeighbours();
      for (size_t ki=0; ki<next.size(); ++ki)
	{
	  RevEngPoint *curr = dynamic_cast<RevEngPoint*>(next[ki]);
	  if (!curr)
	    continue;  // Should not happen
	  if (curr->hasRegion())
	    continue;  // Already belonging to a segment
	  if (curr->isEdge())
	    continue;  // An edge point
	  int type2 = curr->surfaceClassification(classification_type_);
	  if (type2 != type)
	    continue;   // Different classification

	  // Add to region
	  curr->setRegion(this);
	  group_points_.push_back(curr);

	  // Continue growing from this point
	  grouped.push_back(curr);
	}
    }
}

//===========================================================================
RevEngPoint* RevEngRegion::seedPoint()
//===========================================================================
{
  int min_next = std::max(10, (int)group_points_.size()/100);
  double rfac = 3;
  double min_out = std::numeric_limits<double>::max();
  int min_ix = -1;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      double local_len = group_points_[ki]->getMeanEdgLen();
      double radius = rfac*local_len;
      vector<RevEngPoint*> nearpts;
      group_points_[ki]->fetchClosePoints2(radius, min_next, nearpts);

      // Count deviant points
      int deviant = 0;
      for (size_t kj=0; kj<nearpts.size(); ++kj)
      {
	if (nearpts[kj]->region() != this)
	  ++deviant;
      }

      if (deviant == 0)
	return group_points_[ki];   // Seed point found

      double curr_dev = (double)deviant/(double)nearpts.size();
      if (curr_dev < min_out)
	{
	  min_out = curr_dev;
	  min_ix = (int)ki;
	}
    }
  return (min_ix >= 0) ? group_points_[min_ix] : 0;
}
