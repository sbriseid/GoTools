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

#include "GoTools/utils/config.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRApproxApp.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;
using std::vector;



int colors[3][3] = {
  {0, 255, 0},
  {255, 255, 255},
  {255, 0, 0},
};



int main(int argc, char *argv[])
{
  if (argc != 6 && argc != 7) {
    std::cout << "Usage: surface in (.g2), point cloud (.g2), points_out.g2, max level, nmb _levels, (use projected distance (0/1))" << std::endl;
    return -1;
  }

  std::ifstream sfin(argv[1]);
  std::ifstream ptsin(argv[2]);
  std::ofstream fileout(argv[3]); 
  
  double max_level = atof(argv[4]);
  int nmb_level = atoi(argv[5]);
  double min_level = -max_level;

  int use_proj = 0;
  if (argc == 7)
    use_proj = atoi(argv[6]);

  ObjectHeader header1;
  header1.read(sfin);
  shared_ptr<LRSplineSurface> sf1(new LRSplineSurface());
  sf1->read(sfin);

  ObjectHeader header2;
  header2.read(ptsin);
  PointCloud3D points;
  points.read(ptsin);

  int nmb_pts = points.numPoints();
  vector<double> data(points.rawData(), points.rawData()+3*nmb_pts);

  int ki;
  vector<double> limits(2*nmb_level+1);
  vector<vector<double> > level_points(2*nmb_level+2);

  // Set distance levels 
  double del = max_level/(double)nmb_level;
  limits[nmb_level] = 0;
  for (ki=1; ki<=nmb_level; ++ki)
    {
      limits[nmb_level-ki] = -ki*del;
      limits[nmb_level+ki] = ki*del;
    }

  double maxdist = 0.0;
  double mindist = 0.0;
  double avdist = 0.0;
  int nmb;
  vector<int> nmb_group;
  LRApproxApp::classifyCloudFromDist(data, sf1, limits, maxdist, mindist,
				     avdist, nmb, level_points, nmb_group,
				     use_proj);

  // Write to file
  for (ki=0; ki<level_points.size(); ++ki)
    {
      std::cout << "Level: " << limits[ki] << ", no. of pts: " << level_points[ki].size()/3 << std::endl;
      if (level_points[ki].size() == 0)
	continue;

      // Make point cloud
      PointCloud3D level_cloud(level_points[ki].begin(), level_points[ki].size()/3);

      double cc[3];
      if (ki <= nmb_level)
      	{
      	  cc[0] = ((nmb_level-ki)*colors[0][0] + ki*colors[1][0])/nmb_level;
      	  cc[1] = ((nmb_level-ki)*colors[0][1] + ki*colors[1][1])/nmb_level;
      	  cc[2] = ((nmb_level-ki)*colors[0][2] + ki*colors[1][2])/nmb_level;
      	}
      else
      	{
      	  cc[0] = ((ki-nmb_level-1)*colors[2][0] + 
      		   (2*nmb_level-ki+1)*colors[1][0])/nmb_level;
      	  cc[1] = ((ki-nmb_level-1)*colors[2][1] + 
      		   (2*nmb_level-ki+1)*colors[1][1])/nmb_level;
      	  cc[2] = ((ki-nmb_level-1)*colors[2][2] + 
      		   (2*nmb_level-ki+1)*colors[1][2])/nmb_level;
      	}

      fileout << "400 1 0 4 " << cc[0] << " " << cc[1];
      fileout << " " << cc[2] << " 255" << std::endl;
      level_cloud.write(fileout);
    }

  std::cout << "Max dist: " << maxdist << "Max dist below: " << mindist;
  std::cout << ", average dist: " << avdist << std::endl;
}

      




 
