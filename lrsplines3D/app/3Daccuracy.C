#include <iostream>
#include <fstream>

#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/GoTools.h"

using namespace std;
using namespace Go;

int colors[3][3] = {
  {0, 255, 0},
  {255, 255, 255},
  {255, 0, 0},
};


int main (int argc, char *argv[]) {

  if (argc != 4) {
    cout << "usage: ./hightAccuracy <input tre pt cloud(.g2)> <input lrspline(.g2)> <dim>" << endl;
    return -1;
  }

  ifstream ifs1(argv[1]);
  ifstream ifs2(argv[2]);
  int dim = atoi(argv[3]);

  int num_pts;
  ifs1 >> num_pts;

  GoTools::init();

  vector<double> pc4d;

  ObjectHeader oh;
  oh.read(ifs2);

  shared_ptr<LRSplineVolume> vol(new LRSplineVolume());
  vol->read(ifs2);

  double maxd = 0.0;
  double avd = 0.0;
  double domain[6];
  domain[0] = domain[2] = domain[4] = 1.0e8;
  domain[1] = domain[3] = domain[5] = -1.0e8;
  for (int ix=0; ix!=num_pts; ++ix)
    {
      double p0, p1, p2, q0;
      vector<double> pos(dim);
      ifs1 >> p0 >> p1 >> p2;
      for (int ka=0; ka<dim; ++ka)
	{
	  ifs1 >> q0;
	  pos[ka] = q0;
	}

      domain[0] = std::min(domain[0], p0);
      domain[1] = std::max(domain[1], p0);
      domain[2] = std::min(domain[2], p1);
      domain[3] = std::max(domain[3], p1);
      domain[4] = std::min(domain[4], p2);
      domain[5] = std::max(domain[5], p2);

      Point volpos;
      vol->point(volpos, p0, p1, p2);
      Point ptpos(pos.begin(), pos.end());
      double dist = volpos.dist(ptpos);
      maxd = std::max(maxd, dist);
      avd += dist/(double)num_pts;
    }
  //avd /= (double)num_pts;
  std::cout << "Max d: " << maxd << ", average: " << avd << std::endl;

  std::cout << "Domain: [" << domain[0] << "," << domain[1] << "]x[" << domain[2];
  std::cout << "," << domain[3] << "]x[" << domain[4] << "," << domain[5] << "]" << std::endl;
}
