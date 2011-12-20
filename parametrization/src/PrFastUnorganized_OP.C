/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1997 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/


/********************************************************************
 FILENAME    : PrFastUnorganized_OP.C
 AUTHOR      : Kai Hormann and Michael Floater, SINTEF
 DATE        : Jan. 2001
 DESCRIPTION : Implementation of methods in the class PrFastUnorganized_OP.
 CHANGE LOG  :
*********************************************************************/

#include "GoTools/parametrization/PrFastUnorganized_OP.h"

using namespace std;

// PRIVATE MEMBER FUNCTIONS

// PUBLIC MEMBER FUNCTIONS

//-----------------------------------------------------------------------------
PrFastUnorganized_OP::PrFastUnorganized_OP(
  int n, int n_int, double* xyz_points, int num_cells)
//-----------------------------------------------------------------------------
{
  cellstruct_.setNumCells(num_cells);
  cellstruct_.attach(n,xyz_points);
  uv_.resize(n);
  int j;
  for(j=0; j<n; j++)
  {
    uv_[j].x() = 0.0;
    uv_[j].y() = 0.0;
  }
  nInt_ = n_int;
  use_k_ = 1;
  knearest_ = 20;
}

//-----------------------------------------------------------------------------
PrFastUnorganized_OP::PrFastUnorganized_OP(int num_cells)
//-----------------------------------------------------------------------------
{
  cellstruct_.setNumCells(num_cells);
  use_k_ = 1;
  knearest_ = 20;
}

//----------------------------------------------------------------------------
void
PrFastUnorganized_OP::findNeighbours(int k, vector<int>& neighbours) const
//-----------------------------------------------------------------------------
//   Return the neighbours of the k-th node in random order.
{
  neighbours.clear();

  // first: get all neighbours
  Vector3D p = cellstruct_.get3dNode(k);
  if(use_k_ == 1) 
  {
    cellstruct_.getKNearest(p,knearest_,neighbours,1);
    //cellstruct_.getBall(p,sqrt(radius2_),neighbours,1);
    //cout << neighbours.size() << endl;
  }
  else 
  {
    cellstruct_.getBall(p,sqrt(radius2_),neighbours,1);
  }

  // if it is a boundary node, put boundary neighbours in the beginning
  // and end
  int i;
  if (isBoundary(k)) 
  {
    i = (k == getNumNodes() - 1 ? nInt_ : k + 1);
    neighbours.insert(neighbours.begin(), i);
    i = (k == nInt_ ? getNumNodes() - 1 : k - 1);
    neighbours.insert(neighbours.end(), i);
      
    // delete multiple entries
    int last = (int)neighbours.size()-1;
    for (i=1; i<last; ++i)
      if ((neighbours[i] == neighbours[0]) ||
	  (neighbours[i] == neighbours[last]))
      {
	neighbours.erase(neighbours.begin()+i);
	--last;
      }
  }
  
  return;
}


//----------------------------------------------------------------------------
void PrFastUnorganized_OP::initNeighbours()
//----------------------------------------------------------------------------
{

  bool too_small = false;
  bool far_too_small = false;

  int n = getNumNodes();
  nbrs.resize(n);

  for (int i=0; i<n; i++) {
    findNeighbours(i, nbrs[i]);
    if (nbrs[i].size() < 3)
      too_small = true;
    if (nbrs[i].size() < 1)
      far_too_small = true;
  }

  if (far_too_small)
    cout << "WARNING! Some vertices have NO neighbours!" << endl;
  else
    if (too_small)
      cout << "WARNING! Some vertices have less than 3 neighbours!" << endl;
}


//----------------------------------------------------------------------------
void
PrFastUnorganized_OP::getNeighbours(int k, vector<int>& neighbours) const
//-----------------------------------------------------------------------------
//   Return the neighbours of the k-th node in random order.
{
  neighbours = nbrs[k];
  return;
}


//-----------------------------------------------------------------------------
void PrFastUnorganized_OP::print(ostream& os)
//-----------------------------------------------------------------------------
{
  os << getNumNodes() << " " << nInt_ << endl;
  int i;
  for(i=0; i<getNumNodes(); i++)
  {
    Vector3D p = cellstruct_.get3dNode(i);
    os << p.x() << " " << p.y() << " " << p.z();
    os << "  " << uv_[i].x() << " " << uv_[i].y() << "\n";
  }
}

/*
//-----------------------------------------------------------------------------
void PrFastUnorganized_OP::scan(istream& is)
//-----------------------------------------------------------------------------
{
  int numpnts;
  is >> numpnts;
  is >> nInt_;
  xyz_.resize(numpnts);
  uv_.resize(numpnts);
  int i;
  for(i=0; i<numpnts; i++)
  {
    is >> xyz_[i].x() >> xyz_[i].y() >> xyz_[i].z();
    is >> uv_[i].x() >> uv_[i].y();
  }
}

//-----------------------------------------------------------------------------
void PrFastUnorganized_OP::printRawData(ostream& os)
//-----------------------------------------------------------------------------
{
  os << xyz_.size() << ' ' << nInt_ << '\n';
  int i;
  for(i=0; i<xyz_.size(); i++)
  {
    os << xyz_[i].x() << " " << xyz_[i].y() << " " << xyz_[i].z();
    os << "\n";
  }
}
*/

//-----------------------------------------------------------------------------
void PrFastUnorganized_OP::scanRawData(istream& is, int num_cells, double noise)
//-----------------------------------------------------------------------------
{
  int numpnts;
  is >> numpnts;
  is >> nInt_;
  uv_.resize(numpnts);
  double *points = new double[3*numpnts];
  int i;
  double x,y,z;
  unsigned int seed = 1;
  srand(seed);
  double random;
  //Alter interior points by up to noise in either direction
  // thus xnew will be between x-noise and x+noise.

  for(i=0; i<nInt_; i++)
  {
    is >> x >> y >> z;
    random = static_cast<double>(rand())/RAND_MAX; // random is a number between 0 and 1
    points[3*i] = x + (2.0 * random - 1.0) * noise;
    random = static_cast<double>(rand())/RAND_MAX; // random is a number between 0 and 1
    points[3*i+1] = y + (2.0 * random - 1.0) * noise;
    random = static_cast<double>(rand())/RAND_MAX; // random is a number between 0 and 1
    points[3*i+2] = z + (2.0 * random - 1.0) * noise;
    uv_[i].x() = 0.0;
    uv_[i].y() = 0.0;
  }
  for(i=nInt_; i<numpnts; i++)
  {
    is >> x >> y >> z;
    points[3*i] = x;
    points[3*i+1] = y;
    points[3*i+2] = z;
    uv_[i].x() = 0.0;
    uv_[i].y() = 0.0;
  }
  cellstruct_.setNumCells(num_cells);
  cellstruct_.attach(numpnts,points);
  delete points;
}

