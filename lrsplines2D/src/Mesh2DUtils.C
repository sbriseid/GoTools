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

#include "GoTools/lrsplines2D/Mesh2DUtils.h"
#include "GoTools/lrsplines2D/LRBSpline2DUtils.h"

#include <map>
#include <array>
#include <math.h>
#include <stdexcept> // @@ for debug purposes only 

using namespace std;

namespace Go {


// =============================================================================
  bool Mesh2DUtils::identify_patch_lower_left(const Mesh2D&m, double u, 
					      double v, int& x_ix, int& y_ix)
// =============================================================================
{
  double tol = 1.0e-8;

  x_ix = last_nonlarger_knotvalue_ix(m, XFIXED, u);
  y_ix = last_nonlarger_knotvalue_ix(m, YFIXED, v);

  // adjustment of index if positioned _exactly_ at upper bound of grid
  if ((x_ix == m.numDistinctKnots(XFIXED) - 1) && (fabs(u-m.maxParam(XFIXED)) < tol))
    --x_ix;
  if ((y_ix == m.numDistinctKnots(YFIXED) - 1) && (fabs(v-m.maxParam(YFIXED)) < tol))
    --y_ix;
  
  // checking if a valid corner was found
  if (x_ix < 0 || x_ix >= m.numDistinctKnots(XFIXED) - 1) return false; // u outside domain
  if (y_ix < 0 || y_ix >= m.numDistinctKnots(YFIXED) - 1) return false; // v outside domain

  // We have now found the largest smaller knot in the u and v direction.  From here we
  // can search downwards to the lower-left corner of the containing mesh element, which
  // defines the sought-for "patch" of the surface.

  x_ix = search_downwards_for_nonzero_multiplicity(m, XFIXED, x_ix, y_ix);
  y_ix = search_downwards_for_nonzero_multiplicity(m, YFIXED, y_ix, x_ix);

  return true;
}

// =============================================================================
  bool Mesh2DUtils::identify_patch_upper_right(const Mesh2D&m, double u, 
					       double v, int& x_ix, int& y_ix)
// =============================================================================
{
  double tol = 1.0e-8;

  x_ix = first_larger_knotvalue_ix(m, XFIXED, u);
  y_ix = first_larger_knotvalue_ix(m, YFIXED, v);

  // We do not need to adjust for a position _exactly_ at lower bound of grid as the value given by the
  // indices are guaranteed to be larger (or equal if at the end).

  // checking if a valid corner was found
  if (x_ix == 0 || x_ix >= m.numDistinctKnots(XFIXED)) return false; // u outside domain
  if (y_ix == 0 || y_ix >= m.numDistinctKnots(YFIXED)) return false; // v outside domain

  // We have now found the largest smaller knot in the u and v direction.  From here we
  // can search downwards to the lower-left corner of the containing mesh element, which
  // defines the sought-for "patch" of the surface.

  x_ix = search_upwards_for_nonzero_multiplicity(m, XFIXED, x_ix, y_ix);
  y_ix = search_upwards_for_nonzero_multiplicity(m, YFIXED, y_ix, x_ix);

  return true;
}

// =============================================================================
int Mesh2DUtils::search_downwards_for_nonzero_multiplicity(const Mesh2D& m, 
							   Direction2D d, 
							   int start_ix, int other_ix)
// =============================================================================
{
  // provided that the mesh is a valid LR mesh, a valid index should always be found.
  int ix;
  for (ix = start_ix; ix != 0 && m.nu(d, ix, other_ix, other_ix + 1) == 0; --ix);
  return ix;
}

//==============================================================================
int Mesh2DUtils::search_upwards_for_nonzero_multiplicity(const Mesh2D& m, 
							 Direction2D d, 
							 int start_ix, int other_ix)
//==============================================================================
{
  // provided that the mesh is a valid LR mesh, a valid index should always be found.
  int ix;
  for (ix = start_ix; ix != m.numDistinctKnots(d) && m.nu(d, ix, other_ix - 1, other_ix) == 0; ++ix);
  return ix;  // @@ not yet tested!
}

// =============================================================================
int Mesh2DUtils::search_downwards_for_nonzero_multiplicity(const Mesh2D& m, 
							   Direction2D d, 
							   int start_ix, 
							   int ix1, int ix2)
// =============================================================================
{
  // provided that the mesh is a valid LR mesh, a valid index should always be found.
  int ix;
  for (ix = start_ix; ix != 0 && m.nu(d, ix, ix1, ix2) == 0; --ix);
  return ix;
}

//==============================================================================
int Mesh2DUtils::search_upwards_for_nonzero_multiplicity(const Mesh2D& m, 
							 Direction2D d, 
							 int start_ix, 
							 int ix1, int ix2)
//==============================================================================
{
  // provided that the mesh is a valid LR mesh, a valid index should always be found.
  int ix;
  for (ix = start_ix; ix != m.numDistinctKnots(d) && m.nu(d, ix, ix1, ix2) == 0; ++ix);
  return ix;  // @@ not yet tested!
}





// =============================================================================
int Mesh2DUtils::last_nonlarger_knotvalue_ix(const Mesh2D&m, Direction2D d, 
					     double par)
// =============================================================================
{
  const double* a = m.knotsBegin(d);
  const double* b = m.knotsEnd(d);
  
  double tol = 0.0; //1.0e-8;
  // if (par < a[0] && par >= a[0]-tol)
  //   par = a[0];
  // if (par > b[-1] && par <= b[-1]+tol)
  //   par = b[-1];

  // searching for last nonlarger knotvalue using bisection
  for (int diff = (b-a)/2; diff != 0; diff = (b-a)/2) {
    if (par < a[0]+tol && par >= a[0]-tol)
      par = a[0];
    if (par > b[-1]-tol && par <= b[-1]+tol)
      par = b[-1];
    const double* mid = a + diff;
    ( (*mid > par) ? b : a) = mid;
  }

  if (b == m.knotsEnd(d))
    --b;
  if (fabs(par-b[0]) < tol && fabs(par-a[0]) > tol)
    return (b - m.knotsBegin(d));
  else
    return (a - m.knotsBegin(d)); // if index becomes negative, it signalizes that 'par' 
                                // is smaller than even the first knot value
}



// // =============================================================================
// int last_nonlarger_knotvalue_ix(const Mesh2D&m, Direction2D d, double par)
// // =============================================================================
// {
//   const int IMAX = m.numDistinctKnots(d);
//   int ix;
//   for (ix = 0; ix != IMAX && m.kval(d, ix) <= par; ++ix); // the first one that is larger

  
//   --ix; // by decrementing index by one, we get to the last nonlarger knot value.  If index becomes
//         // negative, it signalizes that 'par' is smaller than even the first knot value.
  

//   return ix;
// }

// =============================================================================
  int Mesh2DUtils::first_larger_knotvalue_ix(const Mesh2D& m, Direction2D d, 
					     double par)
// =============================================================================
{
  int ix = last_nonlarger_knotvalue_ix(m, d, par);
  if (ix < m.numDistinctKnots(d)-1)
    ix++;
  return  ix; 
}

// =============================================================================
  void Mesh2DUtils::derive_Bspline_mesh(const Mesh2D& m, 
					int degree1, int degree2,
					int start1, int start2,
					vector<vector<int> >& kvec_u, 
					vector<vector<int> >& kvec_v)
// =============================================================================
{
  kvec_u.resize(0);
  kvec_v.resize(0);
  int end1 = m.numDistinctKnots(XFIXED);
  int end2 = m.numDistinctKnots(YFIXED);
  int ix1 = start1; //search_upwards_for_nonzero_multiplicity(m, XFIXED, start1, start2);
  int ix2 = start2; //search_upwards_for_nonzero_multiplicity(m, YFIXED, start2, start1);
  if (ix1 != start1 || ix2 != start2)
    return;   // No valid meshline starting in given point

  int ki, kj;
  int del1=1, del2=1;

  // Until all valid set of knot vectors are found
  vector<int> mkvec_u_prev, mkvec_v_prev;
  while (true)
    {
      vector<int> mkvec_u, mkvec_v;
      int incr_red = 0;

      // While no valid set of knot vectors is found and there are still 
      // possibilities of finding one
      while (true)
	{
	  // Get a valid knot vector in the first parameter direction
	  size_t incr;
	  vector<int> mkvec2_u, mkvec2_v;
	  for (; del1<end1-degree1; ++del1)
	    {
	      mkvec2_u =
		LRBSpline2DUtils::derive_knots(m, XFIXED, ix1, 
					       std::min(ix1+degree1+del1,end1-1), 
					       ix2, std::min(ix2+degree2+del2,end2-1));
	      mkvec_u = mkvec2_u;

	      // Check if the knot vector is long enough
	      for (ki=1; ki<(int)mkvec_u.size() && mkvec_u[ki] == mkvec_u[ki-1]; ++ki);
	      if ((int)mkvec_u.size()-ki+1 >= degree1+2)
		break;
	    }
      
	  // Get a valid knot vector in the second parameter direction
	  incr = 0;
	  for (incr=0; del2<end2-degree2; del2++)
	    {
	      mkvec2_v = 
		LRBSpline2DUtils::derive_knots(m, YFIXED, ix2, 
					       std::min(ix2+degree2+del2,end2-1),
					       ix1, std::min(ix1+degree1+del1,end1-1));
	      if (mkvec_v.size() > 0 && mkvec2_v.size() > mkvec_v.size())
		{
		  // Check if the first knot vector is still valid
		  size_t kr;
		  for (kr=0; kr<mkvec_u.size(); ++kr)
		    {
		      int st = m.nu(XFIXED, mkvec_u[kr], 
				    ix2, std::min(ix2+degree2+del2,end2-1));
		      if (st == 0)
			break;
		    }
		  if (kr == mkvec_u.size())
		    incr = 0;  // Still valid
		  else
		    incr += (mkvec2_v.size() - mkvec_v.size());  // First knot 
		  // vector must be recomputed
		}

	      mkvec_v = mkvec2_v;
	      for (kj=1; ki<(int)mkvec_v.size() && mkvec_v[kj] == mkvec_v[kj-1]; ++kj);
	      if ((int)mkvec_v.size()-kj+1 >= degree2+2)
		break;
	    }
	  if (incr > 0 && (int)mkvec_v.size()-incr >= degree2+2)
	    {
	      mkvec_v.erase(mkvec_v.end()-incr, mkvec_v.end());
	      incr_red = 1;
	      incr = 0;
	    }
	  if (incr == 0)
	    break;
	}

      if (mkvec_u_prev.size() > 0 && mkvec_u.size() > 0 &&
	  mkvec_v_prev.size() > 0 && mkvec_v.size() > 0)
	{
	  // Check if the knots differ
	  size_t size1 = std::min(mkvec_u_prev.size(), mkvec_u.size());
	  size_t size2 = std::min(mkvec_v_prev.size(), mkvec_v.size());
	  for (ki=0; ki<size1; ++ki)
	    if (mkvec_u_prev[ki] != mkvec_u[ki])
	      break;
	  for (kj=0; kj<size2; ++kj)
	    if (mkvec_v_prev[kj] != mkvec_v[kj])
	      break;
	  if ((ki == size1 || kj == size2) && mkvec_u.size() == size1 && 
	      mkvec_v.size() == size2)
	    break;
	}

       // Check completeness of result
      if (mkvec_u.size() < degree1+2 || mkvec_v.size() < degree1+2)
	return;
      if (mkvec_u[0] != ix1 || mkvec_v[0] != ix2)
	return;

      // Collect B-splines
      for (kj=0; kj<(int)mkvec_v.size()-degree2-1; ++kj)
	{
	  if (mkvec_v[kj] != ix2)
	    break;
	  if (mkvec_v_prev.size() > 0)
	    {
	      int kj2;
	      int sizej = std::min(kj+degree2+2,(int)mkvec_v_prev.size());
	      for (kj2=kj; kj2<sizej; ++kj2)
		if (mkvec_v_prev[kj2] != mkvec_v[kj2])
		  break;
	      if (kj2 == kj+degree2+2)
		continue;
	    }
	  vector<int> curr_v(mkvec_v.begin()+kj, mkvec_v.begin()+kj+degree2+2);
	  for (ki=0; ki<(int)mkvec_u.size()-degree1-1; ++ki)
	    {
	      if (mkvec_u[ki] != ix1)
		break;

	      if (mkvec_u_prev.size() > 0)
		{
		  int ki2;
		  int sizei = std::min(ki+degree1+2,(int)mkvec_u_prev.size());
		  for (ki2=ki; ki2<sizei; ++ki2)
		    if (mkvec_u_prev[ki2] != mkvec_u[ki2])
		      break;
		  if (ki2 == ki+degree1+2)
		    continue;
		}

	      vector<int> curr_u(mkvec_u.begin()+ki, mkvec_u.begin()+ki+degree1+2);
	      kvec_u.push_back(curr_u);
	      kvec_v.push_back(curr_v);
	    }
	  if (ki<(int)mkvec_u.size()-degree1-1)
	    continue;
	}
      if (kj<(int)mkvec_v.size()-degree2-1)
	break;
      del1++;
      del2 += (1 - incr_red);
      mkvec_u_prev = mkvec_u;
      mkvec_v_prev = mkvec_v;
    }
  return;
}

// =============================================================================
}; // end namespace Go
// =============================================================================




