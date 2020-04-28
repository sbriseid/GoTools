/* (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
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

#include "GoTools/creators/CutCellQuad.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/geometry/LoopUtils.h"
#include "GoTools/geometry/BoundedUtils.h"
#include <iostream>
#include <fstream>

using std::vector;
using std::pair;
using std::make_pair;
using namespace Go;

//==============================================================================
CutCellQuad::CutCellQuad(vector<shared_ptr<ParamCurve> >& bd_curves,
			 double tol)
//==============================================================================
  : tol_(tol), angtol_(1.0e-3)
{
  //vector<shared_ptr<CurveLoop> > loops;
  sortBoundary(bd_curves, loops_);

  //domain_ = shared_ptr<CurveBoundedDomain>(new CurveBoundedDomain(loops));
}

//==============================================================================
int CutCellQuad::cellStat(const Point& ll, const Point& ur)
//==============================================================================
{
  // Create cell boundaries
  vector<shared_ptr<SplineCurve> > cell_bd(4);
  cell_bd[0] = shared_ptr<SplineCurve>(new SplineCurve(ll, Point(ur[0],ll[1])));
  cell_bd[1] = shared_ptr<SplineCurve>(new SplineCurve(Point(ur[0],ll[1]), ur));
  cell_bd[2] = shared_ptr<SplineCurve>(new SplineCurve(ur, Point(ll[0],ur[1])));
  cell_bd[3] = shared_ptr<SplineCurve>(new SplineCurve(Point(ll[0],ur[1]), ll));

  // Check for intersections between the boundary curves and the curve model
  double eps = 1.0e-10;
  shared_ptr<CurveBoundedDomain> domain(new CurveBoundedDomain(loops_));
  int nmb_outside = 0;
  for (size_t ki=0; ki<4; ++ki)
    {
      vector<double> start_end_par;
      domain->findPcurveInsideSegments(*cell_bd[ki], tol_, start_end_par, false);
      if (start_end_par.size() > 2)
	return 2;

      if (start_end_par.size() == 0)
	nmb_outside++;

      double t1 = cell_bd[ki]->startparam();
      double t2 = cell_bd[ki]->endparam();
      for (size_t kj=0; kj<start_end_par.size(); ++kj)
	if (start_end_par[kj] > t1 + eps && start_end_par[kj] < t2 - eps)
	  return 2;
    }
  if (nmb_outside > 0 && nmb_outside < 4)
    return 2;

  // Inside/outside test
  Vector2D mid(0.5*(ll[0]+ur[0]), 0.5*(ll[1]+ur[1]));
  bool inside = domain->isInDomain(mid, tol_);
  return (inside) ? 1 : 0;
}

//==============================================================================
void
CutCellQuad::quadrature(const Point& ll, const Point& ur,
			vector<vector<double> >& quadraturepoints,
			vector<vector<shared_ptr<ParamCurve> > >& unresolved_cells,
			int stat)
//==============================================================================
{
  // Create cell boundaries
  vector<shared_ptr<SplineCurve> > cell_bd(4);
  cell_bd[0] = shared_ptr<SplineCurve>(new SplineCurve(ll, Point(ur[0],ll[1])));
  cell_bd[1] = shared_ptr<SplineCurve>(new SplineCurve(Point(ur[0],ll[1]), ur));
  cell_bd[2] = shared_ptr<SplineCurve>(new SplineCurve(ur, Point(ll[0],ur[1])));
  cell_bd[3] = shared_ptr<SplineCurve>(new SplineCurve(Point(ll[0],ur[1]), ll));

  std::ofstream ofbd("cell_bd.g2");
  for (int kk=0; kk<4; ++kk)
    {
      cell_bd[kk]->writeStandardHeader(ofbd);
      cell_bd[kk]->write(ofbd);
    }
  
  shared_ptr<CurveBoundedDomain> domain(new CurveBoundedDomain(loops_));
  if (stat < 0 || stat > 2)
    {
      double eps = 1.0e-10;
      int nmb_outside = 0;
      for (size_t ki=0; ki<4; ++ki)
	{
	  vector<double> start_end_par;
	  domain->findPcurveInsideSegments(*cell_bd[ki], tol_, start_end_par, false);
	  if (start_end_par.size() > 2)
	    {
	      stat = 2;
	      break;
	    }

	  if (start_end_par.size() == 0)
	    nmb_outside++;

	  double t1 = cell_bd[ki]->startparam();
	  double t2 = cell_bd[ki]->endparam();
	  for (size_t kj=0; kj<start_end_par.size(); ++kj)
	    if (start_end_par[kj] > t1 + eps && start_end_par[kj] < t2 - eps)
	      {
		stat = 2;
		break;
	      }
	  if (stat == 2)
	    break;
	}
      if (nmb_outside > 0 && nmb_outside < 4)
	stat = 2;
      
      if (stat != 2)
	{
	  Vector2D mid(0.5*(ll[0]+ur[0]), 0.5*(ll[1]+ur[1]));
	  bool inside = domain->isInDomain(mid, tol_);
	  stat = (inside) ? 1 : 0;
	}
    }

  if (stat == 0)
    return;  // Outside, no quadrature points

  else if (stat == 1)
    {
      // Not a cut cell. Adapt quadrature points to cell size
      vector<double> quadpar1(quadpar_.size());
      vector<double> quadpar2(quadpar_.size());
      double del1 = ur[0] - ll[0];
      double del2 = ur[1] - ll[1];
      for (size_t ki=0; ki<quadpar_.size(); ++ki)
	quadpar1[ki] = ll[0] + quadpar_[ki]*del1;
      for (size_t ki=0; ki<quadpar_.size(); ++ki)
	quadpar2[ki] = ll[1] + quadpar_[ki]*del2;
      quadraturepoints.resize(1);
      quadraturepoints[0].resize(2*quadpar_.size()*quadpar_.size());
      for (size_t kj=0; kj<quadpar_.size(); ++kj)
	for (size_t ki=0; ki<quadpar_.size(); ++ki)
	  {
	    quadraturepoints[0][2*(kj*quadpar_.size()+ki)] = quadpar1[ki];
	    quadraturepoints[0][2*(kj*quadpar_.size()+ki)+1] = quadpar2[kj];
	  }
    }
  else
    {
      // Cut cell
      std::cout << "Cut cell" << std::endl;

      // Create cut cell representation
      vector<vector<shared_ptr<CurveLoop> > > cut_bd;
      createCutCell(loops_, cell_bd, cut_bd, true);

      for (size_t ki=0; ki<cut_bd.size(); ++ki)
	{
	  quadraturePoints(cut_bd[ki], quadraturepoints, unresolved_cells);
	}
    }

  std::ofstream ofq("quadpts.g2");
  for (size_t ki=0; ki<quadraturepoints.size(); ++ki)
    {
      ofq << "400 1 0 4 100 100 55 255" << std::endl;
      ofq << quadpar_.size()*quadpar_.size() << std::endl;
      for (size_t kj=0; kj<quadraturepoints[ki].size(); kj+=2)
	{
	  Point pt(quadraturepoints[ki][kj], quadraturepoints[ki][kj+1], 0.0);
	  ofq << pt << std::endl;
	}
    }
  int stop_break = 1;
}

//==============================================================================
void CutCellQuad::sortBoundary(vector<shared_ptr<ParamCurve> >& bd_curves,
			       vector<shared_ptr<CurveLoop> >& bd_loops)
//==============================================================================
{
  if (bd_curves.size() == 0)
    return;
  
  // Compute end points of curves and ensure that the curves are in 2D
  vector<Point> endpt(2*bd_curves.size());
  for (size_t ki=0; ki<bd_curves.size(); ++ki)
    {
      if (bd_curves[ki]->dimension() != 2)
	THROW("Curve dimension different from 2");

      endpt[2*ki] = bd_curves[ki]->point(bd_curves[ki]->startparam());
      endpt[2*ki+1] = bd_curves[ki]->point(bd_curves[ki]->endparam());
      Point mid = bd_curves[ki]->point(0.5*(bd_curves[ki]->startparam()+
					    bd_curves[ki]->endparam()));
      if (endpt[2+ki].dist(mid) + endpt[2*ki+1].dist(mid) <= tol_)
	THROW("Curve too short");
    }

  // Sort curves into loops and turn direction of curves if necessary.
  // Two curves should meet at every end point (possibly the two ends of
  // a closed curve).
  size_t ix = 0;
  vector<vector<shared_ptr<ParamCurve> > > loop_cvs;
  vector<Point> endloop;
  for (size_t ki=0; ki<bd_curves.size(); ++ki)
    {
      size_t kj;
      if (endpt[2*ki].dist(endpt[2*ki+1]) <= tol_)
	kj = loop_cvs.size();  // Closed curve
      else
	{
	  for (kj=0; kj<loop_cvs.size(); ++kj)
	    {
	      if (endloop[2*kj].dist(endloop[2*kj+1]) <= tol_)
		continue;  // Loop already closed
	  
	      double d1 = endloop[2*kj].dist(endpt[2*ki]);
	      double d2 = endloop[2*kj].dist(endpt[2*ki+1]);
	      double d3 = endloop[2*kj+1].dist(endpt[2*ki]);
	      double d4 = endloop[2*kj+1].dist(endpt[2*ki+1]);
	      if (d1 <= tol_ || d2 <= tol_)
		{
		  if (d1 < d2)
		    {
		      bd_curves[ki]->reverseParameterDirection();
		      endloop[2*kj] = endpt[2*ki+1];
		    }
		  else
		    endloop[2*kj] = endpt[2*ki];
		  loop_cvs[kj].insert(loop_cvs[kj].begin(), bd_curves[ki]);
		  break;
		}
	      else if (d3 <= tol_ || d4 <= tol_)
		{
		  if (d4 < d3)
		    {
		      bd_curves[ki]->reverseParameterDirection();
		      endloop[2*kj+1] = endpt[2*ki];
		    }
		  else
		    endloop[2*kj+1] = endpt[2*ki+1];
		  loop_cvs[kj].push_back(bd_curves[ki]);
		  break;
		}
	    }
		
	}
      if (kj == loop_cvs.size())
	{
	  vector<shared_ptr<ParamCurve> > curr_loop;
	  curr_loop.push_back(bd_curves[ki]);
	  loop_cvs.push_back(curr_loop);
	  endloop.push_back(endpt[2*ki]);
	  endloop.push_back(endpt[2*ki+1]);
	}
    }

  // Check if all loops are closed
  for (size_t kj=0; kj<loop_cvs.size(); ++kj)
    {
      if (endloop[2*kj].dist(endloop[2*kj+1]) > tol_)
	THROW("Boundary loop not closed");
    }
 
  // Ensure correct orientation of loops. First check orientation
  // Create curve loops and domains
  double int_tol = tol_; //std::max(0.1*tol_, std::min(1.0e-4, tol_));
  vector<bool> CCW(loop_cvs.size());
  bd_loops.resize(loop_cvs.size());
  vector<shared_ptr<CurveBoundedDomain> > cvdom(loop_cvs.size());
  vector<Vector2D> internalPt(loop_cvs.size());
  for (size_t kj=0; kj<loop_cvs.size(); ++kj)
    {
      CCW[kj] = LoopUtils::loopIsCCW(loop_cvs[kj], tol_, int_tol);
      bd_loops[kj] = shared_ptr<CurveLoop>(new CurveLoop(loop_cvs[kj], tol_));
      cvdom[kj] =
	shared_ptr<CurveBoundedDomain>(new CurveBoundedDomain(bd_loops[kj]));

      double u, v;
      cvdom[kj]->getInternalPoint(u, v);
      internalPt[kj] = Vector2D(u, v);
    }

  // Find outer loop
  int ixout = 0;
  for (size_t ki=1; ki<loop_cvs.size(); ++ki)
    {
      bool in1 = cvdom[ixout]->isInDomain(internalPt[ki], int_tol);
      bool in2 = cvdom[ki]->isInDomain(internalPt[ixout], int_tol);
      if (in2 && (!in1))
	ixout = (int)ki;
    }

  // Make sure that outer loop has a counter clockwise orientation while
  // the remaining have a clockwise orientation
  // NB! This will fail if the input curves has several levels of curve
  // inside curve
  if (!CCW[ixout])
    bd_loops[ixout]->turnOrientation();
  for (size_t ki=0; ki<loop_cvs.size(); ++ki)
    {
      if ((int)ki == ixout)
	continue;
  
      if (CCW[ki])
	bd_loops[ki]->turnOrientation();
    }

  if (ixout != 0)
    std::swap(bd_loops[ixout], bd_loops[0]);  // Largest loop in first
}

//==============================================================================
void
CutCellQuad::createCutCell(vector<shared_ptr<CurveLoop> >& bd_loops,
			   vector<shared_ptr<SplineCurve> >& cell_bd,
			   vector<vector<shared_ptr<CurveLoop> > >& trim_loops,
			   bool test_inside)
//==============================================================================
{
  // First intersect fetch the part of the cell boundares lying inside the curve model
  shared_ptr<CurveBoundedDomain> domain(new CurveBoundedDomain(bd_loops));
  vector<shared_ptr<ParamCurve> > split_cvs;
  for (size_t ki=0; ki<cell_bd.size(); ++ki)
    {
      vector<double> par_start_end;
      domain->findPcurveInsideSegments(*cell_bd[ki], tol_, par_start_end, true, true);
      for (size_t kj=0; kj<par_start_end.size(); kj+=2)
	split_cvs.push_back(shared_ptr<ParamCurve>(cell_bd[ki]->subCurve(par_start_end[kj],
									 par_start_end[kj+1])));
    }
  std::ofstream ofcv("split_cvs.g2");
  for (size_t ki=0; ki<split_cvs.size(); ++ki)
    {
      split_cvs[ki]->writeStandardHeader(ofcv);
      split_cvs[ki]->write(ofcv);
    }

  if (split_cvs.size() == 0)
    {
      // Touch situation
      return;
    }
      
  // Prepare for boundary loop construction of cut cell(s)
  size_t nmb_cvs = split_cvs.size();
  for (size_t ki=0; ki<nmb_cvs; ++ki)
    {
      shared_ptr<ParamCurve> tmp_crv = 
	shared_ptr<ParamCurve>(split_cvs[ki]->clone());
      tmp_crv->reverseParameterDirection();
      split_cvs.push_back(tmp_crv);
    }

  // Define tolerances
  double min_loop_tol = 1.1*loops_[0]->getMaxCurveDist();
  min_loop_tol = std::max(min_loop_tol, tol_);
  double scaling = 0.5;
  Point par_eps(scaling*min_loop_tol, scaling*min_loop_tol);
  
  CurveBoundedDomain cvdom(bd_loops);
  RectDomain dom = cvdom.containingDomain();

  // Compute trimmed loops
  vector<vector<shared_ptr<ParamCurve> > > trim_cvs =
    BoundedUtils::getBoundaryLoops(bd_loops, split_cvs, tol_, min_loop_tol,
				   min_loop_tol, par_eps);

  std::ofstream of("cut_bd.g2");
  for (size_t ki=0; ki<trim_cvs.size(); ++ki)
    for (size_t kj=0; kj<trim_cvs[ki].size(); ++kj)
      {
	trim_cvs[ki][kj]->writeStandardHeader(of);
	trim_cvs[ki][kj]->write(of);
      }

  vector<shared_ptr<ParamCurve> > cell_bd2(cell_bd.begin(), cell_bd.end());
  sortLoops(trim_cvs, cell_bd2, trim_loops, test_inside);
  std::cout << "Number of cell pieces: " << trim_loops.size() << std::endl;

  
  std::ofstream of2("cut_bd2.g2");
  for (size_t ki=0; ki<trim_loops.size(); ++ki)
    for (size_t kj=0; kj<trim_loops[ki].size(); ++kj)
      {
	int nmb = trim_loops[ki][kj]->size();
	for (int ka=0; ka<nmb; ++ka)
	  {
	    (*trim_loops[ki][kj])[ka]->writeStandardHeader(of2);
	    (*trim_loops[ki][kj])[ka]->write(of2);
	  }
      }

  int stop_break = 1;
}


//==============================================================================
void
CutCellQuad::sortLoops(vector<vector<shared_ptr<ParamCurve> > >& loop_cvs,
		       vector<shared_ptr<ParamCurve> >& cell_cvs,
		       vector<vector<shared_ptr<CurveLoop> > >& trim_loops,
		       bool test_inside)
//==============================================================================
{
  // Collect the loops lying inside the cell
  double int_tol = tol_; //std::max(1.0e-6, 0.1*tol_);
  shared_ptr<CurveLoop> cell_loop;
  shared_ptr<CurveBoundedDomain> cell_domain;
  if (test_inside)
    {
      cell_loop = shared_ptr<CurveLoop>(new CurveLoop(cell_cvs, tol_));
      cell_domain = shared_ptr<CurveBoundedDomain>(new CurveBoundedDomain(cell_loop));
    }
  vector<shared_ptr<CurveLoop> > in_loops_ccw;
  vector<shared_ptr<CurveLoop> > in_loops_cw;
  for (size_t ki=0; ki<loop_cvs.size(); ++ki)
    {
      shared_ptr<CurveLoop> curr_loop(new CurveLoop(loop_cvs[ki], tol_));
      bool inside = true;
      if (test_inside)
	{
	  shared_ptr<CurveBoundedDomain> trim_domain(new CurveBoundedDomain(curr_loop));
	  double upar, vpar;
	  trim_domain->getInternalPoint(upar, vpar);
	  inside = cell_domain->isInDomain(Vector2D(upar,vpar), tol_);
	  if (inside)
	    {
	      // Test also a boundary points
	      for (size_t kj=0; kj<loop_cvs[ki].size(); ++kj)
		{
		  double t1 = 0.5*(loop_cvs[ki][kj]->startparam() +
				   loop_cvs[ki][kj]->endparam());
		  Point bd_pt = loop_cvs[ki][kj]->point(t1);
		  bool bd_in = cell_domain->isInDomain(Vector2D(bd_pt[0], bd_pt[1]), tol_);
		  if (!bd_in)
		    {
		      inside = false;
		      break;
		    }
		}
	    }
	}
      if (inside)
	{
	  bool ccw = LoopUtils::loopIsCCW(loop_cvs[ki], tol_, int_tol);
	  if (ccw)
	    in_loops_ccw.push_back(curr_loop);
	  else
	    in_loops_cw.push_back(curr_loop);
	}
    }

  // Combine loops
  for (size_t ki=0; ki<in_loops_ccw.size(); ++ki)
    {
      vector<shared_ptr<CurveLoop> > outer_loop;
      outer_loop.push_back(in_loops_ccw[ki]);
      trim_loops.push_back(outer_loop);
    }

  // Distribute inner loops
  for (size_t kj=0; kj<in_loops_cw.size(); ++kj)
    {
      shared_ptr<CurveBoundedDomain> innerdom(new CurveBoundedDomain(in_loops_cw[kj]));
      double upar, vpar;
      innerdom->getInternalPoint(upar, vpar);
      size_t ki;
      for (ki=0; ki<trim_loops.size(); ++ki)
	{
	  shared_ptr<CurveBoundedDomain> outdom(new CurveBoundedDomain(trim_loops[ki][0]));
	  bool inside = outdom->isInDomain(Vector2D(upar, vpar), tol_);
	  if (inside)
	    {
	      trim_loops[ki].push_back(in_loops_cw[kj]);
	      break;
	    }
	}
    }
  int stop_break = 1;
}

//==============================================================================
void
CutCellQuad::quadraturePoints(vector<shared_ptr<CurveLoop> >& cell_loops,
			      vector<vector<double> >& quadraturepoints,
			      vector<vector<shared_ptr<ParamCurve> > >& unresolved_cells)
//==============================================================================
{
  // Cell domain
  shared_ptr<CurveBoundedDomain> cvdom(new CurveBoundedDomain(cell_loops));
  RectDomain domain = cvdom->containingDomain();
  if (domain.umax() - domain.umin() < min_cell_size_ ||
      domain.vmax() - domain.vmin() < min_cell_size_)
    {
      vector<shared_ptr<ParamCurve> > unresolved_cvs;
      for (size_t ki=0; ki<cell_loops.size(); ++ki)
	{
	  int nmb;
	  for (int ka=0; ka<nmb; ++ka)
	    unresolved_cvs.push_back((*cell_loops[ki])[ka]);
	}
      return;   // Cell not handled
    }

  // Collect split points and select parameter direction
  vector<double> splitpar;
  int basedir = -1;
  splitPars(cell_loops, cvdom, domain, splitpar, basedir);

  if (splitpar.size() > 0)
    {
      // Split
      vector<shared_ptr<SplineCurve> > div_lines;
      for (size_t ki=0; ki<splitpar.size(); ++ki)
	{
	  Point pt1, pt2;
	  if (basedir == 1)
	    {
	      pt1 = Point(splitpar[ki], domain.vmin());
	      pt2 = Point(splitpar[ki], domain.vmax());
	    }
	  else
	    {
	      pt1 = Point(domain.umin(), splitpar[ki]);
	      pt2 = Point(domain.umax(), splitpar[ki]);
	    }
	  shared_ptr<SplineCurve> splitcv(new SplineCurve(pt1,pt2));
	  div_lines.push_back(splitcv);
	}

      vector<vector<shared_ptr<CurveLoop> > > trim_loops;
      createCutCell(cell_loops, div_lines, trim_loops, false);
      for (size_t ki=0; ki<trim_loops.size(); ++ki)
	{
	  quadraturePoints(trim_loops[ki], quadraturepoints, unresolved_cells);
	  int stop_break = 1;
	}
    }
  else
    {
      // Find height direction, 1 = runs in u, 2 = runs in v
      // There should be only the outer loop left now
      if (basedir < 0)
	basedir = heightDirection(cell_loops[0], domain);
 
      // Compute quadrature points
      vector<double> quadpts;
      computeQuadraturePoints(cvdom, domain, basedir, quadpts);
      quadraturepoints.push_back(quadpts);
    }
}

//==============================================================================
int CutCellQuad::heightDirection(shared_ptr<CurveLoop> cell_loop,
				 const RectDomain& domain)
//==============================================================================
{
  int basedir = -1;
  int nmbcvs = cell_loop->size();
  double umaxlen = 0.0, vmaxlen = 0.0;
  for (int ka=0; ka<nmbcvs; ++ka)
    {
      shared_ptr<ParamCurve> cv = (*cell_loop)[ka];
      Point dir;
      if (cv->isLinear(dir, angtol_))
  	{
  	  Point pt1 = cv->point(cv->startparam());
  	  Point pt2 = cv->point(cv->endparam());
  	  int bd_ix = domain.whichBoundary(Vector2D(pt1[0],pt1[1]),
  					   Vector2D(pt2[0],pt2[1]), tol_);
  	  if (bd_ix == 0 || bd_ix == 1)
  	    {
	      vmaxlen = std::max(vmaxlen, fabs(pt1[1]-pt2[1]));
  	      // if (fabs(std::min(pt1[1],pt2[1])-domain.vmin()) < tol_ &&
  	      // 	  fabs(std::max(pt1[1],pt2[1])-domain.vmax()) < tol_)
  	      // 	{
  	      // 	  basedir = 2;
  	      // 	  break;
  	      // 	}
  	    }
  	  else if (bd_ix == 2 || bd_ix == 3)
  	    {
	      umaxlen = std::max(umaxlen, fabs(pt1[0]-pt2[0]));
  	      // if (fabs(std::min(pt1[0],pt2[0])-domain.umin()) < tol_ &&
  	      // 	  fabs(std::max(pt1[0],pt2[0])-domain.umax()) < tol_)
  	      // 	{
  	      // 	  basedir = 1;
  	      // 	  break;
  	      // 	}
  	    }
  	}
    }

  basedir = (umaxlen >= vmaxlen) ? 1 : 2;
  if (basedir < 0)
    {
      // For the time being
      THROW("Unresolved cell");
    }
  return basedir;
}

//==============================================================================
void
CutCellQuad::computeQuadraturePoints(shared_ptr<CurveBoundedDomain> cvdom,
				     const RectDomain& domain,
				     int dir,
				     vector<double>& quadraturepoints)
//==============================================================================
{
  // Adapt quadrature values
  int dir2 = 3 - dir;
  double min1 = (dir == 1) ? domain.umin() : domain.vmin();
  double max1 = (dir == 1) ? domain.umax() : domain.vmax();
  double min2 = (dir == 1) ? domain.vmin() : domain.umin();
  double max2 = (dir == 1) ? domain.vmax() : domain.umax();
  double del1 = max1 - min1;
  vector<double> quadpar1(quadpar_.size());
  for (size_t ki=0; ki<quadpar_.size(); ++ki)
    quadpar1[ki] = min1 + quadpar_[ki]*del1;

  // Fetch constant parameter curves
  quadraturepoints.resize(2*quadpar_.size()*quadpar_.size());
  for (size_t ki=0; ki<quadpar1.size(); ++ki)
    {
      Point pt1, pt2;
      if (dir == 1)
	{
	  pt1 = Point(quadpar1[ki], min2);
	  pt2 = Point(quadpar1[ki], max2);
	}
      else
	{
	  pt1 = Point(min2, quadpar1[ki]);
	  pt2 = Point(max2, quadpar1[ki]);
	}
      shared_ptr<SplineCurve> crv(new SplineCurve(pt1, pt2));
      vector<double> inside;
      cvdom->findPcurveInsideSegments(*crv, tol_, inside, true, true);
      if (inside.size() != 2)
	THROW("Height function not properly defined");

      // Compute quadrature
      Point pt3 = crv->ParamCurve::point(inside[0]);
      Point pt4 = crv->ParamCurve::point(inside[1]);
      int ix = (dir == 1) ? 1 : 0;
      double del2 = pt4[ix] - pt3[ix];
      for (size_t kj=0; kj<quadpar_.size(); ++kj)
	{
	  double quad2 = pt3[ix] + quadpar_[kj]*del2;
	  if (dir == 1)
	    {
	      quadraturepoints[2*(kj*quadpar_.size()+ki)] = quadpar1[ki];
	      quadraturepoints[2*(kj*quadpar_.size()+ki)+1] = quad2;
	    }
	  else
	    {
	      quadraturepoints[2*(ki*quadpar_.size()+kj)] = quad2;
	      quadraturepoints[2*(ki*quadpar_.size()+kj)+1] = quadpar1[ki];
	    }
	}
    }
}

//==============================================================================
void CutCellQuad::splitPars(vector<shared_ptr<CurveLoop> >& cell_loops,
			    shared_ptr<CurveBoundedDomain> cvdom,
			    const RectDomain& domain, vector<double>& splitpar,
			    int& dir)
//==============================================================================
{
  // Fetch loop corners including derivatives
  vector<Point> corner;
  for (size_t ki=0; ki<cell_loops.size(); ++ki)
    {
      int nmb = cell_loops[ki]->size();
      int kj, kk;
      for (kj=0, kk=1; kj<nmb; ++kj, kk=(kk+1)%nmb)
	{
	  vector<Point> der1(2), der2(2);
	  (*cell_loops[ki])[kj]->point(der1, (*cell_loops[ki])[kj]->endparam(), 1);
	  (*cell_loops[ki])[kk]->point(der2, (*cell_loops[ki])[kk]->startparam(), 1);
	  bool incorner = domain.isOnCorner(Vector2D(der1[0][0],der1[0][1]), tol_);
	  if ((!incorner) && der1[1].angle(der2[1]) > angtol_)
	    {
	      corner.insert(corner.end(), der1.begin(), der1.end());
	      corner.push_back(der2[1]);
	    }
	}
    }

  if (corner.size() > 0)
    {
      std::ofstream cp("corner.g2");
      cp << "400 1 0 4 255 0 0 255" << std::endl;
      cp << corner.size() << std::endl;
      for (size_t ki=0; ki<corner.size(); ki+=3)
	cp << corner[ki] << "  0.0" << std::endl;
    }

  vector<double> candpar;
  defineSplits1(corner, cvdom, domain, splitpar, dir, candpar);

  if (splitpar.size() == 0)
    {
      // Check for closed inner loops
      vector<double> cell_split1, cell_split2;
      int cell_dir = -1;
      for (size_t ki=1; ki<cell_loops.size(); ++ki)
	{
	  // Check if any candidate parameter opens the loop
	  shared_ptr<CurveBoundedDomain> hole_domain(new CurveBoundedDomain(cell_loops[ki]));
	  RectDomain hole_box = hole_domain->containingDomain();
	  size_t kj;
	  for (kj=0; kj<candpar.size(); ++kj)
	    {
	      Point pt1, pt2;
	      if (dir == 1)
		{
		  pt1 = Point(hole_box.umin(), candpar[kj]);
		  pt2 = Point(hole_box.umax(), candpar[kj]);
		}
	      else
		{
		  pt1 = Point(candpar[kj], hole_box.vmin());
		  pt2 = Point(candpar[kj], hole_box.vmax());
		}
	      shared_ptr<SplineCurve> cv(new SplineCurve(pt1, pt2));
	      if (hole_domain->doIntersect(*cv, tol_))
		{
		  cell_split1.push_back(candpar[kj]);
		  cell_dir = 3 - dir;
		  break;
		}
	    }
	  if (kj = candpar.size())
	    {
	      double upar, vpar;
	      hole_domain->getInternalPoint(upar, vpar);
	      if (cell_dir == 1 || dir == 2)
		{
		  cell_split1.push_back(upar);
		  cell_split2.push_back(vpar);
		}
	      else
		{
		  cell_split1.push_back(vpar);
		  cell_split2.push_back(upar);
		  if (cell_dir < 0)
		    cell_dir = 1;
		}
	    }
	}
      if (cell_split1.size() > cell_split2.size())
	{
	  splitpar = cell_split1;
	  dir = cell_dir;
	}
    }

  if (splitpar.size() == 0)
    {
      // Search for non-monotone parameter direction of trimming curves
      vector<Point> turnpts;
      for (size_t ki=0; ki<cell_loops.size(); ++ki)
	{
	  int nmb = cell_loops[ki]->size();
	  for (int ka=0; ka<nmb; ++ka)
	    {
	      shared_ptr<ParamCurve> cv = (*cell_loops[ki])[ka];
	      Point dir;
	      if (!cv->isLinear(dir, tol_))
		{
		  fetchTurningPoints(cv, turnpts);
		}
	    }
	}
      if (turnpts.size() > 0)
	defineSplits2(turnpts, cvdom, domain, splitpar, dir, candpar);
    }
  int stop_break = 1;
}

//==============================================================================
void CutCellQuad::defineSplits1(const vector<Point>& corner,
				shared_ptr<CurveBoundedDomain> cvdom,
				const RectDomain& domain,
				vector<double>& splitpar, int& dir,
				vector<double>& candpar)
//==============================================================================
{
  double umin = domain.umin();
  double umax = domain.umax();
  double vmin = domain.vmin();
  double vmax = domain.vmax();
  Point udir(1.0, 0.0);
  Point vdir(0.0, 1.0);
  vector<double> split_u1, split_u2, split_uin, split_uout;
  vector<double> split_v1, split_v2, split_vin, split_vout;

  // Sort corner information
  double eps = 1.0e-8;
  for (size_t ki=0; ki<corner.size(); ki+=3)
    {
      int on_v = (corner[ki][0]-umin < tol_) ? 1 : ((umax-corner[ki][0] < tol_) ? 2 : 0);
      int on_u = (corner[ki][1]-vmin < tol_) ? 1 : ((vmax-corner[ki][1] < tol_) ? 2 : 0);
      double uang1 = corner[ki+1].angle(udir);
      uang1 = std::min(uang1, M_PI-uang1);
      double uang2 = corner[ki+2].angle(udir);
      uang2 = std::min(uang2, M_PI-uang2);
      double vang1 = corner[ki+1].angle(vdir);
      vang1 = std::min(vang1, M_PI-vang1);
      double vang2 = corner[ki+2].angle(vdir);
      vang2 = std::min(vang2, M_PI-vang2);
      if (on_u == 1 && (uang1<angtol_ || uang2<angtol_))
	split_u1.push_back(corner[ki][0]);
      else if (on_u == 2 && (uang1<angtol_ || uang2<angtol_))
	split_u2.push_back(corner[ki][0]);
      else if (on_v == 1 && (vang1<angtol_ || vang2<angtol_))
	split_v1.push_back(corner[ki][1]);
      else if (on_v == 2 && (vang1<angtol_ || vang2<angtol_))
	split_v2.push_back(corner[ki][1]);
      else
	{
	  Point pt1 = corner[ki] + 2*tol_*udir;
	  Point pt2 = corner[ki] - 2*tol_*udir;
	  Point pt3 = corner[ki] + 2*tol_*vdir;
	  Point pt4 = corner[ki] - 2*tol_*vdir;
	  bool in1 = cvdom->isInDomain(Vector2D(pt1[0], pt1[1]), eps);
	  bool in2 = cvdom->isInDomain(Vector2D(pt2[0], pt2[1]), eps);
	  bool in3 = cvdom->isInDomain(Vector2D(pt3[0], pt3[1]), eps);
	  bool in4 = cvdom->isInDomain(Vector2D(pt4[0], pt4[1]), eps);
	  if (in1 || in2)
	    split_vin.push_back(corner[ki][1]);
	  else
	    split_vout.push_back(corner[ki][1]);
	  if (in3 || in4)
	    split_uin.push_back(corner[ki][0]);
	  else
	    split_uout.push_back(corner[ki][0]);
	}
    }


  size_t num_u = split_u1.size() + split_u2.size() + split_uin.size();
  size_t num_v = split_v1.size() + split_v2.size() + split_vin.size();
  if (num_u == 0 && num_v > 0)
    {
      if (split_v1.size() < 2 && split_v2.size() < 2 && split_vin.size() == 0 &&
	  split_vout.size() == 0)
	{
	  dir = 1;
	  candpar.insert(candpar.end(), split_v1.begin(), split_v1.end());
	  candpar.insert(candpar.end(), split_v2.begin(), split_v2.end());
	  std::sort(candpar.begin(), candpar.end());
	}
      else 
	{
	  dir = 2;
	  splitpar.insert(splitpar.end(), split_v1.begin(), split_v1.end());
	  splitpar.insert(splitpar.end(), split_v2.begin(), split_v2.end());
	  if (split_uin.size() > 0)
	    splitpar.insert(splitpar.end(), split_vin.begin(), split_vin.end());
	  std::sort(splitpar.begin(), splitpar.end());

	  // Check for close split parameters
	  bool OK = checkSplits(splitpar);
	  if (!OK)
	    {
	      splitpar.clear();
	      splitpar.push_back(0.5*(umin+umax));
	      dir = 1;
	    }
	}
    }
  else if (num_v == 0 && num_u > 0)
    {
      if (split_u1.size() < 2 && split_u2.size() < 2 && split_uin.size() == 0 &&
	  split_vout.size() == 0)
	{
	  dir = 2;
	  candpar.insert(candpar.end(), split_u1.begin(), split_u1.end());
	  candpar.insert(candpar.end(), split_u2.begin(), split_u2.end());
	  std::sort(candpar.begin(), candpar.end());
	}
      else 
	{
	  dir = 1;
	  splitpar.insert(splitpar.end(), split_u1.begin(), split_u1.end());
	  splitpar.insert(splitpar.end(), split_u2.begin(), split_u2.end());
	  if (split_uin.size() > 0)
	    splitpar.insert(splitpar.end(), split_uin.begin(), split_uin.end());
	  std::sort(splitpar.begin(), splitpar.end());

	  // Check for close split parameters
	  bool OK = checkSplits(splitpar);
	  if (!OK)
	    {
	      splitpar.clear();
	      splitpar.push_back(0.5*(vmin+vmax));
	      dir = 2;
	    }
	}
    }
  else if (num_u > 0 && num_v > 0)
    {
      vector<double> tmpsplit1, tmpsplit2;
      tmpsplit1.insert(tmpsplit1.end(), split_u1.begin(), split_u1.end());
      tmpsplit1.insert(tmpsplit1.end(), split_u2.begin(), split_u2.end());
      if (split_uin.size() > 0)
	tmpsplit1.insert(tmpsplit1.end(), split_uin.begin(), split_uin.end());
      std::sort(tmpsplit1.begin(), tmpsplit1.end());
      tmpsplit2.insert(tmpsplit2.end(), split_v1.begin(), split_v1.end());
      tmpsplit2.insert(tmpsplit2.end(), split_v2.begin(), split_v2.end());
      if (split_vin.size() > 0)
	tmpsplit2.insert(tmpsplit2.end(), split_vin.begin(), split_vin.end());
      std::sort(tmpsplit2.begin(), tmpsplit2.end());

      double u1 = umin;
      double minu = umax - umin;
      for (size_t ki=0; ki<tmpsplit1.size(); ++ki)
	{
	  minu = std::min(minu, tmpsplit1[ki]-u1);
	  u1 = tmpsplit1[ki];
	}
      minu = std::min(minu, umax-u1);
			  
      double v1 = vmin;
      double minv = vmax - vmin;
      for (size_t ki=0; ki<tmpsplit2.size(); ++ki)
	{
	  minv = std::min(minv, tmpsplit2[ki]-v1);
	  v1 = tmpsplit2[ki];
	}
      minv = std::min(minv, vmax-v1);
      double ufrac = minu/(umax - umin);
      double vfrac = minv/(vmax - vmin);
      if (fabs(ufrac-vfrac) < 0.1)
	ufrac = (umax-umin > vmax-vmin) ? 2*vfrac : 0.5*vfrac;
      
      if (split_uin.size() > 0 ||
	  ((num_u > num_v || (num_u == num_v && ufrac > vfrac)) && split_vin.size() == 0))
	{
	  dir = 1;
	  splitpar = tmpsplit1;

	  // Check for close split parameters
	  bool OK = checkSplits(splitpar);
	  if (!OK)
	    {
	      // Try the other direction
	      dir = 2;
	      splitpar = tmpsplit2;

	      // Check for close split parameters
	      bool OK2 = checkSplits(splitpar);
	      if (!OK2)
		{
		  std::cout << "Must select a good split parameter" << std::endl;
		}
	    }
	}
      else
	{
	  dir = 2;
	  splitpar = tmpsplit2;

	  // Check for close split parameters
	  bool OK = checkSplits(splitpar);
	  if (!OK)
	    {
	      // Try the other direction
	      dir = 1;
	      splitpar = tmpsplit1;

	      // Check for close split parameters
	      bool OK2 = checkSplits(splitpar);
	      if (!OK2)
		{
		  std::cout << "Must select a good split parameter" << std::endl;
		}
	    }
	}
    }
  else if (split_uout.size() > 0 || split_vout.size() > 0)
    {
      if (split_vout.size() > split_uout.size())
	{
	  dir = 2;
	  splitpar = split_uout;
	}
      else
	{
	  dir = 1;
	  splitpar = split_vout;
	}
    }

  // Check candidate parameters
  bool OKcand = checkSplits(candpar);
  if (!OKcand)
    candpar.clear();
}

//==============================================================================
void CutCellQuad::defineSplits2(const vector<Point>& turnpts,
				shared_ptr<CurveBoundedDomain> cvdom,
				const RectDomain& domain,
				vector<double>& splitpar, int& dir,
				vector<double>& candpar)
//==============================================================================
{
}
 
//==============================================================================
bool CutCellQuad::checkSplits(vector<double>& splitpar)
//==============================================================================
{
  double eps = 1.0e-10;
  size_t kj;
  for (kj=1; kj<splitpar.size();)
    {
      if (splitpar[kj] - splitpar[kj-1] <= eps)
	{
	  splitpar[kj-1] = 0.5*(splitpar[kj-1]+splitpar[kj]);
	  splitpar.erase(splitpar.begin()+kj);
	}
      else if (splitpar[kj] - splitpar[kj-1] < min_cell_size_)
	break;
      else
	++kj;
    }
  return (kj < splitpar.size()) ? false : true;
}
//==============================================================================
void CutCellQuad::fetchTurningPoints(shared_ptr<ParamCurve> cv,
				     vector<Point>& turnpts)
//==============================================================================
{
  shared_ptr<SplineCurve> splcv(cv->geometryCurve());
  for (int pardir=0; pardir<2; ++pardir)
    {
      int prev_sgn = 0;
      int sgn;
      vector<double>::iterator c1 = splcv->coefs_begin();
      vector<double>::iterator c2 = c1 + 2;
      int numcf = splcv->numCoefs();
      for (int ki=1; ki<numcf; ++ki, c1+=2, c2+=2)
	{
	  double dd = *(c2+pardir) - *(c1+pardir);
	  sgn = (dd < 0) ? -1 : 1;
	  if (sgn*prev_sgn < 0)
	    {
	      double t1 = splcv->basis().grevilleParameter(ki-1);
	      double t2 = splcv->basis().grevilleParameter(ki);
	      Point pos = splcv->ParamCurve::point(0.5*(t1+t2));
	      turnpts.push_back(pos);
	    }
	}
    }
}
