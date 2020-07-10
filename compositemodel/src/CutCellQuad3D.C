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

#include "GoTools/compositemodel/CutCellQuad3D.h"
#include "GoTools/intersections/Identity.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/BoundedUtils.h"
#include <fstream>

using std::vector;
using std::pair;
using std::make_pair;
using namespace Go;

//==============================================================================
CutCellQuad3D::CutCellQuad3D(shared_ptr<SurfaceModel> sfmodel, double tol)
//==============================================================================
  : tol_(tol)
{
  body_ = shared_ptr<Body>(new Body(sfmodel));
}

//==============================================================================
int CutCellQuad3D::cellStat(const Point& ll, const Point& ur, int& coinc)
//==============================================================================
{
  coinc = 0;
  vector<shared_ptr<SplineSurface> > cell_sfs;
  createCellSfs(ll, ur, cell_sfs);

  vector<BoundingBox> sidebb(cell_sfs.size());
  for (size_t kr=0; kr<cell_sfs.size(); ++kr)
    sidebb[kr] = cell_sfs[kr]->boundingBox();
  
  tpTolerances top = body_->getTolerances();
  vector<shared_ptr<ParamSurface> > cell_sfs2(cell_sfs.begin(), cell_sfs.end());
  shared_ptr<SurfaceModel> cell_mod(new SurfaceModel(top.gap, top.gap, top.neighbour,
						     top.kink, top.bend, cell_sfs2));
  shared_ptr<Body> cell_body(new Body(cell_mod));
  BoundingBox cellbb = cell_mod->boundingBox();

  // Check for intersection and coincidence
  double eps = 1.0e-6; 
  int nmbshell = body_->nmbOfShells();
  int ki, kj;
  Identity ident;
  int intern = -1; 
  for (ki=0; ki<nmbshell; ++ki)
    {
      shared_ptr<SurfaceModel> shell = body_->getShell(ki);
      int nmbface = shell->nmbEntities();
      for (kj=0; kj<nmbface; ++kj)
	{
	  shared_ptr<ParamSurface> sf = shell->getSurface(kj);
	  double u, v;
	  Point sfmid = sf->getInternalPoint(u, v);
	  if (((fabs(ll[0]-sfmid[0]) < tol_ || fabs(ur[0]-sfmid[0]) < tol_) &&
	       sfmid[1] >= ll[1] && sfmid[1] <= ur[1] && sfmid[2] >= ll[2] && sfmid[2] <= ur[2]) ||
	      ((fabs(ll[1]-sfmid[1]) < tol_ || fabs(ur[1]-sfmid[1]) < tol_) &&
		sfmid[0] >= ll[0] && sfmid[0] <= ur[0] && sfmid[2] >= ll[2] && sfmid[2] <= ur[2]) ||
	      ((fabs(ll[2]-sfmid[2]) < tol_ || fabs(ur[2]-sfmid[2]) < tol_) &&
		sfmid[0] >= ll[0] && sfmid[0] <= ur[0] && sfmid[1] >= ll[1] && sfmid[1] <= ur[1]))
	    {
	      // Potentially coincident surfaces
	      ;
	    }
	  else
	    {
	      bool sfin = cell_body->isInside(sfmid);
	      if (intern == -1)
		intern = (sfin) ? 1 : 0;
	      else if ((sfin && (!intern)) || ((!sfin) && intern))
		return 2;
	    }
	  BoundingBox sfbb = sf->boundingBox();
	  if (!sfbb.overlaps(cellbb))
	    continue;
	  
	  for (size_t kr=0; kr<cell_sfs.size(); ++kr)
	    {
	      if (!sfbb.overlaps(sidebb[kr]))
		continue;

	      std::ofstream of("curr_sfs.g2");
	      sf->writeStandardHeader(of);
	      sf->write(of);
	      cell_sfs[kr]->writeStandardHeader(of);
	      cell_sfs[kr]->write(of);

	      // Possible intersection. Check for coincidence
	      int stat = ident.identicalSfs(sf, cell_sfs[kr], tol_);
	      if (stat != 0)
		{
		  coinc = 1;
		  continue;   // Avoid performing surface-surface intersection
		}

	      // Intersect
	      shared_ptr<BoundedSurface> bd1, bd2;
	      vector<shared_ptr<CurveOnSurface> > int_cv1, int_cv2;
	      BoundedUtils::getSurfaceIntersections(sf, cell_sfs[kr], eps,
						    int_cv1, bd1,
						    int_cv2, bd2, true);
	      if (int_cv1.size() > 0 || int_cv2.size() > 0)
		return 2;
	    }
	}
    }
  // Not a cut cell (except possible in connection with coincident boundary surfaces)
  // Check if cell midpoint is inside the shells
  Point mid = 0.5*(ll + ur);
  bool inside = body_->isInside(mid);
  return (inside) ? 1 : 0;
}

//==============================================================================
void
CutCellQuad3D::quadrature(const Point& ll, const Point& ur,
			  vector<vector<double> >& quadraturepoints,
			  vector<vector<double> >& pointsweights,
			  vector<vector<shared_ptr<ParamSurface> > >& unresolved_cells,
			  vector<vector<double> >& surfquads,
			  vector<vector<double> >& sfptweights,
			  vector<vector<shared_ptr<ParamSurface> > >& small_sfs,
			  int stat, int coinc)
//==============================================================================
{
  if (stat < 0 || stat > 2)
    stat = cellStat(ll, ur, coinc);

    if (stat == 0)
    return;  // Outside, no quadrature points

  else if (stat == 1)
    {
      // Not a cut cell. Adapt quadrature points to cell size
      vector<double> quadpar1(quadpar_.size());
      vector<double> quadpar2(quadpar_.size());
      vector<double> quadpar3(quadpar_.size());
      vector<double> wgt1(quadpar_.size());
      vector<double> wgt2(quadpar_.size());
      vector<double> wgt3(quadpar_.size());
      double del1 = ur[0] - ll[0];
      double del2 = ur[1] - ll[1];
      double del3 = ur[2] - ll[2];
      for (size_t ki=0; ki<quadpar_.size(); ++ki)
	{
	  quadpar1[ki] = ll[0] + quadpar_[ki]*del1;
	  wgt1[ki] = weights_[ki]*del1;
	}
      for (size_t ki=0; ki<quadpar_.size(); ++ki)
	{
	  quadpar2[ki] = ll[1] + quadpar_[ki]*del2;
	  wgt2[ki] = weights_[ki]*del2;
	}
      for (size_t ki=0; ki<quadpar_.size(); ++ki)
	{
	  quadpar3[ki] = ll[2] + quadpar_[ki]*del3;
	  wgt3[ki] = weights_[ki]*del3;
	}
      quadraturepoints.resize(1);
      quadraturepoints[0].resize(3*quadpar_.size()*quadpar_.size()*quadpar_.size());
      pointsweights.resize(1);
      pointsweights[0].resize(quadpar_.size()*quadpar_.size()*quadpar_.size());
      for (size_t kk=0; kk<quadpar_.size(); ++kk)
	for (size_t kj=0; kj<quadpar_.size(); ++kj)
	  for (size_t ki=0; ki<quadpar_.size(); ++ki)
	    {
	      quadraturepoints[0][3*((kk*quadpar_.size()+kj)*quadpar_.size()+ki)] = quadpar1[ki];
	      quadraturepoints[0][3*((kk*quadpar_.size()+kj)*quadpar_.size()+ki)+1] = quadpar2[kj];
	      quadraturepoints[0][3*((kk*quadpar_.size()+kj)*quadpar_.size()+ki)+2] = quadpar3[kk];
	      pointsweights[0][(kk*quadpar_.size()+kj)*quadpar_.size()+ki] = wgt1[ki]*wgt2[kj]*wgt3[kk];
	    }

      // Check for coincidence between a piece of the boundary loop and the
      // cell boundaries
      if (coinc != 0)
	{
	  // Coincidence possible or identified
	}
    }
  else
    {
      // Cut cell
      std::cout << "Cut cell" << std::endl;
      shared_ptr<Body> cutcell;
      createCutCell(ll, ur, cutcell);

      int stop_break0 = 1;
    }

    int stop_break = 1; 
}

//==============================================================================
void CutCellQuad3D::createCutCell(const Point& ll, const Point& ur,
				  shared_ptr<Body>& cutcell)
//==============================================================================
{
  vector<shared_ptr<SplineSurface> > cell_sfs;
  createCellSfs(ll, ur, cell_sfs);

  // Make surface model
  tpTolerances top = body_->getTolerances();
  vector<shared_ptr<ParamSurface> > cell_sfs2(cell_sfs.begin(), cell_sfs.end());
  shared_ptr<SurfaceModel> cell_mod(new SurfaceModel(top.gap, top.gap, top.neighbour,
						     top.kink, top.bend, cell_sfs2));
  shared_ptr<Body> cell_body(new Body(cell_mod));

  // Intersect CAD model (exclusive voids) and element and sort pieces according to
  // inside and outside
  shared_ptr<SurfaceModel> outer = body_->getShell(0);
  vector<shared_ptr<SurfaceModel> > sub_mod = cell_mod->splitSurfaceModels(outer);

  // Set code for surface origin
  int nmb1 = sub_mod[0]->nmbEntities();
  for (int ki=0; ki<nmb1; ++ki)
    {
      shared_ptr<ParamSurface> sf = sub_mod[0]->getSurface(ki);
      if (sf->getFlagCode() != 3)
	sf->setFlagCode(1);
    }
  // Collect faces being internal to both surface models
  Identity identity;
  if (sub_mod[2]->nmbEntities() > 0)
    {
      vector<shared_ptr<ftSurface> > addfaces;
      int nmb2 = sub_mod[2]->nmbEntities();
      for (int ki=0; ki<nmb2; ++ki)
	{
	  shared_ptr<ParamSurface> sf = sub_mod[2]->getSurface(ki);

	  // Check for coincidence
	  bool coinc = false;
	  for (int kj=0; kj<nmb1; ++kj)
	    {
	      shared_ptr<ParamSurface> sf1 = sub_mod[0]->getSurface(kj);
	      int stat = identity.identicalSfs(sf1, sf, tol_);
	      if (stat > 0)
		{
		  sf1->setFlagCode(3);
		  coinc = true;
		  break;
		}
	    }
	  if (!coinc)
	    {
	      sf->setFlagCode(2);
	      addfaces.push_back(sub_mod[2]->getFace(ki));
	    }
	}

      if (addfaces.size() > 0)
	{
	  for (size_t kr=0; kr<addfaces.size(); ++kr)
	    sub_mod[2]->removeFace(addfaces[kr]);
	  sub_mod[0]->append(addfaces);
	}
    }

  // Voids

  std::ofstream of("cutcell.g2");
  int nmb = sub_mod[0]->nmbEntities();
  for (int ki=0; ki<nmb; ++ki)
    {
      shared_ptr<ParamSurface> sf = sub_mod[0]->getSurface(ki);
      sf->writeStandardHeader(of);
      sf->write(of);
    }
  
  int nmbbd = sub_mod[0]->nmbBoundaries();
  if (nmbbd > 0)
    {
      vector<shared_ptr<ftEdge> > bdedg = sub_mod[0]->getBoundaryEdges();
      std::ofstream of2("ededg.g2");
      for (size_t kr=0; kr<bdedg.size(); ++kr)
	{
	  shared_ptr<ParamCurve> cv = bdedg[kr]->geomCurve();
	  shared_ptr<SplineCurve> splcv(cv->geometryCurve());
	  splcv->writeStandardHeader(of2);
	  splcv->write(of2);
	}
      
      THROW("Cut cell not closed");
    }
  cutcell = shared_ptr<Body>(new Body(sub_mod[0]));

  int stop_break = 1;
}


//==============================================================================
void CutCellQuad3D::createCellSfs(const Point& ll, const Point& ur,
				  vector<shared_ptr<SplineSurface> >& cell_sfs)
//==============================================================================
{
  cell_sfs.resize(6);

  double coefs1[] = {ll[0], ll[1], ll[2], ll[0], ur[1], ll[2], ur[0], ll[1], ll[2], ur[0], ur[1], ll[2]};
  double coefs2[] = {ll[0], ll[1], ur[2], ur[0], ll[1], ur[2], ll[0], ur[1], ur[2], ur[0], ur[1], ur[2]};
  double coefs3[] = {ll[0], ll[1], ll[2], ll[0], ll[1], ur[2], ll[0], ur[1], ll[2], ll[0], ur[1], ur[2]};
  double coefs4[] = {ur[0], ll[1], ll[2], ur[0], ur[1], ll[2], ur[0], ll[1], ur[2], ur[0], ur[1], ur[2]};
  double coefs5[] = {ll[0], ll[1], ll[2], ur[0], ll[1], ll[2], ll[0], ll[1], ur[2], ur[0], ll[1], ur[2]};
  double coefs6[] = {ll[0], ur[1], ll[2], ll[0], ur[1], ur[2], ur[0], ur[1], ll[2], ur[0], ur[1], ur[2]};

  double knots1[] = {ll[0], ll[0], ur[0], ur[0]};
  double knots2[] = {ll[1], ll[1], ur[1], ur[1]};
  double knots3[] = {ll[2], ll[2], ur[2], ur[2]};

  cell_sfs[0] = shared_ptr<SplineSurface>(new SplineSurface(2, 2, 2, 2, &knots1[0], &knots2[0],
							    &coefs1[0], 3, false));
  cell_sfs[1] = shared_ptr<SplineSurface>(new SplineSurface(2, 2, 2, 2, &knots1[0], &knots2[0],
							    &coefs2[0], 3, false));
  cell_sfs[2] = shared_ptr<SplineSurface>(new SplineSurface(2, 2, 2, 2, &knots2[0], &knots3[0],
							    &coefs3[0], 3, false));
  cell_sfs[3] = shared_ptr<SplineSurface>(new SplineSurface(2, 2, 2, 2, &knots2[0], &knots3[0],
							    &coefs4[0], 3, false));
  cell_sfs[4] = shared_ptr<SplineSurface>(new SplineSurface(2, 2, 2, 2, &knots1[0], &knots3[0],
							    &coefs5[0], 3, false));
  cell_sfs[5] = shared_ptr<SplineSurface>(new SplineSurface(2, 2, 2, 2, &knots1[0], &knots3[0],
							    &coefs6[0], 3, false));

  std::ofstream of("cell_sfs.g2");
  for (size_t ki=0; ki<cell_sfs.size(); ++ki)
    {
      cell_sfs[ki]->writeStandardHeader(of);
      cell_sfs[ki]->write(of);
    }
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << "1" << std::endl;
  of << ll << std::endl;
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << "1" << std::endl;
  of << ur << std::endl;

  int stop_break = 1;
}
