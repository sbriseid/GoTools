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
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/compositemodel/ftPlane.h"
#include "GoTools/compositemodel/ftCurve.h"
#include "GoTools/compositemodel/ModifyFaceSet.h"
#include "GoTools/creators/CutCellQuad.h"
#include <fstream>

using std::vector;
using std::set;
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
      vector<shared_ptr<Body> > cutcells;
      createCutCell(ll, ur, cutcells);

      for (size_t ki=0; ki<cutcells.size(); ++ki)
	{
	  quadraturePoints(ll, ur, cutcells[ki], quadraturepoints, pointsweights, unresolved_cells,
			   surfquads, sfptweights, small_sfs);
	}

      int stop_break0 = 1;
    }

    int stop_break = 1; 
}

//==============================================================================
void CutCellQuad3D::splitCell(shared_ptr<Body> body,
			      shared_ptr<ParamSurface> splitsf,
			      vector<shared_ptr<Body> >& subcell)
//==============================================================================
{
  tpTolerances top = body->getTolerances();
  vector<shared_ptr<ParamSurface> > split_sfs;
  split_sfs.push_back(splitsf);
  shared_ptr<SurfaceModel> other_mod(new SurfaceModel(top.gap, top.gap, top.neighbour,
						     top.kink, top.bend, split_sfs));

  // Intersect CAD model (exclusive voids) and element and sort pieces according to
  // inside and outside
  shared_ptr<SurfaceModel> outer = body->getShell(0);
  vector<shared_ptr<SurfaceModel> > split_mod = outer->splitSurfaceModels(other_mod);
  
  // Check for identical faces between the surfaces of the initial face
  // internal to the initial volume that are coincident with a surface in any
  // of the new volume pieces. These surfaces must be removed
  removeCoincFaces(split_mod[0], split_mod[1], split_mod[2], top.gap);

  int nmb = split_mod[2]->nmbEntities();
  for (int ki=0; ki<nmb; ++ki)
    {
      shared_ptr<ftSurface> face1 = split_mod[2]->getFace(ki);
      shared_ptr<ParamSurface> surf1 = split_mod[2]->getSurface(ki);
      shared_ptr<ParamSurface> surf2(surf1->clone());
      surf2->swapParameterDirection();
      
      shared_ptr<ftSurface> face1_2(new ftSurface(surf1, -1));
      split_mod[0]->append(face1_2, true, false, true);
      shared_ptr<ftSurface> face2(new ftSurface(surf2, -1));
      split_mod[1]->append(face2, true, false, true);
    }

  if (body->nmbOfShells() > 1)
    MESSAGE("Voids not handled");

  std::ofstream of1("submod1.g2");
  std::ofstream of2("submod2.g2");
  int nmb1 = split_mod[0]->nmbEntities();
  int nmb2 = split_mod[1]->nmbEntities();
  for (int ki=0; ki<nmb1; ++ki)
    {
      shared_ptr<ParamSurface> sf = split_mod[0]->getSurface(ki);
      sf->writeStandardHeader(of1);
      sf->write(of1);
    }
  for (int ki=0; ki<nmb2; ++ki)
    {
      shared_ptr<ParamSurface> sf = split_mod[1]->getSurface(ki);
      sf->writeStandardHeader(of2);
      sf->write(of2);
    }
  
  for (int ki=0; ki<2; ++ki)
    {
        int nmbbd = split_mod[ki]->nmbBoundaries();
	if (nmbbd > 0)
	  {
	    vector<shared_ptr<ftEdge> > bdedg = split_mod[ki]->getBoundaryEdges();
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

	  vector<shared_ptr<SurfaceModel> > connected_mod = split_mod[ki]->getConnectedModels();

	  for (size_t kj=0; kj<connected_mod.size(); ++kj)
	    subcell.push_back(shared_ptr<Body>(new Body(connected_mod[kj])));


    }
  int stop_break = 1;
}

//==============================================================================
void CutCellQuad3D::createCutCell(const Point& ll, const Point& ur,
				  vector<shared_ptr<Body> >& cutcell)
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

  vector<shared_ptr<SurfaceModel> > connected_mod = sub_mod[0]->getConnectedModels();

  for (size_t kj=0; kj<connected_mod.size(); ++kj)
    cutcell.push_back(shared_ptr<Body>(new Body(connected_mod[kj])));

  int stop_break = 1;
}

//===========================================================================
// 
// 
void
CutCellQuad3D::removeCoincFaces(shared_ptr<SurfaceModel>& mod1,
				shared_ptr<SurfaceModel>& mod2,
				shared_ptr<SurfaceModel>& mod3,
				double tol)

//===========================================================================
{
  // Check for surfaces in mod3 that are identical with surfaces in
  // mod1 or mod2. Remove those surfaces
  int nmb = (!mod3.get()) ? 0 : mod3->nmbEntities();
  for (int kj=0; kj<nmb; ++kj)
    {
      shared_ptr<ParamSurface> surf1 = mod3->getSurface(kj);
      int nmb2 = (!mod1.get()) ? 0 : mod1->nmbEntities();
      int kr;
      for (kr=0; kr<nmb2; ++kr)
	{
	  shared_ptr<ParamSurface> surf2 = mod1->getSurface(kr);
	  Identity ident;
	  int coinc = ident.identicalSfs(surf1, surf2, tol);
	  if (coinc > 0)
	    break;
	}
      if (kr < nmb2)
	{
	  mod3->removeFace(mod3->getFace(kj));
	  --kj;
	  --nmb;
	}
      else
	{
	  nmb2 = (!mod2.get()) ? 0 : mod2->nmbEntities();
	  for (kr=0; kr<nmb2; ++kr)
	    {
	      shared_ptr<ParamSurface> surf2 = mod2->getSurface(kr);
	      Identity ident;
	      int coinc = ident.identicalSfs(surf1, surf2, tol);
	      if (coinc > 0)
		break;
	    }
	  if (kr < nmb2)
	    {
	      mod3->removeFace(mod3->getFace(kj));
	      --kj;
	      --nmb;
	    }
	}
    }
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

//==============================================================================
void
CutCellQuad3D::quadraturePoints(const Point& ll, const Point& ur,
				shared_ptr<Body> body,
				vector<vector<double> >& quadraturepoints,
				vector<vector<double> >& pointsweights,
				vector<vector<shared_ptr<ParamSurface> > >& unresolved_cells,
				vector<vector<double> >& surfquads,
				vector<vector<double> >& sfptweights,
				vector<vector<shared_ptr<ParamSurface> > >& small_sfs)
//==============================================================================
{
  int splitdir;
  double splitval;
  BoundingBox bbox = body->boundingBox();
  Point ll3 = ll;
  Point ur3 = ur;
  shared_ptr<SplineCurve> rulecv;
  Point ruledir;
  int basedir = selectBaseDir(ll3, ur3, body, splitdir, splitval, rulecv, ruledir);
  if (basedir < 0)
    {
      if (splitdir >= 0)
	{
	  // Split cut cell model and treat sub models recursively
	  Point pos = 0.5*(bbox.low() + bbox.high());
	  pos[splitdir] = splitval;
	  Point norm(0.0, 0.0, 0.0);
	  norm[splitdir] = 1.0;

	  shared_ptr<Plane> splitplane(new Plane(pos, norm));
	  int udir = (splitdir == 0) ? 1 : 0;
	  int vdir = (splitdir == 2) ? 1 : 2;
	  double udel = bbox.high()[udir] - bbox.low()[udir];
	  double vdel = bbox.high()[vdir] - bbox.low()[vdir];
	  splitplane->setParameterBounds(-udel, -vdel, udel, vdel);

	  vector<shared_ptr<Body> > subcell;
	  splitCell(body, splitplane, subcell);

	  for (size_t ki=0; ki<subcell.size(); ++ki)
	    quadraturePoints(ll, ur, subcell[ki], quadraturepoints, pointsweights,
			     unresolved_cells, surfquads, sfptweights, small_sfs);

	  int stop_break = 1;
	}
      else if (rulecv.get())
	{
	  // Linearily extend curve to make sure it covers the entire model
	  double len = bbox.low().dist(bbox.high());
	  rulecv->enlarge(0.5*len, true);
	  rulecv->enlarge(0.5*len, false);

	  // Create ruled surface
	  shared_ptr<SplineCurve> cv2(rulecv->clone());
	  Point vec = len*ruledir;
	  cv2->translateCurve(vec);

	  vector<double> coefs(rulecv->coefs_begin(), rulecv->coefs_end());
	  coefs.insert(coefs.end(), cv2->coefs_begin(), cv2->coefs_end());
	  vector<double> knots2(4);
	  knots2[0] = knots2[1] = 0.0;
	  knots2[2] = knots2[3] = len;
	  shared_ptr<SplineSurface> ruledsf(new SplineSurface(rulecv->numCoefs(), 2,
							      rulecv->order(), 2,
							      rulecv->basis().begin(),
							      &knots2[0], &coefs[0],
							      rulecv->dimension(), rulecv->rational()));
	  std::ofstream ofrule("ruledsf.g2");
	  ruledsf->writeStandardHeader(ofrule);
	  ruledsf->write(ofrule);
	  
	  vector<shared_ptr<Body> > subcell;
	  splitCell(body, ruledsf, subcell);

	  for (size_t ki=0; ki<subcell.size(); ++ki)
	    quadraturePoints(ll, ur, subcell[ki], quadraturepoints, pointsweights,
			     unresolved_cells, surfquads, sfptweights, small_sfs);

	  int stop_rule = 1;
							      
	}
      else
	std::cout << "No suitable base direction found" << std::endl;
    }
  else
    {
      std::ofstream ofp("currquads.g2");
      
      // 2D cell
      Point ll2(2), ur2(2);
      int ix = 0;
      for (int kc=0; kc<3; ++kc)
	{
	  if (kc == basedir)
	    continue;
	  ll2[ix] = ll[kc];
	  ur2[ix++] = ur[kc];
	}
      
      double tmin = ll3[basedir];
      double tmax = ur3[basedir];
      double del = tmax - tmin;

      // Create 2D models
      double quadpar;
      double wgt;
      Point pos = 0.5*(bbox.low() + bbox.high());
      Point norm(0.0, 0.0, 0.0);
      norm[basedir] = 1.0;
      for (size_t ki=0; ki<quadpar_.size(); ++ki)
	{
	  quadpar = tmin + quadpar_[ki]*del;
	  wgt = weights_[ki]*del;

	  // Intersect with plane through quadrature parameter
	  pos[basedir] = quadpar;
	  ftPlane plane(norm, pos);

	  std::ofstream of("intcvs.g2");
	  
	  int numshell = body->nmbOfShells();
	  vector<shared_ptr<ParamCurve> > bd_curves;
	  for (int kd=0; kd<numshell; ++kd)
	    {
	      ftCurve intres = body->getShell(kd)->intersect(plane);

	      // Collect intersection curves and remove basedir coordinate
	      int num = intres.numSegments();
	      for (int ka=0; ka<num; ++ka)
		{
		  shared_ptr<ParamCurve> cv = intres.segment(ka).spaceCurve();
		  shared_ptr<CurveOnSurface> sfcv =
		    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv);
		  if (sfcv)
		    cv = sfcv->spaceCurve();
		  if (cv.get())
		    {
		      cv->writeStandardHeader(of);
		      cv->write(of);
		    }
		  
		  shared_ptr<SplineCurve> spl2;
		  SplineCurve *spl = cv->getSplineCurve();
		  if (!spl)
		    {
		      spl2 = shared_ptr<SplineCurve>(cv->geometryCurve());
		      spl = spl2.get();
		    }

		  int dim = spl->dimension();
		  bool rat = spl->rational();
		  int numc = spl->numCoefs();
		  vector<double> coefs;
		  vector<double>::iterator cf = spl->coefs_begin();
		  for (int kb=0; kb<numc; ++kb)
		    {
		      for (int kc=0; kc<dim; ++kc, cf++)
			{
			  if (kc == basedir)
			    continue;
			  coefs.push_back(*cf);
			  if (rat)
			    {
			      cf++;
			      coefs.push_back(*cf);
			    }
			}
		    }
		  shared_ptr<SplineCurve> bdspl(new SplineCurve(spl->basis(), &coefs[0],
								dim-1, rat));
		  bd_curves.push_back(bdspl);
		}
	    }

	  // Compute 2D quadrature points
	  // NB! What about disjunct curve loops? Are they handled in 2D code?
	  CutCellQuad quad2D(bd_curves, tol_);
	  quad2D.setQuadratureInfo(quadpar_, weights_, min_cell_size_);
	  vector<vector<double> > quadpts;
	  vector<vector<double> > wgts;
	  vector<vector<shared_ptr<ParamCurve> > > unresolved;
	  vector<vector<double> > cvquads;
	  vector<vector<double> > cvwgts;
	  vector<vector<shared_ptr<ParamCurve> > > shortcvs;
	  quad2D.quadrature(ll2, ur2, quadpts, wgts, unresolved, cvquads, cvwgts, shortcvs);

	  // Expand quadrature results with basedir information
	  for (size_t kj=0; kj<quadpts.size(); ++kj)
	    {
	      vector<double> currquads;
	      vector<double> currwgts;
	      size_t kr, kh;
	      for (kr=0, kh=0; kr<quadpts[kj].size(); kr+=2, ++kh)
		{
		  currwgts.push_back(wgts[kj][kh]*wgt);
		  size_t kw = 0;
		  for (int ka=0; ka<3; ++ka)
		    {
		      if (ka == basedir)
			currquads.push_back(quadpar);
		      else
			{
			  currquads.push_back(quadpts[kj][kr+kw]);
			  ++kw;
			}
		    }
		}
	      quadraturepoints.push_back(currquads);
	      pointsweights.push_back(currwgts);

	      ofp << "400 1 0 4 255 0 0 255" << std::endl;
	      ofp << currwgts.size() << std::endl;
	      for (kr=0; kr<currquads.size(); kr+=3)
		{
		  Point pt(currquads[kr], currquads[kr+1], currquads[kr+2]);
		  ofp << pt << std::endl;
		}
	      int stop_break = 1;
	    }
	  
	}
      int stop_break2 = 1;
    }
}

struct edgeInfo
{
  int splitdir_;
  double splitval_;
  double angle_;
  int sfix1_;
  int sfix2_;
  bool convex_;

  edgeInfo()
  {
    splitdir_ = -1;
    sfix1_ = sfix2_ = -1;
    convex_ = false;
  }

  edgeInfo(int splitdir, double splitval, double angle, int sfix1, int sfix2, bool convex)
  {
    splitdir_ = splitdir;
    splitval_ = splitval;
    angle_ = angle;
    sfix1_ = sfix1;
    sfix2_ = sfix2;
    convex_ = convex;
  }
  
};

struct ruledInfo
{
  int ruleddir_;
  int sgn_;
  int sfix1_;
  int sfix2_;
  shared_ptr<SplineCurve> crv_;

  ruledInfo()
  {
    ruleddir_ = -1;
    sgn_ = 1;
    sfix1_ = sfix2_ = -1;
  }
  
  ruledInfo(int ruleddir, int sgn, int sfix1, int sfix2, shared_ptr<SplineCurve> crv)
  {
    ruleddir_ = ruleddir;
    sgn_ = sgn;
    sfix1_ = sfix1;
    sfix2_ = sfix2;
    crv_ = crv;
  }

};
  
//==============================================================================
int
CutCellQuad3D::selectBaseDir(Point& ll, Point& ur,
			     shared_ptr<Body> body,
			     int& splitdir, double& splitval,
			     shared_ptr<SplineCurve>& rulecv,
			     Point& ruledir)
//==============================================================================
{
  splitdir = -1;
  splitval = 0.0;
  
  // Collect surfaces and generate normal cones corresponding to all surfaces
  double angtol = 0.001;
  double anglim = M_PI/6.0;
  int numshells = body->nmbOfShells();
  int numfaces = body->nmbOfFaces();
  vector<shared_ptr<ParamSurface> > sfs(numfaces);
  vector<DirectionCone> ncones(numfaces);
  vector<shared_ptr<DirectionCone> > orthcone[3];
  vector<shared_ptr<DirectionCone> > alongcone[3];
  int curved[3];
  for (int kc=0; kc<3; ++kc)
    {
      orthcone[kc].resize(numfaces);
      alongcone[kc].resize(numfaces);
      curved[kc] = 0;
    }
  int ix = 0;
  for (int ka=0; ka<numshells; ++ka)
    {
      shared_ptr<SurfaceModel> shell = body->getShell(ka);
      int num = shell->nmbEntities();
      for (int kb=0; kb<num; ++kb)
	{
	  sfs[ix] = shell->getSurface(kb);
	  ncones[ix] = sfs[ix]->normalCone();
	  shared_ptr<DirectionCone> tmp1[3], tmp2[3];
	  //sfs[ix]->normalCones(orthcone[0][ix], orthcone[1][ix], orthcone[2][ix]);
	  sfs[ix]->normalCones(tmp1, tmp2);
	  for (int kc=0; kc<3; ++kc)
	    {
	      if (!tmp1[kc].get())
		tmp1[kc] = shared_ptr<DirectionCone>(new DirectionCone(Point(0.0, 0.0, 0.0)));
	      if (!tmp2[kc].get())
		tmp2[kc] = shared_ptr<DirectionCone>(new DirectionCone(Point(0.0, 0.0, 0.0)));
	      orthcone[kc][ix] = tmp1[kc];
	      alongcone[kc][ix] = tmp2[kc];
	      if (orthcone[kc][ix]->angle() > angtol)
		curved[kc]++;
	    }
	  ++ix;
	}
    }

  std::ofstream of("selsfs.g2");
  for (size_t ki=0; ki<sfs.size(); ++ki)
    {
      sfs[ki]->writeStandardHeader(of);
      sfs[ki]->write(of);
    }
  
  // Check the coordinate directions for a suitable base direction
  vector<size_t> cand[3];
  int nonlin[3];
  for (ix=0; ix<3; ++ix)
    {
      nonlin[ix] = 0;
      Point dir(0.0, 0.0, 0.0);
      dir[ix] = 1.0;
      Point dir2 = -dir;

      for (size_t ki=0; ki<ncones.size(); ++ki)
	{
	  // Check normal cone overlap
	  bool overlap1 = alongcone[ix][ki]->containsDirection(dir);
	  bool overlap2 = alongcone[ix][ki]->containsDirection(dir2);
	  if (overlap1 || overlap2)
	    {
	      cand[ix].push_back(ki);
	      if (ncones[ki].angle() > angtol)
		nonlin[ix]++;
	    }
	}

    }

 BoundingBox bbox = body->boundingBox();
 for (int ix=0; ix<3; ++ix)
   {
     ll[ix] = std::max(ll[ix], bbox.low()[ix]);
     ur[ix] = std::min(ur[ix], bbox.high()[ix]);
   }

  // Fetch sharp concave edges
  vector<ftEdge*> convex;
  vector<pair<int,int> > convex_sfs;
  std::ofstream ofconv("convex.g2");
  for (int ka=0; ka<numshells; ++ka)
    {
       shared_ptr<SurfaceModel> shell = body->getShell(ka);
       ModifyFaceSet fset(shell);
       vector<ftEdge*> convex2 = fset.fetchSharpEdges();
       if (convex2.size() > 0)
	 {
	   convex.insert(convex.end(), convex2.begin(), convex2.end());
	   for (size_t ki=0; ki<convex2.size(); ++ki)
	     {
	       shared_ptr<ParamCurve> cv = convex2[ki]->geomCurve();
	       SplineCurve *spl = cv->getSplineCurve();
	       if (spl)
		 {
		   spl->writeStandardHeader(ofconv);
		   spl->write(ofconv);
		 }
	       shared_ptr<ParamSurface> sf1 =
		 convex2[ki]->face()->asFtSurface()->surface();
	       shared_ptr<ParamSurface> sf2 =
		 convex2[ki]->twin()->face()->asFtSurface()->surface();
	       int ix1 = -1, ix2 = -1;
	       for (size_t kj=0; kj<sfs.size(); ++kj)
		 {
		   if (sfs[kj].get() == sf1.get())
		     ix1 = (int)kj;
		   if (sfs[kj].get() == sf2.get())
		     ix2 = (int)kj;
		 }
	       if (ix2 < ix1)
		 std::swap(ix1, ix2);
	       convex_sfs.push_back(make_pair(ix1, ix2));
	     }
	 }
    }
  
  if (nonlin[0] == 0 && nonlin[1] == 0 && nonlin[2] == 0 & convex.size() == 0)
   {
     for (int ix=0; ix<3; ++ix)
       {
	 Point dir(0.0, 0.0, 0.0);
	 dir[ix] = 1.0;
	 Point dir2 = -dir;
	 int nmb_dir = 0;
	 for (size_t ki=0; ki<cand[ix].size(); ++ki)
	   {
	     if (ncones[cand[ix][ki]].containsDirection(dir) ||
		 ncones[cand[ix][ki]].containsDirection(dir2))
	       nmb_dir++;
	   }
       if (nmb_dir == 2)
	 return ix;
       }
   }

  // // Search for simple cases
  // if (nonlin[0] == 0 && nonlin[1] == 0)
  //   {
  //     int d1 = (curved[0] < curved[1]) ? 0 : 1;
  //     int d2 = (curved[0] < curved[1]) ? 1 : 0;
  //     if (curved[d1] == 0 && cand[d2].size() == 2)
  // 	return d2;
  //   }
    
  // if (nonlin[0] == 0 && nonlin[2] == 0)
  //   {
  //     int d1 = (curved[0] < curved[2]) ? 0 : 2;
  //     int d2 = (curved[0] < curved[2]) ? 2 : 0;
  //     if (curved[d1] == 0 && cand[d2].size() == 2)
  // 	return d2;
  //   }
    
  // if (nonlin[1] == 0 && nonlin[2] == 0)
  //   {
  //     int d1 = (curved[1] < curved[2]) ? 1 : 2;
  //     int d2 = (curved[1] < curved[2]) ? 2 : 1;
  //     if (curved[d1] == 0 && cand[d2].size() == 2)
  // 	return d2;
  //   }

  // Count number of connected groups of candidate surfaces in each parameter
  // direction, and check for sharp edges and significant curvature in each group
  tpTolerances tptol = body->getTolerances();
  int numzero = 0;
  for (ix=0; ix<3; ++ix)
    if (nonlin[ix] == 0)
      numzero++;
  if (numzero > 0)
    {
      for (ix=0; ix<3; ++ix)
	{
	  Point dir(0.0, 0.0, 0.0);
	  dir[ix] = 1.0;
	  Point dir2 = -dir;
	  if (nonlin[ix] == 0)
	    {
	      size_t ki;
	      // Check if the side surface coincides with the coordinate directions
	      for (ki=0; ki<cand[ix].size(); ++ki)
		{
		  if (dir.angle(ncones[cand[ix][ki]].centre()) > angtol &&
		      dir2.angle(ncones[cand[ix][ki]].centre()) > angtol)
		    break;
		}
	      if (ki < cand[ix].size())
		continue;
	      
	      vector<shared_ptr<ParamSurface> > cand_sfs;
	      for (size_t ki=0; ki<cand[ix].size(); ++ki)
		cand_sfs.push_back(sfs[cand[ix][ki]]);
	      shared_ptr<SurfaceModel> sfmod(new SurfaceModel(tptol.gap, tptol.gap, tptol.neighbour,
							      tptol.kink, tptol.bend, cand_sfs));
	      vector<shared_ptr<SurfaceModel> > connmod = sfmod->getConnectedModels();
	      if (connmod.size() > 2)
		continue;
	      
	      vector<ftEdge* > edgs;
	      for (size_t ki=0; ki<connmod.size(); ++ki)
		connmod[ki]->getCorners(edgs);
	      if (edgs.size() > 0)
		continue;

	      for (ki=0; ki<sfs.size(); ++ki)
		{
		  size_t kj;
		  for (kj=0; kj<cand[ix].size(); ++kj)
		    {
		      if (cand[ix][kj] == ki)
			break;
		    }
		  if (kj == cand[ix].size())
		    {
		      bool overlap1 = alongcone[ix][ki]->containsDirection(dir);
		      bool overlap2 = alongcone[ix][ki]->containsDirection(dir2);
		      if (overlap1 || overlap2)
			break;
		    }
		}
	      if (ki == sfs.size())
		return ix;   // Should also check direction of edges
	    }
	}
    }

 vector<edgeInfo> split_info;
 vector<ruledInfo> ruled_info;
 std::ofstream edgof("planar_edgs.g2");
 for (ix=0; ix<3; ++ix)
    {
      vector<shared_ptr<ParamSurface> > cand_sfs;
      for (size_t ki=0; ki<cand[ix].size(); ++ki)
	cand_sfs.push_back(sfs[cand[ix][ki]]);
      shared_ptr<SurfaceModel> sfmod(new SurfaceModel(tptol.gap, tptol.gap, tptol.neighbour,
						      tptol.kink, tptol.bend, cand_sfs));
      vector<shared_ptr<SurfaceModel> > connmod = sfmod->getConnectedModels();
      for (size_t ki=0; ki<connmod.size(); ++ki)
	{
	  int numsf = connmod[ki]->nmbEntities();
	  vector<ftEdge* > edgs;
	  connmod[ki]->getCorners(edgs);

	  size_t nedgs = edgs.size();
	  if (convex.size() > 0)
	    edgs.insert(edgs.end(), convex.begin(), convex.end());
	  
	   for (size_t kr=0; kr<edgs.size(); ++kr)
	     {
	       ftEdgeBase *twin = edgs[kr]->twin();
	       if (!twin)
		 continue;

	       // Fetch associated surfaces
	       shared_ptr<ParamSurface> surf1 = edgs[kr]->face()->asFtSurface()->surface();
	       shared_ptr<ParamSurface> surf2 = twin->face()->asFtSurface()->surface();
	       int sfix1=-1, sfix2=-1;
	       for (size_t kj=0; kj<sfs.size(); ++kj)
		 {
		   if (sfs[kj].get() == surf1.get())
		     sfix1 = (int)kj;
		   else if (sfs[kj].get() == surf2.get())
		     sfix2 = (int)kj;
		 }

	       // Compute normal information
	       double t1 = edgs[kr]->tMin();
	       double t2 = edgs[kr]->tMax();
	       int nn = 3;
	       double tdel = (t2 - t1)/(double)(nn);
	       double tpar = t1 + 0.5*tdel;
	       double tang =0.0;
	       for (int ka=0; ka<nn; ++ka, tpar+=tdel)
		 {
		   Point pos1 = edgs[kr]->point(tpar);
		   Point norm1 =edgs[kr]->normal(tpar);
		   double tpar2, dist2;
		   Point pos2;
		   twin->closestPoint(pos1, tpar2, pos2, dist2);
		   Point norm2 = twin->normal(tpar2);
		   double ang = norm1.angle(norm2);
		   if (norm1*norm2 < 0.0)
		     ang = M_PI - ang;
		   tang += ang;
		 }
	       tang /= (double)nn;

	       // Check if current edge is feasible for model split
	       shared_ptr<ParamCurve> crv = edgs[kr]->geomCurve();
	       if (crv.get())
		 {
		   BoundingBox edgebb = crv->boundingBox();
		   BoundingBox sfbox1 = sfs[sfix1]->boundingBox();
		   BoundingBox sfbox2 = sfs[sfix2]->boundingBox();
		   bool planar = false;
		   for (int ix2=0; ix2<3; ++ix2)
		     {
		       double elow = edgebb.low()[ix2];
		       double ehigh = edgebb.high()[ix2];
		       if (ehigh-elow < tol_ &&
			   ur[ix2]-ehigh > tol_ && elow-ll[ix2] > tol_ /*&&
			   (sfbox1.high()[ix2]-ehigh > tol_ || elow-sfbox1.low()[ix2] > tol_) &&
			   (sfbox2.high()[ix2]-ehigh > tol_ || elow-sfbox2.low()[ix2] > tol_)*/)
			 {
			   edgeInfo edginfo(ix2, 0.5*(elow+ehigh),
					    tang, sfix1, sfix2, kr>=nedgs);
			   split_info.push_back(edginfo);
			   planar = true;
			   shared_ptr<SplineCurve> crv2(crv->geometryCurve());
			   crv2->writeStandardHeader(edgof);
			   crv2->write(edgof);
			 }
		     }
		   if (kr >= nedgs && !planar)
		     {
		       std::cout << "Convex non-planar curve identified";
		       // Find ruled direction
		       double minang = M_PI;
		       int minix = -1;
		       int sgn = 0;
		       for (int ix2=0; ix2<3; ++ix2)
			 {
			   Point dir(0.0, 0.0, 0.0);
			   dir[ix2] = 1.0;
			   double ang1 = dir.angle(ncones[sfix1].centre());
			   double ang2 = dir.angle(ncones[sfix2].centre());
			   if (ang1 < minang || M_PI-ang1 < minang)
			     {
			       minang = ang1;
			       minix = ix2;
			       sgn = (ang1 < M_PI-ang1) ? -1 : 1;
			     }
			   if (ang2 < minang || M_PI-ang2 < minang)
			     {
			       minang = ang2;
			       minix = ix2;
			       sgn = (ang2 < M_PI-ang2) ? -1 : 1;
			     }
			 }
		       std::cout << "Ruled direction: " << minix << std::endl;
		       shared_ptr<SplineCurve> rulecv(crv->geometryCurve());
		       shared_ptr<SplineCurve> subrule(rulecv->subCurve(edgs[kr]->tMin(),
									edgs[kr]->tMax()));
		       ruledInfo rinfo(minix, sgn, sfix1, sfix2, subrule);
		       ruled_info.push_back(rinfo);
		     }
		 }
	     }

	  DirectionCone sidecone = connmod[ki]->getSurface(0)->normalCone();
	  for (int ka=1; ka<numsf; ++ka)
	    sidecone.addUnionWith(connmod[ki]->getSurface(ka)->normalCone());

	  int stop_break = 1;
	}
      
    }

 if (ruled_info.size() > 0)
   {
     // Should merge curves if possible or select best curve if more than one exists.
     // Simple first approach
     rulecv = ruled_info[0].crv_;
     ruledir = Point(0.0, 0.0, 0.0);
     ruledir[ruled_info[0].ruleddir_] = 1.0;
     ruledir *= ruled_info[0].sgn_;
     return -1;
   }
 
 // Identify surfaces that occur in more than one direction
 vector<size_t> multiple;
 for (size_t ki=0; ki<sfs.size(); ++ki)
   {
     int nmb = 0;
     for (int ix=0; ix<3; ++ix)
       {
	 for (size_t kj=0; kj<cand[ix].size(); ++kj)
	   if (cand[ix][kj] == ki)
	     nmb++;
       }
     if (nmb > 1)
       multiple.push_back(ki);
   }

 // Check if multiply occuring surfaces have too large curvature
 size_t ki;
 for (ki=0; ki<multiple.size(); ++ki)
   if (ncones[ki].angle() > anglim)
     break;
 if (ki == multiple.size())
   {
     // Associate the multiple surfaces to one direction only
     vector<size_t> cand2[3];
     for (ix=0; ix<3; ++ix)
       cand2[ix].insert(cand2[ix].end(), cand[ix].begin(), cand[ix].end());
     for (ki=0; ki<multiple.size(); ++ki)
       {
	 double minang = std::numeric_limits<double>::max();
	 int minix = -1;
	 for (ix=0; ix<3; ++ix)
	   {
	     Point dir(0.0, 0.0, 0.0);
	     dir[ix] = 1.0;
	     size_t kj;
	     for (kj=0; kj<cand2[ix].size(); ++kj)
	       if (cand2[ix][kj] == multiple[ki])
		 break;

	     if (kj < cand2[ix].size())
	       {
		 double ang = dir.angle(ncones[cand2[ix][kj]].centre());
		 ang = std::min(ang, M_PI-ang);
		 if (ang < minang)
		   {
		     minang = ang;
		     minix = ix;
		   }
	       }
	   }

	 // Remove the multiple surface from the other directions
	 for (ix=0; ix<3; ++ix)
	   {
	     if (ix == minix)
	       continue;
	     for (size_t kj=0; kj<cand2[ix].size(); ++kj)
	       if (cand2[ix][kj] == multiple[ki])
		 {
		   cand2[ix].erase(cand2[ix].begin()+kj);
		   break;
		 }
	   }
       }
 
     // There is potential for a legal height direction. Check.
     // Should also check for edges. Not done
     int ix2 = -1;
     size_t nmb_conn = 0;
     double bbdiff = std::numeric_limits<double>::max();
     Point ll2 = ll;
     Point ur2 = ur;

     for (ix=0; ix<3; ++ix)
       {
	 size_t kj;
	 vector<shared_ptr<ParamSurface> > cand_sfs;
	 double bmin = std::numeric_limits<double>::max();
	 double bmax = std::numeric_limits<double>::lowest();
	 for ( kj=0; kj<cand2[ix].size(); ++kj)
	   {
	     cand_sfs.push_back(sfs[cand2[ix][kj]]);
	     BoundingBox sfbox = sfs[cand2[ix][kj]]->boundingBox();
	     bmin = std::min(bmin, sfbox.low()[ix]);
	     bmax = std::max(bmax, sfbox.high()[ix]);
	   }
	 ll2[ix] = bmin;
	 ur2[ix] = bmax;
	 
	 shared_ptr<SurfaceModel> sfmod(new SurfaceModel(tptol.gap, tptol.gap, tptol.neighbour,
							 tptol.kink, tptol.bend, cand_sfs));
	 vector<shared_ptr<SurfaceModel> > connmod = sfmod->getConnectedModels();
	 if (connmod.size() > 2)
	   continue;  // Not a basedir candidate

	 // Check for equality with basedir direction
	 Point dir(0.0, 0.0, 0.0);
	 dir[ix] = 1.0;
	 for (kj=0; kj<connmod.size(); ++kj)
	   {
	     int nmbsfs = connmod[kj]->nmbEntities();
	     int ka;
	     for (ka=0; ka<nmbsfs; ++ka)
	       {
		 shared_ptr<ParamSurface> surf = connmod[kj]->getSurface(ka);
		 DirectionCone sfcone = surf->normalCone();
		 if (sfcone.angle() > angtol)
		   break;
		 double ang = dir.angle(sfcone.centre());
		 ang = std::min(ang, M_PI-ang);
		 if (ang > angtol)
		   break;
	       }
	     if (ka < nmbsfs)
	       break;
	   }
	 if (kj < connmod.size())
	   continue;
	 
	 double lendiff = 0.0;
	 for (kj=0; kj<connmod.size(); ++kj)
	   {
	     BoundingBox bbmod = sfmod->boundingBox();
	     double bblen = bbmod.low().dist(bbmod.high());
	     lendiff = (kj == 0) ? bblen : fabs(bblen - lendiff);
	   }
		 
	 if (connmod.size() > nmb_conn ||
	     (connmod.size() == nmb_conn && lendiff < bbdiff))
	   {
	     ix2 = ix;
	     nmb_conn = connmod.size();
	     bbdiff = lendiff;
	   }
       }
     if (ix2 >= 0)
       {
	 // Check for identified edges in remaining directions
	 vector<edgeInfo> split_info2;
	 int nmbdir = 0;
	 for (size_t kr=0; kr<split_info.size(); ++kr)
	   {
	     if (split_info[kr].splitdir_ == ix2)
	       nmbdir++;
	     // int ix3 = split_info[kr].splitdir_;
	     for (int ix3=0; ix3<3; ++ix3)
	       {
	     	 if (ix3 == ix2)
	     	   continue;
		 int found = 0;
		 for (size_t kj=0; kj<cand[ix3].size(); ++kj)
		   {
		     if ((int)cand[ix3][kj] == split_info[kr].sfix1_ ||
			 (int)cand[ix3][kj] == split_info[kr].sfix2_)
		       found++;
		   }
		 if (found > 1)
		   split_info2.push_back(split_info[kr]);
	       }
	     if (split_info[kr].convex_)
	       split_info2.push_back(split_info[kr]);
	   }
	 if (split_info2.size() == 0 || nmbdir == 0)
	   {
	     if (cand2[ix2].size() > 1)
	       {
		 ll[ix2] = ll2[ix2];
		 ur[ix2] = ur2[ix2];
	       }
	     return ix2;
	   }
	 split_info = split_info2;
       }
   }
  
 if (split_info.size() > 0)
   {
     // Select the candidate with the largest opening angle
     double minang = split_info[0].angle_;
     size_t minix = 0;
     for (size_t ki=1; ki<split_info.size(); ++ki)
       if (split_info[ki].angle_ < minang || split_info[ki].convex_)
	 {
	   minang = split_info[ki].angle_;
	   minix = ki;
	 }
     splitdir = split_info[minix].splitdir_;
     splitval = split_info[minix].splitval_;
   }

  return -1;  // No suitable direction found
}
