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

#include "GoTools/compositemodel/Body.h"
#include "GoTools/geometry/SplineCurve.h"
#include <fstream>

//#define DEBUG

using std::vector;

namespace Go
{

//---------------------------------------------------------------------------
Body::Body()
  : material_id_(-1), toptol_(0.0, 0.0, 0.0, 0.0)
//---------------------------------------------------------------------------
{
}

// //---------------------------------------------------------------------------
// Body::Body(const CoordinateSystem<3> xyz)
//     : coordinate_(xyz), toptol_(0.0, 0.0, 0.0, 0.0)
// //---------------------------------------------------------------------------
// {
// }

// //---------------------------------------------------------------------------
// Body::Body(const CoordinateSystem<3> xyz, vector<shared_ptr<SurfaceModel> >& shells)
//     : coordinate_(xyz), shells_(shells), toptol_(0.0, 0.0, 0.0, 0.0)
// //---------------------------------------------------------------------------
// {
//     double gap=0.0, neighbour=0.0, kink=0.0, bend=0.0;
//     for (size_t ki=0; ki<shells.size(); ki++)
//     {
// 	tpTolerances toptol = shells[ki]->getTolerances();
// 	gap = std::max(toptol.gap, gap);
// 	neighbour = std::max(toptol.neighbour, neighbour);
// 	kink = std::max(toptol.kink, kink);
// 	bend = std::max(toptol.bend, bend);
//     }
//     toptol_ = tpTolerances(gap, neighbour, kink, bend);

//     addBodyPointers();
// }

// //---------------------------------------------------------------------------
// Body::Body(const CoordinateSystem<3> xyz, shared_ptr<SurfaceModel>  shell)
//     : coordinate_(xyz), toptol_(shell->getTolerances())
// //---------------------------------------------------------------------------
// {
//     shells_.push_back(shell);
//     addBodyPointers();
// }

  Body::Body(vector<shared_ptr<SurfaceModel> >& shells, int material_id)
    : shells_(shells), material_id_(material_id), toptol_(0.0, 0.0, 0.0, 0.0)
//---------------------------------------------------------------------------
{
    double gap=0.0, neighbour=0.0, kink=0.0, bend=0.0;
    for (size_t ki=0; ki<shells.size(); ki++)
    {
	tpTolerances toptol = shells[ki]->getTolerances();
	gap = std::max(toptol.gap, gap);
	neighbour = std::max(toptol.neighbour, neighbour);
	kink = std::max(toptol.kink, kink);
	bend = std::max(toptol.bend, bend);
    }
    toptol_ = tpTolerances(gap, neighbour, kink, bend);

    addBodyPointers();
}

//---------------------------------------------------------------------------
  Body::Body(shared_ptr<SurfaceModel>  shell, int material_id)
    : material_id_(material_id), toptol_(shell->getTolerances())
//---------------------------------------------------------------------------
{
    shells_.push_back(shell);
    addBodyPointers();
}

  Body::Body(shared_ptr<Body> body)
    : shells_(body->getAllShells()), material_id_(body->getMaterial()), 
      toptol_(0.0, 0.0, 0.0, 0.0)
//---------------------------------------------------------------------------
{
    double gap=0.0, neighbour=0.0, kink=0.0, bend=0.0;
    for (size_t ki=0; ki<shells_.size(); ki++)
    {
	tpTolerances toptol = shells_[ki]->getTolerances();
	gap = std::max(toptol.gap, gap);
	neighbour = std::max(toptol.neighbour, neighbour);
	kink = std::max(toptol.kink, kink);
	bend = std::max(toptol.bend, bend);
    }
    toptol_ = tpTolerances(gap, neighbour, kink, bend);

    addBodyPointers();
}

//---------------------------------------------------------------------------
Body::~Body()
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void Body::addshell(shared_ptr<SurfaceModel> shell)
//---------------------------------------------------------------------------
{
    shells_.push_back(shell);
    if (shells_.size() == 1)
    {
	toptol_ = shell->getTolerances();
    }
    else
    {
	double gap=0.0, neighbour=0.0, kink=0.0, bend=0.0;
	tpTolerances toptol = shell->getTolerances();
	gap = std::max(toptol_.gap, toptol.gap);
	neighbour = std::max(toptol_.neighbour, toptol.neighbour);
	kink = std::max(toptol_.kink, toptol.kink);
	bend = std::max(toptol_.bend, toptol.bend);

	toptol_ = tpTolerances(gap, neighbour, kink, bend);
    }
    addBodyPointers();
}

//---------------------------------------------------------------------------
BoundingBox Body::boundingBox() const
//---------------------------------------------------------------------------
{
  if (shells_.size() > 0)
    return shells_[0]->boundingBox();
  else
    {
      BoundingBox dummy;
      return dummy;
    }
}

//---------------------------------------------------------------------------
int Body::nmbOfFaces() const
//---------------------------------------------------------------------------
{
  int nmb = 0;
  for (size_t ki=0; ki<shells_.size(); ++ki)
    {
      nmb += shells_[ki]->nmbEntities();
    }
  return nmb;
}

//---------------------------------------------------------------------------
vector<shared_ptr<Vertex> > Body::vertices() const
//---------------------------------------------------------------------------
{
  // Collect all vertices
 vector<shared_ptr<Vertex> > vertices;
 std::set<shared_ptr<Vertex> > all_vertices;  // All vertices in the model represented once

  for (size_t ki=0; ki<shells_.size(); ++ki)
    {
      vector<shared_ptr<Vertex> > curr_vertices;
      shells_[ki]->getAllVertices(curr_vertices);
      all_vertices.insert(curr_vertices.begin(), curr_vertices.end());
    }

  vertices.insert(vertices.end(), all_vertices.begin(), all_vertices.end());
  return vertices;
}


//---------------------------------------------------------------------------
bool Body::areNeighbours(Body *other, shared_ptr<ftSurface>& bd_face1,
			 shared_ptr<ftSurface>& bd_face2, int adj_idx) const
//---------------------------------------------------------------------------
{
  for (size_t ki=0; ki<shells_.size(); ++ki)
    {
      int nmb_faces = shells_[ki]->nmbEntities();
      for (int kj=0; kj<nmb_faces; ++kj)
	{
	  shared_ptr<ftSurface> curr = shells_[ki]->getFace(kj);
	  if (!curr->twin())
	    continue;

	  for (size_t kr=0; kr<other->shells_.size(); ++kr)
	    {
	      int nmb_faces2 = other->shells_[kr]->nmbEntities();
	      for (int kh=0; kh<nmb_faces2; ++kh)
		{
		  shared_ptr<ftSurface> curr2 = other->shells_[kr]->getFace(kh);
		  if (!curr2->twin())
		    continue;

		  if (curr->twin() == curr2.get() && curr2->twin() == curr.get())
		    {
		      // Decrement local adjecency index, return when desired index is found
		      if (adj_idx-- > 0)
			continue;
		      bd_face1 = curr;
		      bd_face2 = curr2;
		      return true;
		    }
		}
	    }
	}
    }
  return false;
}

//---------------------------------------------------------------------------
void Body::eraseBodyAdjacency()
//---------------------------------------------------------------------------
{
  for (size_t ki=0; ki<shells_.size(); ++ki)
    {
      int nmb_faces = shells_[ki]->nmbEntities();
      for (int kj=0; kj<nmb_faces; ++kj)
	{
	  shared_ptr<ftSurface> curr = shells_[ki]->getFace(kj);
	  if (curr->hasTwin())
	    curr->disconnectTwin();
	}
    }
}

//---------------------------------------------------------------------------
  bool Body::isInside(const Point& pnt) const
//---------------------------------------------------------------------------
{
  double tol = 1.0e-9;  // 1.0e-10; This tolerance is very significant regarding
  // whether or not the correct result is returned. Should think about a more robus
  // solution

  // Candidate points for curve creation
  vector<Point> cand;

  // Fetch the midpoint of the bounding box
  BoundingBox box = boundingBox();
  Point mid = 0.5*(box.low() + box.high());
  if (mid.dist(pnt) > toptol_.neighbour)
    cand.push_back(mid);
  cand.push_back(box.low());
  cand.push_back(box.high());

  // Make a curve through the input point and this midpoint
  for (size_t kr=0; kr<cand.size(); ++kr)
    {
      Point vec = pnt - cand[kr];
      vec.normalize();
      Point vec2 = box.high() - box.low();
      double len = vec2.length();

      Point start = pnt - 2.0*len*vec;
      Point end = pnt + 2.0*len*vec;
      shared_ptr<SplineCurve> crv = 
	shared_ptr<SplineCurve>(new SplineCurve(start, -2.0*len, end, 2.0*len));

      // Intersect all boundary shells with the curve
      vector<bool> segment;
      vector<pair<ftPoint, double> > int_pts;
      size_t ki, kj;
      for (ki=0; ki<shells_.size(); ++ki)
	{
#ifdef DEBUG
	  std::ofstream of("insidetest.g2");
	  int nmb = shells_[ki]->nmbEntities();
	  for (kj=0; kj<nmb; ++kj)
	    {
	      shared_ptr<ParamSurface> sf = shells_[ki]->getSurface(kj);
	      sf->writeStandardHeader(of);
	      sf->write(of);
	    }
	  crv->writeStandardHeader(of);
	  crv->write(of);
	  of << "400 1 0 4 255 0 0 255" << std::endl;
	  of << "1" << std::endl;
	  of << pnt << std::endl;
#endif
	  vector<bool> seg0;
	  vector<pair<ftPoint, double> > int_pts0 = 
	    shells_[ki]->intersect(crv, seg0, false);
	  int_pts.insert(int_pts.end(), int_pts0.begin(), int_pts0.end());
	  segment.insert(segment.end(), seg0.begin(), seg0.end());
	}

      // Remove touch points
      for (ki=0; ki<int_pts.size();)
	{
	  if (segment[ki])
	    {
	      int_pts.erase(int_pts.begin()+ki);
	      segment.erase(segment.begin()+ki);
	    }
	  else
	    ki++;
	}

      // Remove duplicates
      for (ki=0; ki<int_pts.size(); ++ki)
	for (kj=ki+1; kj<int_pts.size(); )
	  {
	    if (fabs(int_pts[ki].second - int_pts[kj].second) < tol)
	      int_pts.erase(int_pts.begin()+kj);
	    else
	      kj++;
	  }


      // Count number of intersection points on either side of the
      // input point
      int nmb1, nmb2;
      for (nmb1=0, nmb2=0, ki=0; ki<int_pts.size(); ++ki)
	{
	  if (int_pts[ki].second < -tol)
	    nmb1++;
	  if (int_pts[ki].second > tol)
	    nmb2++;
	}

      if (nmb1+nmb2 < (int)int_pts.size())
	return true;  // On boundary

      if (nmb1 % 2 == 0 && nmb2 % 2 == 0)
	return false;  // Outside
      else if (nmb1 % 2 == 1 && nmb2 % 2 == 1)
	return true;
    }
  return false;
}

//---------------------------------------------------------------------------
  void Body::addBodyPointers()
//---------------------------------------------------------------------------
{
  for (size_t ki=0; ki<shells_.size(); ++ki)
    {
      int nmb = shells_[ki]->nmbEntities();
      for (int kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ftSurface> face = shells_[ki]->getFace(kj);
	  face->setBody(this);
	}
    }
}

//---------------------------------------------------------------------------
  bool Body::isInside(const Point& pnt, double& dist, double& ang) const
//---------------------------------------------------------------------------
{
  bool inside = isInside(pnt);
  dist = std::numeric_limits<double>::max();
  ang = M_PI;
  for (size_t ki=0; ki<shells_.size(); ++ki)
    {
      Point pnt1 = pnt;
      Point clo_pnt;
      int idx;
      double par[2];
      double dist2;
      shells_[ki]->closestPoint(pnt1, clo_pnt, idx, par, dist2);
      if (dist2 < dist)
	{
	  dist = dist2;
	  // if (dist < 1.0e-6)
	  //   ang = 0.0;
	  // else
	    // {
	      Point vec = pnt - clo_pnt;
	      Point norm = shells_[ki]->getFace(idx)->normal(par[0],par[1]);
	      ang = vec.angle(norm);
	    // }
	}
    }

  return inside;
}

//---------------------------------------------------------------------------
  shared_ptr<SurfaceModel> Body::getShell(ftSurface* face) const
//---------------------------------------------------------------------------
{
  shared_ptr<SurfaceModel> model;
  for (size_t ki=0; ki<shells_.size(); ++ki)
    {
      int nmb = shells_[ki]->nmbEntities();
      int kj;
      for (kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ftSurface> curr = shells_[ki]->getFace(kj);
	  if (curr.get() == face)
	    {
	      model = shells_[ki];
	      break;
	    }
	}
      if (kj < nmb)
	break;
    }
  return model;
}

} // namespace Go
