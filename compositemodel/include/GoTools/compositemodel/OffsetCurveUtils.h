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

#ifndef _OFFSETCURVEUTILS_H
#define _OFFSETCURVEUTILS_H


#include "GoTools/geometry/SplineCurve.h"


namespace Go
{

namespace OffsetCurveUtils
{
    /// Top function which returns a smooth offset curve without loops.  param_curve_ assumed to be
    /// without self-intersections.  We expect that input is a loop.  Will usually only return one offset
    /// curve, but might return several (the island scenario).
    std::vector<shared_ptr<SplineCurve> >
    createSmoothOffsetCurves(const SplineCurve& param_curve,
			     const SplineCurve& dist_function,
			     const SplineCurve& tol_function,
			     int nmb_check_pts, bool& last_loop, bool loop = false);


    /// Routine assumes param_loop forms a simple loop. Returns simple offset loops.
    /// Will usually only return one offset curve, but might actually
    /// return several (the island scenario)
    std::vector<shared_ptr<SplineCurve> >
    createSimpleOffsetLoops(const SplineCurve& param_loop,
			    const SplineCurve& dist_function,
			    const SplineCurve& tol_function,
			    std::vector<std::pair<double, double> >& segments_param,
			    int nmb_check_pts,
			    bool sampled_based_offsetting = true);


    /// Returns a spline that is (hopefully) a smoother version of the input offset curve
    /// Might also return a pointer to NULL, if the function was unable to create
    /// a smooth curve of the input curve, in this case OK is set to false.
    SplineCurve*
    createSmoothedOutOffsetCurve(const SplineCurve& offset_curve,
				 const SplineCurve& reference_curve,
				 const SplineCurve& dist_function,
				 const SplineCurve& tol_function,
				 const int nmb_check_pts,
				 const std::vector<std::pair<double, double> >& segments_param,
                                 const double tol, bool& OK);

    // @@sbr201802 Hide all functions that do not need to be public.

    /// Given a 2D-curve param_curve, preferably with no self-intersections, compute
    /// the offset-curve (direction = (0,0,1) X tangent), where the offset-distance
    /// is allowed to vary continuously.
    /// We use the method of Tiller & Hanson: "Offset of two-dimensional profiles",
    /// IEEE Computer Graphics and Applications, Vol. 4, No 9, Sept. 1984, pp.36-46,
    /// with the natural adjustion to our case of varying offset distance (instead of
    /// translating each leg of the control polygon a distance d, we translate the
    /// end points of each leg the corresponding distances).
    /// We demand an exactness within tolerance when sampling the offset curve.
    /// Number of samples is as given by nmb_check_pts.
    /// If curve is a loop, we make sure the same applies for the offset curve.
    // This version of offsetting a curve is based on offsetting the control points
    // of param_curve.
    SplineCurve* createCtrlPntBasedOffsetCurve(const SplineCurve& param_curve,
						 const SplineCurve& dist_function,
						 bool loop = false);

    // This version of offsetting a curve is based on offsetting sampled points from
    // param_curve (# samples given by ctrl pts in dist_function).
    SplineCurve* createSamplePntBasedOffsetCurve(const SplineCurve& param_curve,
						   const SplineCurve& dist_function,
						   bool loop = false);


    /// We offset control point number index in param_curve a distance dist,
    /// in direction given by the two control legs in given point.
    /// If param_curve is a loop, end control points get special treatment.  Instead of
    /// offsetting the legs, we offset in median direction given by ctrl legs.
    Go::Point offsetCtrlPoint(const SplineCurve& param_curve,
				int index, double dist, bool loop = false);

    Go::Point offsetPoint(const SplineCurve& param_curve,
			double tpar, double dist, bool loop);

    /// Approximate the input curve using SISL s1940.
    SplineCurve*
    createSISLApproximatedCurve(const SplineCurve& curve,
				const double epsge);

    // Slightly altered, handling chaotic pattern better, hopefully. Replace above when stable.
    std::vector<std::vector<shared_ptr<Go::SplineCurve> > >
      getSimpleLoops2(const SplineCurve& loop, double loop_tol);

    /// Function erases all loops which are not:
    /// - oriented CCW.
    /// - totally inside bnd_loop.
    /// - within distance given by dist_function and tol_function.
    void removeInvalidLoops(std::vector<std::vector<shared_ptr<Go::SplineCurve> > >&
			    simple_offset_loops,
			    const SplineCurve& bnd_loop,
			    const SplineCurve& dist_function,
			    const SplineCurve& tol_function,
			    double epsge);

} // namespace OffsetCurveUtils

} // namespace Go


#endif // _OFFSETCURVEUTILS_H

