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


#include "GoTools/compositemodel/OffsetCurveUtils.h"
#include "GoTools/implicitization/ImplicitizeCurveAlgo.h"
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/geometry/CurveBoundedDomain.h"
#include "GoTools/geometry/LoopUtils.h"
#include "GoTools/utils/errormacros.h"
#include "GoTools/creators/CurveCreators.h"
#include "GoTools/geometry/GoIntersections.h"
#include "sislP.h"

#if 0
#include "SelfIntersectCurveAlgo.h"
#endif

#include <vector>
#include <utility>
#include <algorithm>

using std::vector;
using std::pair;
using std::make_pair;
using std::max;
using std::min;

namespace Go
{


//===========================================================================
vector<shared_ptr<SplineCurve> >
OffsetCurveUtils::createSmoothOffsetCurves(const SplineCurve& param_curve,
					  const SplineCurve& dist_function,
					  const SplineCurve& tol_function,
					  int nmb_check_pts, bool& last_loop,
					  bool loop)
//===========================================================================
{
#ifndef NDEBUG
    {
	double epspar = 1e-10;
	DEBUG_ERROR_IF((fabs(param_curve.startparam() - dist_function.startparam())> epspar)
		    ||
		    (fabs(param_curve.endparam() - dist_function.endparam()) > epspar) ||
		    (fabs(param_curve.startparam() - tol_function.startparam()) > epspar)
		    ||
		    (fabs(param_curve.endparam() - tol_function.endparam()) > epspar),
		    "HeroOffsetUtils::createSmoothOffsetCurves(): Parameter values NOK");
    }
#endif // NDEBUG

    // As we intend to make a smooth offset curve, we demand order >= 4 (cubic).
    // As curve is to be constant, we make a copy.
    SplineCurve pcurve = param_curve;
    int order = pcurve.order();
    if (order < 4)
	pcurve.raiseOrder(4 - order);

#ifndef NDEBUG
    {
	double epsgeo = 1e-05;
	vector<shared_ptr<SplineCurve> > cvs(1);
	cvs[0] = shared_ptr<SplineCurve>(pcurve.clone());
	bool ccw = LoopUtils::loopIsCCW(cvs, epsgeo);
	if (!ccw)
	    MESSAGE("Input loop is not ccw!");
    }
#endif // NDEBUG

    // vector containing parameter intervals over which the remaining curve segments
    // are defined after cutting off loops resulting from the offset strategy
    vector<pair<double, double> > segments_param;
    // We offset the curve. We make sure the resulting curve is simple (and closed
    // if loop).
    vector<shared_ptr<SplineCurve> > offset_curves;
    try {
	offset_curves = createSimpleOffsetLoops(pcurve, dist_function, tol_function,
						segments_param, nmb_check_pts, true);
	if (offset_curves.size() == 0)
	    last_loop = true;
    } catch (...) {
	// Seems the first approach did not succeed, trying the old way.
	offset_curves = createSimpleOffsetLoops(pcurve, dist_function, tol_function,
						segments_param, nmb_check_pts, false);
	if (offset_curves.size() == 0)
	    last_loop = true;
    }

#ifndef NDEBUG
    {
	std::ofstream curveout("data/bnd_loop.g2");
	writeSpaceParamCurve(param_curve, curveout);
	std::ofstream debug("data/debug.g2");
	for (size_t i = 0; i < offset_curves.size(); ++i)
	    writeSpaceParamCurve(*offset_curves[i], debug);
    }
#endif // NDEBUG

    // for the approximation of the curve we use a fraction of the smallest value
    // from the tol_function to be sure to be within limits.
    double min_tol = 1e+37;    // some very large value
    for (vector<double>::const_iterator it=tol_function.coefs_begin();
	 it!=tol_function.coefs_end(); ++it) {
	if (*it<min_tol) min_tol = *it;
    }
    //const HeroParameters& hero_parameters = HeroParameters::Instance();
    // "playground" of new curve.
    const double smoothing_ratio = 0.15;//hero_parameters.getSmoothingRatio();
    min_tol *= smoothing_ratio;
    for (size_t i=0; i<offset_curves.size(); ++i) {
	shared_ptr<SplineCurve> smooth_offset_curve;
	bool smoothing_OK = false;
	if (!smoothing_OK) {
	    //    const double max_curvature = 1.0;    // example value
	    try {
		smooth_offset_curve = shared_ptr<SplineCurve>
		    (createSmoothedOutOffsetCurve(*offset_curves[i],
						  param_curve, dist_function,
						  tol_function, nmb_check_pts,
						  segments_param, min_tol,
						  smoothing_OK));
	    } catch (...) {
		MESSAGE("Something went wrong when trying to smooth out "
			   "offset curve. This should not happen!");
		continue; // We move on to next loop.
	    }
	}
	if (smoothing_OK) {
	    assert(smooth_offset_curve.get()!=NULL);
	    offset_curves[i] = smooth_offset_curve;
	}
    }

    return offset_curves;
}



//===========================================================================
vector<shared_ptr<SplineCurve> >
OffsetCurveUtils::createSimpleOffsetLoops(const SplineCurve& param_loop,
					 const SplineCurve& dist_function,
					 const SplineCurve& tol_function,
					 vector<pair<double, double> >&
					 segments_param,
					 int nmb_check_pts,
					 bool sampled_based_offsetting)
//===========================================================================
{
//   const HeroTolerances& hero_tol = HeroTolerances::Instance();
//	setNumericalTolerance(1e-9);
//	setGeometricalTolerance(1e-5);
    double num_tol = 1.0e-09;//hero_tol.getNumericalTolerance();
    double epsgeo = 1.0e-05;//hero_tol.getGeometricalTolerance();

    DEBUG_ERROR_IF((fabs(param_loop.startparam() - dist_function.startparam()) >
		 num_tol) ||
		(fabs(param_loop.endparam() - dist_function.endparam()) >
		 num_tol) ||
		(fabs(tol_function.startparam() - dist_function.startparam()) >
		 num_tol) ||
		(fabs(tol_function.endparam() - dist_function.endparam()) > num_tol),
		"Curves' parameter intervals do not match!");
    DEBUG_ERROR_IF(param_loop.dimension() != 2,
		"Curve assumed to be 2-dimensional.");
    DEBUG_ERROR_IF((dist_function.dimension() != 1) || (tol_function.dimension() != 1),
		"'dist_function' & 'tol_function' assumed to be 1-dimensional.");

    // As we most likely will have to insert knots into the curve, we make a copy.
    shared_ptr<SplineCurve> pcurve(new SplineCurve(param_loop));
    vector<vector<shared_ptr<SplineCurve> > > simple_offset_loops;
//     vector<shared_ptr<CurveLoop> > prev_offset_loops;

#ifndef NDEBUG
    {
	vector<shared_ptr<SplineCurve> > cvs(1);
	cvs[0] = pcurve;
	bool ccw = LoopUtils::loopIsCCW(cvs, epsgeo);
	if (!ccw)
	    MESSAGE("Input loop is not ccw!");
    }
#endif // NDEBUG

    int max_iter = 15; // If no success, caller should try with smaller offset distance.
    int iter = 0;
    bool tolerance_achieved = false;
    bool loop = true;
    size_t i, j;

    shared_ptr<SplineCurve> offset_loop;
    if (sampled_based_offsetting)
	offset_loop = // We sample points of pcurve, offset, and return linear interpolation.
	    shared_ptr<SplineCurve>(createSamplePntBasedOffsetCurve(*pcurve, dist_function, loop));
    else
	offset_loop = // We offset the control points of previous loop.
	    shared_ptr<SplineCurve>(createCtrlPntBasedOffsetCurve(*pcurve, dist_function, loop));

#ifndef NDEBUG
    {
	std::ofstream curveout2("data/appr_offset_loop.g2");
	writeSpaceParamCurve(*pcurve, curveout2);
	writeSpaceParamCurve(*offset_loop, curveout2);
    }
#endif // NDEBUG

    // Given closed offset curve, extract simple loops.
    double loop_tol = epsgeo; //hero_tol.getIGESFileTolerance();
    try {
	simple_offset_loops = getSimpleLoops2(*offset_loop, loop_tol);
    } catch (...) {
	// When knot spacing of loop gets dense, the self-intersection routine
	// may become unstable. This may cause multiples of a parameter, which
	// should not occur.
	THROW("Failed dividing offset curve into loops.");
    }

#ifndef NDEBUG
    {
	std::ofstream curveout4("data/simple_offset_loops.g2");
	for (i = 0; i < simple_offset_loops.size(); ++i)
	    for (j = 0; j < simple_offset_loops[i].size(); ++j)
		writeSpaceParamCurve(*(simple_offset_loops[i][j]),
				     curveout4);
    }
#endif // NDEBUG

    // If loop is CCW, inside boundary and fulfill dist demand, it's a keeper!
    removeInvalidLoops(simple_offset_loops, param_loop,
		       dist_function, tol_function, epsgeo);

#ifndef NDEBUG
    {
	std::ofstream curveout3("data/cut_offset_loop.g2");
	for (i = 0; i < simple_offset_loops.size(); ++i)
	    for (j = 0; j < simple_offset_loops[i].size(); ++j)
		writeSpaceParamCurve(*(simple_offset_loops[i][j]),
				     curveout3);
    }
#endif // NDEBUG

    if (sampled_based_offsetting) { // As offset curve is piecewise linear,we must try to approximate.
	// We approximate each segment by a cubic spline curve.
	for (i = 0; i < simple_offset_loops.size(); ++i) {
	    // In the future we probably should approximate more than one segment,
	    // given appr 180 angle.
	    for (j = 0; j < simple_offset_loops[i].size(); ++j) {
		double tmin = simple_offset_loops[i][j]->startparam();
		double tmax = simple_offset_loops[i][j]->endparam();
		int min_ind =
		    tol_function.basis().knotIntervalFuzzy(tmin) + tol_function.order() - 1;
		int max_ind =
		    tol_function.basis().knotIntervalFuzzy(tmax) + tol_function.order() - 1;
		vector<double>::const_iterator min_elem =
		    min_element(tol_function.coefs_begin(), tol_function.coefs_end());
		double min_tol = *min_elem;
		try {
		    int max_iter = 10; // @@sbr Hardcoded!
		    vector<Point> start_pt = simple_offset_loops[i][j]->
			ParamCurve::point(simple_offset_loops[i][j]->startparam(), 1);
		    start_pt[1].normalize();
		    vector<Point> end_pt = simple_offset_loops[i][j]->
			ParamCurve::point(simple_offset_loops[i][j]->endparam(), 1);
		    end_pt[1].normalize();
		    start_pt.clear(); // @@sbr We should include end pts.
		    end_pt.clear();
		    double appr_length =
			simple_offset_loops[i][j]->ParamCurve::estimatedCurveLength(100);
		    if (appr_length < min_tol) {
			MESSAGE("Segments smaller than tolerance, using linear interpolant.");
			continue;
		    }
// 		    double tol = min(min_tol, appr_length/10.0); // To great tolerance leads to trouble.
		    double max_dist;
                    vector<shared_ptr<ParamCurve> > param_curves(simple_offset_loops[i].begin() + j,
                                                                  simple_offset_loops[i].begin() + j + 1);

		    shared_ptr<SplineCurve> appr_cv
			(CurveCreators::approxCurves(&param_curves[0],
						      &param_curves[param_curves.size()],
                                                     start_pt, end_pt, min_tol, max_dist, max_iter));
		    if (max_dist > min_tol) {
			MESSAGE("Failed approximating within tolerance.");
		    }
#ifndef NDEBUG
		    {
			std::ofstream debug("data/debug.g2");
			writeSpaceParamCurve(*simple_offset_loops[i][j], debug);
			writeSpaceParamCurve(*appr_cv, debug);
		    }
#endif // NDEBUG

		    simple_offset_loops[i][j] = appr_cv;
		    // @@sbr Check for self intersections!
		} catch (...) {
// 		    ERROR("Failed approximating offset curve.", UnknownError());
		    // @@sbr This is not that ideal, but since curve is to be sampled we do it this way.
		    MESSAGE("Using linear interpolat for segment.");
		}
	    }
	}
    }

#ifndef NDEBUG
    {
	std::ofstream curveout32("data/cut_offset_loop2.g2");
	for (i = 0; i < simple_offset_loops.size(); ++i)
	    for (j = 0; j < simple_offset_loops[i].size(); ++j)
		writeSpaceParamCurve(*(simple_offset_loops[i][j]),
				     curveout32);
    }
#endif // NDEBUG

    // For each loop, we make a continuous spline curve.
    vector<shared_ptr<SplineCurve> > return_loops;
    double dist;
    for (i = 0; i < simple_offset_loops.size(); ++i) {
	shared_ptr<SplineCurve> curve = simple_offset_loops[i][0];
	segments_param.push_back
	    (pair<double,double>(curve->startparam(), curve->endparam()));
	for (j = 1; j < simple_offset_loops[i].size(); ++j) {
	    // remember which parts of the offset curve that have been removed
	    segments_param.push_back
		(pair<double,double>(simple_offset_loops[i][j]->startparam(),
				     simple_offset_loops[i][j]->endparam()));
	    curve->appendCurve((simple_offset_loops[i][j]).get(), 0, dist);
	}
	return_loops.push_back(curve);
    }


#ifndef NDEBUG
    {
	std::ofstream curveout5("data/approved_offset_loop.g2");
	for (i = 0; i < return_loops.size(); ++i)
	    writeSpaceParamCurve(*(return_loops[i]), curveout5);
    }
#endif // NDEBUG

    return return_loops;
}



//===========================================================================
SplineCurve*
OffsetCurveUtils::createSmoothedOutOffsetCurve(const SplineCurve& offset_curve,
					      const SplineCurve& reference_curve,
					      const SplineCurve& dist_function,
					      const SplineCurve& tol_function,
					      const int nmb_check_pts,
					      const vector<pair<double, double> >& segments_param,
					      const double tol, bool& OK)
//===========================================================================
    // @ers: TODO: max_curvature not used
    // returns NULL or a curve
{
#ifndef NDEBUG
    {
	cout << "s1940: tolerance = " << tol << endl;
    }
#endif // NDEBUG

    // @@sbr Assumed to be OK, as offset routine was based on sampled points.
//     ERROR_IF(!isCurveWithinTolerance(offset_curve, reference_curve, dist_function,
// 					tol_function, nmb_check_pts, segments_param),
// 		"Input offset curve is not within tolerances!",
// 		InputError());

    OK = true;
    shared_ptr<SplineCurve> current_curve(offset_curve.clone());

    const int old_no_coefs = current_curve->numCoefs();   // @ers: debug only

    shared_ptr<SplineCurve> approx
	(createSISLApproximatedCurve(*current_curve, tol));

    vector<double> int_par;
//     vector<int> orders;

#if 0
    ImplicitizeCurveAlgo self_int_obj;
    self_int_obj.useSplineCurve(*approx);
    self_int_obj.perform();
    self_int_obj.getResultData(int_par); //, orders);
//     find_self_intersections(*approx, int_par);
#else
    MESSAGE("Removed call to ImplicitizeCurveAlgo!");
#endif    

    if (int_par.size()!=0) {
// 	// must either remove self intersections, or not use the approximated curve
// 	MESSAGE("Curve smoothing: smoothed curve had self intersections, "
// 		   "not used.");
	OK = false;
	return NULL;
    }
	
//     // Now check that we are still within given tolerance
//     OK = isCurveWithinTolerance(*approx, reference_curve, dist_function,
// 				tol_function, nmb_check_pts, segments_param);
//     if (!OK){
// 	MESSAGE("Curve smoothing: smoothed curve not within tolerance, "
// 		   "not used.");
// 	// OK is already false
// 	return NULL;
//     }

    current_curve = approx;

#ifndef NDEBUG
    {
	const int new_no_coefs = current_curve->numCoefs();   // @ers: debug only
	cout << "Curve smoothing: Removed " << old_no_coefs-new_no_coefs <<
	    " coefficients";
    }
#endif // NDEBUG
    /*
      // @ers: here i am, remove knots

    // ensure that knots are max (d-2)-tupple, except at start/end
    shared_ptr<SplineCurve>
	C2_curve(createC2ContinuousLoopCurve(offset_curve, reference_curve,
					     // C2_curve(createC2ContinuousLoopCurve(approx.get(), reference_curve,
					     dist_function,
					     tol_function, nmb_check_pts,
					     segments_param));

    */



    /*
      // use max_curvature to smooth out curve, not implemented
    // get curvature everywhere on curve
    const double start_param = offset_curve.startparam();
    const double end_param = offset_curve.endparam();
    const double param_step = (end_param - start_param) / (nmb_check_pts - 1);

    vector<double> params;
    vector<double> curvature;
    vector<double> dist;
    vector<double> tol;

    for (int i=0; i<nmb_check_pts; ++i) {
	const double tpar = start_param + i*param_step;
	params.push_back(tpar);
	Point dist_pt;
	Point tol_pt;
	
	curvature.push_back(HeroUtils::curvature(offset_curve, tpar));

	dist_function.point(dist_pt, tpar);
	dist.push_back(dist_pt[0]);
	tol_function.point(tol_pt, tpar);
	dist.push_back(tol_pt[0]);
    }
    */


    return current_curve->clone();    // curve goes out of scope
}



//===========================================================================
SplineCurve*
OffsetCurveUtils::createCtrlPntBasedOffsetCurve(const SplineCurve& pcurve,
					       const SplineCurve& dist_function,
					       bool loop)
//===========================================================================
{
    // @@sbr An idea could be to insert knots near both sides of corner (if not already
    //       present). Will make the offset curve near the corner more 'true'.
    //       Not necessary in an area of high curvature as local ctrl pts already are dense.

    Point new_ctrl_pt;

    int dim = pcurve.dimension();
    int num_coefs = pcurve.numCoefs();
    vector<double> new_coefs(num_coefs*dim);
    Point dist;
    for (int i = 0; i < num_coefs; ++i) {
	double tpar = pcurve.basis().grevilleParameter(i);
	dist_function.point(dist, tpar);
	new_ctrl_pt = offsetCtrlPoint(pcurve, i, dist[0], false); //loop);
	for (int j = 0; j < dim; ++j)
	    new_coefs[i*dim+j] = new_ctrl_pt[j];
    }
    if (loop) // The pcurve may be a loop, while being not totally closed.
	for (int j = 0; j < dim; ++j) // We make sure the loop is totally closed.
	    new_coefs[j] = new_coefs[(num_coefs-1)*dim+j]
		= (new_coefs[j] + new_coefs[(num_coefs-1)*dim+j])/2;

    SplineCurve* offset_curve =
	new SplineCurve(pcurve.numCoefs(), pcurve.order(), pcurve.basis().begin(),
			  new_coefs.begin(), pcurve.dimension());

    return offset_curve;
}


//===========================================================================
SplineCurve*
OffsetCurveUtils::createSamplePntBasedOffsetCurve(const SplineCurve& pcurve,
						 const SplineCurve& dist_function,
						 bool loop)
//===========================================================================
{
    DEBUG_ERROR_IF(dist_function.dimension() != 1,
		"Input dist function must be 1-dimensional.");

    int dim = pcurve.dimension();
    int nmb_samples = dist_function.numCoefs();
    Point dist, offset_pt;
    // For each Greville parameter in dist_function we offset a point.
    vector<double> offset_params, offset_points;
    offset_params.push_back(dist_function.basis().grevilleParameter(0));
    for (int i = 0; i < nmb_samples; ++i) {
	double tpar = dist_function.basis().grevilleParameter(i);
	dist_function.point(dist, tpar);
	offset_pt = offsetPoint(pcurve, tpar, dist[0], false); //loop);
	offset_params.push_back(tpar);
	offset_points.insert(offset_points.end(), offset_pt.begin(), offset_pt.end());
    }
    offset_params.push_back(dist_function.basis().grevilleParameter(nmb_samples));
    if (loop) // The pcurve may be a loop, while being not totally closed.
	for (int j = 0; j < dim; ++j) // We make sure the loop is totally closed.
	    offset_points[j] = offset_points[(nmb_samples-1)*dim+j]
		= 0.5*(offset_points[j] + offset_points[(nmb_samples-1)*dim+j]);

    int offset_num_coefs = dist_function.numCoefs();
    int offset_order = 2; // Linear interpolation.
    int offset_dim = pcurve.dimension();
    SplineCurve* offset_curve =
	new SplineCurve(offset_num_coefs, offset_order, offset_params.begin(),
			  offset_points.begin(), offset_dim);

    return offset_curve;
}


//===========================================================================
Point
OffsetCurveUtils::offsetCtrlPoint(const SplineCurve& param_curve,
				 int index, double dist, bool loop)
//===========================================================================
{
    int num_coefs = param_curve.numCoefs();
    // When using normal to curve, evaluated in Greville parameter, bad parametrization
    // will, to a certain extent, be fixed.
    // Currently direction is computed from the curve itself, not the neighbour control legs.
    double tpar = param_curve.basis().grevilleParameter(index);
    vector<Point> eval = param_curve.ParamCurve::point(tpar, 1);
    double tpar_left = (loop && index == 0) ? param_curve.basis().grevilleParameter(num_coefs - 1) :
	param_curve.basis().grevilleParameter(index);
    double tpar_right = (loop && index == num_coefs - 1) ? param_curve.basis().grevilleParameter(0) :
	param_curve.basis().grevilleParameter(index);
    vector<Point> eval_left = param_curve.ParamCurve::point(tpar_left, 1, false);
    vector<Point> eval_right = param_curve.ParamCurve::point(tpar_right, 1, true);
    Point normal_left(-eval_left[1][1], eval_left[1][0]);
    normal_left.normalize();
    Point normal_right(-eval_right[1][1], eval_right[1][0]);
    normal_right.normalize();
//     Point normal(-eval[1][1], eval[1][0]);
//     normal.normalize();
    // Average tends to cause trouble, better to let user insert knots.
//     average_normal = normal;
    Point average_normal = normal_left + normal_right;
    average_normal.normalize();

//     Point return_point = center + dist*average_normal;
    Point return_point = eval_right[0] + dist*average_normal;
    return return_point;
}




//===========================================================================
Point OffsetCurveUtils::offsetPoint(const SplineCurve& param_curve,
				     double tpar, double dist, bool loop)
//===========================================================================
{
    // When using normal to curve, bad parametrization will (to a certain extent) be fixed.
    vector<Point> eval = param_curve.ParamCurve::point(tpar, 1);
    double tpar_left = (tpar == param_curve.startparam()) ? param_curve.endparam() : tpar;
    double tpar_right = (tpar == param_curve.endparam()) ? param_curve.startparam() : tpar;
    vector<Point> eval_left = param_curve.ParamCurve::point(tpar_left, 1, false);
    vector<Point> eval_right = param_curve.ParamCurve::point(tpar_right, 1, true);
    Point normal_left(-eval_left[1][1], eval_left[1][0]);
    normal_left.normalize();
    Point normal_right(-eval_right[1][1], eval_right[1][0]);
    normal_right.normalize();
    Point normal(-eval[1][1], eval[1][0]);
    normal.normalize();
    // Average tends to cause trouble, better to let user insert knots.
    Point average_normal = normal_left + normal_right;
    average_normal.normalize();
//     average_normal = normal;

    Point return_point;
    if (loop)
	return_point = eval[0] + dist*normal_right; //average_normal;
    else
	return_point = eval[0] + dist*normal;
    return return_point;
}


//===========================================================================
SplineCurve*
OffsetCurveUtils::createSISLApproximatedCurve(const SplineCurve& curve,
				       const double epsge)
//===========================================================================
// @ers: not tested
{
    // We would like the curve to be of order 3.
    SplineCurve in_crv = *(curve.clone());
    if (in_crv.order() < 4)
	in_crv.raiseOrder(4 - in_crv.order());

    DEBUG_ERROR_IF(epsge < 0.0,
		"Geometric tolerance must be larger than zero.");

    const shared_ptr<SISLCurve> sisl_curve(Curve2SISL(in_crv));
    
    vector<double> eps;
    eps.push_back(epsge);
    eps.push_back(epsge);

    SISLCurve* new_curve = NULL;
    double maxerr[2];
    int stat;

    s1940(sisl_curve.get(), &eps[0], in_crv.order(), in_crv.order(),
	  0, 10, &new_curve, maxerr, &stat);
    DEBUG_ERROR_IF(stat!=0,
		"Error in s1940, aborting.");

    SplineCurve* approx = SISLCurve2Go(new_curve);
    if (new_curve) freeCurve(new_curve);
    //    if (maxerr)    free(maxerr);    // is set in any case

    return approx;
}



//===========================================================================
vector<vector<shared_ptr<SplineCurve> > >
OffsetCurveUtils::getSimpleLoops2(const SplineCurve& loop, double loop_tol)
//===========================================================================
{
    // We find all self-intersections.
    vector<double> loop_params;
//     vector<int> orders;

#if 0
    SelfIntersectCurveAlgo self_int_obj;
    self_int_obj.useSplineCurve(loop);
    self_int_obj.perform();
    self_int_obj.getResultData(loop_params); //, orders);
#else
MESSAGE("Missing the call to SelfIntersectCurveAlgo!");
#endif

    // We convert to a vector of pairs.
    size_t i, j;
    vector<pair<double, double> > int_pairs;
    for (i = 0; i < loop_params.size() / 2; ++i)
	if (loop_params[2*i] != loop_params[2*i+1]) // We do not care about cusps.
	    int_pairs.push_back(make_pair(loop_params[2*i], loop_params[2*i+1]));
    sort(int_pairs.begin(), int_pairs.end());

    // If intersecting end parameters are missing, these are added.
    double startparam = loop.startparam();
    double endparam = loop.endparam();
    if ((int_pairs.size() == 0) ||
	(fabs(int_pairs[0].first - startparam) > loop_tol) ||
	(fabs(int_pairs[0].second - endparam) > loop_tol))
	int_pairs.insert(int_pairs.begin(), make_pair(startparam, endparam));

    // We transfer all intersections to a vector, and extract all loop segments.
    loop_params.clear();
    for (i = 0; i < int_pairs.size(); ++i) {
	loop_params.push_back(int_pairs[i].first);
	loop_params.push_back(int_pairs[i].second);	
    }
    sort(loop_params.begin(), loop_params.end());

    vector<vector<shared_ptr<SplineCurve> > > offset_segments;
    double dist;
    for (i = 0; i < loop_params.size() - 1; ++i) {
	if (loop_params[i] == loop_params[i+1]) {
	    // Suspect self_intersection routine gets unstable w/dense knot spacing.
	    MESSAGE("Intersection params are equal, this should never happen.");
	    continue;
	} else {
	    vector<shared_ptr<SplineCurve> > dummy_vec;
	    dummy_vec.push_back(shared_ptr<SplineCurve>
				 (loop.subCurve(loop_params[i],
						loop_params[i+1])));
	    offset_segments.push_back(dummy_vec);
	}
    }

    size_t ki, kj;
#ifndef NDEBUG
    {
	std::ofstream debug("data/debug.g2");
	for (ki = 0; ki < offset_segments.size(); ++ki)
	    for (kj = 0; kj < offset_segments[ki].size(); ++kj) {
		writeSpaceParamCurve(*(offset_segments[ki][kj]), debug);
	    }
    }
#endif // NDEBUG

    // Method is somewhat as the one described in the above version, but we will start by extracting
    // loops which consist of one segment only. We will then go on extracting consecutive segments
    // which form a loop (i.e. no other segment start/ends in interior of loop).

    vector<vector<shared_ptr<SplineCurve> > > loops; // Simple loops are extracted to loops.
    int next_element = -1;
    // We count number of dead ends. If all edgbes used, abort! Should never happen.
    int nmb_dead_ends = 0;
    bool finished = false;
    bool leftmost_cv = true; // We will pick curve the most to the left from curr_crv.

    // Here we start routine which extracts loops.
    while ((offset_segments.size() != 0) && (finished == false)) {
	vector<int> available(offset_segments.size(), 1);
	// We remove existing loops.
	while (true) {
	    for (ki = 0; ki < offset_segments.size(); ++ki) {
		vector<int> chosen;
		chosen.push_back(ki);
		vector<int> next_segments =
		    nextSegments(offset_segments, chosen, loop_tol);
		if ((next_segments.size() == 1) && (next_segments[0] != ki) &&
		    (available[next_segments[0]])) {
		    offset_segments[ki].insert(offset_segments[ki].end(),
					       offset_segments[next_segments[0]].begin(),
					       offset_segments[next_segments[0]].end());
		    offset_segments.erase(offset_segments.begin() + next_segments[0]);
		    available.erase(available.begin() + next_segments[0]);
		    if (next_segments[0] < ki)
			--ki;
		}
	    }
	    // We look for loops.
	    int nmb_new_loops = 0;
	    for (ki = 0; ki < offset_segments.size(); ++ki) {
		Point start_pt = offset_segments[ki].front()->ParamCurve::point
		    (offset_segments[ki].front()->startparam());
		Point end_pt = offset_segments[ki].back()->ParamCurve::point
		    (offset_segments[ki].back()->endparam());
		if (start_pt.dist(end_pt) < loop_tol) {
		    loops.push_back(offset_segments[ki]);
		    offset_segments.erase(offset_segments.begin() + ki);
		    available.erase(available.begin() + ki);
		    --ki;
		    ++nmb_new_loops;
		}
	    }
	    if (nmb_new_loops == 0)
		break;
	}

	if (offset_segments.size() == 0) {
	    break;
	}

	// We look upon a vector in consecutive segments as one segment.
	vector<int> chosen_segments; // We store indices of chosen segments.
	if (next_element == -1) {
	    double max_length = -1.0;
	    int max_index = -1;
	    for (ki = 0; ki < offset_segments.size(); ++ki) {
		double length = 0.0;
		for (kj = 0; kj < offset_segments[ki].size(); ++kj)
		    length += offset_segments[ki][kj]->estimatedCurveLength(1000);
		if (length > max_length) {
		    max_length = length;
		    max_index = ki;
		}
	    }
	    available[max_index] = 0;
	    chosen_segments.push_back(max_index);
	} else {
	    MESSAGE("Should not happen!");
	    available[next_element] = 0;
	    chosen_segments.push_back(next_element);
	    next_element = -1;
	}
	vector<shared_ptr<SplineCurve> > curr_segment = offset_segments[chosen_segments.back()];
#ifndef NDEBUG
	{
	    std::ofstream debug("data/debug.g2");
	    writeSpaceParamCurve(*(curr_segment.back()), debug);
	}
#endif // NDEBUG
	Point total_start_pt =
	    curr_segment.front()->ParamCurve::point(curr_segment.front()->startparam());
	Point end_pt =
	    curr_segment.back()->ParamCurve::point(curr_segment.back()->endparam());
	while ((end_pt - total_start_pt).length() > loop_tol) {
	    int index = getNextSegment(offset_segments, chosen_segments, available, loop_tol, leftmost_cv);
	    if (index != -1) {
		available[index] = 0;
		chosen_segments.push_back(index);
		curr_segment = offset_segments[chosen_segments.back()];
#ifndef NDEBUG
		{
		    std::ofstream debug("data/debug.g2");
		    writeSpaceParamCurve(*(curr_segment.back()), debug);
		}
#endif // NDEBUG
		end_pt = curr_segment.back()->ParamCurve::point(curr_segment.back()->endparam());
	    } else {
		chosen_segments.erase(chosen_segments.end() - 1);
		if (chosen_segments.size() == 0)
		    THROW("Should not happen, improve method!");
// 		++nmb_dead_ends;
// 		if (nmb_dead_ends == offset_segments.size()) {
// 		    // This may happen when a parameter is marked as a double
// 		    // intersection point. Seems strange that match would be exact.
// 		    MESSAGE("All available edges lead to a dead end! Should never "
// 			       "happen. Suspecting double intersection points.");
// 		    finished = true;
// 		} else {
// 		    // It appers that our selected start edge was part of a CW loop.
// 		    next_element = (chosen_segments[0] + 1) % offset_segments.size();
// 		    break;
// 		}
	    }
	}

// 	if ((next_element == -1) || (offset_segments.size() == 1)) {
	// We must extract the loop.
	vector<shared_ptr<SplineCurve> > loop_crvs;
	for (i = 0; i < chosen_segments.size(); ++i)
	    loop_crvs.insert(loop_crvs.end(),
			     offset_segments[chosen_segments[i]].begin(),
			     offset_segments[chosen_segments[i]].end());
        const double int_tol = 1.0e-08;
	if (!LoopUtils::loopIsCCW(loop_crvs, loop_tol, int_tol) && leftmost_cv == true) {
	    leftmost_cv = false; // We select cvs once more, now picking the rightmost piece.
	} else {
	    loops.push_back(loop_crvs);
	    // We must remove the extracted segments.
	    sort(chosen_segments.begin(), chosen_segments.end());
	    for (i = chosen_segments.size() - 1; i > -1; --i)
		offset_segments.erase(offset_segments.begin() + chosen_segments[i]);
	    nmb_dead_ends = 0; // We reset the boolean dead end.
	    leftmost_cv = true; // It may already be true, just making sure.
	}
    }

    return loops;
}


//===========================================================================
void
OffsetCurveUtils::removeInvalidLoops(vector<vector<shared_ptr<Go::SplineCurve> > >&
				    simple_offset_loops,
				    const SplineCurve& bnd_loop,
				    const SplineCurve& dist_function,
				    const SplineCurve& tol_function,
				    double epsge)
//===========================================================================
{
    size_t i, j, k, l;
    // Remove all CW loops.
    const double int_tol = 1.0e-08;
    for (i = 0; i < simple_offset_loops.size(); ++i) {
	bool ccw = LoopUtils::loopIsCCW(simple_offset_loops[i], epsge, int_tol);
	if (!ccw) {
	    simple_offset_loops.erase(simple_offset_loops.begin() + i,
			       simple_offset_loops.begin() + i + 1);
	    --i;
	}
    }

    // Remove all loops intersecting bnd_loop.
    for (i = 0; i < simple_offset_loops.size(); ++i) {
	for (j = 0; j < simple_offset_loops[i].size(); ++j) {
	    vector<pair<double, double> > intersections;
	    vector<int> pre_top;
            std::vector<std::pair<std::pair<double,double>, 
                                  std::pair<double,double> > > int_crvs;
            intersect2Dcurves(simple_offset_loops[i][j].get(),
			      bnd_loop.clone(),
			      epsge, intersections, pre_top, int_crvs);
	    if (intersections.size() != 0)
		break;
	}
	if (j < simple_offset_loops[i].size()) {
	    simple_offset_loops.erase(simple_offset_loops.begin() + i,
			       simple_offset_loops.begin() + i + 1);
	    --i;
	}
    }

    // @@sbr Should test for more sample points. Then check that intersection is empty.
    // Remove all loops lying outside bnd_loop.
    if (true) {
	shared_ptr<ParamCurve> crv(bnd_loop.clone());
	vector<shared_ptr<ParamCurve> > crv_vec;
	crv_vec.push_back(crv);
        double loop_tol = epsge; //hero_tol.getIGESFileTolerance();
	// const HeroTolerances& hero_tolerances = HeroTolerances::Instance();
	// double loop_tol = hero_tolerances.getIGESFileTolerance();
	shared_ptr<CurveLoop> loop(new CurveLoop(crv_vec, loop_tol));
	CurveBoundedDomain domain(loop);
	for (i = 0; i < simple_offset_loops.size(); ++i) {
	    // We pick a random point.
	    Point pnt = simple_offset_loops[i][0]->
		ParamCurve::point(simple_offset_loops[i][0]->startparam());
	    if (!(domain.isInDomain(Vector2D(pnt[0], pnt[1]), epsge))) {
		simple_offset_loops.erase(simple_offset_loops.begin() + i,
					  simple_offset_loops.begin() + i + 1);
		--i;
	    }
	}
    }

    // Remove all loops which, in an inner point, given by a weighted average
    // (convex combinztion) of sampled bnd points, fail distance test.
    // @@sbr Another method would be to use end parameter values of all segments.
    // Use parameter values of the surrounding loop.
    size_t nmb_loop_samples = 20;
    vector<double> loop_t_values;
    for (i = 0; i < nmb_loop_samples; ++i) {
	double tmin = bnd_loop.startparam();
	double tmax = bnd_loop.endparam();
	double step = (tmax - tmin) / nmb_loop_samples;
	double t_val = tmin + i*step;
	loop_t_values.push_back(t_val);
    }
    for (i = 0; i < simple_offset_loops.size(); ++i) {
	vector<double> params;
	int nmb_samples = 5; // nmb_samples for each segment of loop.
	Point conv_comb_pt(2);
	conv_comb_pt[0] = conv_comb_pt[1] = 0.0;
	size_t nmb_segments = simple_offset_loops[i].size();
	for (j = 0; j < nmb_segments; ++j) {
	    // We sample 5 points on each segment of loop.	    
	    size_t nmb_samples = 10;
	    double tmin = simple_offset_loops[i][j]->startparam();
	    double tmax = simple_offset_loops[i][j]->endparam();
	    double step = (tmax - tmin) / nmb_samples;
	    for (k = 0; k < nmb_samples; ++k) {
 		double tpar = tmin + k*step;
 		tpar = min(tmax, max(tmin, tpar)); // We move it inside domain.
		Point pt = simple_offset_loops[i][j]->ParamCurve::point(tpar);
// 		params.push_back(tpar);
		conv_comb_pt += pt / (nmb_samples * nmb_segments);
	    }
	}
	// We check whether conv_comb_pt is inside domain.
	vector<shared_ptr<ParamCurve> > cvs(simple_offset_loops[i].begin(),
					      simple_offset_loops[i].end());
	shared_ptr<CurveLoop> loop(new CurveLoop(cvs, epsge));
	CurveBoundedDomain domain(loop);
	Vector2D vec(conv_comb_pt[0], conv_comb_pt[1]);
	if (!domain.isInDomain(vec, epsge)) {
	    SplineCurve ploop = *(simple_offset_loops[i][0]->clone());
	    double tpar = 0.5*(ploop.startparam() + ploop.endparam());
	    double cont = 0;
	    double dist;
	    for (j = 1; j < simple_offset_loops[i].size(); ++j) {
		SplineCurve next_cv(*simple_offset_loops[i][j]);
		ploop.appendCurve(&next_cv, cont, dist);
	    }
	    if (ploop.basis().knotMultiplicity(tpar) != 0) {
		int knot_int = ploop.basis().knotIntervalFuzzy(tpar);
		vector<double>::const_iterator knot_iter = ploop.basis().begin() + knot_int;
		while (knot_iter[0] == knot_iter[1])
		    ++knot_iter;
		tpar = 0.5*(knot_iter[0] + knot_iter[1]);
	    }
	    double diag_length = 1e05; // @@sbr Not exactly, but usually large enough...
	    conv_comb_pt = findInnerPoint(ploop, tpar, diag_length);
	}

	// For conv_comb_pt we run a closest point tests against bnd_loop.
	for (j = 0; j < loop_t_values.size(); ++j) {
	    double clo_dist, clo_t;
	    Point clo_pt;
	    bnd_loop.closestPoint(conv_comb_pt,
				  bnd_loop.startparam(), bnd_loop.endparam(),
				  clo_t, clo_pt, clo_dist, &loop_t_values[j]);
	    Point dist_val = dist_function.ParamCurve::point(clo_t);
// 	    Point tol_val = tol_function.ParamCurve::point(clo_t);
	    if (clo_dist < dist_val[0]) {
// 		MESSAGE("Clipping offset loop too close to boundary loop.");
		simple_offset_loops.erase(simple_offset_loops.begin() + i,
					  simple_offset_loops.begin() + i + 1);
		--i;
		break;
	    }
	}
    }

    // @@sbr This seems to be somwhat ad hoc...
    // Finally we check each pair of loops for intersections (not counting touching)
    // and, if they intersect, remove the shortest of the two.
    // We will only be searching for intersections in the interior of segments.
    for (i = 0; i < simple_offset_loops.size(); ++i) {
	vector<shared_ptr<ParamCurve> > cvs1(simple_offset_loops[i].begin(),
					      simple_offset_loops[i].end());
	CurveLoop loop1(cvs1, epsge);
	for (j = i + 1; j < simple_offset_loops.size(); ++j) {
	    vector<shared_ptr<ParamCurve> > cvs2(simple_offset_loops[j].begin(),
						   simple_offset_loops[j].end());
	    CurveLoop loop2(cvs2, epsge);
	    if (loopsIntersect(loop1, loop2))
		// We remove the shortest loop...
		if (estimatedLoopLength(loop1) < estimatedLoopLength(loop2)) {
		    simple_offset_loops.erase(simple_offset_loops.begin() + i);
		    --i;
		    --j;
		} else {
		    simple_offset_loops.erase(simple_offset_loops.begin() + j);
		    --j;
		}
	}
    }
}


//===========================================================================
vector<int> OffsetCurveUtils::nextSegments(const vector<vector<shared_ptr<SplineCurve> > >& segments,
					  vector<int> chosen_segments, double loop_tol)
//===========================================================================
{
    Point curr_end_pt = segments[chosen_segments.back()].back()->ParamCurve::point
	(segments[chosen_segments.back()].back()->endparam());
    vector<int> next_segments; // At most two candidates.
    for (size_t ki = 0; ki < segments.size(); ++ki) {
// 	if (available[ki]) {
	    Point start_pt =
		segments[ki][0]->ParamCurve::point(segments[ki][0]->startparam());
// 	    if (((end_pt - start_pt).length() < loop_tol) ||
// 		((end_pt - curr_end_pt[0]).length() < loop_tol))
// 		available[i] = 0;
	    if (start_pt.dist(curr_end_pt) < loop_tol) {
		// We then make sure that segment does not end in end of a chosen segment (in the inner).
		size_t kj;
		for (kj = 0; kj < chosen_segments.size() - 1; ++kj) {
		    Point end_pt =
			segments[ki].back()->ParamCurve::point(segments[ki].back()->endparam());
		    Point local_end_pt = segments[chosen_segments[kj]].back()->
			ParamCurve::point(segments[chosen_segments[kj]].back()->endparam());
		    if (end_pt.dist(local_end_pt) < loop_tol)
			break;
		}
		if (kj == chosen_segments.size() - 1)
		    next_segments.push_back(ki);
// 		available[i] = 0;
	    }
	}

    return next_segments;
}



//===========================================================================
int OffsetCurveUtils::getNextSegment(const vector<vector<shared_ptr<SplineCurve> > >& segments,
				    vector<int> chosen_segments, vector<int> available,
				    double loop_tol, bool leftmost_cv)
//===========================================================================
{
    vector<int> candidates = nextSegments(segments, chosen_segments, loop_tol);
    // Remove unavailable
    for (size_t ki = candidates.size() - 1; ki > -1; --ki) {
	if (!available[candidates[ki]])
	    candidates.erase(candidates.begin() + ki);
    }

    if (candidates.size() > 2)
	THROW("More than two curves start in curve's end point. "
              "This should never happen!");
    else if (candidates.size() == 0)
	return -1;
    else if (candidates.size() == 1)
	return candidates[0];
    else { // 2 elements.
	vector<Point> curr_end_pt = segments[chosen_segments.back()].back()->
	    ParamCurve::point(segments[chosen_segments.back()].back()->endparam(), 1);
	// To choose, we must compute angle (in (-180, 180)) between end tangents.
	// We choose the curve with the highest angle, if leftmost_cv == true.
	vector<Point> first_start_pt = segments[candidates[0]][0]->
	    ParamCurve::point(segments[candidates[0]][0]->startparam(), 1);
	vector<Point> second_start_pt = segments[candidates[1]][0]->
	    ParamCurve::point(segments[candidates[1]][0]->startparam(), 1);
	Point curr_end_tg = curr_end_pt[1];
	curr_end_tg.normalize();
	Point first_start_tg = first_start_pt[1];
	first_start_tg.normalize();
	Point second_start_tg = second_start_pt[1];
	second_start_tg.normalize();
	double cos1 = curr_end_tg * first_start_tg;
	double sin1 =
	    curr_end_tg[0]*first_start_tg[1] - curr_end_tg[1]*first_start_tg[0];
	double cos2 = curr_end_tg * second_start_tg;
	double sin2 =
	    curr_end_tg[0]*second_start_tg[1] - curr_end_tg[1]*second_start_tg[0];
	if (sin1 * sin2 > 0) {
	    if (sin1 > 0) // i.e. both above curve_segment defined by curr_end_tg.
		if (cos1 < cos2)
		    return (leftmost_cv ? candidates[0] : candidates[1]);
		else
		    return (leftmost_cv ? candidates[1] : candidates[0]);
	    else // both under line.
		if (cos1 < cos2)
		    return (leftmost_cv ? candidates[1] : candidates[0]);
		else
		    return (leftmost_cv ? candidates[0] : candidates[1]);
	} else // They lie on opposite sides of line.
	    if (sin1 > 0)
		return (leftmost_cv ? candidates[0] : candidates[1]);
	    else
		return (leftmost_cv ? candidates[1] : candidates[0]);
    }
}


//===========================================================================
double OffsetCurveUtils::estimatedLoopLength(const CurveLoop& loop)
//===========================================================================
{
    double length = 0.0;
    for (size_t i = 0; i < loop.size(); ++i)
	length += loop[i]->estimatedCurveLength(1000);

    return length;
}


//===========================================================================
Point OffsetCurveUtils::findInnerPoint(SplineCurve& p_loop, double tpar,
                                       double domain_diagonal_length)
//===========================================================================
{
    double tang_tol = 1e-04;
    vector<Point> pts_left(2), pts_right(2);
    p_loop.point(pts_left, tpar, 1, false);
    p_loop.point(pts_right, tpar, 1, true);
    DEBUG_ERROR_IF(pts_left[1].dist(pts_right[1]) > tang_tol,
                   "Assuming curve is c1 in input tpar.");

    // We construct a linear pline curve going through pts_left[0], in direction
    // normal to pts_left[1] (choosing left direction as loop is ccw).
    Point dir(-pts_left[1][1], pts_left[1][0]);
    Point to_pt = pts_left[0] + (domain_diagonal_length + 1.0)*dir;
    SplineCurve linear_cv(pts_left[0], to_pt);
    double epsge = 1e-05;
    vector<pair<double, double> > intersections;
    vector<int> pretop;
    std::vector<std::pair<std::pair<double,double>, 
                          std::pair<double,double> > > int_crvs;
    intersect2Dcurves(&linear_cv, &p_loop, epsge, intersections, pretop, int_crvs);

    double tmin = linear_cv.startparam();
    sort(intersections.begin(), intersections.end());
    double bd_tpar;
    if (intersections.size() == 1)
	bd_tpar = intersections[0].first;
    else if (intersections.size() > 1)
	if (intersections[0].first == tmin)
	    bd_tpar = intersections[1].first;
	else
	    bd_tpar = intersections[0].first;
    else
	THROW("Failed finding inner pt of loop.");

    double inner_tpar = 0.5*(tmin + bd_tpar);
    Point inner_pt = linear_cv.ParamCurve::point(inner_tpar);
    return inner_pt;
}


//===========================================================================
bool OffsetCurveUtils::loopsIntersect(const CurveLoop& loop1, const CurveLoop& loop2,
				     bool test_end_pts)
//===========================================================================
{
    int i, j;
    for (i = 0; i < loop1.size(); ++i) {
	shared_ptr<SplineCurve> cv1 =
	    dynamic_pointer_cast<SplineCurve>(loop1[i]);
	DEBUG_ERROR_IF(cv1.get() == 0,
		    "Unexpected curve type.");
	for (j = 0; j < loop2.size(); ++j) {
	    shared_ptr<SplineCurve> cv2 =
		dynamic_pointer_cast<SplineCurve>(loop2[j]);
	    DEBUG_ERROR_IF(cv2.get() == 0,
			"Unexpected curve type.");

#if 0
	    IntersectCurveAlgo int_obj;
	    int_obj.useFirstSplineCurve(*cv1);
	    int_obj.useSecondSplineCurve(*cv2);
	    int_obj.perform();
	    vector<double> params;
// 	    vector<int> orders;
	    int_obj.getResultData(params); //, orders);
// 	    vector<pair<double, double> > params;
// 	    find_intersections(*cv1, *cv2, params);
	    if (params.size() != 0)
		break;
#else
            MESSAGE("Missing IntersectCurveAlgo!");
#endif

	}
	if (j < loop2.size())
	    break;
    }

    if (test_end_pts) {
	double tol = min(loop1.getSpaceEpsilon(), loop2.getSpaceEpsilon());
	Point start_pt1 = loop1[0]->point(loop1[0]->startparam());
	Point start_pt2 = loop2[0]->point(loop2[0]->startparam());
	if (start_pt1.dist(start_pt2) < tol)
	    i = 0; // Intersection in first segment.
    }
    return (i < loop1.size() ? true : false);
}


} // namespace Go
