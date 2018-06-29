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


#define BOOST_TEST_MODULE gotools-core/CurveOnSurfaceTest
#include <boost/test/included/unit_test.hpp>

#include <fstream>
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/BoundedUtils.h"

using namespace std;
using namespace Go;


#if 0
struct Config {
public:
    Config()
    {

        const std::string datadir = "data/"; // Relative to build/gotools-core

        infiles.push_back(datadir + "spline_cylinder.g2");
        GoTools::init();
    }

public:
    ObjectHeader header;
    vector<string> infiles;
};
#endif


BOOST_AUTO_TEST_CASE(sameOrientation)
{
    // We read a cylinder.
    const std::string infile = "data/spline_cylinder.g2"; // Relative to build/gotools-core

    ifstream in(infile.c_str());
    BOOST_CHECK_MESSAGE(in.good(), "Input file not found or file corrupt");
    ObjectHeader header;
    header.read(in);
    shared_ptr<SplineSurface> cyl(new SplineSurface());
    cyl->read(in);

    // We then extract an iso-circle around the cylinder and create the corresponding
    // parameter curve, store is as a CurveOnSurface.
    const RectDomain& rect_dom = cyl->containingDomain();
    const bool pardir_is_u = true;
    const double iso_par = (cyl->isBounded()) ? rect_dom.vmin() : 0.0;
    std::vector<shared_ptr<ParamCurve> > iso_cvs = cyl->constParamCurves(iso_par, pardir_is_u);

    Point start(rect_dom.umin(), iso_par);
    Point end(rect_dom.umax(), iso_par);

    shared_ptr<Line> par_cv(new Line(start, end, rect_dom.umin(), rect_dom.umax()));

    const bool pref_par = false;
    CurveOnSurface cv_on_sf(cyl, par_cv, iso_cvs[0], pref_par);

    bool same_orientation = cv_on_sf.sameOrientation();
    BOOST_CHECK_EQUAL(same_orientation, true);

    // We then reverse the parameter curve. The sameOrientation() function should now return the opposite
    // value.
    par_cv->reverseParameterDirection();
    same_orientation = cv_on_sf.sameOrientation();
    BOOST_CHECK_EQUAL(same_orientation, false);

}



struct ensureParCrvExistenceConfig {
public:
    ensureParCrvExistenceConfig()
    {

        const std::string datadir = "data/"; // Relative to build/gotools-core

        infiles.push_back(datadir + "cone_and_iso_circle.g2");
        success.push_back(true);

#if 1
        infiles.push_back(datadir + "cylinder_seam_cv.g2");
        success.push_back(false); // This case should fail as the curve lies on the seam.
#endif

        GoTools::init();
    }

public:
    ObjectHeader header;
    vector<string> infiles;
    vector<bool> success;
};


BOOST_FIXTURE_TEST_CASE(ensureParCrvExistence, ensureParCrvExistenceConfig)
{
    // We read a cone and a circular iso-curve lying on the cone.
    //const std::string infile = "data/cone_and_iso_circle.g2"; // Relative to build/gotools-core

    int cntr = -1;
    for (auto infile : infiles)
    {
        ++cntr;
        ifstream in(infile.c_str());
        BOOST_CHECK_MESSAGE(in.good(), "Input file not found or file corrupt");
        header.read(in);
        shared_ptr<ParamSurface> sf;
        if (header.classType() == Class_Cone)
        {
            sf = shared_ptr<Cone>(new Cone());
        }
        else if (header.classType() == Class_Cylinder)
        {
            sf = shared_ptr<Cylinder>(new Cylinder());
        }
        else
        {
            THROW("Unexpected surface type!");
        }
        sf->read(in);

        header.read(in);
        shared_ptr<ParamCurve> space_cv;
        if (header.classType() == Class_Circle)
        {
            space_cv = shared_ptr<Circle>(new Circle);
        }
        else if (header.classType() == Class_Line)
        {
            space_cv = shared_ptr<Line>(new Line);
        }
        else
        {
            THROW("Unexpected surface type!");
        }
        space_cv->read(in);

        const bool par_pref = false;
        CurveOnSurface cv_on_sf(sf, space_cv, par_pref);
    
        const double epsgeo = 1.0e-04;
        const bool cv_proj = cv_on_sf.ensureParCrvExistence(epsgeo);

        BOOST_CHECK_EQUAL(cv_proj, success[cntr]);

        // We also verify the end points of the projected curve lies in the domain of the surface. This test
        // was added to trigger a bug in the ensureParCrvExistence() function.
        shared_ptr<ParamCurve> pcv = cv_on_sf.parameterCurve();
        if (pcv)
        {
            Point start_pt = pcv->point(pcv->startparam());
            Point end_pt = pcv->point(pcv->endparam());
            double min_u = std::min(start_pt[0], end_pt[0]);
            double max_u = std::max(start_pt[0], end_pt[0]);
            double min_v = std::min(start_pt[1], end_pt[1]);
            double max_v = std::max(start_pt[1], end_pt[1]);
            const RectDomain& dom = sf->containingDomain();
            BOOST_CHECK(min_u >= dom.umin());
            BOOST_CHECK(max_u <= dom.umax());
            BOOST_CHECK(min_v >= dom.vmin());
            BOOST_CHECK(max_v <= dom.vmax());
        }
    }
}
