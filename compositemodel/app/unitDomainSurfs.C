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

#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/igeslib/IGESconverter.h"
#include "GoTools/utils/errormacros.h"
#include "GoTools/creators/CurveCreators.h"
#include <fstream>

using std::vector;
using std::cout;
using std::endl;
using namespace Go;
using std::vector;

int main( int argc, char* argv[] )
{
    if (argc != 3)
    {
        cout << "Input parameters : Input file on g2 format, output file" << std::endl;
        exit(-1);
    }

    // Read input arguments
    std::ifstream file1(argv[1]);
    ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

    std::ofstream file2(argv[2]);
    ALWAYS_ERROR_IF(file2.bad(), "Bad or no output filename");

    IGESconverter conv;
    conv.readgo(file1);
    vector<shared_ptr<GeomObject> > gogeom = conv.getGoGeom();
    int nmbgeom = (int)gogeom.size();

    vector<shared_ptr<ParamCurve> > all_3d_cvs;
    std::ofstream file4("tmp/all_3d_cvs.g2");

    for (int i=0; i<nmbgeom; i++)
    {
        if (gogeom[i].get() == 0)
            continue;
        shared_ptr<GeomObject> lg = gogeom[i];

        shared_ptr<ParamSurface> sf =
            dynamic_pointer_cast<ParamSurface, GeomObject>(lg);

        const double umin = 0.0;
        const double umax = 1.0;
        const double vmin = 0.0;
        const double vmax = 1.0;

        sf->setParameterDomain(umin, umax, vmin, vmax);

        sf->writeStandardHeader(file2);
        sf->write(file2);

        // We extract all curve on surface, write to file the lifted 2d curves and the 3d curves.
        if (sf->instanceType() == Class_BoundedSurface)
        {
            shared_ptr<BoundedSurface> bd_sf = dynamic_pointer_cast<BoundedSurface>(sf);
            shared_ptr<ParamSurface> under_sf = bd_sf->underlyingSurface();
            cout << "class type: " << under_sf->instanceType() << endl;
            const double epsgeo = bd_sf->getEpsGeo();
            std::vector<CurveLoop> bd_loops = bd_sf->absolutelyAllBoundaryLoops();
            for (auto loop : bd_loops)
            {
                std::vector<shared_ptr<ParamCurve> > cvs = loop.getCurves();
                for (auto cv : cvs)
                {
                    if (cv->instanceType() == Class_CurveOnSurface)
                    {
                        shared_ptr<CurveOnSurface> cv_on_sf = dynamic_pointer_cast<CurveOnSurface>(cv);
                        shared_ptr<ParamCurve> par_cv = cv_on_sf->parameterCurve();
                        if (par_cv)
                        {
                            shared_ptr<SplineCurve> lifted_cv(CurveCreators::liftParameterCurve(par_cv, under_sf, epsgeo));
                            if (lifted_cv)
                            {
                                all_3d_cvs.push_back(lifted_cv);
                            }
                            else
                            {
                                cout << "Failed lifting the 2d curve!" << endl;
                            }
                        }

                        shared_ptr<ParamCurve> space_cv = cv_on_sf->spaceCurve();
                        if (space_cv)
                        {
                            all_3d_cvs.push_back(space_cv);
                            space_cv->writeStandardHeader(file4);
                            space_cv->write(file4);
                        }
                        else
                        {
                            cout << "Missing 3d curve!" << endl;
                        }
                    }
                }
            }
        }
    }

    std::ofstream file3("tmp/all_lifted_2d_and_3d_cvs.g2");
    for (auto cv : all_3d_cvs)
    {
        cv->writeStandardHeader(file3);
        cv->write(file3);
    }

}



