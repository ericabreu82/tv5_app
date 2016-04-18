/*  Copyright (C) 2001-2009 National Institute For Space Research (INPE) - Brazil.

    This file is part of the TerraLib - a Framework for building GIS enabled applications.

    TerraLib is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License,
    or (at your option) any later version.

    TerraLib is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with TerraLib. See COPYING. If not, write to
    TerraLib Team at <terralib-team@terralib.org>.
 */

/*!
  \file terralib/app/core/Utils.h
 
  \brief Utils functions.
*/

#ifndef __TERRALIB_APP_INTERNAL_UTILS_H
#define __TERRALIB_APP_INTERNAL_UTILS_H

// Terralib
#include "../Config.h"

// STL
#include <map>
#include <memory>
#include <vector>

namespace te
{
  //forward declarations
  namespace gm  { class Geometry; class LineString; }
  namespace rst { class Raster; }

  namespace app
  {
    /*!
      \brief Function used to create a vectorial representation from a raster.
      \param raster       Input raster
      \param band         Band from input raster to be vectorized.
      \param dataSetName  Name of the output data set.
      \param dsType       Output data source type.
      \param connInfo     The necessary information to create the output data set.
    */
    APPPLUGINDLLEXPORT void Raster2Vector(te::rst::Raster* raster, int band, std::string dataSetName, std::string dsType, std::map<std::string, std::string> connInfo);

    /*!
      \brief Function used to export a vectorial representation to a data source.
      \param geoms        Input geometries
      \param dataSetName  Name of the output data set.
      \param dsType       Output data source type.
      \param connInfo     The necessary information to create the output data set.
    */
    APPPLUGINDLLEXPORT void ExportVector(std::vector<te::gm::Geometry*>& geomVec, std::string dataSetName, std::string dsType, std::map<std::string, std::string> connInfo);

    APPPLUGINDLLEXPORT void ExportVectorLine(std::vector<te::gm::LineString*>& geomVec, std::string dataSetName, std::string dsType, std::map<std::string, std::string> connInfo);

    /* 
      \brief Function used to get single geometry element from collection. 
      \param g      Geometry objetc to get the single geometry
      \param geoms  vector with single geometries
    */
    APPPLUGINDLLEXPORT void Multi2Single(te::gm::Geometry* g, std::vector<te::gm::Geometry*>& geoms);

    APPPLUGINDLLEXPORT std::auto_ptr<te::rst::Raster> CalculateSlope(te::rst::Raster const* inputRst);

  } // end namespace app
}   // end namespace te

#endif //__TERRALIB_APP_INTERNAL_UTILS_H

