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
  \file terralib/app/core/Utils.cpp
 
  \brief Utils functions.
*/

#include "Utils.h"

//TerraLib Includes
#include <terralib/common/progress/TaskProgress.h>
#include <terralib/dataaccess/dataset/DataSetType.h>
#include <terralib/dataaccess/dataset/PrimaryKey.h>
#include <terralib/dataaccess/datasource/DataSource.h>
#include <terralib/dataaccess/datasource/DataSourceFactory.h>
#include <terralib/datatype/SimpleProperty.h>
#include <terralib/geometry/Geometry.h>
#include <terralib/geometry/GeometryCollection.h>
#include <terralib/geometry/GeometryProperty.h>
#include <terralib/geometry/LineString.h>
#include <terralib/geometry/Utils.h>
#include <terralib/memory/DataSet.h>
#include <terralib/memory/DataSetItem.h>
#include <terralib/raster/Band.h>
#include <terralib/raster/BandProperty.h>
#include <terralib/raster/Grid.h>
#include <terralib/raster/Raster.h>
#include <terralib/raster/RasterFactory.h>
#include <terralib/srs/Datum.h>
#include <terralib/srs/Ellipsoid.h>
#include <terralib/srs/GeographicCoordinateSystem.h>
#include <terralib/srs/SpatialReferenceSystemManager.h>

//STL Includes
#include <cassert>


#define M_PI       3.14159265358979323846

void te::app::Raster2Vector(te::rst::Raster* raster, int band, std::string dataSetName, std::string dsType, std::map<std::string, std::string> connInfo)
{
  assert(raster);

//create vectors
  std::vector<te::gm::Geometry*> geomVec;

  raster->vectorize(geomVec, band, 0);

//export vectors
  ExportVector(geomVec, dataSetName, dsType, connInfo);
}

void te::app::ExportVector(std::vector<te::gm::Geometry*>& geomVec, std::string dataSetName, std::string dsType, std::map<std::string, std::string> connInfo)
{
  assert(!geomVec.empty());

  //create dataset type
  std::auto_ptr<te::da::DataSetType> dataSetType(new te::da::DataSetType(dataSetName));

  //create id property-
  te::dt::SimpleProperty* idProperty = new te::dt::SimpleProperty("id", te::dt::INT32_TYPE);
  dataSetType->add(idProperty);

  //create geometry property
  te::gm::GeometryProperty* geomProperty = new te::gm::GeometryProperty("geom", geomVec[0]->getSRID(), te::gm::PolygonType);
  dataSetType->add(geomProperty);

  //create primary key
  std::string pkName = "pk_id";
              pkName+= "_" + dataSetName;
  te::da::PrimaryKey* pk = new te::da::PrimaryKey(pkName, dataSetType.get());
  pk->add(idProperty);

//create data set
  std::auto_ptr<te::mem::DataSet> dataSetMem(new te::mem::DataSet(dataSetType.get()));

  for(std::size_t t = 0; t < geomVec.size(); ++t)
  {
    //create dataset item
    te::mem::DataSetItem* item = new te::mem::DataSetItem(dataSetMem.get());

    //set id
    item->setInt32("id", (int)t);

    //set geometry
    item->setGeometry("geom", geomVec[t]);

    dataSetMem->add(item);
  }

  dataSetMem->moveBeforeFirst();

//save data set
  std::auto_ptr<te::da::DataSource> dataSource = te::da::DataSourceFactory::make(dsType);
  dataSource->setConnectionInfo(connInfo);
  dataSource->open();

  std::map<std::string, std::string> options;
  dataSource->createDataSet(dataSetType.get(), options);
  dataSource->add(dataSetName, dataSetMem.get(), options);
}

void te::app::ExportVectorLine(std::vector<te::gm::LineString*>& geomVec, std::string dataSetName, std::string dsType, std::map<std::string, std::string> connInfo)
{
  assert(!geomVec.empty());

  //create dataset type
  std::auto_ptr<te::da::DataSetType> dataSetType(new te::da::DataSetType(dataSetName));

  //create id property-
  te::dt::SimpleProperty* idProperty = new te::dt::SimpleProperty("id", te::dt::INT32_TYPE);
  dataSetType->add(idProperty);

  //create geometry property
  te::gm::GeometryProperty* geomProperty = new te::gm::GeometryProperty("geom", geomVec[0]->getSRID(), te::gm::LineStringType);
  dataSetType->add(geomProperty);

  //create primary key
  std::string pkName = "pk_id";
  pkName += "_" + dataSetName;
  te::da::PrimaryKey* pk = new te::da::PrimaryKey(pkName, dataSetType.get());
  pk->add(idProperty);

  //create data set
  std::auto_ptr<te::mem::DataSet> dataSetMem(new te::mem::DataSet(dataSetType.get()));

  for (std::size_t t = 0; t < geomVec.size(); ++t)
  {
    //create dataset item
    te::mem::DataSetItem* item = new te::mem::DataSetItem(dataSetMem.get());

    //set id
    item->setInt32("id", (int)t);

    //set geometry
    item->setGeometry("geom", geomVec[t]);

    dataSetMem->add(item);
  }

  dataSetMem->moveBeforeFirst();

  //save data set
  std::auto_ptr<te::da::DataSource> dataSource = te::da::DataSourceFactory::make(dsType);
  dataSource->setConnectionInfo(connInfo);
  dataSource->open();

  std::map<std::string, std::string> options;
  dataSource->createDataSet(dataSetType.get(), options);
  dataSource->add(dataSetName, dataSetMem.get(), options);
}

void te::app::Multi2Single(te::gm::Geometry* g, std::vector<te::gm::Geometry*>& geoms)
{
  te::gm::GeometryCollection* gc = dynamic_cast<te::gm::GeometryCollection*>(g);

  if (gc)
  {
    for (std::size_t i = 0; i < gc->getNumGeometries(); ++i)
    {
      te::app::Multi2Single(gc->getGeometryN(i), geoms);
    }
  }
  else
  {
    geoms.push_back(g);
  }
}

std::auto_ptr<te::rst::Raster> te::app::CalculateSlope(te::rst::Raster const* inputRst)
{
  //create slope raster
  std::vector< te::rst::BandProperty* > bandsProperties;
  te::rst::BandProperty* bandProp = new te::rst::BandProperty(0, te::dt::DOUBLE_TYPE);
  bandProp->m_noDataValue = -9999;
  bandsProperties.push_back(bandProp);

  te::rst::Grid* grid = new te::rst::Grid(*(inputRst->getGrid()));

  std::map<std::string, std::string> rInfo;

  te::rst::Raster* outRaster = te::rst::RasterFactory::make("MEM", grid, bandsProperties, rInfo);

  //initialize
  for (unsigned int i = 0; i < outRaster->getNumberOfRows(); ++i)
  {
    for (unsigned int j = 0; j < outRaster->getNumberOfColumns(); ++j)
    {
      outRaster->setValue(j, i, -9999, 0);
    }
  }

  unsigned int nlines = outRaster->getNumberOfRows();
  unsigned int ncolumns = outRaster->getNumberOfColumns();
  double rasterDummy = inputRst->getBand(0)->getProperty()->m_noDataValue;

  double resx = outRaster->getGrid()->getResolutionX(); //Gd
  double resy = outRaster->getGrid()->getResolutionY(); //Ge
  double A, F, E2, DY, DX;

  std::string authName = "EPSG"; // Now: So far it is the only one supported by TerraLib 5. Future: Review this line!
  bool isGeographic = te::srs::SpatialReferenceSystemManager::getInstance().isGeographic(outRaster->getGrid()->getSRID(), authName);
  te::srs::SpatialReferenceSystem* srs = te::srs::SpatialReferenceSystemManager::getInstance().getSpatialReferenceSystem(outRaster->getGrid()->getSRID()).release();

  te::srs::GeographicCoordinateSystem* gcs = dynamic_cast<te::srs::GeographicCoordinateSystem*>(srs);

  if (isGeographic && gcs)
  {
    A = gcs->getDatum()->getEllipsoid()->getRadium() / 1000.0;
    F = gcs->getDatum()->getEllipsoid()->getInverseFlattening();
    E2 = 2 * F - F * F; //!QUADRADO DA EXCENTRICIDADE
    DY = (A * sin(resy * M_PI / 180)) * 1000;

    resx = resx * 111133;
    resy = resy * 111133;
  }

  // Execute the operation
  te::common::TaskProgress task("Generating Slope");
  task.setTotalSteps(inputRst->getNumberOfRows());

  double rx, ry;

  //calculate slope
  for (unsigned int line = 1; line<nlines - 1; line++)
  {
    if (isGeographic)
    {
      //// A área varia somente na lina / y / latitude
      //te::gm::Coord2D coord = inputRst->getGrid()->gridToGeo(0, line);

      //double YLAT = coord.getY();

      //double RN = A / sqrt((1 - E2 * sin(YLAT * M_PI / 180)*sin(YLAT * M_PI / 180))); //!RAIO DE CURVATURA DA TERRA NA LATITUDE

      //// =B8*SEN(3/3600*PI()/180)*COS((A8+1.5/3600)*PI()/180)
      //DX = (RN * sin(resx * M_PI / 180) * cos(YLAT * M_PI / 180)) * 1000; // em metros

      rx = resx;
      ry = resy;
    }
    else
    {
      rx = resx;
      ry = resy;
    }

    for (unsigned int column = 1; column<ncolumns - 1; column++)
    {
      double value = rasterDummy;

      inputRst->getValue(column, line, value, 0);

      if (value != rasterDummy)
      {
        double z1 = 0, z2 = 0, z3 = 0, z4 = 0, z6 = 0, z7 = 0, z8 = 0, z9 = 0;

        inputRst->getValue(column - 1, line - 1, value, 0);
        if (value != rasterDummy)
        {
          z1 = value;
        }

        inputRst->getValue(column, line - 1, value, 0);
        if (value != rasterDummy)
        {
          z2 = value;
        }

        inputRst->getValue(column + 1, line - 1, value, 0);
        if (value != rasterDummy)
        {
          z3 = value;
        }

        inputRst->getValue(column - 1, line, value, 0);
        if (value != rasterDummy)
        {
          z4 = value;
        }

        inputRst->getValue(column + 1, line, value, 0);
        if (value != rasterDummy)
        {
          z6 = value;
        }

        inputRst->getValue(column - 1, line + 1, value, 0);
        if (value != rasterDummy)
        {
          z7 = value;
        }

        inputRst->getValue(column, line + 1, value, 0);
        if (value != rasterDummy)
        {
          z8 = value;
        }

        inputRst->getValue(column + 1, line + 1, value, 0);
        if (value != rasterDummy)
        {
          z9 = value;
        }

        double d = (z3 + z6 + z9 - z1 - z4 - z7) / (6 * rx);

        double e = (z1 + z2 + z3 - z7 - z8 - z9) / (6 * ry);

        outRaster->setValue(column, line, (float)sqrt(d*d + e*e), 0);
      }
    }

    task.pulse();
  }

  //set output raster into auto_ptr
  std::auto_ptr<te::rst::Raster> outputRst(outRaster);

  return outputRst;
}
