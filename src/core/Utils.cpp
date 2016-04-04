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
#include <terralib/raster/Raster.h>


//STL Includes
#include <cassert>

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
