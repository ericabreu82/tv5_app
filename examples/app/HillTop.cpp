#include "APPExamples.h"

// TerraLib
#include <terralib/common/progress/ConsoleProgressViewer.h>
#include <terralib/common/progress/ProgressManager.h>
#include <terralib/common/STLUtils.h>
#include <terralib/dataaccess.h>
#include <terralib/dataaccess/dataset/DataSet.h>
#include <terralib/geometry/GeometryProperty.h>
#include <terralib/geometry/MultiLineString.h>
#include <terralib/geometry/MultiPolygon.h>
#include <terralib/raster/Grid.h>
#include <terralib/raster/RasterFactory.h>
#include <terralib/app/core/HillTop.h>

// STL
#include <string>
#include <map>
#include <iostream>

te::gm::MultiLineString* GetContourLines2();

te::gm::MultiLineString* GetDrainageLines2();

te::gm::MultiPolygon* GetWatershedsPolygons2();

te::gm::MultiPolygon* GetWatershedsBufferPolygons2();

void HillTop()
{
  te::common::ConsoleProgressViewer* cpViewer = new te::common::ConsoleProgressViewer();
  int id = te::common::ProgressManager::getInstance().addViewer(cpViewer);

  try
  {
    std::cout << "HillTop example using APP module." << std::endl << std::endl;

    // open input raster
    std::map<std::string, std::string> rinfo;
    rinfo["URI"] = TERRALIB_DATA_DIR"/app/Alt_FazMtesClars_grade.tif";

    te::rst::Raster* rasterSRTM = te::rst::RasterFactory::open(rinfo);
    rasterSRTM->getGrid()->setSRID(22523);

    // open input raster
    rinfo["URI"] = TERRALIB_DATA_DIR"/app/Alt_FazMtesClars_grade_slope.tif";

    te::rst::Raster* rasterSlope = te::rst::RasterFactory::open(rinfo);
    rasterSlope->getGrid()->setSRID(22523);

    // get input geometries
    te::gm::MultiLineString* contourLines = GetContourLines2();

    te::gm::MultiLineString* drainageLines = GetDrainageLines2();

    te::gm::MultiPolygon* watershedsPolygosn = GetWatershedsPolygons2();

    te::gm::MultiPolygon* watershedsBufferPolygosn = GetWatershedsBufferPolygons2();

    //execute
    bool executeok = false;
    bool initok = false;
    {
      std::cout << "Executing HillTop APP Operation" << std::endl;

      // create output data
      std::map<std::string, std::string> orinfo;
      orinfo["URI"] = TERRALIB_DATA_DIR"/app/Alt_FazMtesClars_hillTop.shp";

      // create Plateaus algorithm parameters
      te::app::HillTop::InputParameters inputParameters;
      inputParameters.m_demRasterPtr = rasterSRTM;
      inputParameters.m_slopeRasterPtr = rasterSlope;
      inputParameters.m_contourLines = contourLines;
      inputParameters.m_drainage = drainageLines;
      inputParameters.m_watersheds = watershedsPolygosn;
      inputParameters.m_watershedsBuffer = watershedsBufferPolygosn;
      inputParameters.m_height = 50;
      inputParameters.m_slope = 17;

      te::app::HillTop::OutputParameters outputParameters;
      outputParameters.m_createdOutDSName = "Alt_FazMtesClars_hillTop";
      outputParameters.m_createdOutInfo = orinfo;
      outputParameters.m_createdOutDSType = "OGR";

      // execute the algorithm
      te::app::HillTop hillTop;

      initok = hillTop.initialize(inputParameters);

      if (initok)
        executeok = hillTop.execute(outputParameters);

      if (!executeok)
        std::cout << "Problems in HillTop operation." << std::endl;
    }

    // clean up
    delete rasterSRTM;
    delete rasterSlope;
    delete contourLines;

    if (executeok)
      std::cout << "Done!" << std::endl << std::endl;
    else
      std::cout << "Problems in HillTop." << std::endl;
  }
  catch (const std::exception& e)
  {
    std::cout << std::endl << "An exception has occurred in HillTop(): " << e.what() << std::endl;
  }
  catch (...)
  {
    std::cout << std::endl << "An unexpected exception has occurred in HillTop()!" << std::endl;
  }

  te::common::ProgressManager::getInstance().removeViewer(id);

  delete cpViewer;
}

te::gm::MultiLineString* GetContourLines2()
{
  std::map<std::string, std::string> connInfo;
  connInfo["URI"] = TERRALIB_DATA_DIR"/app/curvas_nivel_lin.shp";

  std::auto_ptr<te::da::DataSource> ds = te::da::DataSourceFactory::make("OGR");
  ds->setConnectionInfo(connInfo);
  ds->open();

  std::auto_ptr<te::da::DataSet> dataSet = ds->getDataSet("curvas_nivel_lin");
  std::auto_ptr<te::da::DataSetType> dataSetType = ds->getDataSetType("curvas_nivel_lin");

  if (!dataSet.get())
  {
    std::cout << std::endl << "Invalid input data in HillTop()." << std::endl;
    return 0;
  }

  std::size_t gpos = te::da::GetFirstPropertyPos(dataSet.get(), te::dt::GEOMETRY_TYPE);
  te::gm::GeometryProperty* geomProp = te::da::GetFirstGeomProperty(dataSetType.get());

  te::gm::MultiLineString* contourLines = new te::gm::MultiLineString(0, te::gm::MultiLineStringType, geomProp->getSRID());

  dataSet->moveBeforeFirst();

  //get geometries
  while (dataSet->moveNext())
  {
    std::auto_ptr<te::gm::Geometry> g(dataSet->getGeometry(gpos));

    if (g->getGeomTypeId() == te::gm::MultiLineStringType)
    {
      te::gm::MultiLineString* mls = dynamic_cast<te::gm::MultiLineString*>(g.release());

      te::gm::LineString* ls = dynamic_cast<te::gm::LineString*>(mls->getGeometryN(0));

      contourLines->add(ls);
    }
    else if (g->getGeomTypeId() == te::gm::LineStringType)
    {
      te::gm::LineString* ls = dynamic_cast<te::gm::LineString*>(g.release());

      contourLines->add(ls);
    }
    else if (g->getGeomTypeId() == te::gm::LineStringZType)
    {
      te::gm::LineString* ls = dynamic_cast<te::gm::LineString*>(g.release());

      contourLines->add(ls);
    }
  }

  return contourLines;
}

te::gm::MultiLineString* GetDrainageLines2()
{
  std::map<std::string, std::string> connInfo;
  connInfo["URI"] = TERRALIB_DATA_DIR"/app/drenagem_lin.shp";

  std::auto_ptr<te::da::DataSource> ds = te::da::DataSourceFactory::make("OGR");
  ds->setConnectionInfo(connInfo);
  ds->open();

  std::auto_ptr<te::da::DataSet> dataSet = ds->getDataSet("drenagem_lin");
  std::auto_ptr<te::da::DataSetType> dataSetType = ds->getDataSetType("drenagem_lin");

  if (!dataSet.get())
  {
    std::cout << std::endl << "Invalid input data in RidgeLines()." << std::endl;
    return 0;
  }

  std::size_t gpos = te::da::GetFirstPropertyPos(dataSet.get(), te::dt::GEOMETRY_TYPE);
  te::gm::GeometryProperty* geomProp = te::da::GetFirstGeomProperty(dataSetType.get());

  te::gm::MultiLineString* drainageLines = new te::gm::MultiLineString(0, te::gm::MultiLineStringType, geomProp->getSRID());

  dataSet->moveBeforeFirst();

  //get geometries
  while (dataSet->moveNext())
  {
    std::auto_ptr<te::gm::Geometry> g(dataSet->getGeometry(gpos));

    if (g->getGeomTypeId() == te::gm::MultiLineStringType)
    {
      te::gm::MultiLineString* mls = dynamic_cast<te::gm::MultiLineString*>(g.release());

      te::gm::LineString* ls = dynamic_cast<te::gm::LineString*>(mls->getGeometryN(0));

      drainageLines->add(ls);
    }
    else if (g->getGeomTypeId() == te::gm::LineStringType)
    {
      te::gm::LineString* ls = dynamic_cast<te::gm::LineString*>(g.release());

      drainageLines->add(ls);
    }
    else if (g->getGeomTypeId() == te::gm::LineStringZType)
    {
      te::gm::LineString* ls = dynamic_cast<te::gm::LineString*>(g.release());

      drainageLines->add(ls);
    }
  }

  return drainageLines;
}

te::gm::MultiPolygon* GetWatershedsPolygons2()
{
  std::map<std::string, std::string> connInfo;
  connInfo["URI"] = TERRALIB_DATA_DIR"/app/subbacias_500.shp";

  std::auto_ptr<te::da::DataSource> ds = te::da::DataSourceFactory::make("OGR");
  ds->setConnectionInfo(connInfo);
  ds->open();

  std::auto_ptr<te::da::DataSet> dataSet = ds->getDataSet("subbacias_500");
  std::auto_ptr<te::da::DataSetType> dataSetType = ds->getDataSetType("subbacias_500");

  if (!dataSet.get())
  {
    std::cout << std::endl << "Invalid input data in RidgeLines()." << std::endl;
    return 0;
  }

  std::size_t gpos = te::da::GetFirstPropertyPos(dataSet.get(), te::dt::GEOMETRY_TYPE);
  te::gm::GeometryProperty* geomProp = te::da::GetFirstGeomProperty(dataSetType.get());

  te::gm::MultiPolygon* watershedsPolygons = new te::gm::MultiPolygon(0, te::gm::MultiPolygonType, geomProp->getSRID());

  dataSet->moveBeforeFirst();

  //get geometries
  while (dataSet->moveNext())
  {
    std::auto_ptr<te::gm::Geometry> g(dataSet->getGeometry(gpos));

    if (g->getGeomTypeId() == te::gm::MultiPolygonType)
    {
      te::gm::MultiPolygon* mp = dynamic_cast<te::gm::MultiPolygon*>(g.release());

      te::gm::Polygon* p = dynamic_cast<te::gm::Polygon*>(mp->getGeometryN(0));

      watershedsPolygons->add(p);
    }
    else if (g->getGeomTypeId() == te::gm::PolygonType)
    {
      te::gm::Polygon* p = dynamic_cast<te::gm::Polygon*>(g.release());

      watershedsPolygons->add(p);
    }
  }

  return watershedsPolygons;
}

te::gm::MultiPolygon* GetWatershedsBufferPolygons2()
{
  std::map<std::string, std::string> connInfo;
  connInfo["URI"] = TERRALIB_DATA_DIR"/app/buffer_subbacias_500.shp";

  std::auto_ptr<te::da::DataSource> ds = te::da::DataSourceFactory::make("OGR");
  ds->setConnectionInfo(connInfo);
  ds->open();

  std::auto_ptr<te::da::DataSet> dataSet = ds->getDataSet("buffer_subbacias_500");
  std::auto_ptr<te::da::DataSetType> dataSetType = ds->getDataSetType("buffer_subbacias_500");

  if (!dataSet.get())
  {
    std::cout << std::endl << "Invalid input data in RidgeLines()." << std::endl;
    return 0;
  }

  std::size_t gpos = te::da::GetFirstPropertyPos(dataSet.get(), te::dt::GEOMETRY_TYPE);
  te::gm::GeometryProperty* geomProp = te::da::GetFirstGeomProperty(dataSetType.get());

  te::gm::MultiPolygon* watershedsBufferPolygons = new te::gm::MultiPolygon(0, te::gm::MultiPolygonType, geomProp->getSRID());

  dataSet->moveBeforeFirst();

  //get geometries
  while (dataSet->moveNext())
  {
    std::auto_ptr<te::gm::Geometry> g(dataSet->getGeometry(gpos));

    if (g->getGeomTypeId() == te::gm::MultiPolygonType)
    {
      te::gm::MultiPolygon* mp = dynamic_cast<te::gm::MultiPolygon*>(g.release());

      te::gm::Polygon* p = dynamic_cast<te::gm::Polygon*>(mp->getGeometryN(0));

      watershedsBufferPolygons->add(p);
    }
    else if (g->getGeomTypeId() == te::gm::PolygonType)
    {
      te::gm::Polygon* p = dynamic_cast<te::gm::Polygon*>(g.release());

      watershedsBufferPolygons->add(p);
    }
  }

  return watershedsBufferPolygons;
}
