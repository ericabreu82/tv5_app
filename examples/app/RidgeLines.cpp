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
#include <terralib/app/core/RidgeLines.h>

// STL
#include <string>
#include <map>
#include <iostream>

te::gm::MultiLineString* GetContourLines(std::string filePath, std::string dataSetName);

te::gm::MultiLineString* GetDrainageLines(std::string filePath, std::string dataSetName);

te::gm::MultiPolygon* GetWatershedsPolygons(std::string filePath, std::string dataSetName);

te::gm::MultiPolygon* GetWatershedsBufferPolygons(std::string filePath, std::string dataSetName);

void RidgeLines()
{
  te::common::ConsoleProgressViewer* cpViewer = new te::common::ConsoleProgressViewer();
  int id = te::common::ProgressManager::getInstance().addViewer(cpViewer);

  try
  {
    std::cout << "RidgeLines example using APP module." << std::endl << std::endl;

    // open input raster DEM
    std::map<std::string, std::string> rinfo;
    rinfo["URI"] = TERRALIB_DATA_DIR"/app/Alt_FazMtesClars_grade.tif";

    te::rst::Raster* rasterSRTM = te::rst::RasterFactory::open(rinfo);
    rasterSRTM->getGrid()->setSRID(22523);

    // open input raster SLOPE
    rinfo["URI"] = TERRALIB_DATA_DIR"/app/Alt_FazMtesClars_grade_slp.tif";

    te::rst::Raster* rasterSlope = te::rst::RasterFactory::open(rinfo);
    rasterSlope->getGrid()->setSRID(22523);

    // open input raster SLOPE MASK
    rinfo["URI"] = TERRALIB_DATA_DIR"/app/slope_mask_10.tif";

    te::rst::Raster* rasterSlopeMask = te::rst::RasterFactory::open(rinfo);
    rasterSlopeMask->getGrid()->setSRID(22523);

    // open input raster PLAIN MASK
    rinfo["URI"] = TERRALIB_DATA_DIR"/app/plains_mask_10_10.tif";

    te::rst::Raster* rasterPlainMask = te::rst::RasterFactory::open(rinfo);
    rasterPlainMask->getGrid()->setSRID(22523);

    // get input geometries
    te::gm::MultiLineString* contourLines = GetContourLines(TERRALIB_DATA_DIR"/app/curvas_nivel_lin.shp", "curvas_nivel_lin");

    te::gm::MultiLineString* drainageLines = GetDrainageLines(TERRALIB_DATA_DIR"/app/Alt_FazMtesClars_grade_vsg_500.shp", "Alt_FazMtesClars_grade_vsg_500");

    te::gm::MultiPolygon* watershedsPolygosn = GetWatershedsPolygons(TERRALIB_DATA_DIR"/app/Alt_FazMtesClars_grade_vsw_500.shp", "Alt_FazMtesClars_grade_vsw_500");

    te::gm::MultiPolygon* watershedsBufferPolygosn = GetWatershedsBufferPolygons(TERRALIB_DATA_DIR"/app/buffer_vsw_500_10.shp", "buffer_vsw_500_10");

    //execute
    bool executeok = false;
    bool initok = false;
    {
      std::cout << "Executing RidgeLines APP Operation" << std::endl;

      // create output data
      std::map<std::string, std::string> orinfo;
      orinfo["URI"] = TERRALIB_DATA_DIR"/app/Alt_FazMtesClars_ridgeLines.shp";

      // create Plateaus algorithm parameters
      te::app::RidgeLines::InputParameters inputParameters;
      inputParameters.m_demRasterPtr = rasterSRTM;
      inputParameters.m_slopeRasterPtr = rasterSlope;
      inputParameters.m_slopeMaskRasterPtr = rasterSlopeMask;
      inputParameters.m_plainMaskRasterPtr = rasterPlainMask;
      inputParameters.m_contourLines = contourLines;
      inputParameters.m_drainage = drainageLines;
      inputParameters.m_watersheds = watershedsPolygosn;
      inputParameters.m_watershedsBuffer = watershedsBufferPolygosn;
      inputParameters.m_height = 50;
      inputParameters.m_slope = 17;

      te::app::RidgeLines::OutputParameters outputParameters;
      outputParameters.m_createdOutDSName = "Alt_FazMtesClars_ridgeLines";
      outputParameters.m_createdOutInfo = orinfo;
      outputParameters.m_createdOutDSType = "OGR";

      // execute the algorithm
      te::app::RidgeLines ridgeLines;

      initok = ridgeLines.initialize(inputParameters);

      if (initok)
        executeok = ridgeLines.execute(outputParameters);

      if (!executeok)
        std::cout << "Problems in RidgeLines operation." << std::endl;
    }

    // clean up
    delete rasterSRTM;
    delete rasterSlope;
    delete contourLines;

    if (executeok)
      std::cout << "Done!" << std::endl << std::endl;
    else
      std::cout << "Problems in ridgeLines." << std::endl;
  }
  catch (const std::exception& e)
  {
    std::cout << std::endl << "An exception has occurred in ridgeLines(): " << e.what() << std::endl;
  }
  catch (...)
  {
    std::cout << std::endl << "An unexpected exception has occurred in ridgeLines()!" << std::endl;
  }

  te::common::ProgressManager::getInstance().removeViewer(id);

  delete cpViewer;
}

te::gm::MultiLineString* GetContourLines(std::string filePath, std::string dataSetName)
{
  std::map<std::string, std::string> connInfo;
  connInfo["URI"] = filePath;

  std::auto_ptr<te::da::DataSource> ds = te::da::DataSourceFactory::make("OGR");
  ds->setConnectionInfo(connInfo);
  ds->open();

  std::auto_ptr<te::da::DataSet> dataSet = ds->getDataSet(dataSetName);
  std::auto_ptr<te::da::DataSetType> dataSetType = ds->getDataSetType(dataSetName);

  if (!dataSet.get())
  {
    std::cout << std::endl << "Invalid input data in RidgeLines()." << std::endl;
    return 0;
  }

  std::size_t gpos = te::da::GetFirstPropertyPos(dataSet.get(), te::dt::GEOMETRY_TYPE);
  te::gm::GeometryProperty* geomProp = te::da::GetFirstGeomProperty(dataSetType.get());

  int srid = 22523;

  te::gm::MultiLineString* contourLines = new te::gm::MultiLineString(0, te::gm::MultiLineStringType, srid);

  dataSet->moveBeforeFirst();

  //get geometries
  while (dataSet->moveNext())
  {
    std::auto_ptr<te::gm::Geometry> g(dataSet->getGeometry(gpos));

    if (g->getGeomTypeId() == te::gm::MultiLineStringType)
    {
      te::gm::MultiLineString* mls = dynamic_cast<te::gm::MultiLineString*>(g.release());

      te::gm::LineString* ls = dynamic_cast<te::gm::LineString*>(mls->getGeometryN(0));

      ls->setSRID(srid);

      contourLines->add(ls);
    }
    else if (g->getGeomTypeId() == te::gm::LineStringType)
    {
      te::gm::LineString* ls = dynamic_cast<te::gm::LineString*>(g.release());

      ls->setSRID(srid);

      contourLines->add(ls);
    }
    else if (g->getGeomTypeId() == te::gm::LineStringZType)
    {
      te::gm::LineString* ls = dynamic_cast<te::gm::LineString*>(g.release());

      ls->setSRID(srid);

      contourLines->add(ls);
    }
  }

  return contourLines;
}

te::gm::MultiLineString* GetDrainageLines(std::string filePath, std::string dataSetName)
{
  std::map<std::string, std::string> connInfo;
  connInfo["URI"] = filePath;

  std::auto_ptr<te::da::DataSource> ds = te::da::DataSourceFactory::make("OGR");
  ds->setConnectionInfo(connInfo);
  ds->open();

  std::auto_ptr<te::da::DataSet> dataSet = ds->getDataSet(dataSetName);
  std::auto_ptr<te::da::DataSetType> dataSetType = ds->getDataSetType(dataSetName);

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

te::gm::MultiPolygon* GetWatershedsPolygons(std::string filePath, std::string dataSetName)
{
  std::map<std::string, std::string> connInfo;
  connInfo["URI"] = filePath;

  std::auto_ptr<te::da::DataSource> ds = te::da::DataSourceFactory::make("OGR");
  ds->setConnectionInfo(connInfo);
  ds->open();

  std::auto_ptr<te::da::DataSet> dataSet = ds->getDataSet(dataSetName);
  std::auto_ptr<te::da::DataSetType> dataSetType = ds->getDataSetType(dataSetName);

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

te::gm::MultiPolygon* GetWatershedsBufferPolygons(std::string filePath, std::string dataSetName)
{
  std::map<std::string, std::string> connInfo;
  connInfo["URI"] = filePath;

  std::auto_ptr<te::da::DataSource> ds = te::da::DataSourceFactory::make("OGR");
  ds->setConnectionInfo(connInfo);
  ds->open();

  std::auto_ptr<te::da::DataSet> dataSet = ds->getDataSet(dataSetName);
  std::auto_ptr<te::da::DataSetType> dataSetType = ds->getDataSetType(dataSetName);

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
