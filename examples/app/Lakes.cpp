#include "APPExamples.h"

// TerraLib
#include <terralib/common/progress/ConsoleProgressViewer.h>
#include <terralib/common/progress/ProgressManager.h>
#include <terralib/dataaccess.h>
#include <terralib/raster.h>

#include <terralib/app/core/Lakes.h>

// STL
#include <string>
#include <map>
#include <iostream>

std::vector<te::gm::Geometry*> GetGeometriesFromDataSource(const std::string& path, const std::string& dataSetName);

void Lakes()
{
  te::common::ConsoleProgressViewer* cpViewer = new te::common::ConsoleProgressViewer();
  int id = te::common::ProgressManager::getInstance().addViewer(cpViewer);
  
  try
  {
    std::cout << "Lakes example using APP module." << std::endl << std::endl;

// open lakes input data
    std::string lakesPath = TERRALIB_DATA_DIR"/app/lakes.shp";
    std::string lakesDataSet = "lakes";

    std::vector<te::gm::Geometry*> lakesGeoms = GetGeometriesFromDataSource(lakesPath, lakesDataSet);

// open urban input data
    std::string lakesUrbanPath = TERRALIB_DATA_DIR"/app/lakes_urban.shp";
    std::string lakesUrbanDataSet = "lakes_urban";

    std::vector<te::gm::Geometry*> lakesUrbanGeoms = GetGeometriesFromDataSource(lakesUrbanPath, lakesUrbanDataSet);

// open rural input data
    std::string lakesRuralPath = TERRALIB_DATA_DIR"/app/lakes_rural.shp";
    std::string lakesRuralDataSet = "lakes_rural";

    std::vector<te::gm::Geometry*> lakesRuralGeoms = GetGeometriesFromDataSource(lakesRuralPath, lakesRuralDataSet);

    bool executeok = false;
    bool initok = false;
    {
      std::cout << "Executing Lakes APP Operation" << std::endl;

// create output data
      std::map<std::string, std::string> orinfo;
      orinfo["URI"] = TERRALIB_DATA_DIR"/app/lakes_app.shp";

// create elevation algorithm parameters
      te::app::Lakes::InputParameters inputParameters;
      inputParameters.m_lakesGeometries = lakesGeoms;
      inputParameters.m_urbanGeometries = lakesUrbanGeoms;
      inputParameters.m_ruralGeometries = lakesRuralGeoms;
      inputParameters.m_urbanBufferDistance = 30;
      inputParameters.m_ruralBigBufferDistance = 100;
      inputParameters.m_ruralSmallBufferDistance = 50;

      te::app::Lakes::OutputParameters outputParameters;
      outputParameters.m_createdOutDSName = "lakes_app";
      outputParameters.m_createdOutInfo = orinfo;
      outputParameters.m_createdOutDSType = "OGR";

// execute the algorithm
      te::app::Lakes lakes;

      initok = lakes.initialize(inputParameters);

      if(initok)
        executeok = lakes.execute(outputParameters);

      if(!executeok)
        std::cout << "Problems in Lakes operation." << std::endl;
    }

// clean up
    te::common::FreeContents(lakesGeoms);
    te::common::FreeContents(lakesUrbanGeoms);
    te::common::FreeContents(lakesRuralGeoms);

    if (executeok)
      std::cout << "Done!" << std::endl << std::endl;
    else
      std::cout << "Problems in Lakes." << std::endl;
  }
  catch(const std::exception& e)
  {
    std::cout << std::endl << "An exception has occurred in Lakes(): " << e.what() << std::endl;
  }
  catch(...)
  {
    std::cout << std::endl << "An unexpected exception has occurred in Lakes()!" << std::endl;
  }

  te::common::ProgressManager::getInstance().removeViewer(id);

  delete cpViewer;
}

std::vector<te::gm::Geometry*> GetGeometriesFromDataSource(const std::string& path, const std::string& dataSetName)
{
  std::vector<te::gm::Geometry*> geomVec;

  std::map<std::string, std::string> connInfo;
  connInfo["URI"] = path;

  std::auto_ptr<te::da::DataSource> ds = te::da::DataSourceFactory::make("OGR");
  ds->setConnectionInfo(connInfo);
  ds->open();

  std::auto_ptr<te::da::DataSet> dataSet = ds->getDataSet(dataSetName);

  if(!dataSet.get())
  {
    std::cout << std::endl << "Invalid input data in Lakes()." << std::endl;
    return geomVec;
  }

  std::size_t gpos = te::da::GetFirstPropertyPos(dataSet.get(), te::dt::GEOMETRY_TYPE);

  dataSet->moveBeforeFirst();

  while(dataSet->moveNext())
  {
    std::auto_ptr<te::gm::Geometry> g(dataSet->getGeometry(gpos));
    geomVec.push_back(g.release());
  }

  return geomVec;
}
