#include "APPExamples.h"

// TerraLib
#include <terralib/common/progress/ConsoleProgressViewer.h>
#include <terralib/common/progress/ProgressManager.h>
#include <terralib/common/STLUtils.h>
#include <terralib/dataaccess.h>
#include <terralib/dataaccess/dataset/DataSet.h>
#include <terralib/app/core/Veredas.h>

// STL
#include <string>
#include <map>
#include <iostream>

void Veredas()
{
  te::common::ConsoleProgressViewer* cpViewer = new te::common::ConsoleProgressViewer();
  int id = te::common::ProgressManager::getInstance().addViewer(cpViewer);
  
  try
  {
    std::cout << "Veredas example using APP module." << std::endl << std::endl;

// open input data
    std::map<std::string, std::string> connInfo;
    connInfo["URI"] = TERRALIB_DATA_DIR"/app/vereda_guariroba_pol.shp";

    std::auto_ptr<te::da::DataSource> ds = te::da::DataSourceFactory::make("OGR");
    ds->setConnectionInfo(connInfo);
    ds->open();

    std::auto_ptr<te::da::DataSet> dataSet = ds->getDataSet("vereda_guariroba_pol");

    if(!dataSet.get())
    {
      std::cout << std::endl << "Invalid input data in Veredas()." << std::endl;
      return;
    }

    std::size_t gpos = te::da::GetFirstPropertyPos(dataSet.get(), te::dt::GEOMETRY_TYPE);

    std::vector<te::gm::Geometry*> geomVec;

    dataSet->moveBeforeFirst();

    while(dataSet->moveNext())
    {
      std::auto_ptr<te::gm::Geometry> g(dataSet->getGeometry(gpos));
      geomVec.push_back(g.release());
    }


    bool executeok = false;
    bool initok = false;
    {
      std::cout << "Executing Veredas APP Operation" << std::endl;

// create output data
      std::map<std::string, std::string> orinfo;
      orinfo["URI"] = TERRALIB_DATA_DIR"/app/vereda_guariroba_result.shp";

// create elevation algorithm parameters
      te::app::Veredas::InputParameters inputParameters;
      inputParameters.m_geometries = geomVec;
      inputParameters.m_bufferDistance = 50.;

      te::app::Veredas::OutputParameters outputParameters;
      outputParameters.m_createdOutDSName = "vereda_guariroba_result";
      outputParameters.m_createdOutInfo = orinfo;
      outputParameters.m_createdOutDSType = "OGR";

// execute the algorithm
      te::app::Veredas veredas;

      initok = veredas.initialize(inputParameters);

      if(initok)
        executeok = veredas.execute(outputParameters);

      if(!executeok)
        std::cout << "Problems in Veredas operation." << std::endl;
    }

    if (executeok)
      std::cout << "Done!" << std::endl << std::endl;
    else
      std::cout << "Problems in Veredas." << std::endl;

    te::common::FreeContents(geomVec);
  }
  catch(const std::exception& e)
  {
    std::cout << std::endl << "An exception has occurred in Veredas(): " << e.what() << std::endl;
  }
  catch(...)
  {
    std::cout << std::endl << "An unexpected exception has occurred in Veredas()!" << std::endl;
  }

  te::common::ProgressManager::getInstance().removeViewer(id);

  delete cpViewer;
}
