#include "APPExamples.h"

// TerraLib
#include <terralib/common/progress/ConsoleProgressViewer.h>
#include <terralib/common/progress/ProgressManager.h>
#include <terralib/common/STLUtils.h>
#include <terralib/dataaccess.h>
#include <terralib/dataaccess/dataset/DataSet.h>
#include <terralib/dataaccess/dataset/DataSetType.h>
#include <terralib/geometry/GeometryCollection.h>
#include <terralib/geometry/GeometryProperty.h>
#include <terralib/geometry/LineString.h>
#include <terralib/geometry/MultiLineString.h>
#include <terralib/app/core/Headwater.h>

// STL
#include <string>
#include <map>
#include <iostream>

void Headwater()
{
  te::common::ConsoleProgressViewer* cpViewer = new te::common::ConsoleProgressViewer();
  int id = te::common::ProgressManager::getInstance().addViewer(cpViewer);
  
  try
  {
    std::cout << "Headwater example using APP module." << std::endl << std::endl;

// open input data
    std::map<std::string, std::string> connInfo;
    connInfo["URI"] = TERRALIB_DATA_DIR"/app/romero_drenagem_lin.shp";

    std::auto_ptr<te::da::DataSource> ds = te::da::DataSourceFactory::make("OGR");
    ds->setConnectionInfo(connInfo);
    ds->open();

    std::auto_ptr<te::da::DataSet> dataSet = ds->getDataSet("romero_drenagem_lin");
    std::auto_ptr<te::da::DataSetType> dataSetType = ds->getDataSetType("romero_drenagem_lin");

    if(!dataSet.get())
    {
      std::cout << std::endl << "Invalid input data in Headwater()." << std::endl;
      return;
    }
    
    std::size_t gpos = te::da::GetFirstPropertyPos(dataSet.get(), te::dt::GEOMETRY_TYPE);
    te::gm::GeometryProperty* geomProp = te::da::GetFirstGeomProperty(dataSetType.get());

    std::vector<te::gm::LineString*> lineVec;

    dataSet->moveBeforeFirst();

    //get geometries
    while(dataSet->moveNext())
    {
      std::auto_ptr<te::gm::Geometry> g(dataSet->getGeometry(gpos));

      if(g->getGeomTypeId() == te::gm::MultiLineStringType)
      {
        te::gm::MultiLineString* mls = dynamic_cast<te::gm::MultiLineString*>(g.release());

        te::gm::LineString* ls = dynamic_cast<te::gm::LineString*>(mls->getGeometryN(0));

        lineVec.push_back(ls);
      }
      else if(g->getGeomTypeId() == te::gm::LineStringType)
      {
        te::gm::LineString* ls = dynamic_cast<te::gm::LineString*>(g.release());

        lineVec.push_back(ls);
      }
    }

    bool executeok = false;
    bool initok = false;
    {
      std::cout << "Executing Headwater APP Operation" << std::endl;

// create output data
      std::map<std::string, std::string> orinfo;
      orinfo["URI"] = TERRALIB_DATA_DIR"/app/romero_drenagem_lin_buffer.shp";

// create elevation algorithm parameters
      te::app::Headwater::InputParameters inputParameters;
      inputParameters.m_riversLines = lineVec;
      inputParameters.m_watercourseBuffer = 30.;
      inputParameters.m_headwaterBuffer = 50.;

      te::app::Headwater::OutputParameters outputParameters;
      outputParameters.m_createdOutDSName = "romero_drenagem_lin_buffer";
      outputParameters.m_createdOutInfo = orinfo;
      outputParameters.m_createdOutDSType = "OGR";

// execute the algorithm
      te::app::Headwater hw;

      initok = hw.initialize(inputParameters);

      if(initok)
        executeok = hw.execute(outputParameters);

      if(!executeok)
        std::cout << "Problems in Headwater operation." << std::endl;
    }

    if (executeok)
      std::cout << "Done!" << std::endl << std::endl;
    else
      std::cout << "Problems in Headwater." << std::endl;

    te::common::FreeContents(lineVec);
  }
  catch(const std::exception& e)
  {
    std::cout << std::endl << "An exception has occurred in Headwater(): " << e.what() << std::endl;
  }
  catch(...)
  {
    std::cout << std::endl << "An unexpected exception has occurred in Headwater()!" << std::endl;
  }

  te::common::ProgressManager::getInstance().removeViewer(id);

  delete cpViewer;
}
