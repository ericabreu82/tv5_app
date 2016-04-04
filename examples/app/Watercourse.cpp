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
#include <terralib/app/core/Utils.h>
#include <terralib/app/core/Watercourse.h>

// STL
#include <string>
#include <map>
#include <iostream>

void Watercourse()
{
  te::common::ConsoleProgressViewer* cpViewer = new te::common::ConsoleProgressViewer();
  int id = te::common::ProgressManager::getInstance().addViewer(cpViewer);
  
  try
  {
    std::cout << "Watercourse example using APP module." << std::endl << std::endl;

// open input data
    std::map<std::string, std::string> connInfo;
    connInfo["URI"] = TERRALIB_DATA_DIR"/app/SaoSebastiao_pol.shp";

    std::auto_ptr<te::da::DataSource> ds = te::da::DataSourceFactory::make("OGR");
    ds->setConnectionInfo(connInfo);
    ds->open();

    std::auto_ptr<te::da::DataSet> dataSet = ds->getDataSet("SaoSebastiao_pol");
    std::auto_ptr<te::da::DataSetType> dataSetType = ds->getDataSetType("SaoSebastiao_pol");

    if(!dataSet.get())
    {
      std::cout << std::endl << "Invalid input data in Watercourse()." << std::endl;
      return;
    }
    
    std::size_t gpos = te::da::GetFirstPropertyPos(dataSet.get(), te::dt::GEOMETRY_TYPE);
    te::gm::GeometryProperty* geomProp = te::da::GetFirstGeomProperty(dataSetType.get());

    std::vector<te::gm::Geometry*> geomVec;

    dataSet->moveBeforeFirst();

    //get geometries
    te::gm::Geometry* geomSeed = 0;

    te::gm::GeometryCollection* geomColl = new te::gm::GeometryCollection(0, te::gm::GeometryCollectionType, geomProp->getSRID());

    while(dataSet->moveNext())
    {
      std::auto_ptr<te::gm::Geometry> g(dataSet->getGeometry(gpos));

      if(geomSeed == 0)
      {
        geomSeed = g.release();
      }
      else
      {
        geomColl->add(g.release());
      }
    }

    //union of all geometries 
    te::gm::Geometry* resultGeom = geomSeed->Union(geomColl);

    //get polygons
    te::app::Multi2Single(resultGeom, geomVec);

    std::vector<te::gm::Polygon*> polyVec;

    for(std::size_t t = 0; t < geomVec.size(); ++t)
    {
      te::gm::Polygon* p = dynamic_cast<te::gm::Polygon*>(geomVec[t]);

      polyVec.push_back(p);
    }

    bool executeok = false;
    bool initok = false;
    {
      std::cout << "Executing Watercourse APP Operation" << std::endl;

// create output data
      std::map<std::string, std::string> orinfo;
      orinfo["URI"] = TERRALIB_DATA_DIR"/app/SaoSebastiao_pol_result.shp";

// create elevation algorithm parameters
      te::app::Watercourse::InputParameters inputParameters;
      inputParameters.m_riversPolygons = polyVec;
      inputParameters.m_fiscalModule = 10000.;

      te::app::Watercourse::OutputParameters outputParameters;
      outputParameters.m_createdOutDSName = "SaoSebastiao_pol_result";
      outputParameters.m_createdOutInfo = orinfo;
      outputParameters.m_createdOutDSType = "OGR";

// execute the algorithm
      te::app::Watercourse wc;

      initok = wc.initialize(inputParameters);

      if(initok)
        executeok = wc.execute(outputParameters);

      if(!executeok)
        std::cout << "Problems in Watercourse operation." << std::endl;
    }

    if (executeok)
      std::cout << "Done!" << std::endl << std::endl;
    else
      std::cout << "Problems in Watercourse." << std::endl;

    te::common::FreeContents(geomVec);
  }
  catch(const std::exception& e)
  {
    std::cout << std::endl << "An exception has occurred in Watercourse(): " << e.what() << std::endl;
  }
  catch(...)
  {
    std::cout << std::endl << "An unexpected exception has occurred in Watercourse()!" << std::endl;
  }

  te::common::ProgressManager::getInstance().removeViewer(id);

  delete cpViewer;
}
