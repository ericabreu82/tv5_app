#include "APPExamples.h"

// TerraLib
#include <terralib/common/progress/ConsoleProgressViewer.h>
#include <terralib/common/progress/ProgressManager.h>
#include <terralib/common/STLUtils.h>
#include <terralib/dataaccess.h>
#include <terralib/dataaccess/dataset/DataSet.h>
#include <terralib/raster.h>
#include <terralib/app/core/Plateaus.h>

// STL
#include <string>
#include <map>
#include <iostream>

void Plateaus()
{
  te::common::ConsoleProgressViewer* cpViewer = new te::common::ConsoleProgressViewer();
  int id = te::common::ProgressManager::getInstance().addViewer(cpViewer);
  
  try
  {
    std::cout << "Plateau example using APP module." << std::endl << std::endl;

// open input raster
    std::map<std::string, std::string> rinfo;
    rinfo["URI"] = TERRALIB_DATA_DIR"/app/srtm_27_15.tif";

    te::rst::Raster* rasterSRTM = te::rst::RasterFactory::open(rinfo);

// open input raster
    rinfo["URI"] = TERRALIB_DATA_DIR"/app/srtm_27_15_slope.tif";

    te::rst::Raster* rasterSlope = te::rst::RasterFactory::open(rinfo);

    bool executeok = false;
    bool initok = false;
    {
      std::cout << "Executing Plateaus APP Operation" << std::endl;

// create output data
      std::map<std::string, std::string> orinfo;
      orinfo["URI"] = TERRALIB_DATA_DIR"/app/srtm_27_15_plateau.shp";

// create Plateaus algorithm parameters
      te::app::Plateaus::InputParameters inputParameters;
      inputParameters.m_bufferDistance = 500.;
      inputParameters.m_demRasterPtr = rasterSRTM;
      inputParameters.m_slopeRasterPtr = rasterSlope;

      te::app::Plateaus::OutputParameters outputParameters;
      outputParameters.m_createdOutDSName = "srtm_27_15_plateau";
      outputParameters.m_createdOutInfo = orinfo;
      outputParameters.m_createdOutDSType = "OGR";

// execute the algorithm
      te::app::Plateaus plateau;

      initok = plateau.initialize(inputParameters);

      if(initok)
        executeok = plateau.execute(outputParameters);

      if(!executeok)
        std::cout << "Problems in Plateaus operation." << std::endl;
    }

// clean up
    delete rasterSRTM;
    delete rasterSlope;

    if (executeok)
      std::cout << "Done!" << std::endl << std::endl;
    else
      std::cout << "Problems in Plateaus." << std::endl;
  }
  catch(const std::exception& e)
  {
    std::cout << std::endl << "An exception has occurred in Plateaus(): " << e.what() << std::endl;
  }
  catch(...)
  {
    std::cout << std::endl << "An unexpected exception has occurred in Plateaus()!" << std::endl;
  }

  te::common::ProgressManager::getInstance().removeViewer(id);

  delete cpViewer;
}
