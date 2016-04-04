#include "APPExamples.h"

// TerraLib
#include <terralib/common/progress/ConsoleProgressViewer.h>
#include <terralib/common/progress/ProgressManager.h>
#include <terralib/dataaccess.h>
#include <terralib/raster.h>

#include <terralib/app/core/Elevation.h>

// STL
#include <string>
#include <map>
#include <iostream>

void Elevation()
{
  te::common::ConsoleProgressViewer* cpViewer = new te::common::ConsoleProgressViewer();
  int id = te::common::ProgressManager::getInstance().addViewer(cpViewer);
  
  try
  {
    std::cout << "Elevation example using APP module." << std::endl << std::endl;

// open input raster
    std::map<std::string, std::string> rinfo;
    rinfo["URI"] = TERRALIB_DATA_DIR"/app/srtm_28_17.tif";

    te::rst::Raster* rin = te::rst::RasterFactory::open(rinfo);

    bool executeok = false;
    bool initok = false;
    {
      std::cout << "Executing Elevation APP Operation" << std::endl;

// create output raster
      std::map<std::string, std::string> orinfo;
      orinfo["URI"] = TERRALIB_DATA_DIR"/app/srtm_28_17_elevation.shp";

// create elevation algorithm parameters
      te::app::Elevation::InputParameters inputParameters;
      inputParameters.m_elevationValue = 1800;
      inputParameters.m_inRasterBand = 0;
      inputParameters.m_inRasterPtr = rin;

      te::app::Elevation::OutputParameters outputParameters;
      outputParameters.m_createdOutDSName = "srtm_28_17_elevation";
      outputParameters.m_createdOutInfo = orinfo;
      outputParameters.m_createdOutDSType = "OGR";

// execute the algorithm
      te::app::Elevation elevation;

      initok = elevation.initialize(inputParameters);

      if(initok)
        executeok = elevation.execute(outputParameters);

      if(!executeok)
        std::cout << "Problems in Elevation operation." << std::endl;
    }

// clean up
    delete rin;

    if (executeok)
      std::cout << "Done!" << std::endl << std::endl;
    else
      std::cout << "Problems in Elevation." << std::endl;
  }
  catch(const std::exception& e)
  {
    std::cout << std::endl << "An exception has occurred in Elevation(): " << e.what() << std::endl;
  }
  catch(...)
  {
    std::cout << std::endl << "An unexpected exception has occurred in Elevation()!" << std::endl;
  }

  te::common::ProgressManager::getInstance().removeViewer(id);

  delete cpViewer;
}
