#include "APPExamples.h"

// TerraLib
#include <terralib/common/progress/ConsoleProgressViewer.h>
#include <terralib/common/progress/ProgressManager.h>
#include <terralib/dataaccess.h>
#include <terralib/raster.h>

#include <terralib/app/core/HighSlope.h>

// STL
#include <string>
#include <map>
#include <iostream>

void HighSlope()
{
  te::common::ConsoleProgressViewer* cpViewer = new te::common::ConsoleProgressViewer();
  int id = te::common::ProgressManager::getInstance().addViewer(cpViewer);
  
  try
  {
    std::cout << "High Slope example using APP module." << std::endl << std::endl;

// open input raster
    std::map<std::string, std::string> rinfo;
    rinfo["URI"] = TERRALIB_DATA_DIR"/app/srtm_28_15_slope.tif";

    te::rst::Raster* rin = te::rst::RasterFactory::open(rinfo);

    bool executeok = false;
    bool initok = false;
    {
      std::cout << "Executing High Slope APP Operation" << std::endl;

// create output raster
      std::map<std::string, std::string> orinfo;
      orinfo["URI"] = TERRALIB_DATA_DIR"/app/srtm_28_15_highSlope.shp";

// create elevation algorithm parameters
      te::app::HighSlope::InputParameters inputParameters;

      inputParameters.m_inRasterBand = 0;
      inputParameters.m_inRasterPtr = rin;
      inputParameters.m_slopeValue = 1.;

      te::app::HighSlope::OutputParameters outputParameters;
      outputParameters.m_createdOutDSName = "srtm_28_15_highSlope";
      outputParameters.m_createdOutInfo = orinfo;
      outputParameters.m_createdOutDSType = "OGR";

// execute the algorithm
      te::app::HighSlope hs;

      initok = hs.initialize(inputParameters);

      if(initok)
        executeok = hs.execute(outputParameters);

      if(!executeok)
        std::cout << "Problems in High Slope operation." << std::endl;
    }

// clean up
    delete rin;

    if (executeok)
      std::cout << "Done!" << std::endl << std::endl;
    else
      std::cout << "Problems in High Slope." << std::endl;
  }
  catch(const std::exception& e)
  {
    std::cout << std::endl << "An exception has occurred in HighSlope(): " << e.what() << std::endl;
  }
  catch(...)
  {
    std::cout << std::endl << "An unexpected exception has occurred in HighSlope()!" << std::endl;
  }

  te::common::ProgressManager::getInstance().removeViewer(id);

  delete cpViewer;
}
