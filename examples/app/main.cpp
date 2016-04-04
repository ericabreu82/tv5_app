// Examples
#include "APPExamples.h"

// TerraLib
#include <terralib/common.h>
#include <terralib/plugin.h>

// STL
#include <cstdlib>
#include <iostream>

int main()
{
  try
  {
    TerraLib::getInstance().initialize();

    LoadModules();

    //Elevation();

    HighSlope();

    //Veredas();

    //Lakes();

    //Plateaus();

    //Headwater();

    //Watercourse();

    //RidgeLines();

    //HillTop();

    te::plugin::PluginManager::getInstance().unloadAll();

    TerraLib::getInstance().finalize();
  }
  catch(const std::exception& e)
  {
    std::cout << std::endl << "An exception has occurred in APP examples: " << e.what() << std::endl;

    return EXIT_FAILURE;
  }
  catch(...)
  {
    std::cout << std::endl << "An unexpected exception has occurred in APP examples!" << std::endl;

    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
