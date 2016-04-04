/*  Copyright (C) 2001-2009 National Institute For Space Research (INPE) - Brazil.

    This file is part of the TerraLib - a Framework for building GIS enabled applications.

    TerraLib is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License,
    or (at your option) any later version.

    TerraLib is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with TerraLib. See COPYING. If not, write to
    TerraLib Team at <terralib-team@terralib.org>.
 */

/*!
  \file terralib/app/Plateus.cpp
 
  \brief Generates the APP information from a plateaus information.
*/

#include "Plateaus.h"
#include "Utils.h"

// TerraLib
#include <terralib/common/UnitOfMeasure.h>
#include <terralib/common/UnitsOfMeasureManager.h>
#include <terralib/common/progress/TaskProgress.h>
#include <terralib/geometry/Geometry.h>
#include <terralib/raster/BandProperty.h>
#include <terralib/raster/Grid.h>
#include <terralib/raster/Raster.h>
#include <terralib/raster/RasterFactory.h>
#include <terralib/raster/Utils.h>
#include <terralib/srs/SpatialReferenceSystemManager.h>


//STL Includes
#include <queue>

#define DEFAULT_SRID_METERS 29193

void exportRaster(te::rst::Raster* rasterIn, std::string fileName);

te::app::Plateaus::InputParameters::InputParameters()
{
  reset();
}

te::app::Plateaus::InputParameters::~InputParameters()
{
  reset();
}

void te::app::Plateaus::InputParameters::reset() throw(te::common::Exception)
{
  m_demRasterPtr = 0;
  m_slopeRasterPtr = 0;
  m_bufferDistance = 500;
}

te::common::AbstractParameters* te::app::Plateaus::InputParameters::clone() const
{
  return new InputParameters(*this);
}

te::app::Plateaus::OutputParameters::OutputParameters()
{
  reset();
}

te::app::Plateaus::OutputParameters::OutputParameters(const OutputParameters& other)
{
  reset();

  operator=(other);
}

te::app::Plateaus::OutputParameters::~OutputParameters()
{
  reset();
}

void te::app::Plateaus::OutputParameters::reset() throw(te::common::Exception)
{
  m_createdOutDSName.clear();
  m_createdOutDSType.clear();
  m_createdOutInfo.clear();

  m_createdOutRasterDSType.clear();
  m_createdOutRasterInfo.clear();
}

const te::app::Plateaus::OutputParameters& te::app::Plateaus::OutputParameters::operator=(const Plateaus::OutputParameters& params)
{
  reset();

  m_createdOutDSName = params.m_createdOutDSName;
  m_createdOutDSType = params.m_createdOutDSType;
  m_createdOutInfo = params.m_createdOutInfo;

  m_createdOutRasterDSType = params.m_createdOutRasterDSType;
  m_createdOutRasterInfo = params.m_createdOutRasterInfo;

  return *this;
}

te::common::AbstractParameters* te::app::Plateaus::OutputParameters::clone() const
{
  return new OutputParameters(*this);
}

te::app::Plateaus::Plateaus()
{
  reset();
}

te::app::Plateaus::~Plateaus()
{
}

bool te::app::Plateaus::execute(AlgorithmOutputParameters& outputParams) throw(te::common::Exception)
{
  if(!m_isInitialized)
  {
    m_logMsg = "Algoritm not initialized.";

    return false;
  }

   m_outputParametersPtr = dynamic_cast<Plateaus::OutputParameters*>(&outputParams);

  if(!m_outputParametersPtr)
  {
    m_logMsg = "Invalid output parameters.";

    return false;
  }

  //first create escarp raster
  std::auto_ptr<te::rst::Raster> escarpRaster = createSlopeMask(1.);

  if(!escarpRaster.get())
  {
    return false;
  }

  exportRaster(escarpRaster.get(), "escarpRaster.tif");

  //second create high elevation low slope mask
  std::auto_ptr<te::rst::Raster> highElevationSlopeRaster = createHighElevationSlopeMask(600., 0.1);

  if(!highElevationSlopeRaster.get())
  {
    return false;
  }

  exportRaster(highElevationSlopeRaster.get(), "highElevationSlopeRaster.tif");

  //create terrace raster mask
  std::auto_ptr<te::rst::Raster> terraceRaster = createTerraceRasterMask(highElevationSlopeRaster.get(), 10.);

  if(!terraceRaster.get())
  {
    return false;
  }

  exportRaster(terraceRaster.get(), "terraceRaster.tif");

  //last step define terrace mask
  if(!defineTerracesfromEscarps(terraceRaster.get(), escarpRaster.get()))
  {
    return false;
  }

  //exportRaster(m_outputParametersPtr->m_createdOutRasterPtr.get(), "plateaus.tif");

  //create vector representation
  te::app::Raster2Vector(m_outputParametersPtr->m_createdOutRasterPtr.get(), 0,  m_outputParametersPtr->m_createdOutDSName, m_outputParametersPtr->m_createdOutDSType,  m_outputParametersPtr->m_createdOutInfo);

  return true;
}

void te::app::Plateaus::reset() throw(te::common::Exception)
{
  m_inputParameters.reset();
  m_outputParametersPtr = 0;
}

bool te::app::Plateaus::initialize(const AlgorithmInputParameters& inputParams) throw(te::common::Exception)
{
  reset();

  Plateaus::InputParameters const* inputParamsPtr = dynamic_cast<Plateaus::InputParameters const*>(&inputParams);

  if(!inputParamsPtr)
  {
    m_logMsg = "Invalid Parameters.";

    return false;
  }

  if(!inputParamsPtr->m_demRasterPtr || inputParamsPtr->m_demRasterPtr->getAccessPolicy() == te::common::NoAccess)
  {
    m_logMsg = "Invalid Raster.";

    return false;
  }

  if(!inputParamsPtr->m_slopeRasterPtr || inputParamsPtr->m_slopeRasterPtr->getAccessPolicy() == te::common::NoAccess)
  {
    m_logMsg = "Invalid Raster.";

    return false;
  }

  m_inputParameters = *inputParamsPtr;

  m_isInitialized = true;

  return true;
}

std::auto_ptr<te::rst::Raster> te::app::Plateaus::createSlopeMask(double threshold_slope)
{
  double slope_undef = -9999.;

  //create escarp raster
  std::vector<te::rst::BandProperty*> bandsProperties;
  te::rst::BandProperty* bandProp = new te::rst::BandProperty(0, te::dt::UCHAR_TYPE);
  bandProp->m_noDataValue = 255.;
  bandsProperties.push_back(bandProp);

  te::rst::Grid* grid = new te::rst::Grid(*( m_inputParameters.m_demRasterPtr->getGrid()));

  std::map<std::string, std::string> rinfo;

  te::rst::Raster* outRaster = te::rst::RasterFactory::make("MEM", grid, bandsProperties, rinfo);

  if(!outRaster)
  {
    throw(te::common::Exception("Terrace raster mask creation error."));
  }

  te::common::TaskProgress task("Creating Slope Mask");
  task.setTotalSteps(outRaster->getNumberOfRows());

  for(unsigned int i = 0; i < outRaster->getNumberOfRows(); ++i)
  {
    for(unsigned int j = 0; j < outRaster->getNumberOfColumns(); ++j)
    {
      double slope_value;

      m_inputParameters.m_slopeRasterPtr->getValue(j, i, slope_value);

      if(slope_value != slope_undef && slope_value >= threshold_slope)
      {
        outRaster->setValue(j, i, 0., 0);
      }
      else
      {
        outRaster->setValue(j, i, 255., 0);
      }
    }

    task.pulse();
  }

  //return raster
  std::auto_ptr<te::rst::Raster> escarpRaster(outRaster);

  return escarpRaster;
}

std::auto_ptr<te::rst::Raster> te::app::Plateaus::createHighElevationSlopeMask(double dem_threshold, double slope_threshold)
{
  //create high elevation low slope raster
  std::vector<te::rst::BandProperty*> bandsProperties;
  te::rst::BandProperty* bandProp = new te::rst::BandProperty(0, te::dt::UCHAR_TYPE);
  bandProp->m_noDataValue = 0.;
  bandsProperties.push_back(bandProp);

  te::rst::Grid* grid = new te::rst::Grid(*( m_inputParameters.m_demRasterPtr->getGrid()));

  std::map<std::string, std::string> rinfo;

  te::rst::Raster* outRaster = te::rst::RasterFactory::make("MEM", grid, bandsProperties, rinfo);

  if(!outRaster)
  {
    throw(te::common::Exception("Terrace raster mask creation error."));
  }

  double slope_undef = -9999.;

  te::common::TaskProgress task("Creating High Elevation Slope Mask");
  task.setTotalSteps(outRaster->getNumberOfRows());

  for(unsigned int i = 0; i < outRaster->getNumberOfRows(); ++i)
  {
    for(unsigned int j = 0; j < outRaster->getNumberOfColumns(); ++j)
    {
      double slope_value, dem_value;

      // read elevation
      m_inputParameters.m_demRasterPtr->getValue(j, i, dem_value);

      // read slope
      m_inputParameters.m_slopeRasterPtr->getValue(j, i, slope_value);

      // if undef, do nothing
      if(slope_value == slope_undef)
      {
        outRaster->setValue(j, i, 0., 0);
        continue;
      }

      // if both thresholds are ok
      if(dem_value >= dem_threshold && slope_value < slope_threshold)
      {
        outRaster->setValue(j, i, 1., 0);
      }
      else
      {
        outRaster->setValue(j, i, 0., 0);
      }
    }

    task.pulse();
  }

  //return raster
  std::auto_ptr<te::rst::Raster> highElevationLowSlopeRaster(outRaster);

  return highElevationLowSlopeRaster;
}

std::auto_ptr<te::rst::Raster> te::app::Plateaus::createTerraceRasterMask(te::rst::Raster* highElevationSlopeMask, double area_threshold)
{
  int dx[8] = {1, 1, 1, 0, -1, -1, -1, 0};
  int dy[8] = {1, 0, -1, -1, -1, 0, 1, 1};

  //create terrace raster
  std::vector<te::rst::BandProperty*> bandsProperties;
  te::rst::BandProperty* bandProp = new te::rst::BandProperty(0, te::dt::UCHAR_TYPE);
  bandProp->m_noDataValue = -1;
  bandsProperties.push_back(bandProp);

  te::rst::Grid* grid = new te::rst::Grid(*( m_inputParameters.m_demRasterPtr->getGrid()));

  std::map<std::string, std::string> rinfo;

  te::rst::Raster* outRaster = te::rst::RasterFactory::make("MEM", grid, bandsProperties, rinfo);

  if(!outRaster)
  {
    throw(te::common::Exception("Terrace raster mask creation error."));
  }

  //fill raster with dummy values
  for(unsigned int i = 0; i < outRaster->getNumberOfRows(); ++i)
  {
    for(unsigned int j = 0; j < outRaster->getNumberOfColumns(); ++j)
    {
      outRaster->setValue(j, i, - 1., 0);
    }
  }

  //create terrace mask boolean matrix
  std::vector< std::vector<bool> > terrace_mask_visited;

  for(unsigned int i = 0;i < m_inputParameters.m_demRasterPtr->getNumberOfRows();i++)
  {
    std::vector<bool> row;
    row.resize(m_inputParameters.m_demRasterPtr->getNumberOfColumns(), false);

    terrace_mask_visited.push_back(row);
  }

  double terrace_undef = 0.;

  //calculate cell area
  te::common::UnitOfMeasurePtr unit = te::srs::SpatialReferenceSystemManager::getInstance().getUnit(m_inputParameters.m_demRasterPtr->getSRID());
  std::string unitSymbol = unit->getSymbol();
  std::string unitDefault = "metre";

  //double conversionFactor = te::common::UnitsOfMeasureManager::getInstance().getConversion(unitSymbol, unitDefault);

  double conversionFactor = 111000.;

  double resX_meter = m_inputParameters.m_demRasterPtr->getGrid()->getResolutionX() * conversionFactor;
  double resY_meter = m_inputParameters.m_demRasterPtr->getGrid()->getResolutionY() * conversionFactor;

  double cell_area = resX_meter * resY_meter;

  te::common::TaskProgress task("Creating Terrace Mask");
  task.setTotalSteps(m_inputParameters.m_demRasterPtr->getNumberOfRows());

  //fill terrace raster
  for(unsigned int i = 0; i < m_inputParameters.m_demRasterPtr->getNumberOfRows(); ++i)
  {
    for(unsigned int j = 0; j < m_inputParameters.m_demRasterPtr->getNumberOfColumns(); ++j)
    {
      if(terrace_mask_visited[i][j] == true)
        continue;

      double terrace_value;

      // read terrace mask
      highElevationSlopeMask->getValue(j, i, terrace_value);

      // if terrace mask is valid
      if(terrace_value != terrace_undef)
      {
        std::queue<te::gm::Coord2D> grid_points_queue;
        std::vector<te::gm::Coord2D> terrace_points;

        te::gm::Coord2D grid_point(j, i);

        terrace_mask_visited[i][j] = true;

        terrace_points.push_back(grid_point);
        grid_points_queue.push(grid_point);
        
        while(grid_points_queue.empty() == false)
        {
          grid_point = grid_points_queue.front();
          grid_points_queue.pop();

          int x = (int)grid_point.getX();
          int y = (int)grid_point.getY();

          for(int k = 0; k < 8; k++)
          {
            if(y+dy[k] == -1 || y+dy[k] == m_inputParameters.m_demRasterPtr->getNumberOfRows() || x+dx[k] == -1 || x+dx[k] == m_inputParameters.m_demRasterPtr->getNumberOfColumns())
              continue;

            if(terrace_mask_visited[y+dy[k]][x+dx[k]] == true)
              continue;

            // read terrace mask
            highElevationSlopeMask->getValue(x+dx[k], y+dy[k], terrace_value);

            if(terrace_value != terrace_undef)
            {
              grid_point.x = x+dx[k];
              grid_point.y = y+dy[k];

              terrace_mask_visited[y+dy[k]][x+dx[k]] = true;

              terrace_points.push_back(grid_point);
              grid_points_queue.push(grid_point);
            }
          }
        }

        // plain area in m^2
        double terrace_area = cell_area * terrace_points.size();
        
        // convert from m^2 to ha
        terrace_area /= 10000.;

        if(terrace_area > area_threshold)
        {
          for(std::size_t k = 0; k < terrace_points.size(); k++)
          {
            int x = (int) terrace_points[k].getX();
            int y = (int) terrace_points[k].getY();

            outRaster->setValue(x, y, 1., 0);
          }
        }
      }
    }

    task.pulse();
  }

  //return raster
  std::auto_ptr<te::rst::Raster> terraceRaster(outRaster);

  return terraceRaster;
}

bool te::app::Plateaus::defineTerracesfromEscarps(te::rst::Raster* terraceRaster, te::rst::Raster* escarpRaster)
{
  // Initializing the output raster
  if(m_outputParametersPtr->m_createdOutRasterDSType.empty())
    m_outputParametersPtr->m_createdOutRasterDSType = "MEM";

  std::vector<te::rst::BandProperty*> bandsProperties;
  te::rst::BandProperty* bandProp = new te::rst::BandProperty(0, te::dt::UCHAR_TYPE);
  bandProp->m_noDataValue = 255.;
  bandsProperties.push_back(bandProp);

  te::rst::Grid* grid = new te::rst::Grid(*( m_inputParameters.m_demRasterPtr->getGrid()));

  te::rst::Raster* outRaster = te::rst::RasterFactory::make(m_outputParametersPtr->m_createdOutRasterDSType, grid, bandsProperties, m_outputParametersPtr->m_createdOutRasterInfo);

  //fill raster with dummy values
  for(unsigned int i = 0; i < outRaster->getNumberOfRows(); ++i)
  {
    for(unsigned int j = 0; j < outRaster->getNumberOfColumns(); ++j)
    {
      outRaster->setValue(j, i, 255., 0);
    }
  }

  m_outputParametersPtr->m_createdOutRasterPtr.reset(outRaster);

  if(!m_outputParametersPtr->m_createdOutRasterPtr.get())
  {
    throw(te::common::Exception("Output raster creation error."));
  }

  // Execute the operation
  te::common::TaskProgress task("Generating Terraces APP");
  task.setTotalSteps(m_inputParameters.m_demRasterPtr->getNumberOfRows());

  //calculate window
  te::common::UnitOfMeasurePtr unit = te::srs::SpatialReferenceSystemManager::getInstance().getUnit(m_inputParameters.m_demRasterPtr->getSRID());
  std::string unitSymbol = unit->getSymbol();
  std::string unitDefault = "metre";

  //double conversionFactor = te::common::UnitsOfMeasureManager::getInstance().getConversion(unitSymbol, unitDefault);
  double conversionFactor = 111000.;

  double resX_meter = m_inputParameters.m_demRasterPtr->getGrid()->getResolutionX() * conversionFactor;
  double resY_meter = m_inputParameters.m_demRasterPtr->getGrid()->getResolutionY() * conversionFactor;

  double res = std::max(resX_meter, resY_meter);

  int wSize = (int) (m_inputParameters.m_bufferDistance / res);
  int bSize = (int) (100./res);

  for(unsigned int i = 0; i < m_inputParameters.m_demRasterPtr->getNumberOfRows(); ++i)
  {
    for(unsigned int j = 0; j < m_inputParameters.m_demRasterPtr->getNumberOfColumns(); ++j)
    {
      double escarp_value;

      escarpRaster->getValue(j, i, escarp_value);

      bool found_terrace_point = false;

      if(escarp_value == 0.)
      {
        found_terrace_point = false;

        for(int lin = -wSize; lin <= wSize; lin++)
        {
          for(int col = -wSize; col <= wSize; col++)
          {
            if(lin == 0 && col == 0)
              continue;

            int x = j + col;
            int y = i + lin;

            if(y < 0 || y >= (int)m_inputParameters.m_demRasterPtr->getNumberOfRows() || x < 0 || x >= (int)m_inputParameters.m_demRasterPtr->getNumberOfColumns())
              continue;

            double terrace_value;

            terraceRaster->getValue(x, y, terrace_value);

            if(terrace_value == 1.)
            {
              found_terrace_point = true;
              break;
            }
          }
        }

        if(found_terrace_point == true)
        {
          double escarp_dem_value;

          m_inputParameters.m_demRasterPtr->getValue(j, i, escarp_dem_value);

          for(int lin = -bSize;lin <= bSize;lin++)
          {
            for(int col = -bSize;col <= bSize;col++)
            {
              if(lin == 0 && col == 0)
                continue;

              int x = j + col;
              int y = i + lin;

              if(y < 0 || y >= (int)m_inputParameters.m_demRasterPtr->getNumberOfRows() || x < 0 || x >= (int)m_inputParameters.m_demRasterPtr->getNumberOfColumns())
                continue;

              escarpRaster->getValue(x, y, escarp_value);
              
              if(escarp_value == 0.)
                continue;

              te::gm::Coord2D c1 = m_inputParameters.m_demRasterPtr->getGrid()->gridToGeo(j, i);
              te::gm::Coord2D c2 = m_inputParameters.m_demRasterPtr->getGrid()->gridToGeo(x, y);
              
              te::gm::Point p1(c1.getX(), c1.getY(), m_inputParameters.m_demRasterPtr->getSRID());
              te::gm::Point p2(c2.getX(), c2.getY(), m_inputParameters.m_demRasterPtr->getSRID());

              double neighbor_dem_value;

              m_inputParameters.m_demRasterPtr->getValue(x, y, neighbor_dem_value);
              
              if(neighbor_dem_value > escarp_dem_value && p1.distance(&p2) <= 100.)
              {
                 m_outputParametersPtr->m_createdOutRasterPtr->setValue(x, y, 1., 0);
              }
            }
          }
        }
      }
    }

    if(!task.isActive())
    {
      m_logMsg = "Operation canceled.";

      return false;
    }

    task.pulse();
  }

  return true;
}

void exportRaster(te::rst::Raster* rasterIn, std::string fileName)
{
  assert(rasterIn);

  std::string uri = TERRALIB_DATA_DIR"/" + fileName;

  te::rst::CreateCopy(*rasterIn, uri);
}