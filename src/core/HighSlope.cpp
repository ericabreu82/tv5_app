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
  \file terralib/app/HighSlope.cpp
 
  \brief Generates the APP information from a high slope raster.
*/

#include "HighSlope.h"
#include "Utils.h"

// TerraLib
#include <terralib/common/progress/TaskProgress.h>
#include <terralib/raster/Band.h>
#include <terralib/raster/BandProperty.h>
#include <terralib/raster/Grid.h>
#include <terralib/raster/Raster.h>
#include <terralib/raster/RasterFactory.h>


te::app::HighSlope::InputParameters::InputParameters()
{
  reset();
}

te::app::HighSlope::InputParameters::~InputParameters()
{
  reset();
}

void te::app::HighSlope::InputParameters::reset() throw(te::common::Exception)
{
  m_inRasterPtr = 0;
  m_inRasterBand = 0;
  m_slopeValue = 1.;
  m_slopeUndef = -9999.;
}

te::common::AbstractParameters* te::app::HighSlope::InputParameters::clone() const
{
  return new InputParameters(*this);
}

te::app::HighSlope::OutputParameters::OutputParameters()
{
  reset();
}

te::app::HighSlope::OutputParameters::OutputParameters(const OutputParameters& other)
{
  reset();

  operator=(other);
}

te::app::HighSlope::OutputParameters::~OutputParameters()
{
  reset();
}

void te::app::HighSlope::OutputParameters::reset() throw(te::common::Exception)
{
  m_createdOutDSName.clear();
  m_createdOutDSType.clear();
  m_createdOutInfo.clear();

  m_createdOutRasterDSType.clear();
  m_createdOutRasterInfo.clear();
}

const te::app::HighSlope::OutputParameters& te::app::HighSlope::OutputParameters::operator=(const HighSlope::OutputParameters& params)
{
  reset();

  m_createdOutDSName = params.m_createdOutDSName;
  m_createdOutDSType = params.m_createdOutDSType;
  m_createdOutInfo = params.m_createdOutInfo;

  m_createdOutRasterDSType = params.m_createdOutRasterDSType;
  m_createdOutRasterInfo = params.m_createdOutRasterInfo;

  return *this;
}

te::common::AbstractParameters* te::app::HighSlope::OutputParameters::clone() const
{
  return new OutputParameters(*this);
}

te::app::HighSlope::HighSlope()
{
  reset();
}

te::app::HighSlope::~HighSlope()
{
}

bool te::app::HighSlope::execute(AlgorithmOutputParameters& outputParams) throw(te::common::Exception)
{
  if(!m_isInitialized)
  {
    m_logMsg = "Algoritm not initialized.";

    return false;
  }

   m_outputParametersPtr = dynamic_cast<HighSlope::OutputParameters*>(&outputParams);

  if(!m_outputParametersPtr)
  {
    m_logMsg = "Invalid output parameters.";

    return false;
  }

  // Initializing the output raster
  if(m_outputParametersPtr->m_createdOutRasterDSType.empty())
    m_outputParametersPtr->m_createdOutRasterDSType = "MEM";

  std::vector< te::rst::BandProperty* > bandsProperties;
  te::rst::BandProperty* bandProp = new te::rst::BandProperty(0, te::dt::UCHAR_TYPE);
  bandProp->m_noDataValue = 255.;
  bandsProperties.push_back(bandProp);

  te::rst::Grid* grid = new te::rst::Grid(*( m_inputParameters.m_inRasterPtr->getGrid()));

  te::rst::Raster* outRaster = te::rst::RasterFactory::make(m_outputParametersPtr->m_createdOutRasterDSType, grid, bandsProperties, m_outputParametersPtr->m_createdOutRasterInfo);

  m_outputParametersPtr->m_createdOutRasterPtr.reset(outRaster);

  if(!m_outputParametersPtr->m_createdOutRasterPtr.get())
  {
    throw(te::common::Exception("Output raster creation error."));
  }

  //calculate slope
  std::auto_ptr<te::rst::Raster> slopeRaster = te::app::CalculateSlope(m_inputParameters.m_inRasterPtr);

  // Execute the operation
  te::common::TaskProgress task("Generating High Slope APP");
  task.setTotalSteps(slopeRaster->getNumberOfRows());

  for (unsigned int i = 0; i < slopeRaster->getNumberOfRows(); ++i)
  {
    for (unsigned int j = 0; j < slopeRaster->getNumberOfColumns(); ++j)
    {
      double value;

      slopeRaster->getValue(j, i, value, 0);

      if(value != m_inputParameters.m_slopeUndef && value > m_inputParameters.m_slopeValue)
      {
        m_outputParametersPtr->m_createdOutRasterPtr->setValue(j, i, 0., 0);
      }
      else
      {
        m_outputParametersPtr->m_createdOutRasterPtr->setValue(j, i, 255., 0);
      }
    }

    if(!task.isActive())
    {
      m_logMsg = "Operation canceled.";

      return false;
    }

    task.pulse();
  }

  //create vector representation
  te::app::Raster2Vector(m_outputParametersPtr->m_createdOutRasterPtr.get(), 0,  m_outputParametersPtr->m_createdOutDSName, m_outputParametersPtr->m_createdOutDSType,  m_outputParametersPtr->m_createdOutInfo);

  return true;
}

void te::app::HighSlope::reset() throw(te::common::Exception)
{
  m_inputParameters.reset();
  m_outputParametersPtr = 0;
}

bool te::app::HighSlope::initialize(const AlgorithmInputParameters& inputParams) throw(te::common::Exception)
{
  reset();

  HighSlope::InputParameters const* inputParamsPtr = dynamic_cast<HighSlope::InputParameters const*>(&inputParams);

  if(!inputParamsPtr)
  {
    m_logMsg = "Invalid Parameters.";

    return false;
  }

  if(!inputParamsPtr->m_inRasterPtr || inputParamsPtr->m_inRasterPtr->getAccessPolicy() == te::common::NoAccess)
  {
    m_logMsg = "Invalid Raster.";

    return false;
  }

  if(inputParamsPtr->m_inRasterBand < 0 || inputParamsPtr->m_inRasterBand >= inputParamsPtr->m_inRasterPtr->getNumberOfBands())
  {
    m_logMsg = "Invalid Band Value.";

    return false;
  }

  m_inputParameters = *inputParamsPtr;

  m_isInitialized = true;

  return true;
}
