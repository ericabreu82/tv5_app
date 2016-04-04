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
  \file terralib/app/Veredas.cpp
 
  \brief Generates the APP information from a veredas information.
*/

#include "Veredas.h"
#include "Utils.h"

// TerraLib
#include <terralib/common/progress/TaskProgress.h>
#include <terralib/geometry/Geometry.h>


#define DEFAULT_SRID_METERS 29193

te::app::Veredas::InputParameters::InputParameters()
{
  reset();
}

te::app::Veredas::InputParameters::~InputParameters()
{
  reset();
}

void te::app::Veredas::InputParameters::reset() throw(te::common::Exception)
{
  m_geometries.clear();
  m_bufferDistance = 50.;
}

const te::app::Veredas::InputParameters& te::app::Veredas::InputParameters::operator=(const te::app::Veredas::InputParameters& params)
{
  reset();

  m_geometries = params.m_geometries;
  m_bufferDistance = params.m_bufferDistance;

  return *this;
}

te::common::AbstractParameters* te::app::Veredas::InputParameters::clone() const
{
  return new InputParameters(*this);
}

te::app::Veredas::OutputParameters::OutputParameters()
{
  reset();
}

te::app::Veredas::OutputParameters::OutputParameters(const OutputParameters& other)
{
  reset();

  operator=(other);
}

te::app::Veredas::OutputParameters::~OutputParameters()
{
  reset();
}

void te::app::Veredas::OutputParameters::reset() throw(te::common::Exception)
{
  m_createdOutDSName.clear();
  m_createdOutDSType.clear();
  m_createdOutInfo.clear();
}

const te::app::Veredas::OutputParameters& te::app::Veredas::OutputParameters::operator=(const Veredas::OutputParameters& params)
{
  reset();

  m_createdOutDSName = params.m_createdOutDSName;
  m_createdOutDSType = params.m_createdOutDSType;
  m_createdOutInfo = params.m_createdOutInfo;

  return *this;
}

te::common::AbstractParameters* te::app::Veredas::OutputParameters::clone() const
{
  return new OutputParameters(*this);
}

te::app::Veredas::Veredas()
{
  reset();
}

te::app::Veredas::~Veredas()
{
}

bool te::app::Veredas::execute(AlgorithmOutputParameters& outputParams) throw(te::common::Exception)
{
  if(!m_isInitialized)
  {
    m_logMsg = "Algoritm not initialized.";

    return false;
  }

   m_outputParametersPtr = dynamic_cast<Veredas::OutputParameters*>(&outputParams);

  if(!m_outputParametersPtr)
  {
    m_logMsg = "Invalid output parameters.";

    return false;
  }

  std::vector<te::gm::Geometry*> geomBufferVec;

  // Execute the operation
  te::common::TaskProgress task("Generating Veredas APP");
  task.setTotalSteps(m_inputParameters.m_geometries.size());

  for(std::size_t t = 0; t < m_inputParameters.m_geometries.size(); ++t)
  {
    te::gm::Geometry* geom = m_inputParameters.m_geometries[t];

    //remap
    bool remap = false;
    int originalSRID = geom->getSRID();

    if(originalSRID != DEFAULT_SRID_METERS)
    {
      geom->transform(DEFAULT_SRID_METERS);
      remap = true;
    }

    //create buffer
    te::gm::Geometry* geomBuffer = geom->buffer(m_inputParameters.m_bufferDistance);

    if(remap)
    {
      geomBuffer->transform(originalSRID);
    }

    geomBufferVec.push_back(geomBuffer);

    if(!task.isActive())
    {
      m_logMsg = "Operation canceled.";

      return false;
    }

    task.pulse();
  }

  //export vector representation
  te::app::ExportVector(geomBufferVec, m_outputParametersPtr->m_createdOutDSName, m_outputParametersPtr->m_createdOutDSType,  m_outputParametersPtr->m_createdOutInfo);

  return true;
}

void te::app::Veredas::reset() throw(te::common::Exception)
{
  m_inputParameters.reset();
  m_outputParametersPtr = 0;
}

bool te::app::Veredas::initialize(const AlgorithmInputParameters& inputParams) throw(te::common::Exception)
{
  reset();

  Veredas::InputParameters const* inputParamsPtr = dynamic_cast<Veredas::InputParameters const*>(&inputParams);

  if(!inputParamsPtr)
  {
    m_logMsg = "Invalid Parameters.";

    return false;
  }

  if(inputParamsPtr->m_geometries.empty())
  {
    m_logMsg = "Input geometries is empty.";

    return false;
  }

  if(inputParamsPtr->m_bufferDistance < 0.)
  {
    m_logMsg = "Invalid buffer distance.";

    return false;
  }

  m_inputParameters = *inputParamsPtr;

  m_isInitialized = true;

  return true;
}
