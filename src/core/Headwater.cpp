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
  \file terralib/app/Headwater.cpp
 
  \brief Generates the APP information from a headwater information.
*/

#include "Headwater.h"
#include "Utils.h"

// TerraLib
#include <terralib/common/progress/TaskProgress.h>
#include <terralib/geometry/LineString.h>
#include <terralib/geometry/Point.h>


te::app::Headwater::InputParameters::InputParameters()
{
  reset();
}

te::app::Headwater::InputParameters::~InputParameters()
{
  reset();
}

void te::app::Headwater::InputParameters::reset() throw(te::common::Exception)
{
  m_riversLines.clear();
  m_watercourseBuffer = 30;
  m_headwaterBuffer = 50.;
}

te::common::AbstractParameters* te::app::Headwater::InputParameters::clone() const
{
  return new InputParameters(*this);
}

te::app::Headwater::OutputParameters::OutputParameters()
{
  reset();
}

te::app::Headwater::OutputParameters::OutputParameters(const OutputParameters& other)
{
  reset();

  operator=(other);
}

te::app::Headwater::OutputParameters::~OutputParameters()
{
  reset();
}

void te::app::Headwater::OutputParameters::reset() throw(te::common::Exception)
{
  m_createdOutDSName.clear();
  m_createdOutDSType.clear();
  m_createdOutInfo.clear();
}

const te::app::Headwater::OutputParameters& te::app::Headwater::OutputParameters::operator=(const Headwater::OutputParameters& params)
{
  reset();

  m_createdOutDSName = params.m_createdOutDSName;
  m_createdOutDSType = params.m_createdOutDSType;
  m_createdOutInfo = params.m_createdOutInfo;

  return *this;
}

te::common::AbstractParameters* te::app::Headwater::OutputParameters::clone() const
{
  return new OutputParameters(*this);
}

te::app::Headwater::Headwater()
{
  reset();
}

te::app::Headwater::~Headwater()
{
}

bool te::app::Headwater::execute(AlgorithmOutputParameters& outputParams) throw(te::common::Exception)
{
  if(!m_isInitialized)
  {
    m_logMsg = "Algoritm not initialized.";

    return false;
  }

   m_outputParametersPtr = dynamic_cast<Headwater::OutputParameters*>(&outputParams);

  if(!m_outputParametersPtr)
  {
    m_logMsg = "Invalid output parameters.";

    return false;
  }

  std::vector<te::gm::Geometry*> geomVec;

  // Execute the operation
  te::common::TaskProgress task("Generating Headwater APP");
  task.setTotalSteps(m_inputParameters.m_riversLines.size());

  for(std::size_t t = 0; t < m_inputParameters.m_riversLines.size(); ++t)
  {
    te::gm::LineString* river = m_inputParameters.m_riversLines[t];

    te::gm::Point* start = river->getStartPoint();
    te::gm::Point* end = river->getEndPoint();

    bool touchStart = false;
    bool touchEnd = false;

    for(std::size_t q = 0; q < m_inputParameters.m_riversLines.size(); ++q)
    {
      if(t == q)
        continue;

      te::gm::LineString* r = m_inputParameters.m_riversLines[q];

      if(!touchStart)
      {
        touchStart = start->touches(r);
      }

      if(!touchEnd)
      {
        touchEnd = end->touches(r);
      }

      if(touchStart && touchEnd)
        break;
    }

    te::gm::Geometry* geomBuffer = river->buffer(m_inputParameters.m_watercourseBuffer, 16, te::gm::CapButtType);

    if(!touchStart)
    {
      te::gm::Geometry* geomStartBuffer = start->buffer(m_inputParameters.m_headwaterBuffer, 32, te::gm::CapButtType);

      te::gm::Geometry* aux = geomStartBuffer->Union(geomBuffer);

      delete geomBuffer;

      geomBuffer = aux;
    }

    if(!touchEnd)
    {
      te::gm::Geometry* geomEndBuffer = end->buffer(m_inputParameters.m_headwaterBuffer, 32, te::gm::CapButtType);

      te::gm::Geometry* aux = geomEndBuffer->Union(geomBuffer);

      delete geomBuffer;

      geomBuffer = aux;
    }

    geomVec.push_back(geomBuffer);


    if(!task.isActive())
    {
      m_logMsg = "Operation canceled.";

      return false;
    }

    task.pulse();
  }

 //union all
  te::gm::Geometry* geomBufferResult = 0;

  for(std::size_t t = 0; t < geomVec.size(); ++t)
  {
     if(geomBufferResult == 0)
      {
        geomBufferResult = geomVec[t];
      }
      else
      {
        te::gm::Geometry* geomAux = geomBufferResult->Union(geomVec[t]);

        delete geomBufferResult;

        geomBufferResult = geomAux;

        delete geomVec[t];
      }
  }

  geomVec.clear();

  geomVec.push_back(geomBufferResult);

  //create vector representation
  te::app::ExportVector(geomVec, m_outputParametersPtr->m_createdOutDSName, m_outputParametersPtr->m_createdOutDSType,  m_outputParametersPtr->m_createdOutInfo);

  return true;
}

void te::app::Headwater::reset() throw(te::common::Exception)
{
  m_inputParameters.reset();
  m_outputParametersPtr = 0;
}

bool te::app::Headwater::initialize(const AlgorithmInputParameters& inputParams) throw(te::common::Exception)
{
  reset();

  Headwater::InputParameters const* inputParamsPtr = dynamic_cast<Headwater::InputParameters const*>(&inputParams);

  if(!inputParamsPtr)
  {
    m_logMsg = "Invalid Parameters.";

    return false;
  }

  if(inputParamsPtr->m_riversLines.empty())
  {
    m_logMsg = "Invalid number of geometries";

    return false;
  }

  m_inputParameters = *inputParamsPtr;

  m_isInitialized = true;

  return true;
}
