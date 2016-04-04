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
  \file terralib/app/Lakes.cpp
 
  \brief Generates the APP information from a lakes information.
*/

#include "Lakes.h"
#include "Utils.h"

// TerraLib
#include <terralib/common/progress/TaskProgress.h>
#include <terralib/geometry/Geometry.h>
#include <terralib/geometry/MultiPolygon.h>
#include <terralib/geometry/Polygon.h>


#define DEFAULT_SRID_METERS 29193

te::app::Lakes::InputParameters::InputParameters()
{
  reset();
}

te::app::Lakes::InputParameters::~InputParameters()
{
  reset();
}

void te::app::Lakes::InputParameters::reset() throw(te::common::Exception)
{
  m_lakesGeometries.clear();
  m_urbanGeometries.clear();
  m_ruralGeometries.clear();

  m_urbanBufferDistance = 30.;
  m_ruralBigBufferDistance = 100.;
  m_ruralSmallBufferDistance = 50.;
}

const te::app::Lakes::InputParameters& te::app::Lakes::InputParameters::operator=(const te::app::Lakes::InputParameters& params)
{
  reset();

  m_lakesGeometries = params.m_lakesGeometries;
  m_urbanGeometries = params.m_urbanGeometries;
  m_ruralGeometries = params.m_ruralGeometries;

  m_urbanBufferDistance = params.m_urbanBufferDistance;
  m_ruralBigBufferDistance = params.m_ruralBigBufferDistance;
  m_ruralSmallBufferDistance = params.m_ruralSmallBufferDistance;

  return *this;
}

te::common::AbstractParameters* te::app::Lakes::InputParameters::clone() const
{
  return new InputParameters(*this);
}

te::app::Lakes::OutputParameters::OutputParameters()
{
  reset();
}

te::app::Lakes::OutputParameters::OutputParameters(const OutputParameters& other)
{
  reset();

  operator=(other);
}

te::app::Lakes::OutputParameters::~OutputParameters()
{
  reset();
}

void te::app::Lakes::OutputParameters::reset() throw(te::common::Exception)
{
  m_createdOutDSName.clear();
  m_createdOutDSType.clear();
  m_createdOutInfo.clear();
}

const te::app::Lakes::OutputParameters& te::app::Lakes::OutputParameters::operator=(const Lakes::OutputParameters& params)
{
  reset();

  m_createdOutDSName = params.m_createdOutDSName;
  m_createdOutDSType = params.m_createdOutDSType;
  m_createdOutInfo = params.m_createdOutInfo;

  return *this;
}

te::common::AbstractParameters* te::app::Lakes::OutputParameters::clone() const
{
  return new OutputParameters(*this);
}

te::app::Lakes::Lakes()
{
  reset();
}

te::app::Lakes::~Lakes()
{
}

bool te::app::Lakes::execute(AlgorithmOutputParameters& outputParams) throw(te::common::Exception)
{
  if(!m_isInitialized)
  {
    m_logMsg = "Algoritm not initialized.";

    return false;
  }

   m_outputParametersPtr = dynamic_cast<Lakes::OutputParameters*>(&outputParams);

  if(!m_outputParametersPtr)
  {
    m_logMsg = "Invalid output parameters.";

    return false;
  }

  std::vector<te::gm::Geometry*> geomBufferVec;

 // Execute the operation
  te::common::TaskProgress task("Generating Veredas APP");
  task.setTotalSteps(m_inputParameters.m_lakesGeometries.size());

  for(std::size_t t = 0; t < m_inputParameters.m_lakesGeometries.size(); ++t)
  {
    te::gm::Geometry* geom = m_inputParameters.m_lakesGeometries[t];

    //remap
    bool remap = false;
    int originalSRID = geom->getSRID();

    if(originalSRID != DEFAULT_SRID_METERS)
    {
      geom->transform(DEFAULT_SRID_METERS);
      remap = true;
    }

    //create buffer
    te::gm::Geometry* geomBuffer = 0;

    //check intersection with urban geometries
    bool urbanIntersec = false;

    for(std::size_t u = 0; u < m_inputParameters.m_urbanGeometries.size(); ++u)
    {
      te::gm::Geometry* urbanGeom = m_inputParameters.m_urbanGeometries[u];

      if(urbanGeom->getSRID() != DEFAULT_SRID_METERS)
        urbanGeom->transform(DEFAULT_SRID_METERS);

      if(geom->intersects(urbanGeom))
      {
        geomBuffer = geom->buffer(m_inputParameters.m_urbanBufferDistance);

        urbanIntersec = true;

        break;
      }
    }

    //if not urban intersection check rural intersection
    if(!urbanIntersec)
    {
      for(std::size_t r = 0; r < m_inputParameters.m_ruralGeometries.size(); ++r)
      {
        te::gm::Geometry* ruralGeom = m_inputParameters.m_ruralGeometries[r];

        if(ruralGeom->getSRID() != DEFAULT_SRID_METERS)
          ruralGeom->transform(DEFAULT_SRID_METERS);

        if(geom->intersects(ruralGeom))
        {
          double lakeArea = calculateArea(geom);

          //convert area from m^2 to hectare, 1 ha = 10000 m^2
          lakeArea /= 10000.;

          if(lakeArea > 20.)
          {
            geomBuffer = geom->buffer(m_inputParameters.m_ruralBigBufferDistance);
          }
          else
          {
            geomBuffer = geom->buffer(m_inputParameters.m_ruralSmallBufferDistance);
          }

          break;
        }
      }
    }


    if(geomBuffer)
    {
      if(remap)
      {
        geomBuffer->transform(originalSRID);
      }

      geomBufferVec.push_back(geomBuffer);
    }

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

void te::app::Lakes::reset() throw(te::common::Exception)
{
  m_inputParameters.reset();
  m_outputParametersPtr = 0;
}

bool te::app::Lakes::initialize(const AlgorithmInputParameters& inputParams) throw(te::common::Exception)
{
  reset();

  Lakes::InputParameters const* inputParamsPtr = dynamic_cast<Lakes::InputParameters const*>(&inputParams);

  if(!inputParamsPtr)
  {
    m_logMsg = "Invalid Parameters.";

    return false;
  }

  if(inputParamsPtr->m_lakesGeometries.empty())
  {
    m_logMsg = "Input Lakes geometries is empty.";

    return false;
  }

  if(inputParamsPtr->m_urbanBufferDistance < 0.)
  {
    m_logMsg = "Invalid urban buffer distance.";

    return false;
  }

  if(inputParamsPtr->m_ruralBigBufferDistance < 0.)
  {
    m_logMsg = "Invalid rural buffer distance.";

    return false;
  }

  if(inputParamsPtr->m_ruralSmallBufferDistance < 0.)
  {
    m_logMsg = "Invalid rural buffer distance.";

    return false;
  }

  m_inputParameters = *inputParamsPtr;

  m_isInitialized = true;

  return true;
}

double te::app::Lakes::calculateArea(te::gm::Geometry* geom)
{
  double area = 0.;

  if(geom->getGeomTypeId() == te::gm::MultiPolygonType)
  {
    te::gm::MultiPolygon* mPolygon = dynamic_cast<te::gm::MultiPolygon*>(geom);

    area = mPolygon->getArea();
  }
  else if(geom->getGeomTypeId() == te::gm::PolygonType)
  {
    te::gm::Polygon* polygon = dynamic_cast<te::gm::Polygon*>(geom);

    area = polygon->getArea();
  }

  return area;
}
