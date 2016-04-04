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
#include <terralib/srs/Datum.h>
#include <terralib/srs/Ellipsoid.h>
#include <terralib/srs/GeographicCoordinateSystem.h>
#include <terralib/srs/SpatialReferenceSystemManager.h>


#define M_PI       3.14159265358979323846

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
  std::auto_ptr<te::rst::Raster> slopeRaster = calculateSlope(m_inputParameters.m_inRasterPtr);

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

std::auto_ptr<te::rst::Raster> te::app::HighSlope::calculateSlope(te::rst::Raster const* inputRst)
{
  //create slope raster
  std::vector< te::rst::BandProperty* > bandsProperties;
  te::rst::BandProperty* bandProp = new te::rst::BandProperty(0, te::dt::DOUBLE_TYPE);
  bandProp->m_noDataValue = -9999;
  bandsProperties.push_back(bandProp);

  te::rst::Grid* grid = new te::rst::Grid(*(m_inputParameters.m_inRasterPtr->getGrid()));

  te::rst::Raster* outRaster = te::rst::RasterFactory::make(m_outputParametersPtr->m_createdOutRasterDSType, grid, bandsProperties, m_outputParametersPtr->m_createdOutRasterInfo);

  //initialize
  for (unsigned int i = 0; i < outRaster->getNumberOfRows(); ++i)
  {
    for (unsigned int j = 0; j < outRaster->getNumberOfColumns(); ++j)
    {
      outRaster->setValue(j, i, -9999, 0);
    }
  }

  unsigned int nlines = outRaster->getNumberOfRows();
  unsigned int ncolumns = outRaster->getNumberOfColumns();
  double rasterDummy = inputRst->getBand(0)->getProperty()->m_noDataValue;

  double resx = outRaster->getGrid()->getResolutionX(); //Gd
  double resy = outRaster->getGrid()->getResolutionY(); //Ge
  double A, F, E2, DY, DX;

  std::string authName = "EPSG"; // Now: So far it is the only one supported by TerraLib 5. Future: Review this line!
  bool isGeographic = te::srs::SpatialReferenceSystemManager::getInstance().isGeographic(outRaster->getGrid()->getSRID(), authName);
  te::srs::SpatialReferenceSystem* srs = te::srs::SpatialReferenceSystemManager::getInstance().getSpatialReferenceSystem(outRaster->getGrid()->getSRID()).release();

  te::srs::GeographicCoordinateSystem* gcs = dynamic_cast<te::srs::GeographicCoordinateSystem*>(srs);
  
  if (isGeographic && gcs)
  {
    A = gcs->getDatum()->getEllipsoid()->getRadium() / 1000.0;
    F = gcs->getDatum()->getEllipsoid()->getInverseFlattening();
    E2 = 2 * F - F * F; //!QUADRADO DA EXCENTRICIDADE
    DY = (A * sin(resy * M_PI / 180)) * 1000;

    resx = resx * 111133;
    resy = resy * 111133;
  }

  // Execute the operation
  te::common::TaskProgress task("Generating Slope");
  task.setTotalSteps(inputRst->getNumberOfRows());

  double rx, ry;

  //calculate slope
  for (unsigned int line = 1; line<nlines - 1; line++)
  {
    if (isGeographic)
    {
      //// A área varia somente na lina / y / latitude
      //te::gm::Coord2D coord = inputRst->getGrid()->gridToGeo(0, line);

      //double YLAT = coord.getY();

      //double RN = A / sqrt((1 - E2 * sin(YLAT * M_PI / 180)*sin(YLAT * M_PI / 180))); //!RAIO DE CURVATURA DA TERRA NA LATITUDE

      //// =B8*SEN(3/3600*PI()/180)*COS((A8+1.5/3600)*PI()/180)
      //DX = (RN * sin(resx * M_PI / 180) * cos(YLAT * M_PI / 180)) * 1000; // em metros

      rx = resx;
      ry = resy;
    }
    else
    {
      rx = resx;
      ry = resy;
    }

    for (unsigned int column = 1; column<ncolumns - 1; column++)
    {
      double value = rasterDummy;

      inputRst->getValue(column, line, value, 0);

      if (value != rasterDummy)
      {
        double z1 = 0, z2 = 0, z3 = 0, z4 = 0, z6 = 0, z7 = 0, z8 = 0, z9 = 0;

        inputRst->getValue(column - 1, line - 1, value, 0);
        if (value != rasterDummy)
        {
          z1 = value;
        }

        inputRst->getValue(column, line - 1, value, 0);
        if (value != rasterDummy)
        {
          z2 = value;
        }

        inputRst->getValue(column + 1, line - 1, value, 0);
        if (value != rasterDummy)
        {
          z3 = value;
        }

        inputRst->getValue(column - 1, line, value, 0);
        if (value != rasterDummy)
        {
          z4 = value;
        }

        inputRst->getValue(column + 1, line, value, 0);
        if (value != rasterDummy)
        {
          z6 = value;
        }

        inputRst->getValue(column - 1, line + 1, value, 0);
        if (value != rasterDummy)
        {
          z7 = value;
        }

        inputRst->getValue(column, line + 1, value, 0);
        if (value != rasterDummy)
        {
          z8 = value;
        }

        inputRst->getValue(column + 1, line + 1, value, 0);
        if (value != rasterDummy)
        {
          z9 = value;
        }

        double d = (z3 + z6 + z9 - z1 - z4 - z7) / (6 * rx);

        double e = (z1 + z2 + z3 - z7 - z8 - z9) / (6 * ry);

        outRaster->setValue(column, line, (float)sqrt(d*d + e*e), 0);
      }
    }

    task.pulse();
  }

  //set output raster into auto_ptr
  std::auto_ptr<te::rst::Raster> outputRst(outRaster);

  return outputRst;
}
