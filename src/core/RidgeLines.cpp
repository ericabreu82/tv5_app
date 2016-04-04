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
  \file terralib/app/RidgeLines.cpp
 
  \brief Generates the APP information from a ridge lines information.
*/

#include "RidgeLines.h"
#include "Utils.h"

// TerraLib
#include <terralib/common/UnitOfMeasure.h>
#include <terralib/common/UnitsOfMeasureManager.h>
#include <terralib/common/progress/TaskProgress.h>
#include <terralib/geometry/Geometry.h>
#include <terralib/geometry/LineString.h>
#include <terralib/geometry/MultiLineString.h>
#include <terralib/geometry/MultiPolygon.h>
#include <terralib/raster/BandProperty.h>
#include <terralib/raster/Grid.h>
#include <terralib/raster/Raster.h>
#include <terralib/raster/RasterFactory.h>
#include <terralib/raster/Utils.h>
#include <terralib/srs/SpatialReferenceSystemManager.h>


//STL Includes
#include <queue>

#define DEFAULT_SRID_METERS 29193

bool cmp(const te::app::RidgeLines::Point_Peak &v1, const te::app::RidgeLines::Point_Peak &v2)
{   
    return(v1.d < v2.d);
}

te::app::RidgeLines::InputParameters::InputParameters()
{
  reset();
}

te::app::RidgeLines::InputParameters::~InputParameters()
{
  reset();
}

void te::app::RidgeLines::InputParameters::reset() throw(te::common::Exception)
{
  m_demRasterPtr = 0;
  m_slopeRasterPtr = 0;
  m_slopeMaskRasterPtr = 0;
  m_plainMaskRasterPtr = 0;
  m_contourLines = 0;
  m_drainage = 0;
  m_watershedsBuffer = 0;
  m_watersheds = 0;
  m_height = 50.;
  m_slope = 17.;
}

te::common::AbstractParameters* te::app::RidgeLines::InputParameters::clone() const
{
  return new InputParameters(*this);
}

te::app::RidgeLines::OutputParameters::OutputParameters()
{
  reset();
}

te::app::RidgeLines::OutputParameters::OutputParameters(const OutputParameters& other)
{
  reset();

  operator=(other);
}

te::app::RidgeLines::OutputParameters::~OutputParameters()
{
  reset();
}

void te::app::RidgeLines::OutputParameters::reset() throw(te::common::Exception)
{
  m_createdOutDSName.clear();
  m_createdOutDSType.clear();
  m_createdOutInfo.clear();
}

const te::app::RidgeLines::OutputParameters& te::app::RidgeLines::OutputParameters::operator=(const RidgeLines::OutputParameters& params)
{
  reset();

  m_createdOutDSName = params.m_createdOutDSName;
  m_createdOutDSType = params.m_createdOutDSType;
  m_createdOutInfo = params.m_createdOutInfo;

  return *this;
}

te::common::AbstractParameters* te::app::RidgeLines::OutputParameters::clone() const
{
  return new OutputParameters(*this);
}

te::app::RidgeLines::RidgeLines()
{
  reset();
}

te::app::RidgeLines::~RidgeLines()
{
}

bool te::app::RidgeLines::execute(AlgorithmOutputParameters& outputParams) throw(te::common::Exception)
{
  if(!m_isInitialized)
  {
    m_logMsg = "Algoritm not initialized.";

    return false;
  }

  m_outputParametersPtr = dynamic_cast<RidgeLines::OutputParameters*>(&outputParams);

  if(!m_outputParametersPtr)
  {
    m_logMsg = "Invalid output parameters.";

    return false;
  }

  ////first create slope mask
  //std::auto_ptr<te::rst::Raster> slopeMask = createSlopeMask(0.1);

  //if(!slopeMask.get())
  //{
  //  return false;
  //}

  ////create plain raster mask
  //std::auto_ptr<te::rst::Raster> plainRasterMask = createPlainsRasterMask(slopeMask.get(), 10.);

  //if(!plainRasterMask.get())
  //{
  //  return false;
  //}

  //create peaks
  std::vector<te::app::RidgeLines::Point_T> peaks = generateHillTopVector(m_inputParameters.m_plainMaskRasterPtr);

  if(peaks.empty())
  {
    return false;
  }

  //create peaks lines
  std::vector< std::vector<Point_T> > peaksLines;
  std::vector<te::gm::LineString*> peaksLinesVec;

  generatePeaksLine(peaks, peaksLines, peaksLinesVec);

  if(peaksLines.empty())
  {
    return false;
  }

  te::app::ExportVectorLine(peaksLinesVec, m_outputParametersPtr->m_createdOutDSName, m_outputParametersPtr->m_createdOutDSType, m_outputParametersPtr->m_createdOutInfo);

  //generate peaks lines app
  std::auto_ptr<te::rst::Raster> rasterPeaksApp = generatePeaksExtendedAPP(peaksLines, m_inputParameters.m_plainMaskRasterPtr);

  if(!rasterPeaksApp.get())
  {
    return false;
  }

  std::string uri = TERRALIB_DATA_DIR"/app/rasterPeaksApp.tiff";

  te::rst::CreateCopy(*rasterPeaksApp.get(), uri);

  //create vector representation
  //te::app::Raster2Vector(rasterPeaksApp.get(), 0,  m_outputParametersPtr->m_createdOutDSName, m_outputParametersPtr->m_createdOutDSType,  m_outputParametersPtr->m_createdOutInfo);

  return true;
}

void te::app::RidgeLines::reset() throw(te::common::Exception)
{
  m_inputParameters.reset();
  m_outputParametersPtr = 0;
}

bool te::app::RidgeLines::initialize(const AlgorithmInputParameters& inputParams) throw(te::common::Exception)
{
  reset();

  RidgeLines::InputParameters const* inputParamsPtr = dynamic_cast<RidgeLines::InputParameters const*>(&inputParams);

  if(!inputParamsPtr)
  {
    m_logMsg = "Invalid Parameters.";

    return false;
  }

  if (!inputParamsPtr->m_demRasterPtr || inputParamsPtr->m_demRasterPtr->getAccessPolicy() == te::common::NoAccess)
  {
    m_logMsg = "Invalid Raster.";

    return false;
  }

  if (!inputParamsPtr->m_slopeRasterPtr || inputParamsPtr->m_slopeRasterPtr->getAccessPolicy() == te::common::NoAccess)
  {
    m_logMsg = "Invalid Raster.";

    return false;
  }

  if (!inputParamsPtr->m_contourLines || inputParamsPtr->m_contourLines->getNumGeometries() == 0)
  {
    m_logMsg = "Invalid Contour Lines.";

    return false;
  }

  m_inputParameters = *inputParamsPtr;

  m_inputParameters.m_slope = m_inputParameters.m_slope * 3.14159265 / 180.;
  m_inputParameters.m_slope = tan(m_inputParameters.m_slope)*100.;

  m_isInitialized = true;

  return true;
}

std::auto_ptr<te::rst::Raster> te::app::RidgeLines::createSlopeMask(double threshold_slope)
{
  double slope_undef = -9999.;

  //create escarp raster
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

  for(unsigned int i = 0; i < outRaster->getNumberOfRows(); ++i)
  {
    for(unsigned int j = 0; j < outRaster->getNumberOfColumns(); ++j)
    {
      double slope_value;

      m_inputParameters.m_slopeRasterPtr->getValue(j, i, slope_value);

      if(slope_value != slope_undef && slope_value < threshold_slope)
      {
        outRaster->setValue(j, i, 1., 0);
      }
      else
      {
        outRaster->setValue(j, i, 0., 0);
      }
    }
  }

  //return raster
  std::auto_ptr<te::rst::Raster> slopeMask(outRaster);

  return slopeMask;
}

std::auto_ptr<te::rst::Raster> te::app::RidgeLines::createPlainsRasterMask(te::rst::Raster* slopeMask, double area_threshold)
{
  int dx[8] = {1, 1, 1, 0, -1, -1, -1, 0};
  int dy[8] = {1, 0, -1, -1, -1, 0, 1, 1};

  //create terrace raster
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

  //fill raster with dummy values
  for(unsigned int i = 0; i < outRaster->getNumberOfRows(); ++i)
  {
    for(unsigned int j = 0; j < outRaster->getNumberOfColumns(); ++j)
    {
      outRaster->setValue(j, i, 0., 0);
    }
  }

  //create slope mask boolean matrix
  std::vector< std::vector<bool> > slope_mask_visited;

  for(unsigned int i = 0;i < m_inputParameters.m_demRasterPtr->getNumberOfRows();i++)
  {
    std::vector<bool> row;
    row.resize(m_inputParameters.m_demRasterPtr->getNumberOfColumns(), false);

    slope_mask_visited.push_back(row);
  }

  double slope_undef = 0.;

  //calculate cell area
  //te::common::UnitOfMeasurePtr unit = te::srs::SpatialReferenceSystemManager::getInstance().getUnit(m_inputParameters.m_demRasterPtr->getSRID());
  //std::string unitSymbol = unit->getSymbol();
  //std::string unitDefault = "metre";

  //double conversionFactor = te::common::UnitsOfMeasureManager::getInstance().getConversion(unitSymbol, unitDefault);

  double conversionFactor = 111000.;

  double resX_meter = m_inputParameters.m_demRasterPtr->getGrid()->getResolutionX() * conversionFactor;
  double resY_meter = m_inputParameters.m_demRasterPtr->getGrid()->getResolutionY() * conversionFactor;

  double cell_area = resX_meter * resY_meter;

  //fill terrace raster
  for(unsigned int i = 0; i < m_inputParameters.m_demRasterPtr->getNumberOfRows(); ++i)
  {
    for(unsigned int j = 0; j < m_inputParameters.m_demRasterPtr->getNumberOfColumns(); ++j)
    {
      if(slope_mask_visited[i][j] == true)
        continue;

      double slope_value;

      // read slope mask
      slopeMask->getValue(j, i, slope_value);

      // if slope mask is valid
      if(slope_value != slope_undef)
      {
        std::queue<te::gm::Coord2D> grid_points_queue;
        std::vector<te::gm::Coord2D> plain_points;

        te::gm::Coord2D grid_point(j, i);

        slope_mask_visited[i][j] = true;

        plain_points.push_back(grid_point);
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

            if(slope_mask_visited[y+dy[k]][x+dx[k]] == true)
              continue;

            // read slope mask
            slopeMask->getValue(x+dx[k], y+dy[k], slope_value);

            if(slope_value != slope_undef)
            {
              grid_point.x = x+dx[k];
              grid_point.y = y+dy[k];

              slope_mask_visited[y+dy[k]][x+dx[k]] = true;

              plain_points.push_back(grid_point);
              grid_points_queue.push(grid_point);
            }
          }
        }

        // plain area in m^2
        double plain_area = cell_area * plain_points.size();
        
        // convert from m^2 to ha
        plain_area /= 10000.;

        if(plain_area > area_threshold)
        {
          for(std::size_t k = 0; k < plain_points.size(); k++)
          {
            int x = (int) plain_points[k].getX();
            int y = (int) plain_points[k].getY();

            outRaster->setValue(x, y, 1., 0);
          }
        }
      }
    }
  }

  //return raster
  std::auto_ptr<te::rst::Raster> plainMaskRaster(outRaster);

  return plainMaskRaster;
}

std::vector<te::app::RidgeLines::Point_T> te::app::RidgeLines::generateHillTopVector(te::rst::Raster* plainRasterMask)
{
  std::vector<te::app::RidgeLines::Point_T> peaks;

  constructFlagMatrix(plainRasterMask);

  std::vector<te::app::RidgeLines::Point_T> selectedPeaks = selectPeaks();

  for (std::size_t i = 0; i < selectedPeaks.size(); i++)
  {
    int x = (int)selectedPeaks[i].x;
    int y = (int)selectedPeaks[i].y;

    if (m_SLgFlag[x][y] != 1 && m_SLgFlag[x][y] != 0)
      continue;

    double maxvalue = selectedPeaks[i].z;

    te::gm::Coord2D pm(x, y);

    resetFlagMatrix(m_inputParameters.m_demRasterPtr->getNumberOfRows(), m_inputParameters.m_demRasterPtr->getNumberOfColumns()); // Retorna matriz de flag na situacao original

    te::app::RidgeLines::Point_T pointMin;

    if (!topToDown(x, y, 2, pointMin))
      continue;

    // guarda no vetor LINHA e COLUNA do ponto maximo e a cota do ponto minimo
    Point_T p3d;
    p3d.x = pm.getX();
    p3d.y = pm.getY();

    double cotacorte = 2 * ((maxvalue - pointMin.z) / 3.) + pointMin.z;

    p3d.z = cotacorte;

    peaks.push_back(p3d); // TopoBase armazena coluna, linha e cota de corte
  }

  selectedPeaks.clear();

  std::auto_ptr<te::rst::Raster> outRaster;

  {
    //create peaks raster
    std::vector<te::rst::BandProperty*> bandsProperties;
    te::rst::BandProperty* bandProp = new te::rst::BandProperty(0, te::dt::UCHAR_TYPE);
    bandProp->m_noDataValue = 255.;
    bandsProperties.push_back(bandProp);

    te::rst::Grid* grid = new te::rst::Grid(*(m_inputParameters.m_demRasterPtr->getGrid()));

    std::map<std::string, std::string> rinfo;
    rinfo["URI"] = "d:/peaks.tif";
    
    outRaster.reset(te::rst::RasterFactory::make("GDAL", grid, bandsProperties, rinfo));

    if (!outRaster.get())
    {
      throw(te::common::Exception("Peaks raster APP creation error."));
    }

    //fill raster with dummy values
    for (unsigned int i = 0; i < outRaster->getNumberOfRows(); ++i)
    {
      for (unsigned int j = 0; j < outRaster->getNumberOfColumns(); ++j)
      {
        outRaster->setValue(j, i, 255., 0);
      }
    }
  }

  for (std::size_t t = 0; t < peaks.size(); t++)
  {
    double x = peaks[t].x;
    double y = peaks[t].y;
    double z = peaks[t].z;

    outRaster->setValue(x, y, 1., 0);
  }

  std::auto_ptr<te::rst::Raster> outRaster2;

  {
    //create peaks raster
    std::vector<te::rst::BandProperty*> bandsProperties;
    te::rst::BandProperty* bandProp = new te::rst::BandProperty(0, te::dt::UCHAR_TYPE);
    bandProp->m_noDataValue = 255.;
    bandsProperties.push_back(bandProp);

    te::rst::Grid* grid = new te::rst::Grid(*(m_inputParameters.m_demRasterPtr->getGrid()));

    std::map<std::string, std::string> rinfo;
    rinfo["URI"] = "d:/peaksApp.tif";

    outRaster2.reset(te::rst::RasterFactory::make("GDAL", grid, bandsProperties, rinfo));

    if (!outRaster2.get())
    {
      throw(te::common::Exception("Peaks raster APP creation error."));
    }

    //fill raster with dummy values
    for (unsigned int i = 0; i < outRaster2->getNumberOfRows(); ++i)
    {
      for (unsigned int j = 0; j < outRaster2->getNumberOfColumns(); ++j)
      {
        outRaster2->setValue(j, i, 255., 0);
      }
    }
  }


  //fill raster out
  for (std::size_t t = 0; t < peaks.size(); t++)
  {
    resetFlagMatrix(m_inputParameters.m_demRasterPtr->getNumberOfRows(), m_inputParameters.m_demRasterPtr->getNumberOfColumns()); // Retorna matriz de flag na situacao original

    int x = (int)peaks[t].x;
    int y = (int)peaks[t].y;

    m_SLgFlag[x][y] = 1;

    double z = peaks[t].z;

    double val;
    outRaster2->getValue((unsigned int)x, (unsigned int)y, val);

    if (val != 255.)
      continue;

    te::app::RidgeLines::Point_T pointMin;

    if (!topToDown(x, y, 5, pointMin))
      continue;

    hillTopDraw(x, y, peaks[t].z, outRaster2.get(), 125.);
  }

  return peaks;
}

void te::app::RidgeLines::generatePeaksLine(std::vector<te::app::RidgeLines::Point_T>& peaksVector, std::vector< std::vector<Point_T> >& peaksLines, std::vector<te::gm::LineString*>& peaksLinesVec)
{
  for(std::size_t t = 0; t < peaksVector.size(); ++t)
  {
    peaksVector[t].visited = false;
  }

  int srid = m_inputParameters.m_demRasterPtr->getSRID();

  //store polygons of closed contour lines
  std::vector<te::gm::Polygon*> polyVec;

  for(std::size_t t = 0; t < m_inputParameters.m_contourLines->getNumGeometries(); ++t)
  {
    te::gm::LineString* ls = dynamic_cast<te::gm::LineString*>(m_inputParameters.m_contourLines->getGeometryN(t));

    if(ls && ls->isClosed())
    {
      te::gm::Polygon* p = new te::gm::Polygon(0, te::gm::PolygonType, ls->getSRID());

      te::gm::LinearRing* ring = new te::gm::LinearRing(ls->getNPoints(), te::gm::LineStringType);

      for(std::size_t i = 0; i < ls->getNPoints(); ++i)
        ring->setPoint(i, ls->getX(i), ls->getY(i));

      p->add(ring);

      polyVec.push_back(p);
    }
  }

  //init vectors
  std::vector< std::vector<te::gm::Coord2D> > peaks_line2Dpoints_vector;
  std::vector<te::app::RidgeLines::Point_Peak> peaks_line_vector;

  peaksLines.clear();
  peaksLines.resize(peaksVector.size());
  peaks_line2Dpoints_vector.resize(peaksVector.size());

  int n2DPeaksLines = 0;

  //select peaks inside closed contour lines connected by drainage divides
  for (std::size_t t = 0; t < peaksVector.size(); ++t)
  {
    if (peaksVector[t].visited == true)
      continue;

    peaksVector[t].visited = true;

    double x = peaksVector[t].x;
    double y = peaksVector[t].y;
    double z = peaksVector[t].z;

    te::gm::Coord2D gp1(x, y);
    te::gm::Coord2D pp1 = m_inputParameters.m_demRasterPtr->getGrid()->gridToGeo(gp1.getX(), gp1.getY());

    std::auto_ptr<te::gm::Point> point(new te::gm::Point(pp1.getX(), pp1.getY(), srid));

    // check if base peak is inside a closed contour line
    bool peak1_inside = false;

    for (std::size_t m = 0; m < polyVec.size(); m++)
    {
      try
      {
        te::gm::Geometry* g = polyVec[m];

        if (g->isValid() && point->isValid() && point->within(g))
        {
          peak1_inside = true;
          break;
        }
      }
      catch (...)
      {
        continue;
      }
    }

    if (peak1_inside == false)
      continue;

    // check if first peak point is "near" a drainage divide considering a distance buffer
    if (m_inputParameters.m_watershedsBuffer)
    {
      bool peak_in_drainage_divide = false;

      for (std::size_t u = 0; u < m_inputParameters.m_watershedsBuffer->getNumGeometries(); ++u)
      {
        te::gm::Geometry* g = m_inputParameters.m_watershedsBuffer->getGeometryN(u);

        if (point->within(g))
        {
          peak_in_drainage_divide = true;
          break;
        }
      }

      if (peak_in_drainage_divide == false)
        continue;
    }

    // insert peak into peaks vector
    // peak in grid coordinates
    te::app::RidgeLines::Point_Peak peak;

    peak.idx = t;
    peak.x = peaksVector[t].x;
    peak.y = peaksVector[t].y;
    peak.z = peaksVector[t].z;
    peak.visited = false;

    peaks_line_vector.push_back(peak);

    // current peak index
    std::size_t current_peak_idx = t;

    // search for peaks of drainage divides
    bool found_near_peak = true;

    while(found_near_peak == true)
    {
      found_near_peak = false;

      // process current line peak
      x = peaksVector[current_peak_idx].x;
      y = peaksVector[current_peak_idx].y;
      z = peaksVector[current_peak_idx].z;

      gp1.x = x;
      gp1.y = y;
      pp1 = m_inputParameters.m_demRasterPtr->getGrid()->gridToGeo(gp1.getX(), gp1.getY());

      std::auto_ptr<te::gm::Point> point_pp1(new te::gm::Point(pp1.getX(), pp1.getY(), srid));

      double current_peak_z;
      m_inputParameters.m_demRasterPtr->getValue(gp1.getX(), gp1.getY(), current_peak_z); //???????????????????

      // sort peaks by distance to current peak
      std::vector<te::app::RidgeLines::Point_Peak> peaks_by_distance;

      for(std::size_t j = 0; j < peaksVector.size(); ++j)
      {
        te::gm::Coord2D pp2 = m_inputParameters.m_demRasterPtr->getGrid()->gridToGeo(peaksVector[j].x, peaksVector[j].y);

        std::auto_ptr<te::gm::Point> point_pp2(new te::gm::Point(pp2.getX(), pp2.getY(), srid));

        peak.idx = j;
        peak.d = point_pp2->distance(point_pp1.get());

        peaks_by_distance.push_back(peak);
      }
      sort(peaks_by_distance.begin(), peaks_by_distance.end(), cmp);

      // check if another peak will be included in the line peak
      for(std::size_t j = 1; j < peaksVector.size(); ++j)
      {
        // index of sorted peaks
        int near_peak_idx = peaks_by_distance[j].idx;

        if(peaksVector[near_peak_idx].visited == true)
          continue;

        te::gm::Coord2D pp2 = m_inputParameters.m_demRasterPtr->getGrid()->gridToGeo(peaksVector[near_peak_idx].x, peaksVector[near_peak_idx].y);

        // check if near peak is inside a closed contour line
        bool peak2_inside = false;

        std::auto_ptr<te::gm::Point> point_pp2(new te::gm::Point(pp2.getX(), pp2.getY(), srid));

        for(std::size_t m = 0 ;m < polyVec.size(); m++)
        {
          try
          {
            te::gm::Polygon* g = polyVec[m];

            if (g->isValid() && point_pp2->isValid() && point_pp2->within(g) == true)
            {
              peak2_inside = true;
              break;
            }
          }
          catch (...)
          {
            continue;
          }
        }

        if(peak2_inside == false)
          continue;

        std::auto_ptr<te::gm::LineString> peaks_line(new te::gm::LineString(2, te::gm::LineStringType, srid));

        peaks_line->setPoint(0, pp1.getX(), pp1.getY());
        peaks_line->setPoint(1, pp2.getX(), pp2.getY());

        // check intersection of peak line segment with drainage segment
        if (m_inputParameters.m_drainage)
        {
          bool drainage_intersection = false;

          for (std::size_t t = 0; t < m_inputParameters.m_drainage->getNumGeometries(); ++t)
          {
            te::gm::Geometry* drainage_line = m_inputParameters.m_drainage->getGeometryN(t);

            if (drainage_line->intersects(peaks_line.get()))
            {
              drainage_intersection = true;
              break;
            }
          }

          if (drainage_intersection == true)
            continue;
        }

        // check if near peak point is "near" a drainage divide considering a distance buffer
        if (m_inputParameters.m_watershedsBuffer)
        {
          bool peak_in_drainage_divide = false;

          for (std::size_t t = 0; t < m_inputParameters.m_watershedsBuffer->getNumGeometries(); ++t)
          {
            te::gm::Geometry* g = m_inputParameters.m_watershedsBuffer->getGeometryN(t);

            if (point_pp2->within(g))
            {
              peak_in_drainage_divide = true;
              break;
            }
          }

          if (peak_in_drainage_divide == false)
            continue;
        }

        double neighbor_peak_z;
        m_inputParameters.m_demRasterPtr->getValue(peaksVector[near_peak_idx].x, peaksVector[near_peak_idx].y, neighbor_peak_z); //???????????????

        // grid points of line between current peak and near peak
        std::vector<te::gm::Coord2D> p2D_vector = bresenhamPoints(gp1, te::gm::Coord2D(peaksVector[near_peak_idx].x, peaksVector[near_peak_idx].y));

        // find watershed polygons that intersect the segment between current and near peaks
        std::vector<te::gm::Polygon*> watershed_polygons;

        for(std::size_t t = 0; t < m_inputParameters.m_watersheds->getNumGeometries(); ++t)
        {
          te::gm::Polygon* poly = dynamic_cast<te::gm::Polygon*>(m_inputParameters.m_watersheds->getGeometryN(t));

          if(peaks_line->intersects(poly))
          {
            watershed_polygons.push_back(poly);
          }
        }

        // copy linear ring points of watershed polygons
        std::vector<te::app::RidgeLines::Point_T> watersheds_LR;

        for(std::size_t t = 0; t < watershed_polygons.size(); ++t)
        {
          te::gm::Polygon* poly = watershed_polygons[t];

          te::gm::LineString* lsCurve = dynamic_cast<te::gm::LineString*>(poly->getRingN(0));

          for(std::size_t q = 0; q < lsCurve->getNPoints(); ++q)
          {
            te::gm::Point* point = lsCurve->getPointN(q);

            te::app::RidgeLines::Point_T watershed_LR_point;

            watershed_LR_point.x = point->getX();
            watershed_LR_point.y = point->getY();
            watershed_LR_point.visited = false;

            bool existingPoint = false;

            for (std::size_t u = 0; u < watersheds_LR.size(); ++u)
            {
              if (watershed_LR_point.x == watersheds_LR[u].x && watershed_LR_point.y == watersheds_LR[u].y)
              {
                existingPoint = true;
                break;
              }
            }

            if (existingPoint == false)
              watersheds_LR.push_back(watershed_LR_point);

            delete point;
          }
        }

        // current and near peaks in projected coordinates
        pp1 = m_inputParameters.m_demRasterPtr->getGrid()->gridToGeo(peaksVector[current_peak_idx].x, peaksVector[current_peak_idx].y);

        point_pp1->setX(pp1.getX());
        point_pp1->setY(pp1.getY());

        pp2 = m_inputParameters.m_demRasterPtr->getGrid()->gridToGeo(peaksVector[near_peak_idx].x, peaksVector[near_peak_idx].y);

        point_pp2->setX(pp2.getX());
        point_pp2->setY(pp2.getY());

        // horizontal plane distance
        double dist = point_pp1->distance(point_pp2.get());

        // if peak is valid, then add as part of the line peak
        double distance_threshold = 500.;
        
        if(dist <= distance_threshold)
        {
          peaksVector[near_peak_idx].visited = true;

          peak.idx = near_peak_idx;
          peak.x = peaksVector[near_peak_idx].x;
          peak.y = peaksVector[near_peak_idx].y;
          peak.z = peaksVector[near_peak_idx].z;
          peak.visited = false;
          peaks_line_vector.push_back(peak);

          // for each bresenham point find the closest unvisited linear ring point
          int inc = 2;

          for(std::size_t m = 0; m < p2D_vector.size(); m += inc)
          {
            pp1 = m_inputParameters.m_demRasterPtr->getGrid()->gridToGeo(p2D_vector[m].x, p2D_vector[m].y);

            point_pp1->setX(pp1.getX());
            point_pp1->setY(pp1.getY());

            double min_dist = 9999999.;
            std::size_t min_idx;

            for(std::size_t n = 0; n < watersheds_LR.size(); n++)
            {
              pp2.x = watersheds_LR[n].x;
              pp2.y = watersheds_LR[n].y;

              point_pp2->setX(pp2.getX());
              point_pp2->setY(pp2.getY());

              dist = point_pp1->distance(point_pp2.get());

              if(dist < min_dist)
              {
                min_idx = n;
                min_dist = dist;
              }
            }

            if(watersheds_LR[min_idx].visited == false)
            {
              te::gm::Coord2D pp;

              pp.x = watersheds_LR[min_idx].x;
              pp.y = watersheds_LR[min_idx].y;

              // add point to coord2D
              peaks_line2Dpoints_vector[n2DPeaksLines].push_back(pp);

              watersheds_LR[min_idx].visited = true;
            }

            if(m + inc >= p2D_vector.size())
              inc = 1;
          }

          te::gm::LineString* ls = new te::gm::LineString(peaks_line2Dpoints_vector[n2DPeaksLines].size(), te::gm::LineStringType, srid);

          for (std::size_t t = 0; t < peaks_line2Dpoints_vector[n2DPeaksLines].size(); ++t)
          {
            ls->setPoint(t, peaks_line2Dpoints_vector[n2DPeaksLines][t].getX(), peaks_line2Dpoints_vector[n2DPeaksLines][t].getY());
          }

          peaksLinesVec.push_back(ls);

          // update current peak
          current_peak_idx = near_peak_idx;

          found_near_peak = true;

          break;
        }
      }
      //break;
    }

    // if new peaks line was traced
    if(peaks_line2Dpoints_vector[n2DPeaksLines].empty() == false)
      n2DPeaksLines++;
  }

  double res = std::max(m_inputParameters.m_demRasterPtr->getResolutionX(), m_inputParameters.m_demRasterPtr->getResolutionY());

  // divide peaks line in segments of 1000 m then associates peaks with each segment
  int n3DSegments = 0;
  double raster_value = 1.;

  std::auto_ptr<te::rst::Raster> outRaster;

  {
    //create peaks raster
    std::vector<te::rst::BandProperty*> bandsProperties;
    te::rst::BandProperty* bandProp = new te::rst::BandProperty(0, te::dt::UCHAR_TYPE);
    bandProp->m_noDataValue = 0.;
    bandsProperties.push_back(bandProp);

    te::rst::Grid* grid = new te::rst::Grid(*(m_inputParameters.m_demRasterPtr->getGrid()));

    std::map<std::string, std::string> rinfo;
    rinfo["URI"] = "d:/peaksLines.tif";

    outRaster.reset(te::rst::RasterFactory::make("GDAL", grid, bandsProperties, rinfo));

    if (!outRaster.get())
    {
      throw(te::common::Exception("Peaks raster APP creation error."));
    }

    //fill raster with dummy values
    for (unsigned int i = 0; i < outRaster->getNumberOfRows(); ++i)
    {
      for (unsigned int j = 0; j < outRaster->getNumberOfColumns(); ++j)
      {
        outRaster->setValue(j, i, 0., 0);
      }
    }
  }

  for(int i = 0;i < n2DPeaksLines; i++)
  {
    double dist_sum = 0.;

    for(std::size_t j = 0; j < peaks_line2Dpoints_vector[i].size()-1; j++)
    {
      double pp1_z;
      te::gm::Coord2D pp1 = peaks_line2Dpoints_vector[i][j];
      te::gm::Coord2D gp1 = m_inputParameters.m_demRasterPtr->getGrid()->geoToGrid(pp1.getX(), pp1.getY());
      m_inputParameters.m_demRasterPtr->getValue((unsigned int)gp1.getX(), (unsigned int)gp1.getY(), pp1_z);

      std::auto_ptr<te::gm::Point> point_pp1(new te::gm::Point(pp1.getX(), pp1.getY(), srid));

      double pp2_z;
      te::gm::Coord2D pp2 = peaks_line2Dpoints_vector[i][j+1];
      te::gm::Coord2D gp2 = m_inputParameters.m_demRasterPtr->getGrid()->geoToGrid(pp2.getX(), pp2.getY());
      m_inputParameters.m_demRasterPtr->getValue((unsigned int)gp2.getX(), (unsigned int)gp2.getY(), pp2_z);

      std::auto_ptr<te::gm::Point> point_pp2(new te::gm::Point(pp2.getX(), pp2.getY(), srid));

      // horizontal plane distance
      double dist = point_pp1->distance(point_pp2.get());

      // actual distance considering elevation
      double z_diff = abs(pp1_z-pp2_z);

      dist = sqrt(dist*dist + z_diff*z_diff);

      if(dist_sum + dist <= 1000.)
      {
        // check if peaks line point is a peak
        for(std::size_t k = 0; k < peaks_line_vector.size(); k++)
        {
          if(peaks_line_vector[k].visited == true)
            continue;

          te::gm::Coord2D pp = m_inputParameters.m_demRasterPtr->getGrid()->gridToGeo(peaks_line_vector[k].x, peaks_line_vector[k].y);

          std::auto_ptr<te::gm::Point> point_pp(new te::gm::Point(pp.getX(), pp.getY(), srid));

          double distPoints = point_pp->distance(point_pp1.get());

          //double distance_threshold = 1000.;
          //if (distPoints <= distance_threshold)
          if (distPoints <= res)
          {
            peaks_line_vector[k].visited = true;

            te::app::RidgeLines::Point_T peak_3Dpoint;

            peak_3Dpoint.x = peaks_line_vector[k].x;
            peak_3Dpoint.y = peaks_line_vector[k].y;
            peak_3Dpoint.z = peaks_line_vector[k].z;

            // save segment peak
            peaksLines[n3DSegments].push_back(peak_3Dpoint);

            outRaster->setValue(peaks_line_vector[k].x, peaks_line_vector[k].y, raster_value, 0);
          }
        }

        // update distance sum
        dist_sum += dist;
      }
      else
      {
        // reinitialize distance sum
        dist_sum = 0.;
        // increment value of raster peaks
        raster_value += 1;
        // increment number of segments
        n3DSegments++;
      }
    }

    // add last peak if not added
    std::size_t j = peaks_line2Dpoints_vector[i].size() - 1;
    te::gm::Coord2D pp1 = peaks_line2Dpoints_vector[i][j];

    std::auto_ptr<te::gm::Point> point_pp1(new te::gm::Point(pp1.getX(), pp1.getY(), srid));

    // check if peaks line point is a peak
    for (std::size_t k = 0; k < peaks_line_vector.size(); k++)
    {
      if (peaks_line_vector[k].visited == true)
        continue;

      te::gm::Coord2D pp = m_inputParameters.m_demRasterPtr->getGrid()->gridToGeo(peaks_line_vector[k].x, peaks_line_vector[k].y);

      std::auto_ptr<te::gm::Point> point_pp(new te::gm::Point(pp.getX(), pp.getY(), srid));

      double distPoints = point_pp->distance(point_pp1.get());

      if (distPoints <= res)
      {
        peaks_line_vector[k].visited = true;

        te::app::RidgeLines::Point_T peak_3Dpoint;

        peak_3Dpoint.x = peaks_line_vector[k].x;
        peak_3Dpoint.y = peaks_line_vector[k].y;
        peak_3Dpoint.z = peaks_line_vector[k].z;

        // save segment peak
        peaksLines[n3DSegments].push_back(peak_3Dpoint);

        outRaster->setValue(peaks_line_vector[k].x, peaks_line_vector[k].y, raster_value, 0);
      }
    }
  }
  
  return;
}

std::auto_ptr<te::rst::Raster> te::app::RidgeLines::generatePeaksExtendedAPP(std::vector< std::vector<Point_T> >& peaksLines, te::rst::Raster* plainRasterMask)
{
  //create peaks raster
  std::vector<te::rst::BandProperty*> bandsProperties;
  te::rst::BandProperty* bandProp = new te::rst::BandProperty(0, te::dt::UCHAR_TYPE);
  bandProp->m_noDataValue = 255.;
  bandsProperties.push_back(bandProp);

  te::rst::Grid* grid = new te::rst::Grid(*( m_inputParameters.m_demRasterPtr->getGrid()));

  std::map<std::string, std::string> rinfo;

  te::rst::Raster* outRaster = te::rst::RasterFactory::make("MEM", grid, bandsProperties, rinfo);

  if(!outRaster)
  {
    throw(te::common::Exception("Peaks raster APP creation error."));
  }

  std::vector<te::app::RidgeLines::Point_T> base_peaks;

  //fill raster with dummy values
  for(unsigned int i = 0; i < outRaster->getNumberOfRows(); ++i)
  {
    for(unsigned int j = 0; j < outRaster->getNumberOfColumns(); ++j)
    {
      outRaster->setValue(j, i, 255., 0);
    }
  }

  //operation
  for(std::size_t i = 0; i < peaksLines.size(); ++i)
  {
    if(peaksLines[i].empty())
      break;

    double min_peak_z = 9999.;

    for(std::size_t j = 0; j < peaksLines[i].size(); ++j)
    {
      if(peaksLines[i][j].z < min_peak_z)
      {
        min_peak_z = peaksLines[i][j].z;
      }
    }

    for(std::size_t j = 0; j < peaksLines[i].size(); ++j)
    {
      te::app::RidgeLines::Point_T peak;

      peak.x = peaksLines[i][j].x;
      peak.y = peaksLines[i][j].y;
      peak.z = min_peak_z;
      peak.b = i;

      base_peaks.push_back(peak);
    }
  }

  //construct flags matrix
  constructFlagMatrix(plainRasterMask);

  //fill raster out
  for(std::size_t i = 0; i < base_peaks.size(); ++i)
  {
    resetFlagMatrix(m_inputParameters.m_demRasterPtr->getNumberOfRows(), m_inputParameters.m_demRasterPtr->getNumberOfColumns()); // Retorna matriz de flag na situacao original

    int x = (int)base_peaks[i].x;
    int y = (int)base_peaks[i].y;

    m_SLgFlag[x][y] = 1;

    double val;
    outRaster->getValue((unsigned int)x, (unsigned int)y, val);

    if(val != 255.)
      continue;

    te::app::RidgeLines::Point_T pointMin;

    if(!topToDown(x, y, 5, pointMin))
      continue;

    hillTopDraw(x, y, base_peaks[i].z, outRaster, base_peaks[i].b);
  }

  //return raster
  std::auto_ptr<te::rst::Raster> peaksAPPRaster(outRaster);

  return peaksAPPRaster;
}

void te::app::RidgeLines::constructFlagMatrix(te::rst::Raster* plainRasterMask)
{
  std::auto_ptr<te::rst::Raster> drainageRaster = createDrainageRaster();

  int numlin =  m_inputParameters.m_demRasterPtr->getNumberOfRows();
  int numcol =  m_inputParameters.m_demRasterPtr->getNumberOfColumns();

//create SLgFlag matirx
  m_SLgFlag.clear();

  for(int i = 0; i < numcol; i++)
  {
    std::vector<int> col;
    col.resize(numlin, 0);

    m_SLgFlag.push_back(col);
  }

  for(unsigned int i = 0; i < drainageRaster->getNumberOfRows(); i++)
  {
    for(unsigned int j = 0; j < drainageRaster->getNumberOfColumns(); j++)
    {
      double value;

      drainageRaster->getValue(j, i, value, 0);

      te::gm::Coord2D pt2 = drainageRaster->getGrid()->gridToGeo(j, i);
      te::gm::Coord2D pt1 = m_inputParameters.m_demRasterPtr->getGrid()->geoToGrid(pt2.getX(), pt2.getY());
      
      int col = (int)pt1.getX();
      int lin = (int)pt1.getY();

      if(col < numcol && lin < numlin && col >= 0. && lin >= 0.)
      {
        if(value == 255.) 
          m_SLgFlag[col][lin] = 3;
        else 
        m_SLgFlag[col][lin] = 1;
      }
    }
  }

  // set plain cells - HNR
  for(int ii = 0; ii < numlin; ii++)
  {
    for(int jj = 0; jj < numcol; jj++)
    {
      double val;

      plainRasterMask->getValue(jj, ii, val);

      if(val == 1.)
        m_SLgFlag[jj][ii] = 3;
    }
  }
}

void te::app::RidgeLines::resetFlagMatrix(int numlin, int numcol)
{
  // Retorna matriz de flag na situacao original
  for(int i = 0; i < numlin; i++)
  {
    for(int j = 0; j < numcol; j++)
    {
      if(m_SLgFlag[j][i] != 3) 
        m_SLgFlag[j][i] = 0;
    }
  }
}

std::auto_ptr<te::rst::Raster> te::app::RidgeLines::createDrainageRaster()
{
  //create drainage raster
  std::vector<te::rst::BandProperty*> bandsProperties;
  te::rst::BandProperty* bandProp = new te::rst::BandProperty(0, te::dt::UCHAR_TYPE);
  bandProp->m_noDataValue = 0.;
  bandsProperties.push_back(bandProp);

  const te::gm::Envelope* env = m_inputParameters.m_demRasterPtr->getExtent();

  int numlin =  m_inputParameters.m_demRasterPtr->getNumberOfRows();
  int numcol =  m_inputParameters.m_demRasterPtr->getNumberOfColumns();

  unsigned int numlinvec = (unsigned int)(env->getWidth() / m_inputParameters.m_demRasterPtr->getResolutionX());
  unsigned int numcolvec = (unsigned int)(env->getHeight() / m_inputParameters.m_demRasterPtr->getResolutionY());

  te::rst::Grid* grid = new te::rst::Grid(numcolvec, numlinvec, new te::gm::Envelope(*env), m_inputParameters.m_demRasterPtr->getSRID());

  std::map<std::string, std::string> rinfo;

  te::rst::Raster* outRaster = te::rst::RasterFactory::make("MEM", grid, bandsProperties, rinfo);

  if(!outRaster)
  {
    throw(te::common::Exception("Drainage raster creation error."));
  }

  //fill raster
  for(unsigned int i = 0; i < outRaster->getNumberOfRows(); ++i)
  {
    for(unsigned int j = 0; j < outRaster->getNumberOfColumns(); ++j)
    {
      outRaster->setValue(j, i, 0., 0);
    }
  }

  //rasterize
  if (m_inputParameters.m_drainage)
  {
    for (std::size_t t = 0; t < m_inputParameters.m_drainage->getNumGeometries(); ++t)
    {
      std::vector<te::gm::Geometry*> vecGeom;
      vecGeom.push_back(m_inputParameters.m_drainage->getGeometryN(t));

      std::vector<double> pixelValue;
      pixelValue.push_back(255.);

      outRaster->rasterize(vecGeom, pixelValue, 0);
    }
  }

  std::auto_ptr<te::rst::Raster> drainageRaster(outRaster);

  return drainageRaster;
}

std::vector<te::app::RidgeLines::Point_T> te::app::RidgeLines::selectPeaks()
{
  std::vector<te::app::RidgeLines::Point_T> peaks;

  short limit = (short) (60./m_inputParameters.m_demRasterPtr->getResolutionX()); // limite para procurar ponto de sela na grade +/- 6 pontos

  //create aux rasters
  std::auto_ptr<te::rst::Raster> imatopo;
  std::auto_ptr<te::rst::Raster> imasela;

  {
    std::vector<te::rst::BandProperty*> bandsProperties;
    te::rst::BandProperty* bandProp = new te::rst::BandProperty(0, te::dt::UCHAR_TYPE);
    bandProp->m_noDataValue = 0.;
    bandsProperties.push_back(bandProp);

    te::rst::Grid* grid = new te::rst::Grid(*( m_inputParameters.m_demRasterPtr->getGrid()));

    std::map<std::string, std::string> rinfo;

    te::rst::Raster* outRaster = te::rst::RasterFactory::make("MEM", grid, bandsProperties, rinfo);

    if(!outRaster)
    {
      throw(te::common::Exception("Imatopo raster creation error."));
    }

    imatopo.reset(outRaster);
  }

  {
    std::vector<te::rst::BandProperty*> bandsProperties;
    te::rst::BandProperty* bandProp = new te::rst::BandProperty(0, te::dt::UCHAR_TYPE);
    bandProp->m_noDataValue = 0.;
    bandsProperties.push_back(bandProp);

    te::rst::Grid* grid = new te::rst::Grid(*( m_inputParameters.m_demRasterPtr->getGrid()));

    std::map<std::string, std::string> rinfo;

    te::rst::Raster* outRaster = te::rst::RasterFactory::make("MEM", grid, bandsProperties, rinfo);

    if(!outRaster)
    {
      throw(te::common::Exception("Imasela raster creation error."));
    }

    imasela.reset(outRaster);
  }

  //fill rasters
  for(unsigned int i = 0; i < m_inputParameters.m_demRasterPtr->getNumberOfRows(); ++i)
  {
    for(unsigned int j = 0; j < m_inputParameters.m_demRasterPtr->getNumberOfColumns(); ++j)
    {
      imatopo->setValue(j, i, 1., 0);
      imasela->setValue(j, i, 1., 0);
    }
  }

  //processing
  for(unsigned int i = limit; i < m_inputParameters.m_demRasterPtr->getNumberOfRows() - limit; i++)
  {
    for(unsigned int j = limit; j < m_inputParameters.m_demRasterPtr->getNumberOfColumns() - limit; j++)
    {
      double top;

      m_inputParameters.m_demRasterPtr->getValue(j, i, top);
      
      if(top > std::numeric_limits<double>::max()) 
        continue;
  
      int numneigh = 0;
      int numvizigual =0;

      for (int ii = -2; ii <= 2; ii++) // Verifica os 8 vizinhos
      {
        for(int jj = -2; jj <= 2; jj++)
        {
          if (ii==0 && jj==0) 
            continue;

          double viz;

          m_inputParameters.m_demRasterPtr->getValue(j+jj, i+ii, viz);

          if(top == viz) 
            numvizigual++;

          if (top > viz)
            numneigh++;
        }
      }

      bool sela = false;

      if (numneigh >= 24 && numvizigual != 24)
      {
        if(i < m_inputParameters.m_demRasterPtr->getNumberOfRows() && j < m_inputParameters.m_demRasterPtr->getNumberOfColumns() && i > 4 && j >4)
        {
          // verifica se pico nao eh ponto de sela
          bool sa, sb, sc, sd, se, sf, sg, sh;

          sa = sb = sc = sd = se = sf = sg = sh = false;
          
          short tt2, tt3;

          for(short tt = 2; tt < limit; tt++)
          {
            double a, b, c, d, e, f, g, h;

            m_inputParameters.m_demRasterPtr->getValue(j-tt, i-tt, a); 
            m_inputParameters.m_demRasterPtr->getValue(j+tt,i+tt,b); 
            m_inputParameters.m_demRasterPtr->getValue(j,i-tt,c); 
            m_inputParameters.m_demRasterPtr->getValue(j,i+tt,d); 
            m_inputParameters.m_demRasterPtr->getValue(j+tt,i-tt,e);
            m_inputParameters.m_demRasterPtr->getValue(j-tt,i+tt,f); 
            m_inputParameters.m_demRasterPtr->getValue(j-tt,i,g);
            m_inputParameters.m_demRasterPtr->getValue(j+tt,i,h);

            tt2 = tt-1;

            double aa,bb,cc,dd,ee,ff,gg,hh;

            m_inputParameters.m_demRasterPtr->getValue(j-tt2,i-tt2,aa); 
            m_inputParameters.m_demRasterPtr->getValue(j+tt2,i+tt2,bb);
            m_inputParameters.m_demRasterPtr->getValue(j,i-tt2,cc); 
            m_inputParameters.m_demRasterPtr->getValue(j,i+tt2,dd); 
            m_inputParameters.m_demRasterPtr->getValue(j+tt2,i-tt2,ee); 
            m_inputParameters.m_demRasterPtr->getValue(j-tt2,i+tt2,ff); 
            m_inputParameters.m_demRasterPtr->getValue(j-tt2,i,gg); 
            m_inputParameters.m_demRasterPtr->getValue(j+tt2,i,hh);

            tt2 = tt;
            tt3 = tt+1;

            // diagonal 
            if(a>aa && !sa ) 
            {
              m_inputParameters.m_demRasterPtr->getValue(j-tt3,i-tt3,a);
              m_inputParameters.m_demRasterPtr->getValue(j-tt2,i-tt2,aa);
              
              if (a>aa)
                sa = true;
            }
  
            if( (b>bb) && !sb) 
            {
              m_inputParameters.m_demRasterPtr->getValue(j+tt3,i+tt3,b);
              m_inputParameters.m_demRasterPtr->getValue(j+tt2,i+tt2,bb);
              
              if(b>bb)
                sb = true;
            }
  
            if (sa && sb)
            {
              sela = true;
              break;
            }
  
            //vertical
            if((c>cc) && !sc) 
            {
              m_inputParameters.m_demRasterPtr->getValue(j,i-tt3,c);
              m_inputParameters.m_demRasterPtr->getValue(j,i-tt2,cc);
  
              if(c>cc) 
                sc = true;
            }

            if(!sd && (d>dd)) 
            {
              m_inputParameters.m_demRasterPtr->getValue(j,i+tt3,d);
              m_inputParameters.m_demRasterPtr->getValue(j,i+tt2,dd);
            
              if(d>dd)
                sd = true;
            }

            if(sc && sd)
            {
              sela = true;
              break;
            }

            //diagonal 
            if(e>ee &&!se) 
            {
              m_inputParameters.m_demRasterPtr->getValue(j+tt3,i-tt3,e); 
              m_inputParameters.m_demRasterPtr->getValue(j+tt2,i-tt2,ee);
  
              if(e>ee)
                se = true;
            }
  
            if(!sf && f>ff)
            {
              m_inputParameters.m_demRasterPtr->getValue(j-tt3,i+tt3,f);
              m_inputParameters.m_demRasterPtr->getValue(j-tt2,i+tt2,ff);
  
              if(f>ff)
                sf = true;
            }

            if(se && sf)
            {
              sela = true;
              break;
            }

            // horizontal
            if((g>gg) && !sg)
            {
              m_inputParameters.m_demRasterPtr->getValue(j-tt3,i,g);
              m_inputParameters.m_demRasterPtr->getValue(j-tt2,i,gg);
  
              if(g>gg) 
                sg = true;
            }

            if(!sh && (h>hh))
            {
              m_inputParameters.m_demRasterPtr->getValue(j+tt3,i,h);
              m_inputParameters.m_demRasterPtr->getValue(j+tt2,i,hh);
  
              if(h>hh)
                sh = true;
            }

            if(sg && sh)
            {
              sela = true;
              break;
            }
          }

          Point_T ptxyz;

          ptxyz.x = j;
          ptxyz.y = i;
          ptxyz.z = top;

          te::gm::Coord2D ptmax;
  
          if(!sela)
          {
            peaks.push_back(ptxyz);
            m_SLgFlag[j][i] = 1;

            imatopo->setValue(j, i, 255., 0); //provisirio
          }
          else
          {
            imasela->setValue(j, i, 255.,0); //provisirio
          }
        }
      }
    }
  }

// Ordena o vetor
  double maxv,elem;
  Point_T pv;

  for(std::size_t i = 0; i < peaks.size()-1; i++)
  {
    maxv = peaks[i].z;

    for(std::size_t j = i+1; j < peaks.size(); j++)
    {
      elem = peaks[j].z;

      if(maxv < elem)
      {
        pv = peaks[i];
        peaks[i] = peaks[j];
        peaks[j] = pv;

        maxv = elem;
      }
    }
  }

  return peaks;
}

bool te::app::RidgeLines::topToDown(int x, int y, int fase, te::app::RidgeLines::Point_T& pointMin)
{
  double seed;

  m_inputParameters.m_demRasterPtr->getValue(x, y, seed, 0);

  double window[3][3];

  m_SLgFlag[x][y] = 1;

  int nivel = 0;

  te::gm::Coord2D ptMin;

  bool flagIncLine = false;

  double minvalue = std::numeric_limits<double>::max();
  double maxvalue = seed;

  while(1)
  {
    short noseedcount = 0;

    for(int i = (y-nivel); i <= (y+nivel); i++)
    {
      if(i <= 0) 
        continue;

      if(i >= (int)m_inputParameters.m_demRasterPtr->getNumberOfRows() - 2) 
        break;

      short increm;

      // a primeira linha e a ultima sao percorridas todoas as colunas
      if(i == (y-nivel) || i == (y+nivel)) 
        increm = 1;
      else 
        increm = 2 * nivel; //somente a primeira coluna e a ultima sao calculadas

      for(int j = (x-nivel); j <= (x+nivel); j = j+increm)
      {
        if(j >= (int)m_inputParameters.m_demRasterPtr->getNumberOfColumns() - 2) 
          break;

        if(j <= 0) 
          continue;

        if(m_SLgFlag[j][i] != 1) 
        {
          noseedcount++;
          continue;
        }

        m_inputParameters.m_demRasterPtr->getValue(j, i, seed, 0);

        m_SLgFlag[j][i] = fase;

        bool foundDre = false;

        for(int i2 = -1; i2 <= 1; i2++)
        {
          for(int j2 = -1; j2 <= 1; j2++)
          {
            if(m_SLgFlag[j+j2][i+i2] == 3) 
              foundDre = true;

            m_inputParameters.m_demRasterPtr->getValue((j+j2), (i+i2), window[1+j2][1+i2], 0);
          }
        }

        // Verifica os 8 vizinhos
        short sela = 0;

        for(int i2 = -1; i2 <= 1; i2++)
        {
          for(int j2 = -1; j2 <= 1; j2++)
          {
            if(i==1 && j==1) 
              continue;

            if(fase == 6 && m_SLgFlag[j+j2][i+i2] == 5) 
              return true;

            if(m_SLgFlag[j+j2][i+i2] != 0 && m_SLgFlag[j+j2][i+i2] != 3) 
              continue;

            double grdvalue = window[j2+1][i2+1];

            if(seed >= grdvalue)
            {
              if(m_SLgFlag[j+j2][i+i2] == 3)
              {
                if(minvalue > grdvalue)
                {
                  minvalue = grdvalue;

                  te::gm::Coord2D c((double)(j+j2), (double)(i+i2));

                  ptMin = m_inputParameters.m_demRasterPtr->getGrid()->gridToGeo(c.getX(), c.getY());
                }
              }
              else if(!foundDre)
              {
                m_SLgFlag[j+j2][i+i2] = 1; // Set point to be seed
              }
            }
            else 
            {
              m_SLgFlag[j+j2][i+i2] = fase; // teste
              sela++;
            }
          }
        }

        //Ponto de sela 
        if(sela == 8)
          m_SLgFlag[j][i] = 4;

        // calculation of the slope
        if(fase == 2)
        {
          if(!flagIncLine)
          {
            double dzy = (window[0][0] - window[0][2]) + 2*(window[1][0] - window[1][2]) + (window[2][0] - window[2][2]);
            
            dzy = dzy / (8.* m_inputParameters.m_demRasterPtr->getResolutionX());

            double dzx = (window[2][0] - window[0][0]) + 2*(window[2][1] -window[0][1]) + (window[2][2] - window[0][2]);

            dzx = dzx / (8.* m_inputParameters.m_demRasterPtr->getResolutionX());

            double decl = (sqrt((double)(dzx * dzx) + (double)(dzy *dzy))* 100.);

            if(decl >= m_inputParameters.m_slope && decl < std::numeric_limits<double>::max())
              flagIncLine = true;
          }
        }
      }
    }

    if(nivel > 0 && noseedcount >= nivel * 8)
      break;

    nivel++;

    // Verifica se ainda tem ponto de semente
    bool foundSeed= false;

    for(int i = (y-nivel); i < (y+nivel); i++)
    {
      if(i >=  (int)m_inputParameters.m_demRasterPtr->getNumberOfRows()) 
        break;

      if(i<0) 
        break;

      short increm;

      if(i == (y-nivel) || i==(y+nivel)) 
        increm = 1; 
      else 
        increm = nivel;

      for(int j = (x-nivel); j <= (x+nivel); j = j+increm)
      {
        if(j >= (int)m_inputParameters.m_demRasterPtr->getNumberOfColumns())
          break;

        if(j < 0)
          continue;

        if(m_SLgFlag[j][i] == 1)
        { 
          foundSeed = true;
          break;
        }
      }

      if(foundSeed) 
        break;
    }

    if(!foundSeed) 
      break;
  }

  bool flagHeigh = false;

  if(fase == 2)
  {
    if(minvalue != std::numeric_limits<double>::max())
    {
      if((maxvalue - minvalue) >= 1800.) // It is a mountain
      {
        flagHeigh = true;
        flagIncLine = true;
      }
      else if(!flagHeigh && (maxvalue - minvalue) >= m_inputParameters.m_height)
      {
        flagHeigh = true;
      }
    }

    if(flagIncLine && flagHeigh) 
    {
      pointMin.x = ptMin.getX();
      pointMin.y = ptMin.getY();
      pointMin.z = minvalue;

      return true;
    }
    else
    {
      return false;
    }
  }

  if(fase == 6) 
    return false;  //NAO ACHOU NEHUM MORRO VIZINHO
  else 
    return true;
}

bool te::app::RidgeLines::hillTopDraw(int indc, int indl, double height, te::rst::Raster* imaout, double SLindex)
{
  if(height > std::numeric_limits<double>::max() )
    return false;

  m_SLgFlag[indc][indl] = 2;

  imaout->setValue(indc, indl, SLindex);

  int nivel = 1;

  while(1)
  {
    int lim = nivel * 8;

    int cont = 0;

    bool achoulin = false;
    bool achoucol = false;
  
    for(int i=(indl-nivel); i<=(indl+nivel); i++)
    {
      if(i<=0) 
        continue;
      
      if(i >= (int)m_inputParameters.m_demRasterPtr->getNumberOfRows() -1) 
        break;
  
      if(i == (indl-nivel) || i== (indl+nivel)) // a primeira linha e a utlma sao percorridas todoas as colunas
      {
        for (int j = (indc-nivel); j <= (indc+nivel); j++)
        {
          if(j >= (int)m_inputParameters.m_demRasterPtr->getNumberOfColumns() - 1) 
            continue;

          if(j<=0) 
            continue;

          if(m_SLgFlag[j][i] == 0 || m_SLgFlag[j][i] == 3)
          {
            cont++;
            continue;
          }

          double val;
          m_inputParameters.m_demRasterPtr->getValue(j, i, val);

          if(val >= height)
          {
            achoulin = true;
            imaout->setValue(j, i, SLindex);
          }
          else 
          {
            imaout->setValue(j, i, 255.);
          }
        }
      }
      else
      {
        for(int k = 0; k < 2; k++)
        {
          int j = indc- ((k*2 -1)*nivel);

          if(j <= 0) 
            continue;

          if(j >= (int)m_inputParameters.m_demRasterPtr->getNumberOfColumns() - 1) 
            continue;
  
          if(m_SLgFlag[j][i] == 0 || m_SLgFlag[j][i] == 3) 
          {
            cont++;
            continue;
          }
  
          double val;
          m_inputParameters.m_demRasterPtr->getValue(j, i, val);

          if(val >= height)
          {
            achoucol = true;
            imaout->setValue(j, i, SLindex);
          }
          else 
          {
            imaout->setValue(j, i, 255.);
          }
        }
      }
    }
    
  if (!achoulin && !achoucol) 
    break;

  nivel++;

  if (cont >= lim) 
    break;
  }

  return true;
}

std::vector<te::gm::Coord2D> te::app::RidgeLines::bresenhamPoints(te::gm::Coord2D gp1, te::gm::Coord2D gp2)
{
  std::vector<te::gm::Coord2D> line_points;

  int x1g, y1g, x2g, y2g, xg, yg, ddx, ddy, sx, sy, err, err2;
  te::gm::Coord2D point;

  x1g = (int)gp1.getX();
  y1g = (int)gp1.getY();
  x2g = (int)gp2.getX();
  y2g = (int)gp2.getY();

  if(x1g == x2g)
  {
    point.x = x1g;
    
    if(y1g <= y2g)
    {
      for(yg = y1g; yg <= y2g; yg++)
      {
        point.y = yg;
        line_points.push_back(point);
      }
    }
    else
    {
      for(yg = y2g;yg <= y1g;yg++)
      {
        point.y = yg;
        line_points.push_back(point);
      }
    }
  }
  else if(y1g == y2g)
  {
    point.y = y1g;
  
    if(x1g <= x2g)
    {
      for(xg = x1g;xg <= x2g;xg++)
      {
        point.x = xg;
        line_points.push_back(point);
      }
    }
    else
    for(xg = x2g; xg <= x1g; xg++)
    {
      point.x = xg;
      line_points.push_back(point);
    }
  }
  else
  {
    ddx = abs(x2g-x1g);
    ddy = abs(y2g-y1g);
    
    if(x1g < x2g)
      sx = 1;
    else
      sx = -1;
  
    if(y1g < y2g)
      sy = 1;
    else
      sy = -1;
  
    if(ddx > ddy)
      err = ddx;
    else
      err = -ddy;
  
    err /= 2;

    while(true)
    {
      point.x = x1g;
      point.y = y1g;

      line_points.push_back(point);

      if(x1g == x2g && y1g == y2g)
        break;

      err2 = err;
  
      if(err2 > -ddx)
      {
        err -= ddy;
        x1g += sx;
      }

      if(err2 < ddy)
      {
        err += ddx; 
        y1g += sy;
      }
    }
  }

  return line_points;
}