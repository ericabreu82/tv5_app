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
  \file terralib/app/Watercourse.cpp
 
  \brief Generates the APP information from a Watercourse information.
*/

#include "Watercourse.h"
#include "Utils.h"

// TerraLib
#include <terralib/common/progress/TaskProgress.h>
#include <terralib/common/STLUtils.h>
#include <terralib/geometry/Envelope.h>
#include <terralib/geometry/LineString.h>
#include <terralib/geometry/MultiPolygon.h>
#include <terralib/geometry/Point.h>
#include <terralib/geometry/Polygon.h>

// STL
#include <cmath>


te::gm::Coord2D TeGetMiddlePoint(te::gm::Point* first, te::gm::Point* last);

bool TeRotate(te::gm::Coord2D pr, te::gm::LineString* l, double angle, te::gm::LineString* lOut);

bool TeAdjustSegment(te::gm::Point* P0, te::gm::Point* P1, double d0, te::gm::Coord2D &P0out, te::gm::Coord2D &P1out);

bool TeSegmentsIntersectPoint(const te::gm::Coord2D& fr0, const te::gm::Coord2D& to0, const te::gm::Coord2D& fr1, const te::gm::Coord2D& to1, te::gm::Coord2D& pi);

te::app::Watercourse::InputParameters::InputParameters()
{
  reset();
}

te::app::Watercourse::InputParameters::~InputParameters()
{
  reset();
}

void te::app::Watercourse::InputParameters::reset() throw(te::common::Exception)
{
  m_riversPolygons.clear();
  m_propPolygons.clear();
  m_roadsPolygons.clear();
}

te::common::AbstractParameters* te::app::Watercourse::InputParameters::clone() const
{
  return new InputParameters(*this);
}

te::app::Watercourse::OutputParameters::OutputParameters()
{
  reset();
}

te::app::Watercourse::OutputParameters::OutputParameters(const OutputParameters& other)
{
  reset();

  operator=(other);
}

te::app::Watercourse::OutputParameters::~OutputParameters()
{
  reset();
}

void te::app::Watercourse::OutputParameters::reset() throw(te::common::Exception)
{
  m_createdOutDSName.clear();
  m_createdOutDSType.clear();
  m_createdOutInfo.clear();
}

const te::app::Watercourse::OutputParameters& te::app::Watercourse::OutputParameters::operator=(const Watercourse::OutputParameters& params)
{
  reset();

  m_createdOutDSName = params.m_createdOutDSName;
  m_createdOutDSType = params.m_createdOutDSType;
  m_createdOutInfo = params.m_createdOutInfo;

  return *this;
}

te::common::AbstractParameters* te::app::Watercourse::OutputParameters::clone() const
{
  return new OutputParameters(*this);
}

te::app::Watercourse::Watercourse()
{
  reset();
}

te::app::Watercourse::~Watercourse()
{
}

bool te::app::Watercourse::execute(AlgorithmOutputParameters& outputParams) throw(te::common::Exception)
{
  if(!m_isInitialized)
  {
    m_logMsg = "Algoritm not initialized.";

    return false;
  }

   m_outputParametersPtr = dynamic_cast<Watercourse::OutputParameters*>(&outputParams);

  if(!m_outputParametersPtr)
  {
    m_logMsg = "Invalid output parameters.";

    return false;
  }

  // Execute the operation
  te::common::TaskProgress task("Generating Watercourse APP");
  task.setTotalSteps(m_inputParameters.m_riversPolygons.size());

  std::vector<te::gm::Geometry*> resultGeomVec;

  for(std::size_t t = 0; t < m_inputParameters.m_riversPolygons.size(); ++t)
  {
    te::gm::Polygon* riverPoly = m_inputParameters.m_riversPolygons[t];

    if(riverPoly->isValid() == false)
      continue;

    te::gm::LineString* lsCurve = dynamic_cast<te::gm::LineString*>(riverPoly->getRingN(0));

    te::gm::Geometry* geomBufferResult = 0;

    for(std::size_t q = 0; q < lsCurve->getNPoints() - 1; ++q)
    {
      std::auto_ptr<te::gm::Point> pointStart(lsCurve->getPointN(q));
      std::auto_ptr<te::gm::Point> pointEnd(lsCurve->getPointN(q + 1));

      std::auto_ptr<te::gm::LineString> segment(new te::gm::LineString(2, te::gm::LineStringType, riverPoly->getSRID()));

      segment->setPoint(0, pointStart->getX(), pointStart->getY());
      segment->setPoint(1, pointEnd->getX(), pointEnd->getY());

      double distance = getDistance(segment.get(), riverPoly);

      std::auto_ptr<te::gm::Geometry> geomBuffer(segment->buffer(distance, 16, te::gm::CapButtType));

      if(geomBufferResult == 0)
      {
        geomBufferResult = geomBuffer.release();
      }
      else
      {
        te::gm::Geometry* geomAux = geomBufferResult->Union(geomBuffer.get());

        delete geomBufferResult;

        geomBufferResult = geomAux;
      }
    }

    te::gm::Geometry* result = riverPoly->Union(geomBufferResult);

    resultGeomVec.push_back(result);

    if(!task.isActive())
    {
      m_logMsg = "Operation canceled.";

      return false;
    }

    task.pulse();
  }

  //union all
  te::gm::Geometry* geomBufferResult = 0;

  for(std::size_t t = 0; t < resultGeomVec.size(); ++t)
  {
     if(geomBufferResult == 0)
      {
        geomBufferResult = resultGeomVec[t];
      }
      else
      {
        te::gm::Geometry* geomAux = geomBufferResult->Union(resultGeomVec[t]);

        delete geomBufferResult;

        geomBufferResult = geomAux;

        delete resultGeomVec[t];
      }
  }

  resultGeomVec.clear();

  resultGeomVec.push_back(geomBufferResult);

  //checkIntersections(resultGeomVec);

  //create vector representation
  te::app::ExportVector(resultGeomVec, m_outputParametersPtr->m_createdOutDSName, m_outputParametersPtr->m_createdOutDSType,  m_outputParametersPtr->m_createdOutInfo);

  return true;
}

void te::app::Watercourse::reset() throw(te::common::Exception)
{
  m_inputParameters.reset();
  m_outputParametersPtr = 0;

  m_lastWidth = 0.;
  m_lastDistance = 0.;
}

bool te::app::Watercourse::initialize(const AlgorithmInputParameters& inputParams) throw(te::common::Exception)
{
  reset();

  Watercourse::InputParameters const* inputParamsPtr = dynamic_cast<Watercourse::InputParameters const*>(&inputParams);

  if(!inputParamsPtr)
  {
    m_logMsg = "Invalid Parameters.";

    return false;
  }

  m_inputParameters = *inputParamsPtr;

  m_isInitialized = true;

  return true;
}

double te::app::Watercourse::getDistance(te::gm::LineString* segment, te::gm::Polygon* poly)
{
  //get center of segment
  te::gm::Coord2D middle = TeGetMiddlePoint(segment->getStartPoint(), segment->getEndPoint());

  //rotate line
  te::gm::LineString* lineRotate = new te::gm::LineString(2, te::gm::LineStringType, segment->getSRID());

  TeRotate(middle, segment, 90., lineRotate);

  double width = 0;

  const te::gm::Envelope* env = poly->getMBR();

  //create big line 
  double diagonal = std::sqrt((env->getWidth() * env->getWidth()) + (env->getHeight() * env->getHeight()));

  std::auto_ptr<te::gm::Point> pStart(lineRotate->getStartPoint());
  std::auto_ptr<te::gm::Point> pEnd(lineRotate->getEndPoint());

  te::gm::Coord2D cStart, cEnd;

  TeAdjustSegment(pStart.get(), pEnd.get(), diagonal, cStart, cEnd);

  te::gm::Coord2D cStart2, cEnd2;

  TeAdjustSegment(pEnd.get(), pStart.get(), diagonal, cStart2, cEnd2);

  //verificar segmentos que interceptam essa nova linha aumentada
  te::gm::LineString* lsCurve = dynamic_cast<te::gm::LineString*>(poly->getRingN(0));

  for(std::size_t q = 0; q < lsCurve->getNPoints() - 1; ++q)
  {
    std::auto_ptr<te::gm::Point> pointStart(lsCurve->getPointN(q));
    std::auto_ptr<te::gm::Point> pointEnd(lsCurve->getPointN(q + 1));

    te::gm::Coord2D c1(pointStart->getX(), pointStart->getY());
    te::gm::Coord2D c2(pointEnd->getX(), pointEnd->getY());

    te::gm::Coord2D cInter;

    if(TeSegmentsIntersectPoint(cEnd, cEnd2, c1, c2, cInter))
    {
      //create line using middle coord and cInter coord to check if this line is covered by polygon
      te::gm::LineString* lsCover = new te::gm::LineString(2, te::gm::LineStringType, segment->getSRID());
      lsCover->setPoint(0, middle.getX(), middle.getY());
      lsCover->setPoint(1, cInter.getX(), cInter.getY());

      if(lsCover->coveredBy(poly))
      {
        //calculate distance
        std::auto_ptr<te::gm::Point> pA(lsCover->getStartPoint());
        std::auto_ptr<te::gm::Point> pB(lsCover->getEndPoint());

        width = pA->distance(pB.get());

        break;
      }
    }
  }

  //get oposite point of center point in polygon

  double distance;

  if(width <=10.)
  {
    distance = 30;
  }
  else if(width > 10. && width <= 50.) 
  {
    distance = 50;
  }
  else if (width >50. && width <= 200.)
  {
    distance = 100;
  }
  else if(width > 200. && width <= 600.)
  {
    distance = 200;
  }
  else
  {
    distance = 500;
  }

  if(m_lastWidth != 0. && m_lastDistance != 0.)
  {
    if(m_lastDistance != distance && (std::abs(width - m_lastWidth) > m_lastWidth * .5))
    {
      width = m_lastWidth;
      distance = m_lastDistance;
    }
  }

  m_lastWidth = width;

  m_lastDistance = distance;

  return distance;
}

void te::app::Watercourse::checkIntersections(std::vector<te::gm::Geometry*>& geomVec)
{
  //create container with all geometries
  te::gm::MultiPolygon* multiPol = new te::gm::MultiPolygon(0, te::gm::PolygonType, geomVec[0]->getSRID());

  te::gm::MultiPolygon* multiPolResult = new te::gm::MultiPolygon(0, te::gm::PolygonType, geomVec[0]->getSRID());

  for(std::size_t t = 0; t < geomVec.size(); ++t)
  {
    multiPol->add(geomVec[t]);
  }

  for(std::size_t t = 0; t < m_inputParameters.m_propPolygons.size(); ++t)
  {
    te::gm::Polygon* geomProp = m_inputParameters.m_propPolygons[t];

    te::gm::Geometry* geomRiver = 0;

    //get river with intersection with property
    for(std::size_t q = 0; q < m_inputParameters.m_riversPolygons.size(); ++q)
    {
      te::gm::Geometry* g = m_inputParameters.m_riversPolygons[q];

      if(geomProp->intersects(g))
      {
        geomRiver = g;
        break;
      }
    }

    if(!geomRiver)
      continue;

    //get road with intersection with property
    double roadArea = 0.;
    for(std::size_t p = 0; p < m_inputParameters.m_roadsPolygons.size(); ++p)
    {
      te::gm::Geometry* geomRoad = m_inputParameters.m_roadsPolygons[p];

      if(geomProp->intersects(geomRoad))
      {
        std::auto_ptr<te::gm::Geometry> geomRoadProp(geomProp->intersection(geomRoad));

        te::gm::Polygon* p = dynamic_cast<te::gm::Polygon*>(geomRoadProp.get());

        if(p)
        {
          roadArea += p->getArea();
        }
      }
    }

    double area = geomProp->getArea() - roadArea;

    double appBuffer;

    if(area <= m_inputParameters.m_fiscalModule)
    {
      appBuffer = 5;
    }
    else if(area > m_inputParameters.m_fiscalModule && area <= 2*m_inputParameters.m_fiscalModule)
    {
      appBuffer = 8;
    }
    else if(area > 2*m_inputParameters.m_fiscalModule && area <= 4*m_inputParameters.m_fiscalModule)
    {
      appBuffer = 15;
    }
    else 
    {
      appBuffer = 30;
    }

    if(geomProp->intersects(geomRiver))
    {
      std::auto_ptr<te::gm::Geometry> geomPropRiver(geomProp->intersection(geomRiver));

      te::gm::Geometry* gBuffer = geomPropRiver->buffer(appBuffer);

      te::gm::Geometry* geomMultiPropDiff = multiPol->difference(geomProp);

      te::gm::Geometry* geomResult = geomMultiPropDiff->Union(gBuffer);

      multiPolResult->add(geomResult);
    }
    else
    {
      te::gm::LineString* lsProp = dynamic_cast<te::gm::LineString*>(geomProp->getRingN(0));

      std::vector< std::pair<double, double> > vecPairXY;

      for(std::size_t p = 0; p < lsProp->getNPoints(); ++p)
      {
        std::auto_ptr<te::gm::Point> point(lsProp->getPointN(p));

        if(point->touches(geomRiver))
        {
          std::pair<double, double> pairXY;
          pairXY.first = point->getX();
          pairXY.second = point->getY();

          vecPairXY.push_back(pairXY);
        }
      }

      std::auto_ptr<te::gm::LineString> lsTouches(new te::gm::LineString(vecPairXY.size(), te::gm::LineStringType, geomProp->getSRID()));

      for(std::size_t p = 0; p < vecPairXY.size(); ++p)
      {
        lsTouches->setPoint(p, vecPairXY[p].first, vecPairXY[p].second);
      }

      te::gm::Geometry* gBuffer = lsTouches->buffer(appBuffer);

      te::gm::Geometry* geomMultiPropDiff = multiPol->difference(geomProp);

      te::gm::Geometry* geomResult = geomMultiPropDiff->Union(gBuffer);

      multiPolResult->add(geomResult);
    }
  }

  //union of multiPolResult
}

te::gm::Coord2D TeGetMiddlePoint(te::gm::Point* first, te::gm::Point* last)
{
  te::gm::Coord2D middle;

  double lenght,parts,curlenght,incx,incy, deltax,deltay,dx,dy;
  short i,nparts;

  lenght = first->distance(last);

  if(lenght == 0.)
  {
    middle.x = first->getX();
    middle.y = first->getY();

    return middle;
  }

  nparts = 2;
  parts = lenght/2.;

  // verify segment orientation
  if(first->getX() < last->getX())
    incx = 1.;
  else
    incx = -1.;

  if(first->getY() < last->getY())
    incy = 1.;
  else
    incy = -1.;

  curlenght = 0.;

  deltax = fabs(first->getX()-last->getX());
  deltay = fabs(first->getY()-last->getY());

  for(i=0;i<(nparts-1);i++)
  {
    curlenght = curlenght + parts;

    // vertical segment
    if(first->getX() == last->getX())
    {
      middle.x = first->getX();
      middle.y = first->getY()+(curlenght*incy);
      continue;
    }

    // horizontal segment
    if(first->getY() == last->getY())
    {
      middle.x = first->getX()+(curlenght*incx);
      middle.y = first->getY();
      continue;
    }

    // inclined segment

    // calculating X coordinate
    dx = curlenght*deltax/lenght;

    // calculating Y coordinate
    dy = curlenght*deltay/lenght;
    middle.x = first->getX()+(dx*incx);
    middle.y = first->getY()+(dy*incy);
  }

  return middle;
}

bool TeRotate(te::gm::Coord2D pr, te::gm::LineString* l, double angle, te::gm::LineString* lOut)
{
  double alfa;
  double dx, dy;
  double x, y;
  double xr,yr;

  try
  {
    if(l->size() < 2)
      return (false);

    alfa = (4.*atan(1.)*angle) / 180.;//transform Degree to Radius

    dx = pr.getX();
    dy = pr.getY();

    for(std::size_t count = 0; count < l->size(); count++)
    {
      std::auto_ptr<te::gm::Point> curPoint(l->getPointN(count));

      x = curPoint->getX() - dx;
      y = curPoint->getY() - dy;

      xr = x * cos(alfa) - y * sin(alfa);
      yr = x * sin(alfa) + y * cos(alfa);

      lOut->setPoint(count, xr+dx, yr+dy);
    }
  }
  catch(...)
  { 
    return false;
  }

  return true;
}

bool TeAdjustSegment(te::gm::Point* P0, te::gm::Point* P1, double d0, te::gm::Coord2D& P0out, te::gm::Coord2D& P1out)
{
  double vL_norm1;
  double vL_norm2;

  try
  {
    if(P0->getX() == P1->getX() && P0->getY() == P1->getY())
      return false;

    te::gm::Coord2D vL1((P1->getX() - P0->getX()), (P1->getY() - P0->getY()));
    te::gm::Coord2D vL2(-1 * (P1->getX() - P0->getX()), -1 * (P1->getY() - P0->getY()));
    vL_norm1 = sqrt(vL1.getX() * vL1.getX() + vL1.getY() * vL1.getY());
    vL_norm2 = sqrt(vL2.getX() * vL2.getX() + vL2.getY() * vL2.getY());

    te::gm::Coord2D uL1((vL1.getX() / vL_norm1 ), ( vL1.getY() / vL_norm1));
    te::gm::Coord2D uL2((vL2.getX() / vL_norm2 ), ( vL2.getY() / vL_norm2));

    te::gm::Point* pFim1 = new te::gm::Point(P0->getX() + uL1.getX() * d0, P0->getY() + uL1.getY() * d0, P0->getSRID());
    te::gm::Point* pFim2 = new te::gm::Point(P0->getX() + uL2.getX() * d0, P0->getY() + uL2.getY() * d0, P1->getSRID());
    
    if(pFim1->distance(P1) <= pFim2->distance(P1))
    {
      P0out.x = P0->getX();
      P0out.y = P0->getY();
      P1out.x = pFim1->getX();
      P1out.y = pFim1->getY();
    }
    else
    {
      P0out.x = P0->getX();
      P0out.y = P0->getY();
      P1out.x = pFim2->getX();
      P1out.y = pFim2->getY();
    }
  }
  catch(...)
  { 
    return false;
  }

  return true;
}

bool TeSegmentsIntersectPoint(const te::gm::Coord2D& fr0, const te::gm::Coord2D& to0, const te::gm::Coord2D& fr1, const te::gm::Coord2D& to1, te::gm::Coord2D& pi)
{
  //Adapted from TWO-DIMENSIONAL CLIPPING: A VECTOR-BASED APPROACH
  //Hans J.W. Spoelder, Fons H. Ullings in:
  //Graphics Gems I, pp.701, 

  double	a, b, c,
  px, py, lx, ly, lpx, lpy,
  s;

  px  = to0.x - fr0.x;
  py  = to0.y - fr0.y;
  lx  = to1.x - fr1.x;
  ly  = to1.y - fr1.y;
  lpx = fr0.x - fr1.x;
  lpy = fr0.y - fr1.y;

  a = py * lx - px * ly;
  b = lpx * ly - lpy * lx;
  c = lpx * py - lpy * px;

  if (a == 0)// Linhas paralelas
  {
    return false;
  }
  else
  {
    if (a > 0)
    {
      if ((b < 0) || (b > a) || (c < 0) || (c > a))
        return false;
    }
    else
    {
      if ((b > 0) || (b < a) || (c > 0) || (c < a))
        return false;
    }

    s = b/a;

    pi.x = fr0.x + (px*s);
    pi.y = fr0.y + (py*s);
  }

  return true;
}
