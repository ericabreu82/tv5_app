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
  \file terralib/app/HillTop.h
 
  \brief Generates the APP information from hill top information.
*/

#ifndef __TERRALIB_APP_INTERNAL_HILLTOP_H
#define __TERRALIB_APP_INTERNAL_HILLTOP_H

// Terralib
#include "Algorithm.h"

// STL
#include <map>
#include <vector>
#include <memory>

namespace te
{
  //forward declarations
  namespace gm { class MultiLineString; class MultiPolygon; }
  namespace rst { class Raster; }

  namespace app
  {
    /*!
      \class HillTop
      \brief APP generator by hill top.
      \details Generates the APP information from a hill top information.
    */
    class APPPLUGINDLLEXPORT HillTop : public Algorithm
    {
      public:

        /*!
          \class InputParameters
          \brief HillTop input parameters
        */
        class APPPLUGINDLLEXPORT InputParameters : public AlgorithmInputParameters
        {
          public:

            InputParameters();

            ~InputParameters();

            //overload
            void reset() throw(te::common::Exception);

            //overload
            AbstractParameters* clone() const;

          public:

            te::rst::Raster const* m_demRasterPtr;      //!< Input raster for dem information

            te::rst::Raster const* m_slopeRasterPtr;    //!< Input raster for slope information

            te::gm::MultiLineString* m_contourLines;    //!< Input level curves

            te::gm::MultiLineString* m_drainage;        //!< Input drainage lines 

            te::gm::MultiPolygon* m_watershedsBuffer;   //!< Input watersheds buffer

            te::gm::MultiPolygon* m_watersheds;         //!< Input watersheds

            double m_slope;

            double m_height;

        };

        /*!
          \class OutputParameters
          \brief HillTop output parameters
        */
        class APPPLUGINDLLEXPORT OutputParameters : public AlgorithmOutputParameters
        {
          public:

            OutputParameters();

            OutputParameters(const OutputParameters&);

            ~OutputParameters();

            //overload
            void reset() throw(te::common::Exception);

            //overload
            const OutputParameters& operator=(const OutputParameters& params);

            //overload
            AbstractParameters* clone() const;

        public:
            std::string m_createdOutDSName;                              //!< Output data set name.

            std::string m_createdOutDSType;                              //!< Output data source type.

            std::map<std::string, std::string> m_createdOutInfo;         //!< The necessary information to create the output data.

        };

        typedef struct 
        {
          double x,y,z,b;
          bool visited;
        } Point_T;

        typedef struct 
        {
          int idx;
          double x,y,z,d;
          bool visited;
        } Point_Peak;
        
        HillTop();
        
        ~HillTop();
       
        //overload
        bool execute(AlgorithmOutputParameters& outputParams) throw(te::common::Exception);
        
        //overload
        void reset() throw(te::common::Exception);
        
        //overload
        bool initialize(const AlgorithmInputParameters& inputParams) throw(te::common::Exception);

      protected:

        std::auto_ptr<te::rst::Raster> createSlopeMask(double threshold_slope);

        std::auto_ptr<te::rst::Raster> createPlainsRasterMask(te::rst::Raster* slopeMask, double area_threshold);

        std::auto_ptr<te::rst::Raster> generateHillTopVector(te::rst::Raster* plainRasterMask);

        void selectIsolatedPeaks(std::vector<Point_T>& peaksVector, std::vector<Point_T>& isolatedPeaks);

//-------------------------------------------------------------------------------
//		Initialize array flag 
//
//		Flag Matrix : 1 - Drenage
//					  0 - others
//-------------------------------------------------------------------------------
        void constructFlagMatrix(te::rst::Raster* plainRasterMask);

//-------------------------------------------------------------------------------
//		Returns array flag in the original situation
//		Input :
//			numlin,numcol:	size of array
//-------------------------------------------------------------------------------
        void resetFlagMatrix(int numlin, int numcol);

        std::auto_ptr<te::rst::Raster> createDrainageRaster();

        std::vector<Point_T> selectPeaks();

        //-------------------------------------------------------------------------------
//		Algorithm to go top to bottom
//		Input :
//			x,y: Coordinates the hilltop
//			phase	= 2 - states: down the hill and up the hill to check
//					= 5 - only delimit the hill
//					= 6 - down the hill to find the hill 5
//-------------------------------------------------------------------------------
        bool topToDown(int x, int y, int fase, Point_T& pointMin);

//-------------------------------------------------------------------------------
//		Delimit APP hilltop
//		Input :
//			indc,indl: Coordinates the hilltop
//			height: Height limit the elevation to be considered APP
//			TRUE: 	Normal terminal
//			FALSE:  APP does not exist
//-------------------------------------------------------------------------------
        bool hillTopDraw(int indc, int indl, double height, te::rst::Raster* imaout);

      protected:

        HillTop::InputParameters m_inputParameters;       //!< Veredas input execution parameters.

        HillTop::OutputParameters* m_outputParametersPtr; //!< Veredas output execution parameters.

        std::vector< std::vector<int> > m_SLgFlag;
    };

  } // end namespace app
}   // end namespace te

#endif //__TERRALIB_APP_INTERNAL_HILLTOP_H
