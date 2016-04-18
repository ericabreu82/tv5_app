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
  \file terralib/app/Plateaus.h
 
  \brief Generates the APP information from plateaus information.
*/

#ifndef __TERRALIB_APP_INTERNAL_PLATEAUS_H
#define __TERRALIB_APP_INTERNAL_PLATEAUS_H

// Terralib
#include "Algorithm.h"

// STL
#include <map>
#include <vector>
#include <memory>

namespace te
{
  //forward declarations
  namespace gm { class Geometry; }
  namespace rst { class Raster; }

  namespace app
  {
    /*!
      \class Plateaus
      \brief APP generator by plateaus.
      \details Generates the APP information from a plateaus information.
    */
    class APPPLUGINDLLEXPORT Plateaus : public Algorithm
    {
      public:

        /*!
          \class InputParameters
          \brief Plateaus input parameters
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

            te::rst::Raster* m_slopeRasterPtr;    //!< Input raster for slope information

            double m_bufferDistance;                    //!< Attribute used to define the buffer for this app operation

        };

        /*!
          \class OutputParameters
          \brief Plateaus output parameters
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

            std::string m_createdOutRasterDSType;                        //!< Output raster data source type (as described in te::raster::RasterFactory ).

            std::map<std::string, std::string> m_createdOutRasterInfo;   //!< The necessary information to create the raster (as described in te::raster::RasterFactory).

            std::auto_ptr<te::rst::Raster> m_createdOutRasterPtr;        //!< A pointer to the created output raster instance.

        };

        Plateaus();
        
        ~Plateaus();
       
        //overload
        bool execute(AlgorithmOutputParameters& outputParams) throw(te::common::Exception);
        
        //overload
        void reset() throw(te::common::Exception);
        
        //overload
        bool initialize(const AlgorithmInputParameters& inputParams) throw(te::common::Exception);

      protected:

        std::auto_ptr<te::rst::Raster> createSlopeMask(double threshold_slope);

        std::auto_ptr<te::rst::Raster> createHighElevationSlopeMask(double dem_threshold, double slope_threshold);

        std::auto_ptr<te::rst::Raster> createTerraceRasterMask(te::rst::Raster* highElevationSlopeMask, double area_threshold);

        bool defineTerracesfromEscarps(te::rst::Raster* terraceRaster, te::rst::Raster* escarpRaster);

      protected:

        Plateaus::InputParameters m_inputParameters;       //!< Veredas input execution parameters.

        Plateaus::OutputParameters* m_outputParametersPtr; //!< Veredas output execution parameters.
    };

  } // end namespace app
}   // end namespace te

#endif //__TERRALIB_APP_INTERNAL_PLATEAUS_H
