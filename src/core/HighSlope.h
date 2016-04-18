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
  \file terralib/app/HighSlope.h
 
  \brief Generates the APP information from a high slope raster.
  \details  A output raster mask is generated.
  - 0 value is used to define the APP area.
  - 255 is the no data value used.
*/

#ifndef __TERRALIB_APP_INTERNAL_HIGHSLOPE_H
#define __TERRALIB_APP_INTERNAL_HIGHSLOPE_H

// Terralib
#include "Algorithm.h"

// STL
#include <map>
#include <memory>

namespace te
{
  //forward declarations
  namespace rst { class Raster; }

  namespace app
  {
    /*!
      \class HighSlope
      \brief APP generator by high slope.
      \details Generates the APP information from a high slope raster.
    */
    class APPPLUGINDLLEXPORT HighSlope : public Algorithm
    {
      public:

        /*!
          \class InputParameters
          \brief Elevation input parameters
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

            te::rst::Raster const* m_inRasterPtr; //!< Input raster.

            unsigned int m_inRasterBand;          //!< Band to be processed from the input raster (default: 0).

            double m_slopeValue;                  //!< Slope Value used to extract the app information.

            double m_slopeUndef;                  //!< Value used to define a dummy slope value (default: -9999.).

        };

        /*!
          \class OutputParameters
          \brief Elevation output parameters
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

        HighSlope();
        
        ~HighSlope();
       
        //overload
        bool execute(AlgorithmOutputParameters& outputParams) throw(te::common::Exception);
        
        //overload
        void reset() throw(te::common::Exception);
        
        //overload
        bool initialize(const AlgorithmInputParameters& inputParams) throw(te::common::Exception);

      protected:

        HighSlope::InputParameters m_inputParameters;       //!< HighSlope input execution parameters.

        HighSlope::OutputParameters* m_outputParametersPtr; //!< HighSlope output execution parameters.
    };

  } // end namespace app
}   // end namespace te

#endif //__TERRALIB_APP_INTERNAL_HIGHSLOPE_H

