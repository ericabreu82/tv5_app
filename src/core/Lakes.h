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
  \file terralib/app/Lakes.h
 
  \brief Generates the APP information from a lake information.

  \note Using intersection with the FIRST geometry that intersects the lake geometry. Once found the urban
  intersection, the rural intersection will not be calculated.
*/

#ifndef __TERRALIB_APP_INTERNAL_LAKES_H
#define __TERRALIB_APP_INTERNAL_LAKES_H

// Terralib
#include "Algorithm.h"

// STL
#include <map>
#include <vector>

namespace te
{
  //forward declarations
  namespace gm { class Geometry; }

  namespace app
  {
    /*!
      \class Elevation
      \brief APP generator by lakes.
      \details Generates the APP information from a lakes information.
    */
    class APPPLUGINDLLEXPORT Lakes : public Algorithm
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
            const InputParameters& operator=( const InputParameters& params );

            //overload
            AbstractParameters* clone() const;

          public:

            std::vector<te::gm::Geometry*> m_lakesGeometries;     //!< Lakes Geometries used to create the buffer information.

            std::vector<te::gm::Geometry*> m_urbanGeometries;     //!< Lakes Geometries used to create the buffer information.

            std::vector<te::gm::Geometry*> m_ruralGeometries;     //!< Lakes Geometries used to create the buffer information.

            double m_urbanBufferDistance;                         //!< Distance used to create buffer from urban geometry (default: 30.).

            double m_ruralBigBufferDistance;                      //!< Distance used to create buffer from rural big geometry(default: 100.).

            double m_ruralSmallBufferDistance;                    //!< Distance used to create buffer from rural small geometry(default: 50.).
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
        };

        Lakes();
        
        ~Lakes();
       
        //overload
        bool execute(AlgorithmOutputParameters& outputParams) throw(te::common::Exception);
        
        //overload
        void reset() throw(te::common::Exception);
        
        //overload
        bool initialize(const AlgorithmInputParameters& inputParams) throw(te::common::Exception);

      protected:

        double calculateArea(te::gm::Geometry* geom);

      protected:

        Lakes::InputParameters m_inputParameters;       //!< Lakes input execution parameters.

        Lakes::OutputParameters* m_outputParametersPtr; //!< Lakes output execution parameters.
    };

  } // end namespace app
}   // end namespace te

#endif //__TERRALIB_APP_INTERNAL_ELEVATION_H

