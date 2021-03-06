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
  \file terralib/app/Veredas.h
 
  \brief Generates the APP information from a veredas information.
*/

#ifndef __TERRALIB_APP_INTERNAL_VEREDAS_H
#define __TERRALIB_APP_INTERNAL_VEREDAS_H

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
      \brief APP generator by vereda.
      \details Generates the APP information from a veredas information.
    */
    class APPPLUGINDLLEXPORT Veredas : public Algorithm
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

            std::vector<te::gm::Geometry*> m_geometries;     //!< Geometries used to create the buffer information.

            double m_bufferDistance;                         //!< Distance used to create buffer (default: 50.).

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

        Veredas();
        
        ~Veredas();
       
        //overload
        bool execute(AlgorithmOutputParameters& outputParams) throw(te::common::Exception);
        
        //overload
        void reset() throw(te::common::Exception);
        
        //overload
        bool initialize(const AlgorithmInputParameters& inputParams) throw(te::common::Exception);

      protected:

        Veredas::InputParameters m_inputParameters;       //!< Veredas input execution parameters.

        Veredas::OutputParameters* m_outputParametersPtr; //!< Veredas output execution parameters.
    };

  } // end namespace app
}   // end namespace te

#endif //__TERRALIB_APP_INTERNAL_VEREDAS_H

