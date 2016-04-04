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
  \file terralib/app/Headwater.h
 
  \brief Generates the APP information from a headwater information.
*/

#ifndef __TERRALIB_APP_INTERNAL_HEADWATER_H
#define __TERRALIB_APP_INTERNAL_HEADWATER_H

// Terralib
#include "Algorithm.h"

// STL
#include <map>
#include <memory>
#include <vector>

namespace te
{
  //forward declarations
  namespace gm { class LineString; }

  namespace app
  {
    /*!
      \class Headwater
      \brief APP generator by headwater.
      \details Generates the APP information from a headwater information.
    */
    class APPPLUGINDLLEXPORT Headwater : public Algorithm
    {
      public:

        /*!
          \class InputParameters
          \brief Headwater input parameters
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

            std::vector<te::gm::LineString*> m_riversLines;     //!< Rivers Geometries used to create the headwater information.

            double m_watercourseBuffer;                         //!< Attribute with buffer distance for watercourse geometries.

            double m_headwaterBuffer;                           //!< Attribute with buffer distance for headerwater geometries.

        };

        /*!
          \class OutputParameters
          \brief Headwater output parameters
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

        Headwater();
        
        ~Headwater();
       
        //overload
        bool execute(AlgorithmOutputParameters& outputParams) throw(te::common::Exception);
        
        //overload
        void reset() throw(te::common::Exception);
        
        //overload
        bool initialize(const AlgorithmInputParameters& inputParams) throw(te::common::Exception);

      protected:

        Headwater::InputParameters m_inputParameters;       //!< Headwater input execution parameters.

        Headwater::OutputParameters* m_outputParametersPtr; //!< Headwater output execution parameters.
    };

  } // end namespace app
}   // end namespace te

#endif //__TERRALIB_APP_INTERNAL_HEADWATER_H

