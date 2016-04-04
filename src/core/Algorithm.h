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
  \file terralib/app/Algorithm.h
  \brief APP algorithm base interface class.
 */

#ifndef __TERRALIB_APP_INTERNAL_ALGORITHM_H
#define __TERRALIB_APP_INTERNAL_ALGORITHM_H

#include <terralib/common/Exception.h>

#include "../Config.h"

#include "AlgorithmInputParameters.h"
#include "AlgorithmOutputParameters.h"

namespace te
{
  namespace app
  {
    /*!
      \class Algorithm
      \brief APP algorithm base interface.
     */
    class APPPLUGINDLLEXPORT Algorithm
    {
      protected:

        Algorithm();

      public:
        
        virtual ~Algorithm();

        /*!
          \brief Initialize the algorithm instance making it ready for execution.
          
          \param inputParams Input parameters.
          
          \return true if OK, false on errors.
         */
        virtual bool initialize(const AlgorithmInputParameters& inputParams) throw(te::common::Exception) = 0;

        /*!
          \brief Returns true if the algorithm instance is initialized and ready for execution.
          
          \return true if the algorithm instance is initialized and ready for execution.
         */        
        bool isInitialized() const;

        /*!
          \brief Returns the log message from the algorithm.
          
          \return String with log information from the algorithm, used in case of errors.
         */        
        std::string logMessage() const;

        /*!
          \brief Executes the algorithm using the supplied parameters.
          
          \param outputParams Output parameters.
          
          \return true if OK, false on errors.
         */
        virtual bool execute(AlgorithmOutputParameters& outputParams) throw(te::common::Exception) = 0;

        /*!
          \brief Clear all internal allocated objects and reset the algorithm to its initial state.
         */
        virtual void reset() throw(te::common::Exception) = 0;

      protected:

        bool m_isInitialized; //!< Tells if this instance is initialized.

        std::string m_logMsg; //!< Log message used to return information from the algorithm.
    };

  } // end namespace app
}   // end namespace te

#endif  // __TERRALIB_APP_INTERNAL_ALGORITHM_H

