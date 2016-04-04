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
  \file terralib/app/AlgorithmInputParameters.h
  \brief APP algorithm input parameters base interface.
 */

#ifndef __TERRALIB_APP_INTERNAL_ALGORITHM_INPUT_PARAMETERS_H
#define __TERRALIB_APP_INTERNAL_ALGORITHM_INPUT_PARAMETERS_H

#include "../Config.h"

#include <terralib/common/AbstractParameters.h>

namespace te
{
  namespace app
  {
    /*!
      \class AlgorithmInputParameters
      \brief APP algorithm input parameters base interface.
    */
    class APPPLUGINDLLEXPORT AlgorithmInputParameters : public te::common::AbstractParameters
    {
      protected:

        AlgorithmInputParameters();

      public:
        
        virtual ~AlgorithmInputParameters();
    };
  } // end namespace app
}   // end namespace te

#endif  // __TERRALIB_APP_INTERNAL_ALGORITHM_INPUT_PARAMETERS_H

