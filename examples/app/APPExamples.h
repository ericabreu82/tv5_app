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
  \file APPExamples.h

  \brief These routines show how to use the APP module.
 */

#include "../Config.h"

#ifndef __TERRALIB_EXAMPLES_INTERNAL_APP_EXAMPLES_H
#define __TERRALIB_EXAMPLES_INTERNAL_APP_EXAMPLES_H

/*! \brief It loads the data source drivers. */
void LoadModules();

/* \brief Elevation example. */
void Elevation();

/* \brief High Slope example. */
void HighSlope();

/* \brief Veredas example. */
void Veredas();

/* \brief Lakes example. */
void Lakes();

/* \brief Plateaus example. */
void Plateaus();

/* \brief Headwater example. */
void Headwater();

/* \brief Ridge Lines example. */
void RidgeLines();

/* \brief Watercourse example. */
void Watercourse();

/* \brief HillTop example. */
void HillTop();

#endif //__TERRALIB_EXAMPLES_INTERNAL_APP_EXAMPLES_H
