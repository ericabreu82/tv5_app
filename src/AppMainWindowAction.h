/*  Copyright (C) 2011-2012 National Institute For Space Research (INPE) - Brazil.
  
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
  \file terralib/qt/plugins/app/AppMainWindowAction.h

  \brief This file defines the class for Main Window Action
*/

#ifndef TE_QT_PLUGINS_APP_APPMAINWINDOWACTION_H
#define TE_QT_PLUGINS_APP_APPMAINWINDOWACTION_H

// Terralib
#include <terralib/maptools/AbstractLayer.h>
#include "Config.h"

// Qt
#include <QObject>
#include <QMenu>
#include <QAction>

namespace te
{
  namespace qt
  {
    namespace af
    {
      namespace evt
      {
        struct Event;
      }
    }

    namespace plugins
    {
      namespace app
      {
        
        class AppMainWindowAction : public QObject
        {
          Q_OBJECT

          public:

            AppMainWindowAction(QMenu* menu);

            ~AppMainWindowAction();

          protected:

            std::list<te::map::AbstractLayerPtr> getLayerList();

            void addNewLayer(te::map::AbstractLayerPtr layer);

          protected slots:

            virtual void onActionActivated(bool checked);

          Q_SIGNALS:

            void triggered(te::qt::af::evt::Event* e);

          protected:

            QMenu* m_menu;          //!< Parent Menu.

          public:

            QAction* m_action;      //!< Action used to call the process.
        };
        
      } // namespace app
    } // namespace plugins
  } // namespace qt
} // namespace te

#endif // TE_QT_PLUGINS_APP_APPMAINWINDOWACTION_H
