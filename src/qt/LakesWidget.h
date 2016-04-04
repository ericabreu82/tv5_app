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
  \file terralib/app/qt/LakesWidget.h

  \brief A widget used to access the lakes method.
*/

#ifndef __TERRALIB_APP_INTERNAL_LAKESWIDGET_H
#define __TERRALIB_APP_INTERNAL_LAKESWIDGET_H

#include "../Config.h"

// TerraLib
#include <terralib/maptools/AbstractLayer.h>

// Qt
#include <QComboBox>
#include <QWidget>

// STL
#include <memory>

namespace Ui { class LakesWidgetForm; }

// Forward declarations

namespace te
{
  namespace app
  {
// Forward declarations

    class APPPLUGINDLLEXPORT LakesWidget : public QWidget
    {
      Q_OBJECT

      public:

        LakesWidget(QWidget* parent = 0);

        ~LakesWidget();

      public:

        /*!
          \brief This method is used to set the list of layers
         
        */
        void setLayerList(std::list<te::map::AbstractLayerPtr>& layerList);

      public slots:

        void onFileToolButtonClicked();

        void execute();

      signals:

        void operationFinished(bool ok);

        void createLayer(const std::string& driverName, const std::map<std::string, std::string>& connInfo);

      private:

        std::auto_ptr<Ui::LakesWidgetForm> m_ui;
        std::vector<te::gm::Geometry*> layerToGeometryVector(QComboBox *layerComboBox, std::string layerTitle);

    };
  }   // end namespace app
}     // end namespace te

#endif  // __TERRALIB_APP_INTERNAL_LAKESWIDGET_H
