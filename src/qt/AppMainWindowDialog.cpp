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
  \file terralib/app/qt/AppMainWindowDialog.cpp

  \brief A dialog used to access the app methods.
*/

// TerraLib
#include <terralib/dataaccess/dataset/DataSetType.h>
#include <terralib/dataaccess/datasource/DataSource.h>
#include <terralib/dataaccess/datasource/DataSourceFactory.h>
#include <terralib/dataaccess/datasource/DataSourceInfoManager.h>
#include <terralib/dataaccess/datasource/DataSourceManager.h>
#include <terralib/qt/widgets/layer/utils/DataSet2Layer.h>

#include "AppMainWindowDialog.h"
#include "ElevationWidget.h"
#include "HighSlopeWidget.h"
#include "VeredasWidget.h"
#include "LakesWidget.h"
#include "PlateausWidget.h"
#include "HeadwaterWidget.h"
#include "WatercourseWidget.h"
#include "RidgeLinesWidget.h"
#include "HillTopWidget.h"
#include "ui_AppMainWindowDialogForm.h"

// Qt
#include <QMessageBox>

// STL
#include <memory>

// Boost
#include <boost/uuid/random_generator.hpp>
#include <boost/uuid/uuid_io.hpp>

te::app::AppMainWindowDialog::AppMainWindowDialog(QWidget* parent, Qt::WindowFlags f)
  : QDialog(parent, f),
    m_ui(new Ui::AppMainWindowDialogForm),
    m_currentWidget(0)
{
// add controls
  m_ui->setupUi(this);

  m_layout = new QGridLayout(m_ui->m_widget);

  m_layout->setContentsMargins(0,0,0,0);

  m_ui->m_widget->setLayout(m_layout);

// add icons
  m_ui->m_imgLabel->setPixmap(QIcon::fromTheme("app-icon").pixmap(112,48));

//connect tool buttons with slot functions
  connect(m_ui->m_elevationToolButton, SIGNAL(clicked()), this, SLOT(onElevationToolButtonClicked()));
  connect(m_ui->m_highSlopeToolButton, SIGNAL(clicked()), this, SLOT(onHighSlopeToolButtonClicked()));
  connect(m_ui->m_veredasToolButton, SIGNAL(clicked()), this, SLOT(onVeredasToolButtonClicked()));
  connect(m_ui->m_lakesToolButton, SIGNAL(clicked()), this, SLOT(onLakesToolButtonClicked()));
  connect(m_ui->m_plateausToolButton, SIGNAL(clicked()), this, SLOT(onPlateausToolButtonClicked()));
  connect(m_ui->m_headwaterToolButton, SIGNAL(clicked()), this, SLOT(onHeadwaterToolButtonClicked()));
  connect(m_ui->m_watercourseToolButton, SIGNAL(clicked()), this, SLOT(onWatercourseToolButtonClicked()));
  connect(m_ui->m_ridgeLinesToolButton, SIGNAL(clicked()), this, SLOT(onRidgeLinesToolButtonClicked()));
  connect(m_ui->m_hillTopToolButton, SIGNAL(clicked()), this, SLOT(onHillTopToolButtonClicked()));
}

te::app::AppMainWindowDialog::~AppMainWindowDialog()
{
  if(m_currentWidget)
    delete m_currentWidget;
}

void te::app::AppMainWindowDialog::setLayerList(std::list<te::map::AbstractLayerPtr>& layerList)
{
  m_layerList = layerList;
}

te::map::AbstractLayerPtr te::app::AppMainWindowDialog::getOutputLayer()
{
  return m_outputLayer;
}

void te::app::AppMainWindowDialog::onElevationToolButtonClicked()
{
  if(m_currentWidget)
    delete m_currentWidget;

  m_ui->m_elevationToolButton->setChecked(true);

  te::app::ElevationWidget* widget = new te::app::ElevationWidget(m_ui->m_widget);

  widget->setLayerList(m_layerList);

  m_currentWidget = widget;

  m_layout->addWidget(widget);

  widget->show();

//connect signal and slots
  connect(m_ui->m_okPushButton, SIGNAL(clicked()), widget, SLOT(execute()));
  connect(widget, SIGNAL(operationFinished(bool)), this, SLOT(onOperationFinished(bool)));
  connect(widget, SIGNAL(createLayer(const std::string&, const std::map<std::string, std::string>&)), this, SLOT(onCreateLayer(const std::string&, const std::map<std::string, std::string>&)));
}

void te::app::AppMainWindowDialog::onHighSlopeToolButtonClicked()
{
  if(m_currentWidget)
    delete m_currentWidget;

  m_ui->m_highSlopeToolButton->setChecked(true);

  te::app::HighSlopeWidget* widget = new te::app::HighSlopeWidget(m_ui->m_widget);

  widget->setLayerList(m_layerList);

  m_currentWidget = widget;

  m_layout->addWidget(widget);

  widget->show();

//connect signal and slots
  connect(m_ui->m_okPushButton, SIGNAL(clicked()), widget, SLOT(execute()));
  connect(widget, SIGNAL(operationFinished(bool)), this, SLOT(onOperationFinished(bool)));
  connect(widget, SIGNAL(createLayer(const std::string&, const std::map<std::string, std::string>&)), this, SLOT(onCreateLayer(const std::string&, const std::map<std::string, std::string>&)));
}

void te::app::AppMainWindowDialog::onVeredasToolButtonClicked()
{
  if(m_currentWidget)
    delete m_currentWidget;

  m_ui->m_veredasToolButton->setChecked(true);

  te::app::VeredasWidget* widget = new te::app::VeredasWidget(m_ui->m_widget);

  widget->setLayerList(m_layerList);

  m_currentWidget = widget;

  m_layout->addWidget(widget);

  widget->show();

//connect signal and slots
  connect(m_ui->m_okPushButton, SIGNAL(clicked()), widget, SLOT(execute()));
  connect(widget, SIGNAL(operationFinished(bool)), this, SLOT(onOperationFinished(bool)));
  connect(widget, SIGNAL(createLayer(const std::string&, const std::map<std::string, std::string>&)), this, SLOT(onCreateLayer(const std::string&, const std::map<std::string, std::string>&)));
}

void te::app::AppMainWindowDialog::onLakesToolButtonClicked()
{
  if(m_currentWidget)
    delete m_currentWidget;

  m_ui->m_lakesToolButton->setChecked(true);

  te::app::LakesWidget* widget = new te::app::LakesWidget(m_ui->m_widget);

  widget->setLayerList(m_layerList);

  m_currentWidget = widget;

  m_layout->addWidget(widget);

  widget->show();

//connect signal and slots
  connect(m_ui->m_okPushButton, SIGNAL(clicked()), widget, SLOT(execute()));
  connect(widget, SIGNAL(operationFinished(bool)), this, SLOT(onOperationFinished(bool)));
  connect(widget, SIGNAL(createLayer(const std::string&, const std::map<std::string, std::string>&)), this, SLOT(onCreateLayer(const std::string&, const std::map<std::string, std::string>&)));
}

void te::app::AppMainWindowDialog::onPlateausToolButtonClicked()
{
  if(m_currentWidget)
    delete m_currentWidget;

  m_ui->m_plateausToolButton->setChecked(true);

  te::app::PlateausWidget* widget = new te::app::PlateausWidget(m_ui->m_widget);

  widget->setLayerList(m_layerList);

  m_currentWidget = widget;

  m_layout->addWidget(widget);

  widget->show();

//connect signal and slots
  connect(m_ui->m_okPushButton, SIGNAL(clicked()), widget, SLOT(execute()));
  connect(widget, SIGNAL(operationFinished(bool)), this, SLOT(onOperationFinished(bool)));
  connect(widget, SIGNAL(createLayer(const std::string&, const std::map<std::string, std::string>&)), this, SLOT(onCreateLayer(const std::string&, const std::map<std::string, std::string>&)));
}

void te::app::AppMainWindowDialog::onHeadwaterToolButtonClicked()
{
  if(m_currentWidget)
    delete m_currentWidget;

  m_ui->m_headwaterToolButton->setChecked(true);

  te::app::HeadwaterWidget* widget = new te::app::HeadwaterWidget(m_ui->m_widget);

  widget->setLayerList(m_layerList);

  m_currentWidget = widget;

  m_layout->addWidget(widget);

  widget->show();

//connect signal and slots
  connect(m_ui->m_okPushButton, SIGNAL(clicked()), widget, SLOT(execute()));
  connect(widget, SIGNAL(operationFinished(bool)), this, SLOT(onOperationFinished(bool)));
  connect(widget, SIGNAL(createLayer(const std::string&, const std::map<std::string, std::string>&)), this, SLOT(onCreateLayer(const std::string&, const std::map<std::string, std::string>&)));
}

void te::app::AppMainWindowDialog::onWatercourseToolButtonClicked()
{
  if(m_currentWidget)
    delete m_currentWidget;

  m_ui->m_watercourseToolButton->setChecked(true);

  te::app::WatercourseWidget* widget = new te::app::WatercourseWidget(m_ui->m_widget);

  widget->setLayerList(m_layerList);

  m_currentWidget = widget;

  m_layout->addWidget(widget);

  widget->show();

//connect signal and slots
  connect(m_ui->m_okPushButton, SIGNAL(clicked()), widget, SLOT(execute()));
  connect(widget, SIGNAL(operationFinished(bool)), this, SLOT(onOperationFinished(bool)));
  connect(widget, SIGNAL(createLayer(const std::string&, const std::map<std::string, std::string>&)), this, SLOT(onCreateLayer(const std::string&, const std::map<std::string, std::string>&)));
}

void te::app::AppMainWindowDialog::onRidgeLinesToolButtonClicked()
{
  if(m_currentWidget)
    delete m_currentWidget;

  m_ui->m_ridgeLinesToolButton->setChecked(true);

  te::app::RidgeLinesWidget* widget = new te::app::RidgeLinesWidget(m_ui->m_widget);

  widget->setLayerList(m_layerList);

  m_currentWidget = widget;

  m_layout->addWidget(widget);

  widget->show();

//connect signal and slots
  connect(m_ui->m_okPushButton, SIGNAL(clicked()), widget, SLOT(execute()));
  connect(widget, SIGNAL(operationFinished(bool)), this, SLOT(onOperationFinished(bool)));
  connect(widget, SIGNAL(createLayer(const std::string&, const std::map<std::string, std::string>&)), this, SLOT(onCreateLayer(const std::string&, const std::map<std::string, std::string>&)));
}

void te::app::AppMainWindowDialog::onHillTopToolButtonClicked()
{
  if(m_currentWidget)
    delete m_currentWidget;

  m_ui->m_hillTopToolButton->setChecked(true);

  te::app::HillTopWidget* widget = new te::app::HillTopWidget(m_ui->m_widget);

  widget->setLayerList(m_layerList);

  m_currentWidget = widget;

  m_layout->addWidget(widget);

  widget->show();

//connect signal and slots
  connect(m_ui->m_okPushButton, SIGNAL(clicked()), widget, SLOT(execute()));
  connect(widget, SIGNAL(operationFinished(bool)), this, SLOT(onOperationFinished(bool)));
  connect(widget, SIGNAL(createLayer(const std::string&, const std::map<std::string, std::string>&)), this, SLOT(onCreateLayer(const std::string&, const std::map<std::string, std::string>&)));
}

void te::app::AppMainWindowDialog::onOperationFinished(bool ok)
{
  if(ok)
    accept();
}

void te::app::AppMainWindowDialog::onCreateLayer(const std::string& driverName, const std::map<std::string, std::string>& connInfo)
{
  static boost::uuids::basic_random_generator<boost::mt19937> gen;

  boost::uuids::uuid valU = gen();
  std::string id = boost::uuids::to_string(valU);

  std::auto_ptr<te::da::DataSource> ds(te::da::DataSourceFactory::make(driverName));
  ds->setConnectionInfo(connInfo);
  ds->open();

  std::vector<std::string> dsNames = ds->getDataSetNames();
  assert(!dsNames.empty());

  //add ds info
  te::da::DataSourceInfoPtr dsInfoPtr(new te::da::DataSourceInfo);
  dsInfoPtr->setConnInfo(connInfo);
  dsInfoPtr->setId(id);
  dsInfoPtr->setTitle(dsNames[0]);
  dsInfoPtr->setAccessDriver(driverName);
  dsInfoPtr->setType(driverName);

  te::da::DataSourceInfoManager::getInstance().add(dsInfoPtr);

  //add ds
  te::da::DataSourcePtr dsPtr(ds.release());
  dsPtr->setId(id);
  
  te::da::DataSourceManager::getInstance().insert(dsPtr);

  //create layer
  te::da::DataSetTypePtr dsType(dsPtr->getDataSetType(dsNames[0]).release());

  te::qt::widgets::DataSet2Layer ds2l(dsPtr->getId());

  m_outputLayer = ds2l(dsType);
}
