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
  \file terralib/qt/plugins/app/AppMainWindowAction.cpp

  \brief This file defines the class for Main Window Action
*/

// Terralib
#include <terralib/qt/af/events/LayerEvents.h>
#include <terralib/qt/af/ApplicationController.h>

#include "qt/AppMainWindowDialog.h"
#include "AppMainWindowAction.h"

//STL
#include <cassert>

te::qt::plugins::app::AppMainWindowAction::AppMainWindowAction(QMenu *menu) : 
  QObject(), 
  m_menu(menu)
{
  assert(m_menu);

  m_action = new QAction(m_menu);

  m_action->setText("APP...");
  m_action->setIcon(QIcon::fromTheme("app-icon"));

  connect(m_action, SIGNAL(triggered(bool)), this, SLOT(onActionActivated(bool)));

  m_menu->addAction(m_action);
}

te::qt::plugins::app::AppMainWindowAction::~AppMainWindowAction()
{
  delete m_action;
}

std::list<te::map::AbstractLayerPtr> te::qt::plugins::app::AppMainWindowAction::getLayerList()
{
  te::qt::af::evt::GetAvailableLayers e;

  emit triggered(&e);

  std::list<te::map::AbstractLayerPtr> allLayers = e.m_layers;

  std::list<te::map::AbstractLayerPtr> layers;

  for (std::list<te::map::AbstractLayerPtr>::iterator it = allLayers.begin(); it != allLayers.end(); ++it)
  {
    if ((*it)->isValid())
      layers.push_back(*it);
  }

  return layers;
}

void te::qt::plugins::app::AppMainWindowAction::addNewLayer(te::map::AbstractLayerPtr layer)
{
  te::qt::af::evt::LayerAdded evt(layer.get());

  emit triggered(&evt);
}

void te::qt::plugins::app::AppMainWindowAction::onActionActivated(bool checked)
{
  //get list of layers
  std::list<te::map::AbstractLayerPtr> list = getLayerList();

  //create main dialog
  te::app::AppMainWindowDialog dlg(te::qt::af::AppCtrlSingleton::getInstance().getMainWindow());

  dlg.setLayerList(list);

  if(dlg.exec() == QDialog::Accepted)
  {
    te::map::AbstractLayerPtr layer = dlg.getOutputLayer();

    if(layer.get())
      addNewLayer(layer);
  }
}
