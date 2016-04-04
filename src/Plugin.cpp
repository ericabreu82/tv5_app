/*  Copyright (C) 2008-2013 National Institute For Space Research (INPE) - Brazil.

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
  \file terralib/qt/plugins/app/Plugin.cpp

  \brief Plugin implementation for the APP Qt Plugin widget.
*/

// TerraLib
#include <terralib/common/Translator.h>
#include <terralib/common/Logger.h>
#include <terralib/qt/af/ApplicationController.h>

#include "AppMainWindowAction.h"
#include "Plugin.h"

// QT
#include <QMenu>
#include <QMenuBar>

te::qt::plugins::app::Plugin::Plugin(const te::plugin::PluginInfo& pluginInfo)
  : QObject(), te::plugin::Plugin(pluginInfo)
{
}

te::qt::plugins::app::Plugin::~Plugin() 
{
}

void te::qt::plugins::app::Plugin::startup()
{
  if(m_initialized)
    return;

  te::qt::af::AppCtrlSingleton::getInstance().addListener(this, te::qt::af::SENDER);

  TE_LOG_TRACE(TE_TR("TerraLib Qt APP Plugin startup!"));

  // add app entry in plugin menu
  QMenu* pluginMenu = te::qt::af::AppCtrlSingleton::getInstance().getMenu("Plugins");

  // Insert action before plugin manager action
  QAction* pluginsSeparator = te::qt::af::AppCtrlSingleton::getInstance().findAction("ManagePluginsSeparator");
  
  m_appMainWindow = new te::qt::plugins::app::AppMainWindowAction(pluginMenu);
  connect(m_appMainWindow, SIGNAL(triggered(te::qt::af::evt::Event*)), SIGNAL(triggered(te::qt::af::evt::Event*)));

  pluginMenu->insertAction(pluginsSeparator, m_appMainWindow->m_action);

  m_initialized = true;
}

void te::qt::plugins::app::Plugin::shutdown()
{
  if(!m_initialized)
    return;

  TE_LOG_TRACE(TE_TR("TerraLib Qt APP Plugin shutdown!"));
  
  delete m_appMainWindow;

  m_initialized = false;
}

PLUGIN_CALL_BACK_IMPL(te::qt::plugins::app::Plugin)
