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
  \file terralib/app/qt/VeredasWidget.cpp

  \brief A widget used to access the veredas method.
*/

#include "VeredasWidget.h"
#include "ui_VeredasWidgetForm.h"

#include "../core/Veredas.h"

// TerraLib
#include <terralib/common/progress/ProgressManager.h>
#include <terralib/common/Exception.h>
#include <terralib/dataaccess/utils/Utils.h>
#include <terralib/qt/widgets/progress/ProgressViewerDialog.h>
#include <terralib/qt/widgets/Utils.h>

// Qt
#include <QFileDialog>
#include <QMessageBox>
#include <QValidator>

Q_DECLARE_METATYPE(te::map::AbstractLayerPtr);

te::app::VeredasWidget::VeredasWidget(QWidget* parent)
  : QWidget(parent),
    m_ui(new Ui::VeredasWidgetForm)
{
// add controls
  m_ui->setupUi(this);

  m_ui->m_bufferLineEdit->setValidator(new QDoubleValidator(this));

  connect(m_ui->m_fileToolButton, SIGNAL(clicked()), this, SLOT(onFileToolButtonClicked()));
}

te::app::VeredasWidget::~VeredasWidget()
{
}

void te::app::VeredasWidget::setLayerList(std::list<te::map::AbstractLayerPtr>& layerList)
{
  m_ui->m_layersComboBox->clear();

  std::list<te::map::AbstractLayerPtr>::iterator it = layerList.begin();

  while(it != layerList.end())
  {
    te::map::AbstractLayerPtr l = *it;

    std::auto_ptr<te::da::DataSetType> dsType = l->getSchema();

    if(dsType->hasGeom())
      m_ui->m_layersComboBox->addItem(l->getTitle().c_str(), QVariant::fromValue(l));

    ++it;
  }

}

void te::app::VeredasWidget::onFileToolButtonClicked()
{
  QString fileName = QFileDialog::getSaveFileName(this, tr("Save data to File"), te::qt::widgets::GetFilePathFromSettings("AppPlugin"), tr("Shape File (*.shp *.SHP)"));

  if (fileName.isEmpty())
    return;

  QFileInfo file(fileName);

  if(file.suffix().isEmpty())
    fileName.append(".shp");

  te::qt::widgets::AddFilePathToSettings(file.absolutePath(), "AppPlugin");

  m_ui->m_fileLineEdit->setText(fileName);
}

void te::app::VeredasWidget::execute()
{
  //get output file name
  if(m_ui->m_fileLineEdit->text().isEmpty())
  {
    QMessageBox::warning(this, tr("Warning"), tr("Output file name not defined."));
    return;
  }

  QFileInfo file( m_ui->m_fileLineEdit->text());

  std::map<std::string, std::string> connInfo;
  connInfo["URI"] = file.absoluteFilePath().toStdString();

  std::string dataSetName = file.baseName().toStdString();

  std::string dataType = "OGR";

  //get buffer distance value
  if(m_ui->m_bufferLineEdit->text().isEmpty())
  {
    QMessageBox::warning(this, tr("Warning"), tr("Invalid buffer distance value."));
    return;
  }

  double bufferValue = m_ui->m_bufferLineEdit->text().toDouble();

  //get vector
  if(m_ui->m_layersComboBox->currentText().isEmpty())
  {
    QMessageBox::warning(this, tr("Warning"), tr("Layer not selected."));
    return;
  }

  int idx = m_ui->m_layersComboBox->currentIndex();
  QVariant varLayer = m_ui->m_layersComboBox->itemData(idx, Qt::UserRole);
  te::map::AbstractLayerPtr layer = varLayer.value<te::map::AbstractLayerPtr>();

  if(!layer.get())
  {
    QMessageBox::warning(this, tr("Warning"), tr("Invalid layer selected."));
    return;
  }

  std::auto_ptr<te::da::DataSet> ds = layer->getData();

  if(!ds.get())
  {
      QMessageBox::warning(this, tr("Warning"), tr("Invalid input layer selected."));
      return;
  }

  std::size_t gpos = te::da::GetFirstPropertyPos(ds.get(), te::dt::GEOMETRY_TYPE);
  std::vector<te::gm::Geometry*> geomVec;

  ds->moveBeforeFirst();

  while(ds->moveNext())
  {
    std::auto_ptr<te::gm::Geometry> g(ds->getGeometry(gpos));
    geomVec.push_back(g.release());
  }


  //run operation
  te::qt::widgets::ProgressViewerDialog* pViewer = new te::qt::widgets::ProgressViewerDialog(this);
  int id = te::common::ProgressManager::getInstance().addViewer(pViewer);

  try
  {
    bool executeok = false;
    bool initok = false;

    // create veredas algorithm parameters
    te::app::Veredas::InputParameters inputParameters;
    inputParameters.m_geometries = geomVec;
    inputParameters.m_bufferDistance = bufferValue;

    te::app::Veredas::OutputParameters outputParameters;
    outputParameters.m_createdOutDSName = dataSetName;
    outputParameters.m_createdOutInfo = connInfo;
    outputParameters.m_createdOutDSType = dataType;

    // execute the algorithm
    te::app::Veredas veredas;

    initok = veredas.initialize(inputParameters);

    if(initok)
      executeok = veredas.execute(outputParameters);

    if(!executeok)
    {
      std::string errorMsg  = "An exception has occurred: ";
                  errorMsg += veredas.logMessage();

      QMessageBox::warning(this, tr("Warning"), errorMsg.c_str());

      te::common::ProgressManager::getInstance().removeViewer(id);

      delete pViewer;

      return;
    }

  }
  catch(te::common::Exception& e)
  {
    std::string errorMsg  = "Problems in Veredas operation. ";
                errorMsg += e.what();

    QMessageBox::warning(this, tr("Warning"), errorMsg.c_str());

    te::common::ProgressManager::getInstance().removeViewer(id);

    delete pViewer;

    return;
  }
  catch(...)
  {
    std::string errorMsg  = "An unexpected exception has occurred in Veredas()!";

    QMessageBox::warning(this, tr("Warning"), errorMsg.c_str());

    te::common::ProgressManager::getInstance().removeViewer(id);

    delete pViewer;

    return;
  }

  te::common::ProgressManager::getInstance().removeViewer(id);

  delete pViewer;

  QMessageBox::information(this, tr("Information"), tr("Veredas operation done."));

  //create output layer
  emit createLayer(dataType, connInfo);

  //close main window
  emit operationFinished(true);
}
