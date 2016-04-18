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
  \file terralib/app/qt/LakesWidget.cpp

  \brief A widget used to access the lakes method.
*/

#include "LakesWidget.h"
#include "ui_LakesWidgetForm.h"

#include "../core/Lakes.h"

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

te::app::LakesWidget::LakesWidget(QWidget* parent)
  : QWidget(parent),
    m_ui(new Ui::LakesWidgetForm)
{
// add controls
  m_ui->setupUi(this);

  m_ui->m_urbanBufferLineEdit->setValidator(new QDoubleValidator(this));
  m_ui->m_ruralSmallBufferLineEdit->setValidator(new QDoubleValidator(this));
  m_ui->m_ruralBigBufferLineEdit->setValidator(new QDoubleValidator(this));

  connect(m_ui->m_fileToolButton, SIGNAL(clicked()), this, SLOT(onFileToolButtonClicked()));
}

te::app::LakesWidget::~LakesWidget()
{
}

void te::app::LakesWidget::setLayerList(std::list<te::map::AbstractLayerPtr>& layerList)
{
  m_ui->m_layersComboBox->clear();

  std::list<te::map::AbstractLayerPtr>::iterator it = layerList.begin();

  while(it != layerList.end())
  {
    te::map::AbstractLayerPtr l = *it;

    std::auto_ptr<te::da::DataSetType> dsType = l->getSchema();

    if(dsType->hasGeom())
      m_ui->m_layersComboBox->addItem(l->getTitle().c_str(), QVariant::fromValue(l));
      m_ui->m_ruralLayersComboBox->addItem(l->getTitle().c_str(), QVariant::fromValue(l));
      m_ui->m_urbanLayersComboBox->addItem(l->getTitle().c_str(), QVariant::fromValue(l));

    ++it;
  }

}

void te::app::LakesWidget::onFileToolButtonClicked()
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

std::vector<te::gm::Geometry*> te::app::LakesWidget::layerToGeometryVector(QComboBox* layerComboBox, std::string layerTitle) {

    std::vector<te::gm::Geometry*> geomVec;

    if(layerComboBox->currentText().isEmpty())
    {

        QMessageBox::warning(this, tr("Warning"), tr((layerTitle + "Layer not selected.").c_str()));
        return geomVec;
    }

    int idx = layerComboBox->currentIndex();
    QVariant varLayer = layerComboBox->itemData(idx, Qt::UserRole);
    te::map::AbstractLayerPtr layer = varLayer.value<te::map::AbstractLayerPtr>();

    if(!layer.get())
    {
        QMessageBox::warning(this, tr("Warning"), tr(("Invalid " + layerTitle + " layer selected.").c_str()));
        return geomVec;
    }

    std::auto_ptr<te::da::DataSet> ds = layer->getData();

    if(!ds.get())
    {
        QMessageBox::warning(this, tr("Warning"), tr(("Invalid " + layerTitle + " input layer selected.").c_str()));
        return geomVec;
    }

    std::size_t gpos = te::da::GetFirstPropertyPos(ds.get(), te::dt::GEOMETRY_TYPE);


    ds->moveBeforeFirst();

    while(ds->moveNext())
    {
        std::auto_ptr<te::gm::Geometry> g(ds->getGeometry(gpos));
        geomVec.push_back(g.release());
    }

    return geomVec;
}

void te::app::LakesWidget::execute()
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
  if(m_ui->m_urbanBufferLineEdit->text().isEmpty())
  {
    QMessageBox::warning(this, tr("Warning"), tr("Invalid urban buffer distance value."));
    return;
  }

  if(m_ui->m_ruralSmallBufferLineEdit->text().isEmpty())
  {
    QMessageBox::warning(this, tr("Warning"), tr("Invalid rural small buffer distance value."));
    return;
  }

  if(m_ui->m_ruralBigBufferLineEdit->text().isEmpty())
  {
    QMessageBox::warning(this, tr("Warning"), tr("Invalid rural big buffer distance value."));
    return;
  }

  double urbanBufferValue = m_ui->m_urbanBufferLineEdit->text().toDouble();
  double ruralSmallBufferValue = m_ui->m_ruralSmallBufferLineEdit->text().toDouble();
  double ruralBigBufferValue = m_ui->m_ruralBigBufferLineEdit->text().toDouble();

  //get vector
  std::vector<te::gm::Geometry*> lakesGeomVec = layerToGeometryVector(m_ui->m_layersComboBox, "Lakes");
  std::vector<te::gm::Geometry*> urbanGeomVec = layerToGeometryVector(m_ui->m_urbanLayersComboBox, "Urban Areas");
  std::vector<te::gm::Geometry*> ruralGeomVec = layerToGeometryVector(m_ui->m_ruralLayersComboBox, "Rural Areas");

  if(lakesGeomVec.empty() || urbanGeomVec.empty() || ruralGeomVec.empty())
      return;

  //run operation
  te::qt::widgets::ProgressViewerDialog* pViewer = new te::qt::widgets::ProgressViewerDialog(this);
  int id = te::common::ProgressManager::getInstance().addViewer(pViewer);

  try
  {
    bool executeok = false;
    bool initok = false;

    // create Lakes algorithm parameters
    te::app::Lakes::InputParameters inputParameters;
    inputParameters.m_lakesGeometries = lakesGeomVec;
    inputParameters.m_urbanGeometries = urbanGeomVec;
    inputParameters.m_ruralGeometries = ruralGeomVec;
    inputParameters.m_urbanBufferDistance = urbanBufferValue;
    inputParameters.m_ruralBigBufferDistance = ruralBigBufferValue;
    inputParameters.m_ruralSmallBufferDistance = ruralSmallBufferValue;

    te::app::Lakes::OutputParameters outputParameters;
    outputParameters.m_createdOutDSName = dataSetName;
    outputParameters.m_createdOutInfo = connInfo;
    outputParameters.m_createdOutDSType = dataType;

    // execute the algorithm
    te::app::Lakes lakes;

    initok = lakes.initialize(inputParameters);

    if(initok)
      executeok = lakes.execute(outputParameters);

    if(!executeok)
    {
      std::string errorMsg  = "An exception has occurred: ";
                  errorMsg += lakes.logMessage();

      QMessageBox::warning(this, tr("Warning"), errorMsg.c_str());

      te::common::ProgressManager::getInstance().removeViewer(id);

      delete pViewer;

      return;
    }

  }
  catch(te::common::Exception& e)
  {
    std::string errorMsg  = "Problems in Lakes operation. ";
                errorMsg += e.what();

    QMessageBox::warning(this, tr("Warning"), errorMsg.c_str());

    te::common::ProgressManager::getInstance().removeViewer(id);

    delete pViewer;

    return;
  }
  catch(...)
  {
    std::string errorMsg  = "An unexpected exception has occurred in Lakes()!";

    QMessageBox::warning(this, tr("Warning"), errorMsg.c_str());

    te::common::ProgressManager::getInstance().removeViewer(id);

    delete pViewer;

    return;
  }

  te::common::ProgressManager::getInstance().removeViewer(id);

  delete pViewer;

  QMessageBox::information(this, tr("Information"), tr("Lakes operation done."));

  //create output layer
  emit createLayer(dataType, connInfo);

  //close main window
  emit operationFinished(true);
}
