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
  \file terralib/app/qt/ElevationWidget.cpp

  \brief A widget used to access the elevation method.
*/


#include "ElevationWidget.h"
#include "ui_ElevationWidgetForm.h"
#include "../core/Elevation.h"

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

te::app::ElevationWidget::ElevationWidget(QWidget* parent)
  : QWidget(parent),
    m_ui(new Ui::ElevationWidgetForm)
{
// add controls
  m_ui->setupUi(this);

  m_ui->m_elevationLineEdit->setValidator(new QDoubleValidator(this));

  connect(m_ui->m_fileToolButton, SIGNAL(clicked()), this, SLOT(onFileToolButtonClicked()));
}

te::app::ElevationWidget::~ElevationWidget()
{
}

void te::app::ElevationWidget::setLayerList(std::list<te::map::AbstractLayerPtr>& layerList)
{
  m_ui->m_layersComboBox->clear();

  std::list<te::map::AbstractLayerPtr>::iterator it = layerList.begin();

  while(it != layerList.end())
  {
    te::map::AbstractLayerPtr l = *it;

    std::auto_ptr<te::da::DataSetType> dsType = l->getSchema();

    if(dsType->hasRaster())
      m_ui->m_layersComboBox->addItem(l->getTitle().c_str(), QVariant::fromValue(l));

    ++it;
  }
}

void te::app::ElevationWidget::onFileToolButtonClicked()
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

void te::app::ElevationWidget::execute()
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

  //get elevation
  if(m_ui->m_elevationLineEdit->text().isEmpty())
  {
    QMessageBox::warning(this, tr("Warning"), tr("Invalid elevation value."));
    return;
  }

  double elevationValue = m_ui->m_elevationLineEdit->text().toDouble();

  //get raster
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
  std::size_t rpos = te::da::GetFirstPropertyPos(ds.get(), te::dt::RASTER_TYPE);
  std::auto_ptr<te::rst::Raster> rst = ds->getRaster(rpos);

  if(!rst.get())
  {
    QMessageBox::warning(this, tr("Warning"), tr("Selected layer does not have raster representation."));
    return;
  }

  //get band
  int band = 0;

  //run operation
  te::qt::widgets::ProgressViewerDialog* pViewer = new te::qt::widgets::ProgressViewerDialog(this);
  int id = te::common::ProgressManager::getInstance().addViewer(pViewer);

  try
  {
    bool executeok = false;
    bool initok = false;

    // create elevation algorithm parameters
    te::app::Elevation::InputParameters inputParameters;
    inputParameters.m_elevationValue = elevationValue;
    inputParameters.m_inRasterBand = (unsigned int)band;
    inputParameters.m_inRasterPtr = rst.get();

    te::app::Elevation::OutputParameters outputParameters;
    outputParameters.m_createdOutDSName = dataSetName;
    outputParameters.m_createdOutInfo = connInfo;
    outputParameters.m_createdOutDSType = dataType;

    if (m_ui->m_outRasterCheckBox->isChecked())
    {
      std::string path = file.absolutePath().toStdString();

      std::string dataSetName = file.baseName().toStdString();

      std::map<std::string, std::string> rasterConnInfo;
      rasterConnInfo["URI"] = path + std::string("/") + dataSetName + std::string(".tif");

      outputParameters.m_createdOutRasterDSType = "GDAL";
      outputParameters.m_createdOutRasterInfo = rasterConnInfo;
    }

    // execute the algorithm
    te::app::Elevation elevation;

    initok = elevation.initialize(inputParameters);

    if(initok)
      executeok = elevation.execute(outputParameters);

    if(!executeok)
    {
      std::string errorMsg  = "An exception has occurred: ";
                  errorMsg += elevation.logMessage();

      QMessageBox::warning(this, tr("Warning"), errorMsg.c_str());

      te::common::ProgressManager::getInstance().removeViewer(id);

      delete pViewer;

      return;
    }

  }
  catch(te::common::Exception& e)
  {
    std::string errorMsg  = "Problems in Elevation operation. ";
                errorMsg += e.what();

    QMessageBox::warning(this, tr("Warning"), errorMsg.c_str());

    te::common::ProgressManager::getInstance().removeViewer(id);

    delete pViewer;

    return;
  }
  catch(...)
  {
    std::string errorMsg  = "An unexpected exception has occurred in Elevation()!";

    QMessageBox::warning(this, tr("Warning"), errorMsg.c_str());

    te::common::ProgressManager::getInstance().removeViewer(id);

    delete pViewer;

    return;
  }

  te::common::ProgressManager::getInstance().removeViewer(id);

  delete pViewer;

  QMessageBox::information(this, tr("Information"), tr("Elevation operation done."));

  //create output layer
  emit createLayer(dataType, connInfo);

  //close main window
  emit operationFinished(true);
}
