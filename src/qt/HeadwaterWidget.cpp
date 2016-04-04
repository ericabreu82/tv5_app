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
  \file terralib/app/qt/HeadwaterWidget.cpp

  \brief A widget used to access the headwater method.
*/

#include "HeadwaterWidget.h"
#include "ui_HeadwaterWidgetForm.h"

#include "../core/Headwater.h"

// TerraLib
#include <terralib/common/progress/ProgressManager.h>
#include <terralib/common/Exception.h>
#include <terralib/dataaccess/utils/Utils.h>
#include <terralib/geometry/MultiLineString.h>
#include <terralib/qt/widgets/progress/ProgressViewerDialog.h>

// Qt
#include <QFileDialog>
#include <QMessageBox>
#include <QValidator>

Q_DECLARE_METATYPE(te::map::AbstractLayerPtr);

te::app::HeadwaterWidget::HeadwaterWidget(QWidget* parent)
  : QWidget(parent),
    m_ui(new Ui::HeadwaterWidgetForm)
{
// add controls
  m_ui->setupUi(this);

  m_ui->m_watercourseBufferLineEdit->setValidator(new QDoubleValidator(this));
  m_ui->m_headwaterBufferLineEdit->setValidator(new QDoubleValidator(this));

  connect(m_ui->m_fileToolButton, SIGNAL(clicked()), this, SLOT(onFileToolButtonClicked()));
}

te::app::HeadwaterWidget::~HeadwaterWidget()
{
}

void te::app::HeadwaterWidget::setLayerList(std::list<te::map::AbstractLayerPtr>& layerList)
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

void te::app::HeadwaterWidget::onFileToolButtonClicked()
{
  QString fileName = QFileDialog::getSaveFileName(this, tr("Save data to File"), "", tr("Shape File (*.shp *.SHP)"));

  if (fileName.isEmpty())
    return;

  QFileInfo file(fileName);

  if(file.suffix().isEmpty())
    fileName.append(".shp");

  m_ui->m_fileLineEdit->setText(fileName);
}

std::vector<te::gm::LineString*> te::app::HeadwaterWidget::layerToGeometryVector(QComboBox* layerComboBox, std::string layerTitle) {

    std::vector<te::gm::LineString*> geomVec;

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

        if(g->getGeomTypeId() == te::gm::MultiLineStringType)
        {
          te::gm::MultiLineString* mls = dynamic_cast<te::gm::MultiLineString*>(g.release());

          te::gm::LineString* ls = dynamic_cast<te::gm::LineString*>(mls->getGeometryN(0));

          geomVec.push_back(ls);
        }
        else if(g->getGeomTypeId() == te::gm::LineStringType)
        {
          te::gm::LineString* ls = dynamic_cast<te::gm::LineString*>(g.release());

          geomVec.push_back(ls);
        }
    }

    return geomVec;
}

void te::app::HeadwaterWidget::execute()
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
  if(m_ui->m_watercourseBufferLineEdit->text().isEmpty())
  {
    QMessageBox::warning(this, tr("Warning"), tr("Invalid watercourse buffer distance value."));
    return;
  }

  if(m_ui->m_headwaterBufferLineEdit->text().isEmpty())
  {
    QMessageBox::warning(this, tr("Warning"), tr("Invalid headwater buffer distance value."));
    return;
  }

  double watercourseBufferValue = m_ui->m_watercourseBufferLineEdit->text().toDouble();
  double headwaterBufferValue = m_ui->m_headwaterBufferLineEdit->text().toDouble();

  //get vector
  std::vector<te::gm::LineString*> riversGeomVec = layerToGeometryVector(m_ui->m_layersComboBox, "Rivers");

  if(riversGeomVec.empty())
      return;

  //run operation
  te::qt::widgets::ProgressViewerDialog* pViewer = new te::qt::widgets::ProgressViewerDialog(this);
  int id = te::common::ProgressManager::getInstance().addViewer(pViewer);

  try
  {
    bool executeok = false;
    bool initok = false;

    // create Headwater algorithm parameters
    te::app::Headwater::InputParameters inputParameters;
    inputParameters.m_riversLines = riversGeomVec;
    inputParameters.m_watercourseBuffer = watercourseBufferValue;
    inputParameters.m_headwaterBuffer = headwaterBufferValue;

    te::app::Headwater::OutputParameters outputParameters;
    outputParameters.m_createdOutDSName = dataSetName;
    outputParameters.m_createdOutInfo = connInfo;
    outputParameters.m_createdOutDSType = dataType;

    // execute the algorithm
    te::app::Headwater headwater;

    initok = headwater.initialize(inputParameters);

    if(initok)
      executeok = headwater.execute(outputParameters);

    if(!executeok)
    {
      std::string errorMsg  = "An exception has occurred: ";
                  errorMsg += headwater.logMessage();

      QMessageBox::warning(this, tr("Warning"), errorMsg.c_str());

      te::common::ProgressManager::getInstance().removeViewer(id);

      delete pViewer;

      return;
    }

  }
  catch(te::common::Exception& e)
  {
    std::string errorMsg  = "Problems in Headwater operation. ";
                errorMsg += e.what();

    QMessageBox::warning(this, tr("Warning"), errorMsg.c_str());

    te::common::ProgressManager::getInstance().removeViewer(id);

    delete pViewer;

    return;
  }
  catch(...)
  {
    std::string errorMsg  = "An unexpected exception has occurred in Headwater()!";

    QMessageBox::warning(this, tr("Warning"), errorMsg.c_str());

    te::common::ProgressManager::getInstance().removeViewer(id);

    delete pViewer;

    return;
  }

  te::common::ProgressManager::getInstance().removeViewer(id);

  delete pViewer;

  QMessageBox::information(this, tr("Information"), tr("Headwater operation done."));

  //create output layer
  emit createLayer(dataType, connInfo);

  //close main window
  emit operationFinished(true);
}
