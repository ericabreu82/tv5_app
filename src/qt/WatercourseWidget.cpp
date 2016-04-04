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
  \file terralib/app/qt/WatercourseWidget.cpp

  \brief A widget used to access the headwater method.
*/

#include "WatercourseWidget.h"
#include "ui_WatercourseWidgetForm.h"

#include "../core/Watercourse.h"
#include "../core/Utils.h"

// TerraLib
#include <terralib/common/progress/ProgressManager.h>
#include <terralib/common/Exception.h>
#include <terralib/dataaccess/utils/Utils.h>
#include <terralib/geometry/GeometryCollection.h>
#include <terralib/geometry/GeometryProperty.h>
#include <terralib/qt/widgets/progress/ProgressViewerDialog.h>

// Qt
#include <QFileDialog>
#include <QMessageBox>
#include <QValidator>

Q_DECLARE_METATYPE(te::map::AbstractLayerPtr);

te::app::WatercourseWidget::WatercourseWidget(QWidget* parent)
  : QWidget(parent),
    m_ui(new Ui::WatercourseWidgetForm)
{
// add controls
  m_ui->setupUi(this);

  m_ui->m_fiscalModuleLineEdit->setValidator(new QDoubleValidator(this));

  connect(m_ui->m_fileToolButton, SIGNAL(clicked()), this, SLOT(onFileToolButtonClicked()));
}

te::app::WatercourseWidget::~WatercourseWidget()
{
}

void te::app::WatercourseWidget::setLayerList(std::list<te::map::AbstractLayerPtr>& layerList)
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

void te::app::WatercourseWidget::onFileToolButtonClicked()
{
  QString fileName = QFileDialog::getSaveFileName(this, tr("Save data to File"), "", tr("Shape File (*.shp *.SHP)"));

  if (fileName.isEmpty())
    return;

  QFileInfo file(fileName);

  if(file.suffix().isEmpty())
    fileName.append(".shp");

  m_ui->m_fileLineEdit->setText(fileName);
}

std::vector<te::gm::Polygon*> te::app::WatercourseWidget::layerToGeometryVector(QComboBox* layerComboBox, std::string layerTitle) {

    std::vector<te::gm::Geometry*> geomVec;
    std::vector<te::gm::Polygon*> polyVec;

    if(layerComboBox->currentText().isEmpty())
    {

        QMessageBox::warning(this, tr("Warning"), tr((layerTitle + "Layer not selected.").c_str()));
        return polyVec;
    }

    int idx = layerComboBox->currentIndex();
    QVariant varLayer = layerComboBox->itemData(idx, Qt::UserRole);
    te::map::AbstractLayerPtr layer = varLayer.value<te::map::AbstractLayerPtr>();

    if(!layer.get())
    {
        QMessageBox::warning(this, tr("Warning"), tr(("Invalid " + layerTitle + " layer selected.").c_str()));
        return polyVec;
    }

    std::auto_ptr<te::da::DataSet> ds = layer->getData();
    std::auto_ptr<te::da::DataSetType> dsType = layer->getSchema();


    if(!ds.get())
    {
        QMessageBox::warning(this, tr("Warning"), tr(("Invalid " + layerTitle + " input layer selected.").c_str()));
        return polyVec;
    }

    std::size_t gpos = te::da::GetFirstPropertyPos(ds.get(), te::dt::GEOMETRY_TYPE);
    te::gm::GeometryProperty* geomProp = te::da::GetFirstGeomProperty(dsType.get());

    //get geometries
    te::gm::Geometry* geomSeed = 0;

    te::gm::GeometryCollection* geomColl = new te::gm::GeometryCollection(0, te::gm::GeometryCollectionType, geomProp->getSRID());

    ds->moveBeforeFirst();

    while(ds->moveNext())
    {
        std::auto_ptr<te::gm::Geometry> g(ds->getGeometry(gpos));

        if(geomSeed == 0)
        {
          geomSeed = g.release();
        }
        else
        {
          geomColl->add(g.release());
        }
    }

    //union of all geometries
    te::gm::Geometry* resultGeom = geomSeed->Union(geomColl);

    //get polygons
    te::app::Multi2Single(resultGeom, geomVec);

    for(std::size_t t = 0; t < geomVec.size(); ++t)
    {
      te::gm::Polygon* p = dynamic_cast<te::gm::Polygon*>(geomVec[t]);

      polyVec.push_back(p);
    }

    return polyVec;
}

void te::app::WatercourseWidget::execute()
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
  if(m_ui->m_fiscalModuleLineEdit->text().isEmpty())
  {
    QMessageBox::warning(this, tr("Warning"), tr("Invalid Fiscal Module value."));
    return;
  }

  //Fiscal Modules converted from hectares to square meters
  double fiscalModuleValue = m_ui->m_fiscalModuleLineEdit->text().toDouble() * 10000.;

  //get vector
  std::vector<te::gm::Polygon*> riversPolyVec = layerToGeometryVector(m_ui->m_layersComboBox, "Rivers");

  if(riversPolyVec.empty())
      return;

  //run operation
  te::qt::widgets::ProgressViewerDialog* pViewer = new te::qt::widgets::ProgressViewerDialog(this);
  int id = te::common::ProgressManager::getInstance().addViewer(pViewer);

  try
  {
    bool executeok = false;
    bool initok = false;

    // create Watercourse algorithm parameters
    te::app::Watercourse::InputParameters inputParameters;
    inputParameters.m_riversPolygons = riversPolyVec;
    inputParameters.m_fiscalModule = fiscalModuleValue;

    te::app::Watercourse::OutputParameters outputParameters;
    outputParameters.m_createdOutDSName = dataSetName;
    outputParameters.m_createdOutInfo = connInfo;
    outputParameters.m_createdOutDSType = dataType;

    // execute the algorithm
    te::app::Watercourse watercourse;

    initok = watercourse.initialize(inputParameters);

    if(initok)
      executeok = watercourse.execute(outputParameters);

    if(!executeok)
    {
      std::string errorMsg  = "An exception has occurred: ";
                  errorMsg += watercourse.logMessage();

      QMessageBox::warning(this, tr("Warning"), errorMsg.c_str());

      te::common::ProgressManager::getInstance().removeViewer(id);

      delete pViewer;

      return;
    }

  }
  catch(te::common::Exception& e)
  {
    std::string errorMsg  = "Problems in Watercourse operation. ";
                errorMsg += e.what();

    QMessageBox::warning(this, tr("Warning"), errorMsg.c_str());

    te::common::ProgressManager::getInstance().removeViewer(id);

    delete pViewer;

    return;
  }
  catch(...)
  {
    std::string errorMsg  = "An unexpected exception has occurred in Watercourse()!";

    QMessageBox::warning(this, tr("Warning"), errorMsg.c_str());

    te::common::ProgressManager::getInstance().removeViewer(id);

    delete pViewer;

    return;
  }

  te::common::ProgressManager::getInstance().removeViewer(id);

  delete pViewer;

  QMessageBox::information(this, tr("Information"), tr("Watercourse operation done."));

  //create output layer
  emit createLayer(dataType, connInfo);

  //close main window
  emit operationFinished(true);
}
