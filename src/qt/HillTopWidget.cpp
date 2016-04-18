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
  \file terralib/app/qt/HillTopWidget.cpp

  \brief A widget used to access the ridge lines method.
*/

#include "HillTopWidget.h"
#include "ui_HillTopWidgetForm.h"

#include "../core/HillTop.h"

// TerraLib
#include <terralib/common/progress/ProgressManager.h>
#include <terralib/common/Exception.h>
#include <terralib/dataaccess/utils/Utils.h>
#include <terralib/geometry/GeometryProperty.h>
#include <terralib/geometry/MultiLineString.h>
#include <terralib/geometry/MultiPolygon.h>
#include <terralib/raster/Grid.h>
#include <terralib/qt/widgets/progress/ProgressViewerDialog.h>
#include <terralib/qt/widgets/Utils.h>

// Qt
#include <QFileDialog>
#include <QMessageBox>
#include <QValidator>

Q_DECLARE_METATYPE(te::map::AbstractLayerPtr);

te::app::HillTopWidget::HillTopWidget(QWidget* parent)
  : QWidget(parent),
    m_ui(new Ui::HillTopWidgetForm)
{
// add controls
  m_ui->setupUi(this);

  m_ui->m_heightLineEdit->setValidator(new QDoubleValidator(this));

  connect(m_ui->m_fileToolButton, SIGNAL(clicked()), this, SLOT(onFileToolButtonClicked()));
}

te::app::HillTopWidget::~HillTopWidget()
{
}

void te::app::HillTopWidget::setLayerList(std::list<te::map::AbstractLayerPtr>& layerList)
{
  m_ui->m_layersComboBox->clear();

  std::list<te::map::AbstractLayerPtr>::iterator it = layerList.begin();

  while(it != layerList.end())
  {
    te::map::AbstractLayerPtr l = *it;

    std::auto_ptr<te::da::DataSetType> dsType = l->getSchema();

    if(dsType->hasRaster()) {
      m_ui->m_layersComboBox->addItem(l->getTitle().c_str(), QVariant::fromValue(l));
      m_ui->m_slopeLayerComboBox->addItem(l->getTitle().c_str(), QVariant::fromValue(l));
    } else if(dsType->hasGeom()) {
      m_ui->m_contourLayerComboBox->addItem(l->getTitle().c_str(), QVariant::fromValue(l));
      m_ui->m_drainageLayerComboBox->addItem(l->getTitle().c_str(), QVariant::fromValue(l));
      m_ui->m_watershedsLayerComboBox->addItem(l->getTitle().c_str(), QVariant::fromValue(l));
      m_ui->m_watershedsBufferLayerComboBox->addItem(l->getTitle().c_str(), QVariant::fromValue(l));
    }
    ++it;
  }
}

void te::app::HillTopWidget::onFileToolButtonClicked()
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

te::gm::MultiLineString* te::app::HillTopWidget::layerToMultiLine(QComboBox* layerComboBox, std::string layerTitle) {

    te::gm::MultiLineString* lines;

    if(layerComboBox->currentText().isEmpty())
    {
        QMessageBox::warning(this, tr("Warning"), tr((layerTitle + "Layer not selected.").c_str()));
        return lines;
    }

    int idx = layerComboBox->currentIndex();
    QVariant varLayer = layerComboBox->itemData(idx, Qt::UserRole);
    te::map::AbstractLayerPtr layer = varLayer.value<te::map::AbstractLayerPtr>();

    if(!layer.get())
    {
        QMessageBox::warning(this, tr("Warning"), tr(("Invalid " + layerTitle + " layer selected.").c_str()));
        return lines;
    }

    std::auto_ptr<te::da::DataSet> ds = layer->getData();
    std::auto_ptr<te::da::DataSetType> dsType = layer->getSchema();

    if(!ds.get())
    {
        QMessageBox::warning(this, tr("Warning"), tr(("Invalid " + layerTitle + " input layer selected.").c_str()));
        return lines;
    }

    std::size_t gpos = te::da::GetFirstPropertyPos(ds.get(), te::dt::GEOMETRY_TYPE);
    te::gm::GeometryProperty* geomProp = te::da::GetFirstGeomProperty(dsType.get());

    lines = new te::gm::MultiLineString(0, te::gm::MultiLineStringType, geomProp->getSRID());

    ds->moveBeforeFirst();

    //get geometries
    while (ds->moveNext())
    {
      std::auto_ptr<te::gm::Geometry> g(ds->getGeometry(gpos));

      if (g->getGeomTypeId() == te::gm::MultiLineStringType)
      {
        te::gm::MultiLineString* mls = dynamic_cast<te::gm::MultiLineString*>(g.release());

        te::gm::LineString* ls = dynamic_cast<te::gm::LineString*>(mls->getGeometryN(0));

        lines->add(ls);
      }
      else if (g->getGeomTypeId() == te::gm::LineStringType)
      {
        te::gm::LineString* ls = dynamic_cast<te::gm::LineString*>(g.release());

        lines->add(ls);
      }
      else if (g->getGeomTypeId() == te::gm::LineStringZType)
      {
        te::gm::LineString* ls = dynamic_cast<te::gm::LineString*>(g.release());

        lines->add(ls);
      }
    }

    return lines;
}

te::gm::MultiPolygon* te::app::HillTopWidget::layerToMultiPolygon(QComboBox* layerComboBox, std::string layerTitle) {

    te::gm::MultiPolygon* polys;

    if(layerComboBox->currentText().isEmpty())
    {
        QMessageBox::warning(this, tr("Warning"), tr((layerTitle + "Layer not selected.").c_str()));
        return polys;
    }

    int idx = layerComboBox->currentIndex();
    QVariant varLayer = layerComboBox->itemData(idx, Qt::UserRole);
    te::map::AbstractLayerPtr layer = varLayer.value<te::map::AbstractLayerPtr>();

    if(!layer.get())
    {
        QMessageBox::warning(this, tr("Warning"), tr(("Invalid " + layerTitle + " layer selected.").c_str()));
        return polys;
    }

    std::auto_ptr<te::da::DataSet> ds = layer->getData();
    std::auto_ptr<te::da::DataSetType> dsType = layer->getSchema();

    if(!ds.get())
    {
        QMessageBox::warning(this, tr("Warning"), tr(("Invalid " + layerTitle + " input layer selected.").c_str()));
        return polys;
    }

    std::size_t gpos = te::da::GetFirstPropertyPos(ds.get(), te::dt::GEOMETRY_TYPE);
    te::gm::GeometryProperty* geomProp = te::da::GetFirstGeomProperty(dsType.get());

    polys = new te::gm::MultiPolygon(0, te::gm::MultiPolygonType, geomProp->getSRID());

    ds->moveBeforeFirst();

    //get geometries
    while (ds->moveNext())
    {
      std::auto_ptr<te::gm::Geometry> g(ds->getGeometry(gpos));

      if (g->getGeomTypeId() == te::gm::MultiPolygonType)
      {
        te::gm::MultiPolygon* mp = dynamic_cast<te::gm::MultiPolygon*>(g.release());

        te::gm::Polygon* p = dynamic_cast<te::gm::Polygon*>(mp->getGeometryN(0));

        polys->add(p);
      }
      else if (g->getGeomTypeId() == te::gm::PolygonType)
      {
        te::gm::Polygon* p = dynamic_cast<te::gm::Polygon*>(g.release());

        polys->add(p);
      }
    }

    return polys;
}


void te::app::HillTopWidget::execute()
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

  //get height
  if(m_ui->m_heightLineEdit->text().isEmpty())
  {
    QMessageBox::warning(this, tr("Warning"), tr("Invalid height value."));
    return;
  }

  double heightValue = m_ui->m_heightLineEdit->text().toDouble();

  //get slope
  if(m_ui->m_slopeLineEdit->text().isEmpty())
  {
    QMessageBox::warning(this, tr("Warning"), tr("Invalid slope value."));
    return;
  }

  double slopeValue = m_ui->m_slopeLineEdit->text().toDouble();


  //get raster
  if(m_ui->m_layersComboBox->currentText().isEmpty())
  {
    QMessageBox::warning(this, tr("Warning"), tr("DEM Layer not selected."));
    return;
  }

  if(m_ui->m_slopeLayerComboBox->currentText().isEmpty())
  {
    QMessageBox::warning(this, tr("Warning"), tr("Slope Layer not selected."));
    return;
  }

  int idx = m_ui->m_layersComboBox->currentIndex();
  QVariant varLayer = m_ui->m_layersComboBox->itemData(idx, Qt::UserRole);
  te::map::AbstractLayerPtr layer = varLayer.value<te::map::AbstractLayerPtr>();

  if(!layer.get())
  {
    QMessageBox::warning(this, tr("Warning"), tr("Invalid DEM layer selected."));
    return;
  }

  std::auto_ptr<te::da::DataSet> ds = layer->getData();
  std::size_t rpos = te::da::GetFirstPropertyPos(ds.get(), te::dt::RASTER_TYPE);
  std::auto_ptr<te::rst::Raster> rst = ds->getRaster(rpos);

  if(!rst.get())
  {
    QMessageBox::warning(this, tr("Warning"), tr("Selected DEM layer does not have raster representation."));
    return;
  }

  te::rst::Raster* rasterSRTM = rst.get();

  if (rasterSRTM->getSRID() == 0) {
      // Gets SRID from layer object
      rasterSRTM->getGrid()->setSRID(layer->getSRID());
  }

  // get input geometries
  te::gm::MultiLineString* contourLines = layerToMultiLine(m_ui->m_contourLayerComboBox, "Contour Lines");
  te::gm::MultiLineString* drainage = layerToMultiLine(m_ui->m_drainageLayerComboBox, "Drainage");

  te::gm::MultiPolygon* watersheds = layerToMultiPolygon(m_ui->m_watershedsLayerComboBox, "Watershed");
  te::gm::MultiPolygon* watershedsBuffer = layerToMultiPolygon(m_ui->m_watershedsBufferLayerComboBox, "Watershed Buffer");

  if(contourLines->isEmpty() || drainage->isEmpty() || watersheds->isEmpty() || watershedsBuffer->isEmpty())
      return;

  //run operation
  te::qt::widgets::ProgressViewerDialog* pViewer = new te::qt::widgets::ProgressViewerDialog(this);
  int id = te::common::ProgressManager::getInstance().addViewer(pViewer);

  try
  {
    bool executeok = false;
    bool initok = false;

    // create ridge lines algorithm parameters
    te::app::HillTop::InputParameters inputParameters;
    inputParameters.m_demRasterPtr = rasterSRTM;
    inputParameters.m_contourLines = contourLines;
    inputParameters.m_drainage = drainage;
    inputParameters.m_watersheds = watersheds;
    inputParameters.m_watershedsBuffer = watershedsBuffer;
    inputParameters.m_height = heightValue;
    inputParameters.m_slope = slopeValue;

    if (m_ui->m_layerSlopeCheckBox->isChecked())
    {
      int idxSlope = m_ui->m_slopeLayerComboBox->currentIndex();
      QVariant varSlopeLayer = m_ui->m_slopeLayerComboBox->itemData(idxSlope, Qt::UserRole);
      te::map::AbstractLayerPtr slopeLayer = varSlopeLayer.value<te::map::AbstractLayerPtr>();

      if (!slopeLayer.get())
      {
        QMessageBox::warning(this, tr("Warning"), tr("Invalid Slope layer selected."));
        return;
      }

      std::auto_ptr<te::da::DataSet> slopeDs = slopeLayer->getData();
      std::size_t slopePos = te::da::GetFirstPropertyPos(slopeDs.get(), te::dt::RASTER_TYPE);
      std::auto_ptr<te::rst::Raster> slopeRst = slopeDs->getRaster(slopePos);

      if (!slopeRst.get())
      {
        QMessageBox::warning(this, tr("Warning"), tr("Selected Slope layer does not have raster representation."));
        return;
      }

      te::rst::Raster* rasterSlope = slopeRst.release();

      if (rasterSlope->getSRID() == 0) 
      {
        // Gets SRID from layer object
        rasterSlope->getGrid()->setSRID(slopeLayer->getSRID());
      }

      inputParameters.m_slopeRasterPtr = rasterSlope;
    }

    te::app::HillTop::OutputParameters outputParameters;
    outputParameters.m_createdOutDSName = dataSetName;
    outputParameters.m_createdOutInfo = connInfo;
    outputParameters.m_createdOutDSType = dataType;

    // execute the algorithm
    te::app::HillTop hillTop;

    initok = hillTop.initialize(inputParameters);

    if(initok)
      executeok = hillTop.execute(outputParameters);

    if(!executeok)
    {
      std::string errorMsg  = "An exception has occurred: ";
                  errorMsg += hillTop.logMessage();

      QMessageBox::warning(this, tr("Warning"), errorMsg.c_str());

      te::common::ProgressManager::getInstance().removeViewer(id);

      delete pViewer;
      delete contourLines;
      delete drainage;
      delete watersheds;
      delete watershedsBuffer;

      return;
    }

  }
  catch(te::common::Exception& e)
  {
    std::string errorMsg  = "Problems in Ridge Lines operation. ";
                errorMsg += e.what();

    QMessageBox::warning(this, tr("Warning"), errorMsg.c_str());

    te::common::ProgressManager::getInstance().removeViewer(id);

    delete pViewer;
    delete contourLines;
    delete drainage;
    delete watersheds;
    delete watershedsBuffer;

    return;
  }
  catch(...)
  {
    std::string errorMsg  = "An unexpected exception has occurred in HillTop()!";

    QMessageBox::warning(this, tr("Warning"), errorMsg.c_str());

    te::common::ProgressManager::getInstance().removeViewer(id);

    delete pViewer;
    delete contourLines;
    delete drainage;
    delete watersheds;
    delete watershedsBuffer;

    return;
  }

  te::common::ProgressManager::getInstance().removeViewer(id);

  delete pViewer;
  delete contourLines;
  delete drainage;
  delete watersheds;
  delete watershedsBuffer;

  QMessageBox::information(this, tr("Information"), tr("HillTop operation done."));

  //create output layer
  emit createLayer(dataType, connInfo);

  //close main window
  emit operationFinished(true);
}
