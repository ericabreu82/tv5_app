/*  Copyright (C) 2001-2009 National Institute For Space Research (INPE) - Brazil.

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
  \file terralib/app/RidgeLines.h
 
  \brief Generates the APP information from ridge lines information.

    - Linhas de Cumeadas: 

    Entrada:
    grade de altimetria (interpolada de curvas de nível no Spring) em UTM,
    grade de declividade gerada no TerraHidro a partir da grade de
    altimetria, drenagem vetorizada gerada no TerraHidro, curvas de nível,
    delimitação de subbacias vetorizada gerada no TerraHidro, buffer da
    delimitação de subbacias gerada pelo TerraView, máscara de declividade
    (tif), máscara de planícies (tif) 

    Saída: máscara de declividade (tif)
    para geração da máscara de planícies, máscara de planícies (tif),
    coordenadas dos picos tem txt (opcional, comentar/descomentar), grade de
    picos (tif), grade de picos isolados (tif) (para detecção de picos de
    morros/montanhas com o programa da Silvia), grade de picos de cumeada
    (tif), grade de APP para picos isolados (tif), segmentos de linhas de
    cumeada (spr) que é convertido no Spring para shp, grade de APP
    extendida para linhas de cumeada (tif) 

    Código: primeiro é gerada a
    máscara de declividade com a função Define_Slope_Mask para gerar a
    máscara de planícies com a função Define_Plains_Mask. Depois, são
    detectados todos os picos na grade com a função Generate_HillTop_vector
    salvando no vetor peaks (os picos podem ser salvos em txt, tem um código
    comentado). A função Define_Isolated_Peaks gera os picos isolados (que
    não são de cumeada) salvando no vetor isolated_peaks para ser utilizado
    na detecção de picos de morros/montanhas (é gerado um tif). Define um
    vetor com os picos de cumeada salvando as linhas de cumeada em um vetor
    (peaks_line2D_vector) e em um lineset (peaks_lineset, para exportar pra
    spr) e também os picos com as bases em um vetor (peaks_line3D_vector),
    sendo também salvo um tif com os picos de cumeada. Em seguida, a função
    Generate_Peaks_Extended_APP gera as novas APPs considerando todos os
    picos de cada linha de cumeada salvando em um tif, cada APP de 1000 m da
    linha de cumeada com um código diferente.
*/

#ifndef __TERRALIB_APP_INTERNAL_RIDGELINES_H
#define __TERRALIB_APP_INTERNAL_RIDGELINES_H

#include "Algorithm.h"

// Terralib
#include <terralib/geometry/Coord2D.h>


// STL
#include <map>
#include <vector>
#include <memory>

namespace te
{
  //forward declarations
  namespace gm { class LineString;  class MultiLineString; class MultiPolygon; }
  namespace rst { class Raster; }

  namespace app
  {
    /*!
      \class RidgeLines
      \brief APP generator by ridge lines.
      \details Generates the APP information from a ridge lines information.
    */
    class APPPLUGINDLLEXPORT RidgeLines : public Algorithm
    {
      public:

        /*!
          \class InputParameters
          \brief RidgeLines input parameters
        */
        class APPPLUGINDLLEXPORT InputParameters : public AlgorithmInputParameters
        {
          public:

            InputParameters();

            ~InputParameters();

            //overload
            void reset() throw(te::common::Exception);

            //overload
            AbstractParameters* clone() const;

          public:

            te::rst::Raster* m_demRasterPtr;      //!< Input raster for dem information

            te::rst::Raster* m_slopeRasterPtr;    //!< Input raster for slope information

            te::rst::Raster* m_slopeMaskRasterPtr;      //!< Input raster for slope mask information

            te::rst::Raster* m_plainMaskRasterPtr;      //!< Input raster for plain mask information

            te::gm::MultiLineString* m_contourLines;    //!< Input level curves

            te::gm::MultiLineString* m_drainage;        //!< Input drainage lines 

            te::gm::MultiPolygon* m_watershedsBuffer;   //!< Input watersheds buffer

            te::gm::MultiPolygon* m_watersheds;         //!< Input watersheds

            double m_slope;

            double m_height;

        };

        /*!
          \class OutputParameters
          \brief RidgeLines output parameters
        */
        class APPPLUGINDLLEXPORT OutputParameters : public AlgorithmOutputParameters
        {
          public:

            OutputParameters();

            OutputParameters(const OutputParameters&);

            ~OutputParameters();

            //overload
            void reset() throw(te::common::Exception);

            //overload
            const OutputParameters& operator=(const OutputParameters& params);

            //overload
            AbstractParameters* clone() const;

        public:
            std::string m_createdOutDSName;                              //!< Output data set name.

            std::string m_createdOutDSType;                              //!< Output data source type.

            std::map<std::string, std::string> m_createdOutInfo;         //!< The necessary information to create the output data.

        };

        typedef struct 
        {
          double x,y,z,b;
          bool visited;
        } Point_T;

        typedef struct 
        {
          int idx;
          double x,y,z,d;
          bool visited;
        } Point_Peak;
        
        RidgeLines();
        
        ~RidgeLines();
       
        //overload
        bool execute(AlgorithmOutputParameters& outputParams) throw(te::common::Exception);
        
        //overload
        void reset() throw(te::common::Exception);
        
        //overload
        bool initialize(const AlgorithmInputParameters& inputParams) throw(te::common::Exception);

      protected:

        std::auto_ptr<te::rst::Raster> createSlopeMask(double threshold_slope);

        std::auto_ptr<te::rst::Raster> createPlainsRasterMask(te::rst::Raster* slopeMask, double area_threshold);

        std::vector<Point_T> generateHillTopVector(te::rst::Raster* plainRasterMask);

        void generatePeaksLine(std::vector<Point_T>& peaksVector, std::vector< std::vector<Point_T> >& peaksLines, std::vector<te::gm::LineString*>& peaksLinesVec);

        std::auto_ptr<te::rst::Raster> generatePeaksExtendedAPP(std::vector< std::vector<Point_T> >& peaksLines, te::rst::Raster* plainRasterMask);

//-------------------------------------------------------------------------------
//		Initialize array flag 
//
//		Flag Matrix : 1 - Drenage
//					  0 - others
//-------------------------------------------------------------------------------
        void constructFlagMatrix(te::rst::Raster* plainRasterMask);

//-------------------------------------------------------------------------------
//		Returns array flag in the original situation
//		Input :
//			numlin,numcol:	size of array
//-------------------------------------------------------------------------------
        void resetFlagMatrix(int numlin, int numcol);

        std::auto_ptr<te::rst::Raster> createDrainageRaster();

        std::vector<Point_T> selectPeaks();

        //-------------------------------------------------------------------------------
//		Algorithm to go top to bottom
//		Input :
//			x,y: Coordinates the hilltop
//			phase	= 2 - states: down the hill and up the hill to check
//					= 5 - only delimit the hill
//					= 6 - down the hill to find the hill 5
//-------------------------------------------------------------------------------
        bool topToDown(int x, int y, int fase, Point_T& pointMin);

//-------------------------------------------------------------------------------
//		Delimit APP hilltop
//		Input :
//			indc,indl: Coordinates the hilltop
//			height: Height limit the elevation to be considered APP
//			TRUE: 	Normal terminal
//			FALSE:  APP does not exist
//-------------------------------------------------------------------------------
        bool hillTopDraw(int indc, int indl, double height, te::rst::Raster* imaout, double SLindex);

        // bresenham
  // http://rosettacode.org/wiki/Bitmap/Bresenham's_line_algorithm#C
        std::vector<te::gm::Coord2D> bresenhamPoints(te::gm::Coord2D gp1, te::gm::Coord2D gp2);

      protected:

        RidgeLines::InputParameters m_inputParameters;       //!< Veredas input execution parameters.

        RidgeLines::OutputParameters* m_outputParametersPtr; //!< Veredas output execution parameters.

        std::vector< std::vector<int> > m_SLgFlag;
    };

  } // end namespace app
}   // end namespace te

#endif //__TERRALIB_APP_INTERNAL_RIDGELINES_H
