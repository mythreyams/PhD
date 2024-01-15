/****************************************************************************/
/*                                                                          */
/* Program:                                                                 */
/*                                                                          */
/* File:      amiraReader.h                                                 */
/*                                                                          */
/* Purpose:   interface class between Amira SpatialGraph, Surface and       */
/*            other files and the internal spatial graph data structure     */
/*                                                                          */
/* Author:    Robert Egger                                                  */
/*            Max-Planck-Florida Institut                                   */
/*            5353 Parkside Drive                                           */
/*            Jupiter, Florida 33458                                        */
/*            USA                                                           */
/*                                                                          */
/* EMail:     Robert.Egger@maxplanckflorida.org                             */
/*                                                                          */
/* History:   22.12.2010                                                    */
/*                                                                          */
/* Remarks:   All rights are reserved by the Max-Planck-Society             */
/*                                                                          */
/****************************************************************************/

#pragma once
#include "typedefs.h"
#include "typedefs_two.h"
#include "basics.h"

#ifndef AMIRAREADER
#define AMIRAREADER


class Reader
{
	public:
		Reader(const char * filename);
		Reader(const char * filename, const char * outputFilename);
		~Reader();
		
		// default method for reading Amira SpatialGraph files
		void readSpatialGraphFile(bool applyTransform);
		// default method for writing Amira SpatialGraph files
		void writeSpatialGraphFile();
                // default method for writing Amira SpatialGraph files with proximity detection labels
                void writeSpatialGraphFile2();
		// alternate method for writing Amira SpatialGraph files
		void writeSpatialGraphFileFromEdges();
                
                void writeSpatialGraphFile(int Axon_ID, int Dendrite_ID);
		
		// default method for reading NeuroConv hoc files
		void readHocFile();
		// method for writing NeuroConv hoc files
		void writeHocFile();
		
		// method for writing NeuroConv hoc files used during
		// registration. writes separate files for neuron
		// morphology and anatomical landmarks
		//void writeSeparateHocFiles();
		
		// default method for writing Amira Surface files
		//void writeAmiraSurfaceFile(PolyDataPointerType triangleData);
		
		// default method for reading Amira Surface files
		//PolyDataPointerType readAmiraSurfaceFile();
		
		// default method for reading Amira Landmark files
		//PointsPointerType readLandmarkFile();
		
		// method for reading Amira Landmark files representing
		// closed contours and converting them into closed contour
		// graphs (i.e., PolyData format)
		//void readLandmarkFile(bool applyTransform);
		
		// default method for writing Amira Landmark files
		//void writeLandmarkFile(PointsPointerType pts);
		
		// default method for reading Amira Vectorfield files
		
		//ImageDataPointerType readVectorField();
		// default method for writing Amira Vectorfield files
		//void writeVectorField(ImageDataPointerType field);
		
		AmiraSpatialGraph * getSpatialGraph() { return inputSpatialGraph; }
		void setSpatialGraph(AmiraSpatialGraph * outputSpatialGraph) { inputSpatialGraph = outputSpatialGraph; }
		
	private:
		const char * inputFilename;
		const char * outputFilename;
		std::string letters;
		std::string numbers;
		std::string signs;
		std::string otherChars;
		std::string whitespace;
		std::list< int > barrelLabels;
		std::map< int, const char * > int2Labels;
		std::list< const char * > hocLabels;
		std::list< const char * > hocNeuronLabels;
		std::list< const char * > hocLandmarkLabels;
		
		AmiraSpatialGraph * inputSpatialGraph;
		
		void initializeConstants();
};



#endif
