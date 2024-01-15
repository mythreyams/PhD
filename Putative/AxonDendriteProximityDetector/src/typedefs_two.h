/****************************************************************************/
/*                                                                          */
/* File:      typedefs.h                                                    */
/*                                                                          */
/* Purpose:   header file for all necessary inculdes and typedefs for files */
/*            associated with the NeuroMorph or CellCount projects          */
/*                                                                          */
/* Author:    Marcel Oberlaender                                            */
/*            Max-Planck-Institute for Neurobiologie                        */
/*            Am Kolpferspitz 18                                            */
/*            D-82152 Martinsried (Munich)                                  */
/*                                                                          */
/* Co-Author: Robert Egger                                                  */
/*            Max-Planck-Institute for Medical Research                     */
/*            Jahnstrasse 19                                                */
/*            D-69120 Heidelberg                                            */
/*                                                                          */
/* EMail:     regger@mpimf-heidelberg.mpg.de                                */
/*                                                                          */
/* History:   17.01.2008                                                    */
/*                                                                          */
/* Remarks:   All rights are reserved by the Max-Planck-Society             */
/*                                                                          */
/****************************************************************************/

#ifndef TYPEDEFTWO
#define TYPEDEFTWO

//#define DEBUG


#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#define PI 3.1415926

#define _CRT_SECURE_NO_DEPRECATE
#define _SCL_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS

//STL includes
#include "string"
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
//#include <map>
#include <cmath>
#include <algorithm>
//#include <complex>
#include <utility>
//#include <ctime>



//ITK includes
#include "itkImage.h"
#include "itkVTKImageImport.h"
#include "itkVTKImageExport.h"
#include "itkNeighborhoodIterator.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSignedDanielssonDistanceMapImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
// // for Voronoi diagram creation & use
// #include "itkPointSet.h"
// #include "itkVoronoiDiagram2D.h"
// #include "itkVoronoiDiagram2DGenerator.h"
// #include "itkMeshSpatialObject.h"
// #include "itkSpatialObjectToImageFilter.h"

//VTK includes
//basics
/*
#include "vtkSmartPointer.h"
#include "vtkMath.h"
#include "vtkPoints.h"
#include "vtkImageData.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkIdTypeArray.h"
#include "vtkTIFFWriter.h"
#include "vtkPolyDataWriter.h"
#include "vtkGenericCell.h"
#include "vtkPolygon.h"
#include "vtkTriangle.h"
#include "vtkTetra.h"
#include "vtkPlane.h"
#include "vtkLine.h"
#include "vtkCylinder.h"
#include "vtkCutter.h"
#include "vtkClipPolyData.h"
#include "vtkTransform.h"
#include "vtkParametricSpline.h"
#include "vtkKochanekSpline.h"
#include "vtkBoundingBox.h"
#include "vtkImageImport.h"
#include "vtkImageExport.h"
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>

//algorithms
#include "vtkContourFilter.h"
#include "vtkMarchingCubes.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkWindowedSincPolyDataFilter.h"
#include "vtkSurfaceReconstructionFilter.h"
#include "vtkCellLocator.h"
#include "vtkSelectPolyData.h"
#include "vtkExtractPolyDataGeometry.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkAppendFilter.h"
#include "vtkAppendPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkHull.h"
#include "vtkPointsProjectedHull.h"
#include "vtkDelaunay2D.h"
#include "vtkDelaunay3D.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkButterflySubdivisionFilter.h"
#include "vtkLoopSubdivisionFilter.h"
#include "vtkSelectEnclosedPoints.h"*/


  typedef float 									CalcPixelType;
  typedef unsigned char									PixelType;
  typedef itk::Image< PixelType, 3 >							ImageType;
  typedef itk::Image< PixelType, 2 >							Image2DType;
  typedef itk::Image< CalcPixelType, 3 >						CalcImageType;
  typedef itk::NeighborhoodIterator< ImageType >					SegNeighborhoodIteratorType;
  typedef itk::NeighborhoodIterator< CalcImageType >					CalcNeighborhoodIteratorType;
  typedef itk::ShapedNeighborhoodIterator< ImageType >					ShapedNeighborhoodIteratorType;
  typedef itk::ShapedNeighborhoodIterator< CalcImageType >				ShapedCalcNeighborhoodIteratorType;
  typedef itk::ImageRegionIterator< ImageType >						IteratorType2;
  typedef itk::ImageRegionConstIterator< ImageType >					ConstIteratorType;
  typedef itk::ImageRegionIterator< Image2DType >					Iterator2DType;
  typedef itk::ImageRegionIterator< CalcImageType >					CalcIteratorType;
  typedef itk::ImageRegionConstIterator< CalcImageType >				ConstCalcIteratorType;
  typedef itk::ImageRegionIteratorWithIndex< ImageType >				IndexIteratorType;
  typedef itk::SignedDanielssonDistanceMapImageFilter< ImageType, CalcImageType >	DistanceMapImageFilterType;
  typedef itk::DanielssonDistanceMapImageFilter< ImageType, CalcImageType >		DistanceMapImageFilterType2;
  typedef itk::VTKImageImport< ImageType >						VTK2ITKImageImportType;
  typedef itk::VTKImageExport< CalcImageType >						ITK2VTKCalcImageExportType;
  typedef std::vector< SegNeighborhoodIteratorType::OffsetType >			NeighborhoodOffsetVectorType;
  typedef std::vector< ShapedNeighborhoodIteratorType::OffsetType >			ShapedNeighborhoodOffsetVectorType;
  
  
//   typedef itk::PointSet< PixelType, 2 >							itkPointSetType;
//   typedef itkPointSetType::PointType							itkPointType;
//   typedef itk::VoronoiDiagram2D< float >						VoronoiDiagramType;
//   typedef itk::VoronoiDiagram2DGenerator< float >					VoronoiDiagramGeneratorType;
//   typedef itk::MeshSpatialObject< VoronoiDiagramType >					VoronoiDiagramSpatialObjectType;
//   typedef itk::SpatialObjectToImageFilter< VoronoiDiagramSpatialObjectType, Image2DType >	VoronoiDiagramToImageFilterType;
/*
  typedef vtkSmartPointer< vtkPoints >		PointsPointerType;
  typedef vtkSmartPointer< vtkDataArray >	DataArrayPointerType;
  typedef vtkSmartPointer< vtkFloatArray >	FloatArrayPointerType;
  typedef vtkSmartPointer< vtkDoubleArray >	DoubleArrayPointerType;
  typedef vtkSmartPointer< vtkIdTypeArray >	IdTypeArrayPointerType;
  typedef vtkSmartPointer< vtkPolygon >		PolygonPointerType;
  typedef vtkSmartPointer< vtkTriangle >	TrianglePointerType;
  typedef vtkSmartPointer< vtkTetra >		TetraPointerType;
  typedef vtkSmartPointer< vtkIdList >		IdListPointerType;
  typedef vtkSmartPointer< vtkPlane >		PlanePointerType;
  typedef vtkSmartPointer< vtkLine >		LinePointerType;
  typedef vtkSmartPointer< vtkCylinder >	CylinderPointerType;
  typedef vtkSmartPointer< vtkCutter >		CutterPointerType;
  typedef vtkSmartPointer< vtkClipPolyData >	ClipPolyDataPointerType;
  typedef vtkSmartPointer< vtkTransform >	TransformPointerType;
  typedef vtkSmartPointer< vtkMatrix4x4 >	HomogeneousMatrixPointerType;
  typedef vtkSmartPointer< vtkParametricSpline >	ParametricSplinePointerType;
  typedef vtkSmartPointer< vtkKochanekSpline >	KochanekSplinePointerType;
  typedef vtkSmartPointer< vtkBoundingBox >	BoundingBoxPointerType;
  typedef vtkSmartPointer< vtkImageData >	ImageDataPointerType;
  typedef vtkSmartPointer< vtkPolyData >	PolyDataPointerType;
  typedef vtkSmartPointer< vtkUnstructuredGrid >	UnstructuredGridPointerType;
  typedef vtkSmartPointer< vtkTIFFWriter >	TiffWriterPointerType;
  typedef vtkSmartPointer< vtkPolyDataWriter >	PolyDataWriterPointerType;
  typedef vtkSmartPointer< vtkGenericCell >	GenericCellPointerType;
  typedef vtkSmartPointer< vtkCell >		CellPointerType;
  typedef vtkSmartPointer< vtkImageImport >	ITK2VTKImageImportPointerType;
  typedef vtkSmartPointer< vtkImageExport >	VTK2ITKImageExportPointerType;
  
  typedef vtkSmartPointer< vtkMarchingCubes >			MarchingCubesPointerType;
  typedef vtkSmartPointer< vtkContourFilter >			ContourFilterPointerType;
  typedef vtkSmartPointer< vtkSmoothPolyDataFilter >		AveragePolyDataFilterType;
  typedef vtkSmartPointer< vtkWindowedSincPolyDataFilter >	LowpassPolyDataFilterType;
  typedef vtkSmartPointer< vtkSurfaceReconstructionFilter >	SurfaceReconstructionFilterType;
  typedef vtkSmartPointer< vtkCellLocator >			CellLocatorPointerType;
  typedef vtkSmartPointer< vtkSelectPolyData >			SelectPolyDataPointerType;
  typedef vtkSmartPointer< vtkExtractPolyDataGeometry >		ExtractPolyDataGeometryPointerType;
  typedef vtkSmartPointer< vtkTransformPolyDataFilter >		TransformFilterType;
  typedef vtkSmartPointer< vtkAppendFilter >			AppendFilterPointerType;
  typedef vtkSmartPointer< vtkAppendPolyData >			AppendPolyDataPointerType;
  typedef vtkSmartPointer< vtkPolyDataNormals >			PolyDataNormalsPointerType;
  typedef vtkSmartPointer< vtkHull >				ConvexHullFilterPointerType;
  typedef vtkSmartPointer< vtkPointsProjectedHull >		ConvexHull2DFilterPointerType;
  typedef vtkSmartPointer< vtkDelaunay2D >			Delaunay2DFilterPointerType;
  typedef vtkSmartPointer< vtkDelaunay3D >			Delaunay3DFilterPointerType;
  typedef vtkSmartPointer< vtkDataSetSurfaceFilter >		DataSetSurfaceFilterPointerType;
  typedef vtkSmartPointer< vtkButterflySubdivisionFilter >	MeshRefinementFilterPointerType;
  typedef vtkSmartPointer< vtkLoopSubdivisionFilter >		MeshRefinementFilter2PointerType;
  typedef vtkSmartPointer< vtkSelectEnclosedPoints >		SelectEnclosedPointsFilterType;*/


  extern float averageSomaRadius;
  extern float zScale;
  extern unsigned long BINSIZE;
  extern unsigned long BOXSIZE;
  
  
  //celltypes for automatic dendrite detection
  #define SUPRA			1
  #define GRAN			2
  #define INFRA			3
  
  //Label IDs for Amira Spatial Graphs
  #define Neuron		2
  #define Dendrite		3
  #define ApicalDendrite	4
  #define BasalDendrite		5
  #define Axon			6
  #define Soma			7
  #define Landmark		8
  #define Pia			9
  #define WhiteMatter		48
  #define Vessel		10
  #define Barrel		11
  #define ZAxis			50
  #define aRow			12
  #define A1			13
  #define A2			14
  #define A3			15
  #define A4			16
  #define bRow			17
  #define B1			18
  #define B2			19
  #define B3			20
  #define B4			21
  #define cRow			22
  #define C1			23
  #define C2			24
  #define C3			25
  #define C4			26
  #define C5			27
  #define C6			28
  #define dRow			29
  #define D1			30
  #define D2			31
  #define D3			32
  #define D4			33
  #define D5			34
  #define D6			35
  #define eRow			36
  #define E1			37
  #define E2			38
  #define E3			39
  #define E4			40
  #define E5			41
  #define E6			42
  #define greekRow		43
  #define Alpha			44
  #define Beta			45
  #define Gamma			46
  #define Delta			47
  #define Septum		0
//168fines
#endif