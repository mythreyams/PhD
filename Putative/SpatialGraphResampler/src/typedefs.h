/****************************************************************************/
/*                                                                          */
/* File:      typedefs.h                                                    */
/*                                                                          */
/* Purpose:   header file for all necessary inculdes and typedefs for files */
/*            associated with the NeuroMorph project                        */
/*                                                                          */
/* Author:    Marcel Oberlaender                                            */
/*            Max-Planck-Institut fuer medizinische Forschung               */
/*            Jahnstrasse 29                                                */
/*            D-69120 Heidelberg                                            */
/*                                                                          */
/*                                                                          */
/* EMail: marcel.oberlaender@mpimf-heidelberg.mpg.de                        */
/*                                                                          */
/* History:   09.02.2006                                                    */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/
#ifndef TYPEDEF
#define TYPEDEF

//#define DEBUG


#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#define PI 3.1415926

#include "string"
#include <sstream>
// #include <string.h>
//#include <direct.h>
#include <stdio.h>
#include "unistd.h"
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <iomanip>
// #include "fann.h"
// #include "floatfann.h"

//#include "errorcode.h"
#include "itkImage.h"
#include <iostream>
#include <vector>
#include <queue>
// for Reading and Writing Images
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
// for simple Filters
#include "itkBinaryThresholdImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "itkInvertIntensityImageFilter.h"
// for Morphology Filters
#include "itkGrayscaleErodeImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkSubtractImageFilter.h"
#include "itkRGBPixel.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkPasteImageFilter.h"
// for Reading and Writing of Stacks
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
#include "itkNumericSeriesFileNames.h"
#include "itkTIFFImageIO.h"
#include "itkBMPImageIO.h"
#include "itkPNGImageIO.h"
#include "itkRGBPixel.h"
#include "itkRegionOfInterestImageFilter.h"
// for Iterators
#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkRescaleIntensityImageFilter.h"
// for Histograms
#include "itkScalarImageToListAdaptor.h"
#include "itkListSampleToHistogramGenerator.h"
#include "itkCastImageFilter.h"
#include "itkJointDomainImageToListAdaptor.h"
#include "itkImageToListAdaptor.h"
#include "itkScalarToArrayCastImageFilter.h"
#include "itkListSample.h"
#include "itkSubsample.h"
#include "itkVector.h"
#include "itkKdTree.h"
#include "itkKdTreeGenerator.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkEuclideanDistance.h"
#include "itkMeanCalculator.h"
#include "itkCovarianceCalculator.h"
#include "itkOtsuThresholdImageCalculator.h"
//-------------------------------------------------------------------------------------------------------------------
//Type Definitions for 3D Grayvalue Images and Iterators
//-------------------------------------------------------------------------------------------------------------------
  
  typedef unsigned char								PixelType;
  typedef unsigned int								ObjectPixelType;
  typedef float 								CalcPixelType;

  typedef itk::Image< PixelType, 3 >						ImageType;
  typedef itk::Image< ObjectPixelType, 3 >					ObjectImageType;  
  typedef itk::Image< PixelType, 2 >						Image2DType;
  typedef itk::Image<CalcPixelType, 2 >						CalcImage2DType;

  typedef itk::ImageFileReader< ImageType >					ReaderType;
  typedef itk::ImageSeriesReader< ImageType >					SeriesReaderType;
  typedef itk::NumericSeriesFileNames						NameGeneratorType;
  typedef itk::ImageFileWriter< ImageType >					WriterType;
  typedef itk::ImageSeriesWriter< ImageType, Image2DType >			Writer2DType;
  typedef itk::ImageFileWriter< Image2DType >					Single2DWriterType;
  typedef itk::RescaleIntensityImageFilter< CalcImage2DType,Image2DType >       RescalerCalcType; 

  
  typedef itk::RelabelComponentImageFilter< ObjectImageType, ObjectImageType >	RelabelType;
  typedef itk::BinaryThresholdImageFilter< ObjectImageType, ImageType >		ObjectFilterType;
  typedef itk::BinaryThresholdImageFilter< ImageType, ImageType >		BinaryFilterType;
  typedef itk::ConnectedComponentImageFilter< ImageType, ObjectImageType >	ConnectedFilterType;			  

  typedef itk::ImageRegionConstIterator< ImageType >                            ConstIteratorType;
  typedef itk::ImageRegionIterator< ObjectImageType >				IteratorType;

  typedef itk::ConstNeighborhoodIterator< ImageType >				ConstNeighborhoodIteratorType;
  typedef itk::ImageRegionIterator< ImageType>					IteratorType2;
  typedef itk::NeighborhoodIterator< ImageType >				SegNeighborhoodIteratorType;
  typedef itk::NeighborhoodIterator< ObjectImageType >				NeighborhoodIteratorType;
  typedef itk::ImageRegionIterator< Image2DType >				Iterator2DType;
  typedef itk::ImageRegionConstIterator< Image2DType >				ConstIterator2DType;
  typedef itk::ImageLinearIteratorWithIndex< ImageType >			LineIteratorType;

  typedef itk::RescaleIntensityImageFilter< ObjectImageType, ImageType > 	RescaleFilterType;
  typedef itk::SigmoidImageFilter< ImageType, ImageType >				SigmoidFilterType;

  typedef itk::BinaryBallStructuringElement< PixelType, 3 >					StructuringElementType;
  typedef itk::GrayscaleErodeImageFilter< ImageType, ImageType, StructuringElementType >  	ErodeFilterType;
  typedef itk::GrayscaleDilateImageFilter< ImageType, ImageType, StructuringElementType > 	DilateFilterType;
  typedef itk::PasteImageFilter< ImageType >							PasteFilterType;
  typedef itk::RegionOfInterestImageFilter< ImageType, ImageType >				InterestFilterType;
//   typedef itk::InvertIntensityImageFilter< Image2DType >					InvertFilterType;
  
  typedef float MeasurementVectorBase;
  typedef itk::Vector< MeasurementVectorBase, 6 >						MeasurementVectorType;	
  typedef itk::Vector< unsigned char, 1 >							MeasurementHistogramType;

  typedef itk::Statistics::ListSample< MeasurementHistogramType >				HistogramSampleType;
  typedef itk::Statistics::MeanCalculator< HistogramSampleType >				MeanAlgorithmType;
  typedef itk::Statistics::CovarianceCalculator< HistogramSampleType >				CovarianceAlgorithmType;
  typedef itk::MinimumMaximumImageCalculator< ImageType >					MinMaxAlgorithmType;
  typedef itk::OtsuThresholdImageCalculator<Image2DType>                                          OtsuThresholdType;

  typedef enum{red, green, blue, yellow, white} ColorTagType;

//stats regarding image resolution
#define CONFOCAL

#ifdef CONFOCAL
  #define XYSAMPLING  	0.092  //confocal
#else
  #define XYSAMPLING  	0.184    //brightfield
#endif

#define ZSAMPLING 	 0.5

//length of the double array used to represent points  
#define ARRAY_LENGTH 32
  
//VoxelFeatures of the Measurement Vector Type
#define X_COORD			0
#define Y_COORD			1
#define Z_COORD			2
#define IDENTIFIER		3
#define SURFACE			4
#define AVERAGE_STD_DEV		5
#define LOCAL_BRIGHTNESS	6
#define THRESHOLD		7
#define LENGTH_UNTIL_HERE	8
#define IS_BOUTON		9
#define IS_TRUE_BOUTON		10
#define LOCAL_SIGMA		11
#define IS_PAIRED		12
#define IS_USER1_BOUTON		13
#define IS_USER2_BOUTON		14
#define IS_USER3_BOUTON		15
#define CONFIDENCE		16
#define AVERAGE_SURFACE         17
  
#define LOCAL_BRIGHTNESS_RAD2	18
#define LOCAL_BRIGHTNESS_RAD3	19
#define LOCAL_SIGMA_RAD3	20
#define AVG_BRIGHT_RANGE3	21
#define AVG_BRIGHT_RANGE50	22
#define AVG_RAD_RANGE3		23
#define AVG_RAD_RANGE4		24
#define AVG_RAD_RANGE5		25
  
#define VALID_FOR_NN		26
#define AVG_SURF                27
#define AVG_BRIGHT              28  
#define MEASURED_DIAMETER       29  
  

  
//Edge directions
#define FORWARD			true
#define BACKWARD		false
  
  
#define AXON_ID_1       15
#define AXON_ID_2       16
#define AXON_ID_3       17
#define UNKNOWN_ID      18
#define DENDRITE_ID_1   8
#define DENDRITE_ID_2   9
#define DENDRITE_ID_3   10  
  
// #define PRINT_TRAINING_DATA	1
// #define TRAIN_NETWORK		

  
typedef struct{
  float coords[3];
  float magnitude;
  
}VECTOR;
  
  
  
  
  

#endif
